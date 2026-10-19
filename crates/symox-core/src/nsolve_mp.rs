//! Arbitrary precision real and complex numerical root finding.
//!
//! Jacobians are differentiated symbolically. Numerical evaluation, Newton
//! steps, pivoting and overdetermined least squares all retain native binary
//! precision; no machine float or external numerical package participates.

use std::cmp::Ordering;

use rustc_hash::FxHashMap;

use crate::bigfloat::BigFloat;
use crate::evalf_mp::{self, MpComplex};
use crate::expr::{Expr, Kind, RelOp};

pub struct Options {
    /// Requested result precision in bits. Iterations use 64 guard bits.
    pub precision: u32,
    /// Positive absolute tolerance. Verification checks squared residuals.
    pub tolerance: BigFloat,
    pub maxsteps: usize,
    pub verify: bool,
}

fn zero() -> MpComplex { MpComplex::real(BigFloat::zero()) }

fn equation(e: &Expr) -> Result<Expr, String> {
    match e.kind() {
        Kind::Rel(RelOp::Eq, left, right) => Ok(crate::ops::sub(left, right)),
        Kind::Rel(..) | Kind::Bool(_) => Err("nsolve requires expressions or equalities".into()),
        _ => Ok(e.clone()),
    }
}

fn check(options: &Options) -> Result<u32, String> {
    if options.precision < 2 || options.precision > 40_000 || options.maxsteps == 0 || options.tolerance.sign() <= 0 {
        return Err("precision, tolerance and maxsteps must be positive and within supported bounds".into());
    }
    Ok(options.precision + 64)
}

/// Exact symbol identities (including assumptions) are respected. Rebuilding
/// the numerical expression also lets existing MP function kernels participate.
fn evaluate(expressions: &[Expr], variables: &[Expr], values: &[MpComplex], precision: u32) -> Result<Vec<MpComplex>, String> {
    let replacements: FxHashMap<Expr, Expr> = variables.iter().cloned()
        .zip(values.iter().cloned().map(|value| evalf_mp::to_expr(value, precision))).collect();
    expressions.iter().map(|expression| {
        let numerical = crate::subs::xreplace(expression, &replacements);
        evalf_mp::eval_mp(&numerical, precision)
            .ok_or_else(|| "Could not numerically evaluate the equation or Jacobian at this point".into())
    }).collect()
}

fn max_norm(values: &[MpComplex], precision: u32) -> BigFloat {
    values.iter().map(|value| value.abs(precision)).max_by(|a, b| a.cmp(b)).unwrap_or_else(BigFloat::zero)
}

fn inner_product(a: &[MpComplex], b: &[MpComplex], precision: u32) -> MpComplex {
    a.iter().zip(b).fold(zero(), |sum, (left, right)| sum.add(&left.conjugate().mul(right, precision), precision))
}

fn back_substitute(a: &[Vec<MpComplex>], b: &[MpComplex], precision: u32) -> Result<Vec<MpComplex>, String> {
    let mut solution = vec![zero(); b.len()];
    for i in (0..b.len()).rev() {
        let mut value = b[i].clone();
        for j in i + 1..b.len() { value = value.sub(&a[i][j].mul(&solution[j], precision), precision); }
        solution[i] = value.div(&a[i][i], precision).ok_or("Jacobian is numerically singular")?;
    }
    Ok(solution)
}

fn square_step(mut a: Vec<Vec<MpComplex>>, mut b: Vec<MpComplex>, precision: u32) -> Result<Vec<MpComplex>, String> {
    for col in 0..b.len() {
        let pivot = (col..b.len()).max_by(|&i, &j| a[i][col].abs(precision).cmp(&a[j][col].abs(precision))).unwrap();
        if a[pivot][col].re.is_zero() && a[pivot][col].im.is_zero() { return Err("Jacobian is numerically singular".into()); }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for i in col + 1..b.len() {
            let factor = a[i][col].div(&a[col][col], precision).ok_or("Jacobian is numerically singular")?;
            a[i][col] = zero();
            for j in col + 1..b.len() { a[i][j] = a[i][j].sub(&factor.mul(&a[col][j], precision), precision); }
            b[i] = b[i].sub(&factor.mul(&b[col], precision), precision);
        }
    }
    back_substitute(&a, &b, precision)
}

/// Modified Gram–Schmidt with reorthogonalization avoids the squared condition
/// number of normal equations, and uses Hermitian products for complex systems.
fn least_squares_step(a: &[Vec<MpComplex>], b: &[MpComplex], precision: u32) -> Result<Vec<MpComplex>, String> {
    let n = a[0].len();
    let mut q: Vec<Vec<MpComplex>> = Vec::new();
    let mut r = vec![vec![zero(); n]; n];
    for j in 0..n {
        let mut column: Vec<_> = a.iter().map(|row| row[j].clone()).collect();
        let initial_norm = inner_product(&column, &column, precision).re.sqrt(precision).ok_or("Jacobian is numerically singular")?;
        for _ in 0..2 {
            for i in 0..j {
                let projection = inner_product(&q[i], &column, precision);
                r[i][j] = r[i][j].add(&projection, precision);
                for k in 0..column.len() { column[k] = column[k].sub(&q[i][k].mul(&projection, precision), precision); }
            }
        }
        let norm = inner_product(&column, &column, precision).re.sqrt(precision).ok_or("Jacobian is numerically singular")?;
        if norm.is_zero() || norm.cmp(&initial_norm.ldexp(-(precision as i64) + 16)) != Ordering::Greater {
            return Err("Jacobian is numerically singular".into());
        }
        let inverse = BigFloat::one().div(&norm, precision).ok_or("Jacobian is numerically singular")?;
        r[j][j] = MpComplex::real(norm);
        q.push(column.iter().map(|value| value.scale(&inverse, precision)).collect());
    }
    let projected: Vec<_> = q.iter().map(|column| inner_product(column, b, precision)).collect();
    back_substitute(&r, &projected, precision)
}

fn verify(expressions: &[Expr], variables: &[Expr], values: &[MpComplex], options: &Options, precision: u32) -> Result<(), String> {
    if options.verify {
        let residual = max_norm(&evaluate(expressions, variables, values, precision)?, precision);
        if residual.mul(&residual, precision).cmp(&options.tolerance) == Ordering::Greater {
            return Err("Could not find root within given tolerance".into());
        }
    }
    Ok(())
}

pub fn newton(equations: &[Expr], variables: &[Expr], starts: &[MpComplex], options: &Options) -> Result<Vec<MpComplex>, String> {
    let minimum = check(options)?;
    // Retain the information already present in high precision coefficients
    // and starting points, even when the requested output precision is lower.
    let coefficient_precision = equations.iter().filter_map(evalf_mp::max_precision).max().unwrap_or(0);
    let start_precision = starts.iter().map(|value| value.re.man.bits().max(value.im.man.bits()) as u32).max().unwrap_or(0);
    let precision = minimum.max(coefficient_precision.max(start_precision) + 64);
    let n = variables.len();
    if n == 0 || starts.len() != n || equations.len() < n {
        return Err("nsolve requires one start per unknown and at least as many equations as unknowns".into());
    }
    let expressions = equations.iter().map(equation).collect::<Result<Vec<_>, _>>()?;
    let jacobian: Vec<Expr> = expressions.iter().flat_map(|expression| variables.iter().map(move |variable| crate::diff::diff(expression, variable))).collect();
    let mut values: Vec<_> = starts.iter().map(|value| value.round_to(precision)).collect();
    for _ in 0..options.maxsteps {
        let residuals = evaluate(&expressions, variables, &values, precision)?;
        let residual = max_norm(&residuals, precision);
        if residual.cmp(&options.tolerance) != Ordering::Greater { break; }
        let evaluated = evaluate(&jacobian, variables, &values, precision)?;
        let rows: Vec<Vec<_>> = evaluated.chunks(n).map(|row| row.to_vec()).collect();
        let rhs: Vec<_> = residuals.iter().map(MpComplex::neg).collect();
        let step = if rows.len() == n { square_step(rows, rhs, precision)? }
            else { least_squares_step(&rows, &rhs, precision)? };
        let step_norm = max_norm(&step, precision);
        let scale = max_norm(&values, precision).add(&BigFloat::one(), precision);
        if step_norm.cmp(&options.tolerance.mul(&scale, precision)) != Ordering::Greater {
            for i in 0..n { values[i] = values[i].add(&step[i], precision); }
            break;
        }
        let mut damping = BigFloat::one();
        let mut accepted = None;
        for _ in 0..24 {
            let candidate: Vec<_> = values.iter().zip(&step).map(|(value, delta)| value.add(&delta.scale(&damping, precision), precision)).collect();
            if let Ok(next) = evaluate(&expressions, variables, &candidate, precision) {
                if max_norm(&next, precision).cmp(&residual) != Ordering::Greater {
                    accepted = Some(candidate);
                    break;
                }
            }
            damping = damping.ldexp(-1);
        }
        match accepted {
            Some(candidate) => values = candidate,
            None if !options.verify => break,
            None => return Err("Could not find root within given tolerance (Newton step did not decrease the residual)".into()),
        }
    }
    verify(&expressions, variables, &values, options, precision)?;
    Ok(values.iter().map(|value| value.round_to(options.precision)).collect())
}

/// Derivative-free scalar secant iteration, including complex starting points.
/// A single point uses the same exact quarter-unit perturbation as the
/// reference solver. No coefficient or iterate is converted through f64.
pub fn secant(equation_expr: &Expr, variable: &Expr, starts: &[MpComplex], options: &Options) -> Result<MpComplex, String> {
    let minimum = check(options)?;
    if starts.is_empty() || starts.len() > 2 {
        return Err("secant requires one or two initial points".into());
    }
    let start_precision = starts.iter().map(|v| v.re.man.bits().max(v.im.man.bits()) as u32).max().unwrap_or(0);
    let precision = minimum.max(evalf_mp::max_precision(equation_expr).unwrap_or(0) + 64).max(start_precision + 64);
    let expressions = vec![equation(equation_expr)?];
    let variables = vec![variable.clone()];
    let at = |point: &MpComplex| -> Result<MpComplex, String> {
        Ok(evaluate(&expressions, &variables, std::slice::from_ref(point), precision)?[0].clone())
    };
    let mut previous = starts[0].round_to(precision);
    let mut current = starts.get(1).map(|v| v.round_to(precision)).unwrap_or_else(||
        previous.add(&MpComplex::real(BigFloat::one().ldexp(-2)), precision));
    let mut previous_value = at(&previous)?;
    if previous_value.re.is_zero() && previous_value.im.is_zero() {
        return Ok(previous.round_to(options.precision));
    }
    let mut current_value = at(&current)?;
    for _ in 0..options.maxsteps {
        if current_value.abs(precision).cmp(&options.tolerance) != Ordering::Greater { break; }
        let denominator = current_value.sub(&previous_value, precision);
        let correction = current_value.mul(&current.sub(&previous, precision), precision)
            .div(&denominator, precision).ok_or("Could not find root: equal residuals in secant iteration")?;
        let next = current.sub(&correction, precision);
        let threshold = options.tolerance.mul(&next.abs(precision).add(&BigFloat::one(), precision), precision);
        previous = current;
        previous_value = current_value;
        current = next;
        current_value = at(&current)?;
        if correction.abs(precision).cmp(&threshold) != Ordering::Greater { break; }
    }
    verify(&expressions, &variables, std::slice::from_ref(&current), options, precision)?;
    Ok(current.round_to(options.precision))
}

pub fn bisect(equation_expr: &Expr, variable: &Expr, a: &MpComplex, b: &MpComplex, options: &Options) -> Result<MpComplex, String> {
    let precision = check(options)?.max(evalf_mp::max_precision(equation_expr).unwrap_or(0) + 64)
        .max(a.re.man.bits().max(b.re.man.bits()) as u32 + 64);
    if !a.is_real() || !b.is_real() || a.re.cmp(&b.re) != Ordering::Less {
        return Err("bisection requires a real bracket with lower endpoint less than upper endpoint".into());
    }
    let expressions = vec![equation(equation_expr)?];
    let variables = vec![variable.clone()];
    let at = |point: &BigFloat| -> Result<BigFloat, String> {
        let values = evaluate(&expressions, &variables, &[MpComplex::real(point.clone())], precision)?;
        if !values[0].is_real() { return Err("bisection requires a real-valued equation".into()); }
        Ok(values[0].re.clone())
    };
    let (mut a, mut b) = (a.re.clone(), b.re.clone());
    let (mut fa, fb) = (at(&a)?, at(&b)?);
    if fa.is_zero() { return Ok(MpComplex::real(a.round_to(options.precision))); }
    if fb.is_zero() { return Ok(MpComplex::real(b.round_to(options.precision))); }
    if fa.sign() == fb.sign() { return Err("bisection requires opposite signs at the bracket endpoints".into()); }
    let mut midpoint = a.clone();
    for _ in 0..options.maxsteps {
        midpoint = a.add(&b, precision).ldexp(-1);
        let fm = at(&midpoint)?;
        let threshold = options.tolerance.mul(&midpoint.abs().add(&BigFloat::one(), precision), precision);
        if fm.is_zero() || b.sub(&a, precision).abs().cmp(&threshold) != Ordering::Greater { break; }
        if fm.sign() == fa.sign() { a = midpoint.clone(); fa = fm; } else { b = midpoint.clone(); }
    }
    let result = MpComplex::real(midpoint);
    verify(&expressions, &variables, std::slice::from_ref(&result), options, precision)?;
    Ok(result.round_to(options.precision))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn expression(source: &str) -> Expr { crate::parser::parse_str(source).unwrap() }

    fn options() -> Options {
        let precision = crate::bigfloat::dps_to_prec(80);
        Options { precision, tolerance: BigFloat::one().ldexp(-(precision as i64) - 8), maxsteps: 100, verify: true }
    }

    #[test]
    fn secant_uses_brackets_and_complex_values_without_derivatives() {
        let x = Expr::symbol("x");
        let options = options();
        let starts = [MpComplex::real(BigFloat::one()), MpComplex::real(BigFloat::from_i64(2))];
        let root = secant(&expression("x**2 - 2"), &x, &starts, &options).unwrap();
        assert!(root.re.to_string_dps(60, false).starts_with("1.41421356237309504880168872420969807856967187537694807317668"));
        let start = MpComplex { re: BigFloat::one(), im: BigFloat::one() };
        let root = secant(&expression("x**2 + 1"), &x, &[start], &options).unwrap();
        assert!(root.re.abs().mag() < -(options.precision as i64));
        assert!(root.im.sub(&BigFloat::one(), options.precision).abs().mag() < -(options.precision as i64));
        assert!(secant(&expression("x**2 + 1"), &x, &starts, &options).is_err());
    }

    #[test]
    fn real_and_complex_newton_retain_arbitrary_precision() {
        let x = Expr::symbol("x");
        let options = options();
        let real = newton(&[expression("x**2 - 2")], &[x.clone()], &[MpComplex::real(BigFloat::one())], &options).unwrap();
        assert!(real[0].re.to_string_dps(70, false).starts_with("1.414213562373095048801688724209698078569671875376948073176679737990732"));
        let initial = MpComplex { re: BigFloat::one(), im: BigFloat::one() };
        let root = newton(&[expression("x**2 + 1")], &[x], &[initial], &options).unwrap();
        assert!(root[0].re.abs().mag() < -(options.precision as i64));
        assert!(root[0].im.sub(&BigFloat::one(), options.precision).abs().mag() < -(options.precision as i64));
        let root = newton(&[expression("exp(x)+1")], &[Expr::symbol("x")], &[MpComplex { re: BigFloat::zero(), im: BigFloat::from_i64(3) }], &options).unwrap();
        assert!(root[0].im.sub(&crate::bigfloat::pi(options.precision), options.precision).abs().mag() < -(options.precision as i64) + 3);
    }

    #[test]
    fn overdetermined_complex_qr_and_bisection() {
        let options = options();
        let vars = vec![Expr::symbol("x"), Expr::symbol("y")];
        let equations = vec![expression("x + I*y - 2"), expression("x - I*y - 4"), expression("2*x - 6")];
        let root = newton(&equations, &vars, &[zero(), zero()], &options).unwrap();
        assert_eq!(root[0].re, BigFloat::from_i64(3));
        assert_eq!(root[1].im, BigFloat::one());
        let options = Options { maxsteps: 400, ..options };
        let root = bisect(&expression("x**2 - 2"), &vars[0], &zero(), &MpComplex::real(BigFloat::from_i64(2)), &options).unwrap();
        assert!(root.re.to_string_dps(60, false).starts_with("1.41421356237309504880168872420969807856967187537694807317668"));
        assert!(bisect(&expression("1/x"), &vars[0], &MpComplex::real(BigFloat::from_i64(-1)), &MpComplex::real(BigFloat::one()), &options).is_err());
    }

    #[test]
    fn nonlinear_systems_and_convergence_controls() {
        let options = options();
        let x = Expr::symbol("x");
        let root = newton(&[expression("cos(x)-x")], &[x.clone()], &[MpComplex::real(BigFloat::one())], &options).unwrap();
        assert!(root[0].re.to_string_dps(60, false).starts_with("0.739085133215160641655312087673873404013411758900757464965681"));
        let vars = vec![x.clone(), Expr::symbol("y")];
        let equations = vec![expression("x*x-2"), expression("y*y-3"), expression("x*x+y*y-5")];
        let roots = newton(&equations, &vars, &[MpComplex::real(BigFloat::one()), MpComplex::real(BigFloat::one())], &options).unwrap();
        assert!(roots[1].re.to_string_dps(50, false).starts_with("1.7320508075688772935274463415058723669428052538104"));
        let mut one_step = Options { maxsteps: 1, ..options };
        let starts = [MpComplex::real(BigFloat::from_i64(100))];
        assert!(newton(&[expression("x*x-2")], &[x.clone()], &starts, &one_step).is_err());
        one_step.verify = false;
        assert!(newton(&[expression("x*x-2")], &[x.clone()], &starts, &one_step).unwrap()[0].re.cmp(&BigFloat::from_i64(40)) == Ordering::Greater);
        one_step.verify = true;
        assert!(newton(&[expression("x-1"), expression("x-2")], &[x], &[zero()], &one_step).is_err());
        assert!(newton(&[expression("x+y-1"), expression("2*x+2*y-2"), expression("3*x+3*y-3")], &vars, &[zero(), zero()], &one_step).is_err());
        let large = expression("Float('1e400', 50)");
        let point = evalf_mp::eval_mp(&large, crate::bigfloat::dps_to_prec(50)).unwrap();
        let options = Options { precision: 53, tolerance: BigFloat::one().ldexp(-61), maxsteps: 20, verify: true };
        let root = newton(&[crate::ops::sub(&vars[0], &large)], &vars[..1], &[point], &options).unwrap();
        assert_eq!(root[0].re, evalf_mp::eval_mp(&large, 53).unwrap().re);
    }
}
