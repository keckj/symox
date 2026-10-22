//! Bounded complex expansion, separate from the ordinary re()/im() evaluators.
//!
//! All identity decisions are structural/exact. Opaque values retain projection
//! nodes; integer/rational powers use the principal complex branch. Physical
//! owners keep shared DAGs and synthesized arguments safe during memoization.
use crate::expand::{expand_with_options, ExpandBudget, ExpandOptions};
use crate::expr::Constant;
use crate::functions as f;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, powi, sub};
use crate::{Expr, Func, Kind, Number};
use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive};
use rustc_hash::FxHashMap;

pub(crate) mod evaluate;

type Parts = (Expr, Expr);
type Result<T> = std::result::Result<T, PartsError>;
#[derive(Debug, Clone)]
pub struct PartsError(pub String);
impl std::fmt::Display for PartsError {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.write_str(&self.0)
    }
}
impl std::error::Error for PartsError {}
// Opt-in guards of the bounded entry points (`as_real_imag_bounded`), used
// only by internal heuristics that have no SymPy counterpart and must stay
// cheap. SymPy's `as_real_imag`/`expand(complex=True)` have no size limit, so
// the public entry points only keep the recursion-depth guard.
const MAX_NODES: usize = 20_000;
const MAX_TERMS: usize = 4096;
const MAX_POWER: u64 = 128;
const MAX_BITS: u64 = 1_000_000;
const MAX_DEPTH: usize = 256;
// Binomial expansion of `(a + I*b)**n` with a non-real base holds n + 1 terms
// whose coefficients reach n bits: past 2**17 the output alone needs more
// than 2 GiB, so it cannot be materialized (SymPy runs out of memory or never
// finishes). Raising keeps the call interruptible instead of hanging in
// native code; nothing partial is ever returned.
const FEASIBLE_POWER: u64 = 1 << 17;
const FEASIBLE_BITS: u64 = 1 << 34;
fn infeasible() -> PartsError {
    PartsError(format!("as_real_imag: result beyond the materialization limit (binomial powers through {FEASIBLE_POWER}, numbers of {FEASIBLE_BITS} bits)"))
}
fn depth_limit() -> PartsError {
    PartsError(format!("as_real_imag: expression nesting deeper than {MAX_DEPTH} levels"))
}
fn limit() -> PartsError {
    PartsError("as_real_imag resource limit exceeded (20000 nodes, 4096 terms, integer powers through 128, 1000000 numeric bits)".into())
}
fn raw(function: Func, value: &Expr) -> Expr {
    Expr::from_kind(Kind::Func(function, vec![value.clone()]))
}
fn is_zero(e: &Expr) -> bool {
    crate::assumptions::is_zero(e) == Some(true)
}
fn is_real(e: &Expr) -> bool {
    crate::assumptions::fact(e, "extended_real") == Some(true)
}
fn bits(number: &Number) -> u64 {
    match number {
        Number::Int(n) => {
            if n.unsigned_abs() <= 1 {
                0
            } else {
                u64::from(64 - n.unsigned_abs().leading_zeros())
            }
        }
        Number::BigInt(n) => n.bits(),
        Number::Rational(n) => n.numer().bits() + n.denom().bits(),
        Number::Float(_) => 64,
        Number::MpFloat(value, p) => {
            u64::from(*p).max(value.man.bits().saturating_add(u64::try_from(value.exp.unsigned_abs()).unwrap_or(u64::MAX)))
        }
    }
}

/// Reserve before any active constructor or assumption query can force a
/// compact held numeric power. The traversal itself never evaluates expressions.
fn inspect(
    e: &Expr,
    memo: &mut FxHashMap<usize, (Expr, (usize, u64))>,
    depth: usize,
) -> Result<(usize, u64)> {
    if let Some((_, value)) = memo.get(&e.dag_id()) {
        return Ok(*value);
    }
    if depth > MAX_DEPTH {
        return Err(depth_limit());
    }
    if memo.len() > MAX_NODES {
        return Err(limit());
    }
    let mut child = |e: &Expr| inspect(e, memo, depth + 1);
    let mut power_bound = |base: &Expr, exponent: &Expr| -> Result<(usize, u64)> {
        let (terms, bits) = child(base)?;
        child(exponent)?;
        if exponent.is_integer_number() {
            let n = exponent
                .as_number()
                .and_then(Number::to_i64)
                .ok_or_else(limit)?
                .unsigned_abs();
            if n > MAX_POWER {
                return Err(limit());
            }
            let mut count = 1u128;
            for k in 1..=(terms.saturating_sub(1) as u64).min(n) {
                count = count.saturating_mul(u128::from(n + terms as u64 - k)) / u128::from(k);
                if count > MAX_TERMS as u128 {
                    return Err(limit());
                }
            }
            Ok((count as usize, bits.saturating_mul(n)))
        } else {
            if let Some((p, _)) = exponent.as_number().and_then(Number::numer_denom) {
                let n = p.abs().to_u64().ok_or_else(limit)?;
                if n > MAX_POWER {
                    return Err(limit());
                }
                return Ok((1, bits.saturating_mul(n.max(1)).saturating_mul(2)));
            }
            Ok((1, bits))
        }
    };
    let value = match e.kind() {
        Kind::Num(number) => {
            if number
                .float_precision()
                .is_some_and(|precision| precision > 4096)
            {
                return Err(limit());
            }
            (1, bits(number))
        }
        Kind::Add(c, terms) => {
            let mut width = usize::from(!c.is_zero());
            let mut height = bits(c);
            for (term, coefficient) in terms {
                let (n, b) = inspect(term, memo, depth + 1)?;
                width = width.saturating_add(n);
                height = height.max(b.saturating_add(bits(coefficient)));
            }
            (
                width.max(1),
                height.saturating_add(usize::BITS as u64 - u64::from(width.max(1).leading_zeros())),
            )
        }
        Kind::Mul(c, factors) => {
            let mut width = 1usize;
            let mut height = bits(c);
            for (base, exponent) in factors {
                let (n, b) = power_bound(base, exponent)?;
                width = width.saturating_mul(n);
                height = height.saturating_add(b);
            }
            (width, height)
        }
        Kind::Pow(base, exponent) => power_bound(base, exponent)?,
        Kind::Func(_, args) => {
            for a in args {
                inspect(a, memo, depth + 1)?;
            }
            (1, 0)
        }
        Kind::Derivative(body, vars) => {
            inspect(body, memo, depth + 1)?;
            for (v, _) in vars {
                inspect(v, memo, depth + 1)?;
            }
            (1, 0)
        }
        Kind::Rel(_, a, b) => {
            inspect(a, memo, depth + 1)?;
            inspect(b, memo, depth + 1)?;
            (1, 0)
        }
        Kind::Piecewise(pairs) => {
            for (a, b) in pairs {
                inspect(a, memo, depth + 1)?;
                inspect(b, memo, depth + 1)?;
            }
            (1, 0)
        }
        _ => (1, 0),
    };
    if value.0 > MAX_TERMS
        || value.1 > MAX_BITS
        || (value.0 as u64).saturating_mul(value.1) > MAX_BITS
    {
        return Err(limit());
    }
    memo.insert(e.dag_id(), (e.clone(), value));
    Ok(value)
}

/// The guards of the unbounded entry points: the recursive split would
/// otherwise overflow the native stack where SymPy raises RecursionError,
/// and a shared numeric DAG (`Mul(v, v)` nested 80 times) would multiply out
/// a number of more than 2 GiB that can never be materialized.
fn height(e: &Expr, memo: &mut FxHashMap<usize, (Expr, (usize, Option<u64>))>, depth: usize) -> Result<(usize, Option<u64>)> {
    if depth > MAX_DEPTH {
        return Err(depth_limit());
    }
    if let Some((_, value)) = memo.get(&e.dag_id()) {
        return Ok(*value);
    }
    let mut value = 0;
    let children = e.args();
    let mut numeric: Vec<Option<u64>> = Vec::with_capacity(children.len());
    for a in &children {
        let (h, n) = height(a, memo, depth + 1)?;
        value = value.max(h + 1);
        numeric.push(n);
    }
    if depth + value > MAX_DEPTH {
        return Err(depth_limit());
    }
    // a size estimate (in bits) of the number an unevaluated numeric
    // product or sum multiplies out to
    let estimate = match e.kind() {
        Kind::Num(number) => Some(bits(number)),
        Kind::Mul(..) => numeric.iter().try_fold(0u64, |acc, n| n.map(|n| acc.saturating_add(n))),
        Kind::Add(..) => numeric
            .iter()
            .try_fold(0u64, |acc, n| n.map(|n| acc.max(n)))
            .map(|n| n.saturating_add(64)),
        // a power of a real number stays whole, so it is never multiplied out
        _ => None,
    };
    if matches!(e.kind(), Kind::Mul(..) | Kind::Add(..)) && estimate.is_some_and(|n| n > FEASIBLE_BITS) {
        return Err(infeasible());
    }
    memo.insert(e.dag_id(), (e.clone(), (value, estimate)));
    Ok((value, estimate))
}
fn precheck(e: &Expr, bounded: bool) -> Result<()> {
    if bounded {
        inspect(e, &mut FxHashMap::default(), 0).map(|_| ())
    } else {
        height(e, &mut FxHashMap::default(), 0).map(|_| ())
    }
}

struct Context {
    expanded: ExpandBudget,
    parts: FxHashMap<(usize, bool), (Expr, Parts)>,
    projections: FxHashMap<usize, (Expr, Parts)>,
    /// The other hints of an enclosing `expand(complex=True, ...)`. SymPy
    /// forwards them only to the outermost `as_real_imag` call of each node;
    /// nested part computations run without hints.
    hints: ExpandOptions,
    hinted: FxHashMap<usize, (Expr, Parts)>,
    complex: FxHashMap<(usize, bool), (Expr, Expr)>,
    work: usize,
    /// Whether the opt-in resource guards apply (internal heuristics only).
    bounded: bool,
}
impl Context {
    fn new(bounded: bool) -> Self {
        Self::with_hints(ExpandOptions::default(), bounded)
    }
    fn with_hints(hints: ExpandOptions, bounded: bool) -> Self {
        Self {
            expanded: if bounded { ExpandBudget::new(MAX_TERMS) } else { ExpandBudget::unlimited() },
            bounded,
            parts: FxHashMap::default(),
            projections: FxHashMap::default(),
            hints: ExpandOptions { deep: true, ..hints },
            hinted: FxHashMap::default(),
            complex: FxHashMap::default(),
            work: 0,
        }
    }
    fn charge(&mut self, n: usize) -> Result<()> {
        self.work = self.work.saturating_add(n);
        if self.bounded && self.work > 200_000 {
            Err(limit())
        } else {
            Ok(())
        }
    }
    fn expand(&mut self, e: &Expr, deep: bool, mul_only: bool) -> Result<Expr> {
        let options = if mul_only {
            ExpandOptions {
                deep,
                multinomial: false,
                power_base: false,
                power_exp: false,
                log: false,
                ..ExpandOptions::default()
            }
        } else {
            ExpandOptions {
                deep,
                ..ExpandOptions::default()
            }
        };
        expand_with_options(e, &options, &mut self.expanded)
            .map_err(|error| PartsError(format!("as_real_imag: {error}")))
    }
    /// SymPy's `value.expand(deep, **hints)` inside `as_real_imag`: the
    /// default expansion, or a nested complex expansion when hinted.
    fn expand_nested(&mut self, e: &Expr, hinted: bool) -> Result<Expr> {
        if hinted {
            self.expand_full(e, true)
        } else {
            self.expand(e, true, false)
        }
    }
    /// SymPy's `hints['complex'] = False; value.expand(deep, **hints)`.
    fn expand_other(&mut self, e: &Expr, hinted: bool, deep: bool) -> Result<Expr> {
        if !hinted {
            return self.expand(e, deep, false);
        }
        let o = ExpandOptions { deep, ..self.hints };
        if !(o.mul || o.multinomial || o.power_base || o.power_exp || o.log || o.trig) {
            return Ok(e.clone());
        }
        expand_with_options(e, &o, &mut self.expanded)
            .map_err(|error| PartsError(format!("expand: {error}")))
    }
    /// `expand(complex=True, **hints)`: SymPy sorts the complex hint before
    /// every other enabled hint, so the complex pass runs first.
    fn expand_full(&mut self, e: &Expr, deep: bool) -> Result<Expr> {
        let value = self.complex_pass(e, deep)?;
        self.expand_other(&value, true, deep)
    }
    /// One `_eval_expand_complex` pass: children first (when deep), then
    /// `re + I*im` of the rebuilt node, kept only when it differs.
    fn complex_pass(&mut self, e: &Expr, deep: bool) -> Result<Expr> {
        if let Some((_, value)) = self.complex.get(&(e.dag_id(), deep)) {
            return Ok(value.clone());
        }
        self.charge(1)?;
        let unsupported = || {
            Err(PartsError(format!(
                "expand(complex=True) of {e} is not implemented"
            )))
        };
        let mut current = e.clone();
        let mut tuple = false;
        if deep {
            let rebuild: fn(&Expr, Vec<Expr>) -> Expr = match e.kind() {
                Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) => |e, _| e.clone(),
                Kind::Add(..) => |_, args| add(&args),
                Kind::Mul(..) => |_, args| mul(&args),
                Kind::Pow(..) => |_, args| pow(&args[0], &args[1]),
                Kind::Func(Func::Undefined(head), _) if head.is_internal("Tuple") => {
                    tuple = true;
                    |e, args| match e.kind() {
                        Kind::Func(head, _) => f::func(head.clone(), args),
                        _ => unreachable!(),
                    }
                }
                Kind::Func(Func::Undefined(head), _)
                    if head.is_internal("UnevaluatedExpr") || !head.internal =>
                {
                    |e, args| match e.kind() {
                        Kind::Func(head, _) => f::func(head.clone(), args),
                        _ => unreachable!(),
                    }
                }
                Kind::Func(
                    Func::Undefined(_)
                    | Func::Order
                    | Func::Subs
                    | Func::Limit
                    | Func::Integral
                    | Func::Sum
                    | Func::Logic(_)
                    | Func::Set(_)
                    | Func::MatrixElement(_)
                    | Func::MatrixProxy(_)
                    | Func::MatrixTrace(_)
                    | Func::MatrixDeterminant(_)
                    | Func::MatrixDotProduct(..),
                    _,
                )
                | Kind::Derivative(..)
                | Kind::Rel(..)
                | Kind::Piecewise(_)
                | Kind::Bool(_) => return unsupported(),
                Kind::Func(Func::Abs, _) => |_, args| sympy_abs(&args[0]),
                Kind::Func(..) => |e, args| match e.kind() {
                    Kind::Func(head, _) => f::func(head.clone(), args),
                    _ => unreachable!(),
                },
            };
            let args = crate::traverse::ordered_args(e);
            let mut mapped = Vec::with_capacity(args.len());
            for arg in &args {
                mapped.push(self.complex_pass(arg, true)?);
            }
            if !mapped.iter().zip(&args).all(|(a, b)| a.same_presentation(b)) {
                current = rebuild(e, mapped);
            }
        } else if matches!(e.kind(), Kind::Func(Func::Undefined(head), _) if head.is_internal("Tuple"))
        {
            tuple = true;
        }
        let value = if tuple || current.has_matrix() {
            if current.has_matrix() {
                return unsupported();
            }
            current
        } else {
            scalar(&current)?;
            let (r, i) = self.split_hinted(&current)?;
            let candidate = add2(&r, &mul2(&Expr::i(), &i));
            if candidate == current {
                current
            } else {
                candidate
            }
        };
        self.complex
            .insert((e.dag_id(), deep), (e.clone(), value.clone()));
        Ok(value)
    }
    /// The outermost `as_real_imag(**hints)` call of a complex expansion.
    fn split_hinted(&mut self, e: &Expr) -> Result<Parts> {
        if let Some((_, p)) = self.hinted.get(&e.dag_id()) {
            return Ok(p.clone());
        }
        let p = match e.kind() {
            Kind::Pow(base, exponent) => self.power(e, base, exponent, true, true)?,
            Kind::Func(function, args) if args.len() == 1 && elementary(function) => {
                self.elementary(e, function, &args[0], true, true)?
            }
            _ => self.split(e, true)?,
        };
        self.hinted.insert(e.dag_id(), (e.clone(), p.clone()));
        Ok(p)
    }
    /// Local projections match complex expansion without changing the global
    /// ordinary re()/im() evaluators. Opaque function products stay held.
    fn project(&mut self, e: &Expr) -> Result<Parts> {
        if let Some((_, p)) = self.projections.get(&e.dag_id()) {
            return Ok(p.clone());
        }
        self.charge(1)?;
        let p = match e.kind() {
            Kind::Num(_)
            | Kind::Const(
                Constant::Pi
                | Constant::E
                | Constant::EulerGamma
                | Constant::Infinity
                | Constant::NegInfinity,
            ) => (e.clone(), Expr::zero()),
            Kind::Const(Constant::NaN | Constant::ComplexInfinity) => (Expr::nan(), Expr::nan()),
            Kind::Const(Constant::I) => (Expr::zero(), Expr::one()),
            Kind::Func(Func::Re | Func::Im | Func::Arg, args)
                if args[0].is_constant(Constant::ComplexInfinity)
                    || (matches!(e.kind(), Kind::Func(Func::Arg, _)) && args[0].is_zero()) =>
            {
                (Expr::nan(), Expr::nan())
            }
            Kind::Func(Func::Re | Func::Im | Func::Abs | Func::Arg, _) => (e.clone(), Expr::zero()),
            Kind::Func(Func::Conjugate, args) => {
                let (r, i) = self.project(&args[0])?;
                (r, neg(&i))
            }
            Kind::Add(c, terms) => {
                let mut real = vec![Expr::num(c.clone())];
                let mut imag = Vec::new();
                for (term, coefficient) in terms {
                    let (r, i) = self.project(term)?;
                    let c = Expr::num(coefficient.clone());
                    real.push(mul2(&c, &r));
                    imag.push(mul2(&c, &i));
                }
                (add(&real), add(&imag))
            }
            _ if e.has_matrix() => (raw(Func::Re, e), raw(Func::Im, e)),
            _ if is_real(e) => (e.clone(), Expr::zero()),
            _ if crate::assumptions::fact(e, "imaginary") == Some(true) => (
                Expr::zero(),
                mul(&[Expr::minus_one(), Expr::i(), e.clone()]),
            ),
            Kind::Mul(..) => self.product(e, true, None)?.unwrap(),
            Kind::Pow(base, exponent)
                if exponent.is_rational_number() || base.is_constant(Constant::E) =>
            {
                self.power(e, base, exponent, true, false)?
            }
            Kind::Func(function, args)
                if args.len() == 1
                    && elementary(function) =>
            {
                self.elementary(e, function, &args[0], true, false)?
            }
            Kind::Func(Func::Sec | Func::Csc, args) if args.len() == 1 => self.reciprocal(e, &args[0])?,
            _ => (raw(Func::Re, e), raw(Func::Im, e)),
        };
        self.projections.insert(e.dag_id(), (e.clone(), p.clone()));
        Ok(p)
    }
    fn split(&mut self, e: &Expr, deep: bool) -> Result<Parts> {
        if let Some((_, p)) = self.parts.get(&(e.dag_id(), deep)) {
            return Ok(p.clone());
        }
        self.charge(1)?;
        let p = match e.kind() {
            Kind::Add(c, terms) => {
                let mut real = vec![Expr::num(c.clone())];
                let mut imag = Vec::new();
                for (term, coefficient) in terms {
                    let (r, i) = self.split(term, deep)?;
                    let c = Expr::num(coefficient.clone());
                    real.push(mul2(&c, &r));
                    imag.push(mul2(&c, &i));
                }
                (add(&real), add(&imag))
            }
            Kind::Mul(..) => self.product(e, true, None)?.unwrap(),
            Kind::Pow(base, exponent) => self.power(e, base, exponent, deep, false)?,
            Kind::Func(function, args)
                if args.len() == 1
                    && elementary(function) =>
            {
                self.elementary(e, function, &args[0], deep, false)?
            }
            Kind::Func(Func::Sec | Func::Csc, args) if args.len() == 1 => self.reciprocal(e, &args[0])?,
            _ => self.project(e)?,
        };
        self.parts
            .insert((e.dag_id(), deep), (e.clone(), p.clone()));
        Ok(p)
    }
    fn product(&mut self, e: &Expr, complex: bool, ignore: Option<&Expr>) -> Result<Option<Parts>> {
        let mut real = Vec::new();
        let mut imaginary = Vec::new();
        let mut other: Vec<Expr> = Vec::new();
        let mut sums = Vec::new();
        for factor in crate::traverse::ordered_args(e) {
            let (r, i) = if complex {
                self.split(&factor, true)?
            } else {
                self.project(&factor)?
            };
            if is_zero(&i) {
                real.push(r);
            } else if is_zero(&r) {
                imaginary.push(i);
            } else {
                self.charge(other.len())?;
                let index = if other.is_empty() {
                    None
                } else {
                    let conjugate = f::conjugate(&factor);
                    other.iter().position(|candidate| *candidate == conjugate)
                };
                if let Some(index) = index {
                    let candidate = other.remove(index);
                    real.push(powi(&f::abs(&candidate), 2));
                } else if complex && factor.is_add() {
                    sums.push(factor);
                } else {
                    other.push(factor);
                }
            }
        }
        let m = mul(&other);
        if ignore == Some(&m) {
            return Ok(None);
        }
        let imco = if imaginary.len() % 2 == 1 {
            imaginary.remove(0)
        } else {
            Expr::zero()
        };
        for pair in imaginary.chunks_exact(2) {
            real.push(neg(&mul2(&pair[0], &pair[1])));
        }
        let coefficient = mul(&real);
        if coefficient.is_constant(Constant::NaN) {
            return Ok(Some((Expr::nan(), Expr::nan())));
        }
        let mut p = if m.is_one() {
            if imco.is_zero() {
                (coefficient, Expr::zero())
            } else {
                (Expr::zero(), mul2(&coefficient, &imco))
            }
        } else {
            let (r, i) = if m.same_presentation(e) {
                (raw(Func::Re, &m), raw(Func::Im, &m))
            } else {
                self.project(&m)?
            };
            let (r, i) = (mul2(&coefficient, &r), mul2(&coefficient, &i));
            if imco.is_zero() {
                (r, i)
            } else {
                (neg(&mul2(&imco, &i)), mul2(&imco, &r))
            }
        };
        if !sums.is_empty() {
            // SymPy distributes only the top-level sums (expand_mul, deep=False).
            let width = sums
                .iter()
                .fold(1usize, |width, sum| width.saturating_mul(sum.args().len()));
            if self.bounded && width > MAX_TERMS {
                return Err(limit());
            }
            self.charge(width)?;
            let value = crate::simp::expand_mul_deep(&mul(&sums), false);
            let (r, i) = self.split(&value, true)?;
            p = (
                sub(&mul2(&p.0, &r), &mul2(&p.1, &i)),
                add2(&mul2(&p.1, &r), &mul2(&p.0, &i)),
            );
        }
        Ok(Some(p))
    }
    /// sec and csc split as SymPy's reciprocal `1/cos` and `1/sin` powers.
    fn reciprocal(&mut self, e: &Expr, argument: &Expr) -> Result<Parts> {
        let base = if matches!(e.kind(), Kind::Func(Func::Sec, _)) {
            f::cos(argument)
        } else {
            f::sin(argument)
        };
        let value = pow(&base, &Expr::minus_one());
        if value.same_presentation(e) {
            return Ok((raw(Func::Re, e), raw(Func::Im, e)));
        }
        self.split(&value, true)
    }
    fn power(
        &mut self,
        e: &Expr,
        base: &Expr,
        exponent: &Expr,
        deep: bool,
        hinted: bool,
    ) -> Result<Parts> {
        if exponent.is_integer_number() {
            let n = exponent.as_i64();
            if self.bounded && n.map_or(true, |n| n.unsigned_abs() > MAX_POWER) {
                return Err(limit());
            }
            let (mut a, mut b) = self.split(base, deep)?;
            if is_zero(&b) {
                return Ok((e.clone(), Expr::zero()));
            }
            let n = match n {
                Some(n) if n.unsigned_abs() <= FEASIBLE_POWER => n,
                _ => return Err(infeasible()),
            };
            if n < 0 {
                let magnitude = add2(&powi(&a, 2), &powi(&b, 2));
                a = div(&a, &magnitude);
                b = div(&neg(&b), &magnitude);
            }
            let mut real = Vec::new();
            let mut imag = Vec::new();
            let mut coefficient = BigInt::one();
            let degree = n.unsigned_abs();
            for k in 0..=degree {
                self.charge(1)?;
                let imaginary_base = if k % 4 == 3 { neg(&b) } else { b.clone() };
                let term = mul(&[
                    Expr::num(Number::from_bigint(coefficient.clone())),
                    powi(&a, (degree - k) as i64),
                    powi(&imaginary_base, k as i64),
                ]);
                match k % 4 {
                    0 => real.push(term),
                    1 | 3 => imag.push(term),
                    _ => real.push(neg(&term)),
                };
                if k < degree {
                    coefficient = coefficient * BigInt::from(degree - k) / BigInt::from(k + 1);
                }
            }
            return Ok((add(&real), add(&imag)));
        }
        if exponent.is_rational_number() {
            let (a, b) = self.split(base, deep)?;
            if exponent == &Expr::rational(1, 2) && is_zero(&b) {
                if crate::assumptions::fact(&a, "extended_nonnegative") == Some(true) {
                    return Ok((e.clone(), Expr::zero()));
                }
                if crate::assumptions::fact(&a, "extended_nonpositive") == Some(true) {
                    return Ok((Expr::zero(), pow(&neg(base), exponent)));
                }
            }
            if is_zero(&a) && is_zero(&b) {
                return Ok(if exponent.as_number().is_some_and(Number::is_positive) {
                    (Expr::zero(), Expr::zero())
                } else {
                    (Expr::nan(), Expr::nan())
                });
            }
            // a and b are real projections, so a**2+b**2 is nonnegative.
            // Combining these powers is valid on the principal branch even
            // when the general-purpose inference cannot prove that fact.
            let magnitude = add2(&powi(&a, 2), &powi(&b, 2));
            let theta = self.angle(&b, &a);
            let phase = mul2(&theta, exponent);
            let scale = pow(&magnitude, &div(exponent, &Expr::int(2)));
            return Ok((mul2(&scale, &f::cos(&phase)), mul2(&scale, &f::sin(&phase))));
        }
        if base.is_constant(Constant::E) {
            return self.elementary(e, &Func::Exp, exponent, deep, hinted);
        }
        let expanded = if deep {
            self.expand_other(e, hinted, true)?
        } else {
            e.clone()
        };
        // SymPy returns the evaluated `(re(expanded), im(expanded))` here.
        Ok((f::re(&expanded), f::im(&expanded)))
    }
    fn angle(&self, y: &Expr, x: &Expr) -> Expr {
        if is_zero(y) {
            if crate::assumptions::fact(x, "extended_positive") == Some(true) {
                return Expr::zero();
            }
            if crate::assumptions::fact(x, "extended_negative") == Some(true) {
                return Expr::pi();
            }
        }
        if crate::assumptions::is_positive(x) == Some(true) {
            return f::atan(&div(y, x));
        }
        f::atan2(y, x)
    }
    fn absolute(&mut self, value: &Expr) -> Result<Expr> {
        if let Kind::Func(Func::Conjugate, args) = value.kind() {
            return self.absolute(&args[0]);
        }
        if value.is_add() {
            let (r, i) = self.project(value)?;
            fn projections(e: &Expr, seen: &mut FxHashMap<usize, Expr>) -> bool {
                if seen.insert(e.dag_id(), e.clone()).is_some() {
                    return false;
                }
                matches!(e.kind(), Kind::Func(Func::Re | Func::Im, _))
                    || e.args().iter().any(|e| projections(e, seen))
            }
            if !projections(&r, &mut FxHashMap::default())
                && !projections(&i, &mut FxHashMap::default())
            {
                return Ok(pow(
                    &add2(&powi(&r, 2), &powi(&i, 2)),
                    &Expr::rational(1, 2),
                ));
            }
        }
        Ok(sympy_abs(value))
    }
    fn argument(&mut self, value: &Expr) -> Result<Expr> {
        if crate::assumptions::fact(value, "extended_positive") == Some(true) {
            return Ok(Expr::zero());
        }
        if crate::assumptions::fact(value, "extended_negative") == Some(true) {
            return Ok(Expr::pi());
        }
        let (r, i) = self.project(value)?;
        if r.is_number() && i.is_number() {
            return Ok(self.angle(&i, &r));
        }
        Ok(f::arg(value))
    }
    fn elementary(
        &mut self,
        e: &Expr,
        function: &Func,
        argument: &Expr,
        deep: bool,
        hinted: bool,
    ) -> Result<Parts> {
        if argument.is_constant(Constant::ComplexInfinity) || argument.is_constant(Constant::NaN) {
            return Ok((Expr::nan(), Expr::nan()));
        }
        if matches!(
            argument.kind(),
            Kind::Const(Constant::Infinity | Constant::NegInfinity)
        ) {
            let negative = argument.is_constant(Constant::NegInfinity);
            let value=match function{
                Func::Sin|Func::Cos|Func::Tan|Func::Cot=>return Err(PartsError("as_real_imag of oscillatory functions at infinity requires unsupported AccumBounds".into())),
                Func::Sinh=>if negative{Expr::neg_oo()}else{Expr::oo()},
                Func::Cosh|Func::Log=>Expr::oo(),
                Func::Tanh|Func::Coth=>Expr::int(if negative{-1}else{1}),
                Func::Exp=>if negative{Expr::zero()}else{Expr::oo()},
                _=>unreachable!(),
            };
            return Ok((
                value,
                if *function == Func::Log && negative {
                    Expr::pi()
                } else {
                    Expr::zero()
                },
            ));
        }
        if *function == Func::Log {
            let value = if deep {
                self.expand_nested(argument, hinted)?
            } else {
                argument.clone()
            };
            let absolute = self.absolute(&value)?;
            if absolute.is_constant(Constant::Infinity) {
                return Ok((Expr::oo(), Expr::zero()));
            }
            if absolute == value {
                return Ok((e.clone(), Expr::zero()));
            }
            let argument = self.argument(&value)?;
            let logarithm = f::log(&absolute);
            if hinted && self.hints.log {
                return Ok((self.expand_other(&logarithm, true, deep)?, argument));
            }
            return Ok((logarithm, argument));
        }
        // Trigonometric and hyperbolic functions test the argument before
        // splitting it; real arguments only receive the non-complex hints.
        if *function != Func::Exp && is_real(argument) {
            let hyperbolic = matches!(
                function,
                Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth
            );
            return Ok(match (hyperbolic, deep) {
                (true, true) => (self.expand_other(e, hinted, true)?, Expr::zero()),
                (true, false) => (e.clone(), Expr::zero()),
                (false, true) => {
                    let value = self.expand_other(argument, hinted, true)?;
                    (f::func(function.clone(), vec![value]), Expr::zero())
                }
                (false, false) => (f::func(function.clone(), vec![argument.clone()]), Expr::zero()),
            });
        }
        let (a, b) = if *function == Func::Exp {
            let (a, b) = self.split(argument, true)?;
            if deep {
                (self.expand_nested(&a, hinted)?, self.expand_nested(&b, hinted)?)
            } else {
                (a, b)
            }
        } else {
            let value = if deep {
                self.expand_nested(argument, hinted)?
            } else {
                argument.clone()
            };
            self.split(&value, true)?
        };
        if a.is_constant(Constant::NaN) || b.is_constant(Constant::NaN) {
            return Ok((Expr::nan(), Expr::nan()));
        }
        if matches!(
            b.kind(),
            Kind::Const(Constant::Infinity | Constant::NegInfinity)
        ) && *function == Func::Exp
        {
            return Ok((Expr::nan(), Expr::nan()));
        }
        if matches!(
            a.kind(),
            Kind::Const(Constant::Infinity | Constant::NegInfinity)
        ) || matches!(
            b.kind(),
            Kind::Const(Constant::Infinity | Constant::NegInfinity)
        ) {
            return Err(PartsError(
                "as_real_imag of a complex nonfinite function argument is not implemented".into(),
            ));
        }
        if is_zero(&b) {
            return Ok((f::func(function.clone(), vec![a]), Expr::zero()));
        }
        Ok(match function {
            Func::Exp => {
                let scale = f::exp(&a);
                (mul2(&scale, &f::cos(&b)), mul2(&scale, &f::sin(&b)))
            }
            Func::Sin => (
                mul2(&f::sin(&a), &f::cosh(&b)),
                mul2(&f::cos(&a), &f::sinh(&b)),
            ),
            Func::Cos => (
                mul2(&f::cos(&a), &f::cosh(&b)),
                neg(&mul2(&f::sin(&a), &f::sinh(&b))),
            ),
            Func::Sinh => (
                mul2(&f::sinh(&a), &f::cos(&b)),
                mul2(&f::cosh(&a), &f::sin(&b)),
            ),
            Func::Cosh => (
                mul2(&f::cosh(&a), &f::cos(&b)),
                mul2(&f::sinh(&a), &f::sin(&b)),
            ),
            Func::Tan | Func::Cot => {
                let twice_a = mul2(&Expr::int(2), &a);
                let twice_b = mul2(&Expr::int(2), &b);
                let denominator = if *function == Func::Tan {
                    add2(&f::cos(&twice_a), &f::cosh(&twice_b))
                } else {
                    sub(&f::cos(&twice_a), &f::cosh(&twice_b))
                };
                let real = div(&f::sin(&twice_a), &denominator);
                (
                    if *function == Func::Tan {
                        real
                    } else {
                        neg(&real)
                    },
                    div(&f::sinh(&twice_b), &denominator),
                )
            }
            Func::Tanh | Func::Coth => {
                let denominator = add2(
                    &powi(&f::sinh(&a), 2),
                    &powi(
                        &if *function == Func::Tanh {
                            f::cos(&b)
                        } else {
                            f::sin(&b)
                        },
                        2,
                    ),
                );
                let imag = div(&mul2(&f::sin(&b), &f::cos(&b)), &denominator);
                (
                    div(&mul2(&f::sinh(&a), &f::cosh(&a)), &denominator),
                    if *function == Func::Tanh {
                        imag
                    } else {
                        neg(&imag)
                    },
                )
            }
            _ => unreachable!(),
        })
    }
}

/// `Abs` with SymPy's final evaluation rule, which the expanded arguments
/// of a complex pass rely on: when conjugation resolves, and an unknown
/// symbol lacks its own conjugate, `Abs(z) = sqrt(expand_mul(z*conjugate(z)))`.
fn sympy_abs(argument: &Expr) -> Expr {
    let value = f::abs(argument);
    if !argument.is_add() || !matches!(value.kind(), Kind::Func(Func::Abs, _)) {
        return value;
    }
    abs_conjugate_rule(argument).unwrap_or(value)
}
/// The final rule of SymPy's `Abs.eval` for a non-real argument.
fn abs_conjugate_rule(argument: &Expr) -> Option<Expr> {
    fn conjugates(e: &Expr, out: &mut Vec<Expr>) {
        if e.is_func(&Func::Conjugate) && !out.contains(e) {
            out.push(e.clone());
        }
        for arg in e.args() {
            conjugates(&arg, out);
        }
    }
    let conjugate = f::conjugate(argument);
    let (mut old, mut new) = (Vec::new(), Vec::new());
    conjugates(argument, &mut old);
    conjugates(&conjugate, &mut new);
    new.retain(|c| !old.contains(c));
    if !new.is_empty()
        && new
            .iter()
            .all(|c| crate::traverse::has(argument, &c.args()[0]))
    {
        return None;
    }
    if conjugate == *argument || conjugate == neg(argument) {
        return None;
    }
    // Symbols under Abs count as real, as SymPy replaces Abs atoms by real dummies.
    fn unknowns(e: &Expr, out: &mut Vec<Expr>) {
        match e.kind() {
            Kind::Func(Func::Abs, _) => {}
            Kind::Sym(_) => {
                if crate::assumptions::fact(e, "extended_real").is_none() && !out.contains(e) {
                    out.push(e.clone());
                }
            }
            _ => {
                for arg in crate::traverse::ordered_args(e) {
                    unknowns(&arg, out);
                }
            }
        }
    }
    let mut unknown = Vec::new();
    unknowns(argument, &mut unknown);
    if unknown.is_empty()
        || !unknown
            .iter()
            .all(|u| crate::traverse::has(&conjugate, &f::conjugate(u)))
    {
        return Some(crate::ops::sqrt(&crate::simp::expand_mul_deep(
            &mul2(argument, &conjugate),
            true,
        )));
    }
    None
}
fn elementary(function: &Func) -> bool {
    matches!(
        function,
        Func::Exp
            | Func::Log
            | Func::Sin
            | Func::Cos
            | Func::Tan
            | Func::Cot
            | Func::Sinh
            | Func::Cosh
            | Func::Tanh
            | Func::Coth
    )
}
fn specialized_dispatch(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Add(..) | Kind::Mul(..))
        || matches!(e.kind(),Kind::Pow(base,exponent) if exponent.is_rational_number()||base.is_constant(Constant::E))
        || matches!(
            e.kind(),
            Kind::Func(
                Func::Re
                    | Func::Im
                    | Func::Exp
                    | Func::Log
                    | Func::Sin
                    | Func::Cos
                    | Func::Tan
                    | Func::Cot
                    | Func::Sinh
                    | Func::Cosh
                    | Func::Tanh
                    | Func::Coth
                    | Func::Sec
                    | Func::Csc,
                _
            )
        )
}
fn scalar(e: &Expr) -> Result<()> {
    if matches!(
        e.kind(),
        Kind::Bool(_)
            | Kind::Rel(..)
            | Kind::Func(Func::Logic(_) | Func::Set(_) | Func::MatrixProxy(_), _)
    ) || matches!(e.kind(),Kind::Func(Func::Undefined(head),_) if head.is_internal("Tuple"))
    {
        Err(PartsError(
            "as_real_imag requires a scalar expression".into(),
        ))
    } else {
        Ok(())
    }
}
/// Split a literal finite numeric complex value without rounding or evaluation.
pub fn literal_parts(e: &Expr) -> Option<Parts> {
    let finite = |n: &Number| !matches!(n, Number::Float(value) if !value.is_finite());
    if !e.is_canonical() { return None; }
    match e.kind() {
        Kind::Num(n) if finite(n) => Some((e.clone(), Expr::zero())),
        Kind::Const(Constant::I) => Some((Expr::zero(), Expr::one())),
        Kind::Mul(c, factors) if finite(c) && factors.len() == 1
            && factors[0].0 == Expr::i() && factors[0].1.is_one() => {
                Some((Expr::zero(), Expr::num(c.clone())))
            },
        Kind::Add(c, terms) if finite(c) && terms.len() == 1
            && terms[0].0 == Expr::i() && finite(&terms[0].1) => {
                Some((Expr::num(c.clone()), Expr::num(terms[0].1.clone())))
            },
        _ => None,
    }
}

/// SymPy's `e.as_real_imag(deep, ignore=ignore)`, without size limits.
pub fn as_real_imag(e: &Expr, deep: bool, ignore: Option<&Expr>) -> Result<Option<Parts>> {
    split_parts(e, deep, ignore, false)
}
/// `as_real_imag` under the opt-in resource guards (20000 nodes, 4096 terms,
/// integer powers through 128): for internal heuristics without a SymPy
/// counterpart, which treat the error as "undecided".
pub fn as_real_imag_bounded(e: &Expr, deep: bool, ignore: Option<&Expr>) -> Result<Option<Parts>> {
    split_parts(e, deep, ignore, true)
}
fn split_parts(e: &Expr, deep: bool, ignore: Option<&Expr>, bounded: bool) -> Result<Option<Parts>> {
    scalar(e)?;
    if let Some(parts) = literal_parts(e) {
        if e.is_atom() { return Ok(if ignore == Some(e) { None } else { Some(parts) }); }
        if ignore.is_none() { return Ok(Some(parts)); }
    }
    precheck(e, bounded)?;
    crate::assumptions::exact_queries(|| {
        let mut context = Context::new(bounded);
        if let Some(ignore) = ignore {
            if e.is_mul() {
                return context.product(e, true, Some(ignore));
            }
            if !specialized_dispatch(e) {
                let expanded = if deep && e.is_pow() {
                    context.expand(e, true, false)?
                } else {
                    e.clone()
                };
                if expanded == *ignore {
                    return Ok(None);
                }
            }
        }
        context.split(e, deep).map(Some)
    })
}
/// SymPy's `expand(e, deep, complex=True, **hints)`. `options` carries the
/// other hints, which run after the complex pass as in SymPy's hint order;
/// `expand_complex` disables all of them.
pub fn expand_complex(e: &Expr, options: &ExpandOptions) -> Result<Expr> {
    scalar(e)?;
    if literal_parts(e).is_some() { return Ok(e.clone()); }
    precheck(e, false)?;
    crate::assumptions::exact_queries(|| {
        Context::with_hints(*options, false).expand_full(e, options.deep)
    })
}
pub fn as_real_imag_many(values: &[Expr], deep: bool) -> Result<Vec<Parts>> {
    let mut seen = FxHashMap::default();
    for value in values {
        scalar(value)?;
        if literal_parts(value).is_none() { height(value, &mut seen, 0)?; }
    }
    crate::assumptions::exact_queries(|| {
        let mut context = Context::new(false);
        values.iter().map(|e| match literal_parts(e) {
            Some(parts) => Ok(parts), None => context.split(e, deep),
        }).collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numeric_literals_split_without_expanding_precision_or_exponents() {
        for exponent in [i64::MIN, i64::MAX] {
            let value = Expr::num(Number::MpFloat(crate::bigfloat::BigFloat {
                man: BigInt::one(), exp: exponent.into(),
            }, 65_536));
            assert_eq!(as_real_imag(&value, true, None).unwrap(), Some((value.clone(), Expr::zero())));
            assert_eq!(as_real_imag(&value, true, Some(&value)).unwrap(), None);
            assert_eq!(expand_complex(&value, &ExpandOptions::default()).unwrap(), value);
        }
        let real = Number::MpFloat(crate::bigfloat::BigFloat::from_ratio(&1.into(), &3.into(), 65_536), 65_536);
        let imaginary = Number::MpFloat(crate::bigfloat::BigFloat::from_ratio(&1.into(), &5.into(), 65_536), 65_536);
        let value = Expr::from_kind(Kind::Add(real.clone(), vec![(Expr::i(), imaginary.clone())]));
        let wanted = (Expr::num(real), Expr::num(imaginary));
        assert_eq!(as_real_imag(&value, true, None).unwrap(), Some(wanted.clone()));
        assert_eq!(as_real_imag_many(&[value.clone()], true).unwrap(), vec![wanted]);
        assert_eq!(expand_complex(&value, &ExpandOptions::default()).unwrap(), value);
    }
    #[test]
    fn exact_polynomial_reconstruction_and_opaque_product_grouping() {
        let r = Expr::symbol_with("r", crate::Assumptions::REAL);
        let t = Expr::symbol_with("t", crate::Assumptions::REAL);
        let z = add2(&r, &mul2(&Expr::i(), &t));
        let value = powi(&z, 7);
        let (a, b) = as_real_imag(&value, true, None).unwrap().unwrap();
        assert_eq!(
            crate::expand::expand(&sub(&add2(&a, &mul2(&Expr::i(), &b)), &value)),
            Expr::zero()
        );
        let product = mul2(&Expr::symbol("x"), &Expr::symbol("y"));
        let (a, b) = as_real_imag(&product, true, None).unwrap().unwrap();
        assert!(matches!(a.kind(),Kind::Func(Func::Re,args) if args==&[product.clone()]));
        assert!(matches!(b.kind(),Kind::Func(Func::Im,args) if args==&[product.clone()]));
    }
    #[test]
    fn projections_preserve_full_function_facts_and_nonfinite_constants() {
        let head = crate::expr::UndefinedFunction::with_facts(
            "f",
            crate::assumption_facts::FactSet::from_kwargs([("imaginary", true)]).unwrap(),
            None,
        )
        .unwrap();
        let value = Expr::from_kind(Kind::Func(Func::Undefined(head), vec![Expr::symbol("x")]));
        let (a, b) = as_real_imag(&value, true, None).unwrap().unwrap();
        assert!(a.is_zero());
        assert_eq!(b, mul(&[Expr::minus_one(), Expr::i(), value]));
        for value in [Expr::nan(), Expr::zoo()] {
            assert_eq!(
                as_real_imag(&value, true, None).unwrap(),
                Some((Expr::nan(), Expr::nan()))
            );
        }
        assert_eq!(
            as_real_imag(&Expr::oo(), true, None).unwrap(),
            Some((Expr::oo(), Expr::zero()))
        );
    }
    #[test]
    fn complex_logarithm_zeros_and_unsupported_oscillation_are_explicit() {
        for value in [
            Expr::zero(),
            Expr::num(Number::Float(-0.0)),
            Expr::num(Number::MpFloat(crate::bigfloat::BigFloat::zero(), 336)),
        ] {
            assert_eq!(f::log(&value), Expr::zoo());
            assert_eq!(
                as_real_imag(&f::log(&value), true, None).unwrap(),
                Some((Expr::nan(), Expr::nan()))
            );
        }
        let held = raw(Func::Log, &Expr::neg_oo());
        assert_eq!(
            as_real_imag(&held, true, None).unwrap(),
            Some((Expr::oo(), Expr::pi()))
        );
        let held = raw(Func::Sin, &Expr::oo());
        assert!(as_real_imag(&held, true, None)
            .unwrap_err()
            .0
            .contains("AccumBounds"));
    }
    #[test]
    fn shared_opaque_dags_are_linear_and_numeric_materialization_is_preflighted() {
        let mut value = Expr::symbol("x");
        for _ in 0..80 {
            value = Expr::from_kind(Kind::Func(
                Func::Undefined("pair".into()),
                vec![value.clone(), value],
            ));
        }
        let (a, b) = as_real_imag(&value, true, None).unwrap().unwrap();
        assert_eq!(a.args(), vec![value.clone()]);
        assert_eq!(b.args(), vec![value]);
        // A compact exponent must not request a gigantic pi for argument
        // reduction: the sine of a real number is its own real part.
        for exponent in [i64::MIN, i64::MAX] {
            let value = Expr::num(Number::MpFloat(
                crate::bigfloat::BigFloat {
                    man: BigInt::one(),
                    exp: exponent.into(),
                },
                53,
            ));
            let (_, b) = as_real_imag(&raw(Func::Sin, &value), true, None).unwrap().unwrap();
            assert!(b.is_zero());
        }
        let huge = Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2), Expr::int(1_000_000_000)));
        // a power of a real number stays whole, as in SymPy
        let (a, b) = as_real_imag(&huge, true, None).unwrap().unwrap();
        assert_eq!((a, b.is_zero()), (huge, true));
        let mut numeric = Expr::int(2);
        for _ in 0..80 {
            numeric = Expr::from_unevaluated_kind(Kind::Mul(
                Number::ONE,
                vec![(numeric.clone(), Expr::one()), (numeric, Expr::one())],
            ));
        }
        assert!(as_real_imag(&numeric, true, None).is_err());
    }
    #[test]
    fn complex_expansion_follows_sympy_hint_order() {
        let x = Expr::symbol("x");
        let (re, im) = (f::re(&x), f::im(&x));
        let only = ExpandOptions {
            basic: false, mul: false, multinomial: false, power_base: false,
            power_exp: false, log: false, ..ExpandOptions::default()
        };
        // expand_complex(x) = re(x) + I*im(x)
        assert_eq!(expand_complex(&x, &only).unwrap(), add2(&re, &mul2(&Expr::i(), &im)));
        // Abs(x) becomes sqrt(re(x)**2 + im(x)**2), as SymPy's Abs.eval does.
        let modulus = pow(&add2(&powi(&re, 2), &powi(&im, 2)), &Expr::rational(1, 2));
        assert_eq!(expand_complex(&f::abs(&x), &only).unwrap(), modulus);
        // With the log hint, log(sqrt(s)) is expanded after the complex pass.
        let value = expand_complex(&f::log(&x), &ExpandOptions::default()).unwrap();
        let half = mul2(&Expr::rational(1, 2), &f::log(&add2(&powi(&re, 2), &powi(&im, 2))));
        let angle = mul2(&Expr::i(), &f::arg(&add2(&re, &mul2(&Expr::i(), &im))));
        assert_eq!(value, add2(&half, &angle));
        // Real arguments only receive the non-complex hints.
        let r = Expr::symbol_with("r", crate::Assumptions::REAL);
        assert_eq!(expand_complex(&f::sin(&r), &only).unwrap(), f::sin(&r));
        // Unsupported structural nodes are rejected instead of guessed.
        let derivative = Expr::from_unevaluated_kind(Kind::Derivative(
            f::sin(&x), vec![(x.clone(), 1)],
        ));
        assert!(expand_complex(&derivative, &only).is_err());
    }
}
