//! Finite Markov chain kernels of `sympy.stats` (SymPy 1.14
//! `DiscreteMarkovChain` / `ContinuousMarkovChain`).
//!
//! Every kernel reproduces the reference algorithm step for step with the
//! native expression arithmetic (structural `!= 0`, the binary `+` of
//! `sum(row) - 1`), or, for the linear algebra, works over exact rationals
//! where the result is unique. A kernel that cannot decide returns `None`
//! and the caller keeps its general path.
use crate::matrix::domain::RatMat;
use crate::matrix::int::Rat;
use crate::matrix::Matrix;
use crate::{ops, Expr, Kind};

/// The communication classes of a chain, in SymPy's order.
#[derive(Clone, Debug, PartialEq)]
pub struct Classes {
    /// State indices of each class (Tarjan order).
    pub classes: Vec<Vec<usize>>,
    /// `sum(row) - 1 != 0` fails for every row of the class submatrix.
    pub recurrent: Vec<bool>,
    pub periods: Vec<u64>,
}

/// Entries are native expressions without host atoms or matrices, so the
/// structural comparisons below are the ones the Python operators perform.
pub fn admitted(t: &Matrix) -> bool {
    t.is_square() && t.data().iter().all(|e| !e.has_host() && !e.has_matrix())
}

fn nonzero(e: &Expr) -> bool {
    // SymPy `T[i, j] != 0`: structural (Float(0.0) != 0 is True).
    *e != Expr::zero()
}

const LINEAR_PROBES: usize = 9;
const PERTURB_SHIFT: u32 = 5;

/// Iteration order of a CPython (3.7+) `set` of non-negative small ints
/// (`hash(i) == i`) built empty and filled by `add` in `items` order: the
/// set comprehension of the reference period search. The reference result
/// depends on this order (its frontier is the successors of the last state
/// iterated).
fn set_order(items: &[usize]) -> Vec<usize> {
    fn probe(table: &mut [Option<usize>], mask: usize, key: usize, check_equal: bool) -> bool {
        let mut perturb = key;
        let mut i = key & mask;
        loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            for j in 0..=probes {
                match table[i + j] {
                    None => {
                        table[i + j] = Some(key);
                        return true;
                    }
                    Some(k) if check_equal && k == key => return false,
                    _ => {}
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb)) & mask;
        }
    }
    let mut mask = 7usize;
    let mut table: Vec<Option<usize>> = vec![None; 8];
    let mut used = 0usize;
    for &key in items {
        if !probe(&mut table, mask, key, true) {
            continue;
        }
        used += 1;
        // No deletions: fill == used.
        if used * 5 >= mask * 3 {
            let minused = if used > 50000 { used * 2 } else { used * 4 };
            let mut size = 8usize;
            while size <= minused {
                size <<= 1;
            }
            let old = std::mem::replace(&mut table, vec![None; size]);
            mask = size - 1;
            for key in old.into_iter().flatten() {
                probe(&mut table, mask, key, false);
            }
        }
    }
    table.into_iter().flatten().collect()
}

fn igcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

/// The period of a strongly connected class: the gcd of
/// `level(u) + 1 - level(v)` over its edges, with breadth-first levels.
fn true_period(t: &Matrix, class: &[usize]) -> u64 {
    let n = t.rows();
    let mut level = vec![usize::MAX; n];
    level[class[0]] = 0;
    let mut queue = std::collections::VecDeque::from([class[0]]);
    while let Some(u) = queue.pop_front() {
        for &v in class {
            if nonzero(t.get(u, v)) && level[v] == usize::MAX {
                level[v] = level[u] + 1;
                queue.push_back(v);
            }
        }
    }
    let mut g = 0u64;
    for &u in class {
        for &v in class {
            if nonzero(t.get(u, v)) {
                g = igcd(g, (level[u] as i64 + 1 - level[v] as i64).unsigned_abs());
            }
        }
    }
    if g == 0 { class.len() as u64 } else { g }
}

/// The reference breadth-first period search of one class. Where the
/// reference loop never terminates (its frontier keeps the successors of
/// its last element only and cycles without reaching the whole class), the
/// period is computed by `true_period` (an intentional correction: SymPy
/// 1.14 does not return).
fn period(t: &Matrix, class: &[usize]) -> Option<u64> {
    let n = t.rows();
    let mut visited = vec![false; n];
    let mut visited_count = 1usize;
    let mut level = vec![0i64; n];
    visited[class[0]] = true;
    let mut frontier: Vec<usize> = vec![class[0]];
    let mut values: Vec<i64> = Vec::new();
    let mut current_level = 0i64;
    let mut done = false;
    let mut stalled = 0usize;
    while !done {
        done = visited_count == class.len();
        current_level += 1;
        let before = visited_count;
        let mut next = frontier.clone();
        for &i in &frontier {
            let successors: Vec<usize> = class.iter().copied().filter(|&j| nonzero(t.get(i, j))).collect();
            let successors = set_order(&successors);
            for &j in &successors {
                if !visited[j] {
                    level[j] = current_level;
                }
            }
            for &j in &successors {
                if visited[j] {
                    let v = level[i] - level[j] + 1;
                    if !values.contains(&v) {
                        values.push(v);
                    }
                }
            }
            for &j in &successors {
                if !visited[j] {
                    visited[j] = true;
                    visited_count += 1;
                }
            }
            next = successors;
        }
        frontier = next;
        if visited_count == before && !done {
            stalled += 1;
            // With `visited` fixed the frontier is the successor set of one
            // state: more than n stalled rounds repeat a frontier forever.
            if stalled > n + 1 {
                return Some(true_period(t, class));
            }
        } else {
            stalled = 0;
        }
    }
    let positive: Vec<u64> = values.iter().filter(|&&v| v > 0).map(|&v| v as u64).collect();
    Some(match positive.len() {
        0 => class.len() as u64,
        1 => positive[0],
        _ => positive.iter().fold(0, |g, &v| igcd(g, v)),
    })
}

/// `DiscreteMarkovChain.communication_classes` (SymPy 1.14): Tarjan's
/// components over the structural nonzero pattern, the recurrence test of
/// each class submatrix and the reference period search.
pub fn communication_classes(t: &Matrix, check: &mut dyn FnMut() -> bool) -> Option<Classes> {
    if !admitted(t) {
        return None;
    }
    let n = t.rows();
    let gmap: Vec<Vec<usize>> = (0..n).map(|i| (0..n).filter(|&j| nonzero(t.get(i, j))).collect()).collect();
    let classes = crate::matrix::linalg::tarjan(&gmap);
    let mut recurrent = Vec::with_capacity(classes.len());
    let mut periods = Vec::with_capacity(classes.len());
    for class in &classes {
        if !check() {
            return None;
        }
        let mut is_recurrent = true;
        for &i in class {
            // Python `sum(row) - 1`: 0 + row[0] + row[1] + ... - 1.
            let mut total = Expr::zero();
            for &j in class {
                total = ops::add2(&total, t.get(i, j));
            }
            if ops::sub(&total, &Expr::one()) != Expr::zero() {
                is_recurrent = false;
                break;
            }
        }
        recurrent.push(is_recurrent);
        periods.push(period(t, class)?);
    }
    Some(Classes { classes, recurrent, periods })
}

/// `decompose`: the recurrent states, then the transient states, and the
/// blocks `A` (recurrent to recurrent), `B` (transient to recurrent) and
/// `C` (transient to transient).
pub fn decompose(t: &Matrix, c: &Classes) -> (Vec<usize>, Matrix, Matrix, Matrix) {
    let mut r = Vec::new();
    let mut tr = Vec::new();
    for (class, &rec) in c.classes.iter().zip(&c.recurrent) {
        if rec {
            r.extend_from_slice(class);
        } else {
            tr.extend_from_slice(class);
        }
    }
    let a = t.extract(&r, &r).expect("class indices are in range");
    let b = t.extract(&tr, &r).expect("class indices are in range");
    let cc = t.extract(&tr, &tr).expect("class indices are in range");
    r.extend_from_slice(&tr);
    (r, a, b, cc)
}

/// `canonical_form`: `BlockMatrix([[A, 0], [B, C]]).as_explicit()`.
pub fn canonical(a: &Matrix, b: &Matrix, c: &Matrix) -> Matrix {
    let nr = a.rows();
    let n = nr + c.rows();
    Matrix::from_fn(n, n, |i, j| match (i < nr, j < nr) {
        (true, true) => a.get(i, j).clone(),
        (true, false) => Expr::zero(),
        (false, true) => b.get(i - nr, j).clone(),
        (false, false) => c.get(i - nr, j - nr).clone(),
    })
}

/// The unique solution of the exact rational system `m x = rhs`, `None`
/// when the system is inconsistent or has free parameters.
fn unique_solution(m: &RatMat, rhs: &[Rat]) -> Option<Vec<Rat>> {
    let (rows, cols) = (m.rows, m.cols);
    let mut data = Vec::with_capacity(rows * (cols + 1));
    for i in 0..rows {
        for j in 0..cols {
            data.push(m.get(i, j).clone());
        }
        data.push(rhs[i].clone());
    }
    let aug = RatMat { rows, cols: cols + 1, data };
    let (red, pivots) = aug.rref();
    if pivots.len() != cols || pivots.iter().enumerate().any(|(k, &p)| p != k) {
        return None;
    }
    Some((0..cols).map(|i| red.get(i, cols).clone()).collect())
}

fn row_matrix(values: Vec<Rat>) -> Matrix {
    Matrix::row_vector(values.iter().map(|v| Expr::num(v.to_number())).collect())
}

/// `DiscreteMarkovChain.stationary_distribution()` for exact rational
/// transition matrices with a unique solution: `(P - I)^T` with its first row
/// replaced by ones, right-hand side `e_0` (the reference linsolve system).
pub fn stationary_rational(t: &Matrix) -> Option<Matrix> {
    if !t.is_square() || t.rows() == 0 {
        return None;
    }
    let p = RatMat::from_matrix(t)?;
    let n = p.rows;
    let mut data = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            data.push(if i == 0 {
                Rat::one()
            } else {
                let v = p.get(j, i).clone();
                if i == j { v.sub(&Rat::one()) } else { v }
            });
        }
    }
    let mut rhs = vec![Rat::zero(); n];
    rhs[0] = Rat::one();
    unique_solution(&RatMat { rows: n, cols: n, data }, &rhs).map(row_matrix)
}

/// `ContinuousMarkovChain.limiting_distribution()` for an exact rational
/// generator with a unique solution of `w G = 0`, `sum(w) = 1`.
pub fn limiting_rational(g: &Matrix) -> Option<Matrix> {
    if !g.is_square() || g.rows() == 0 {
        return None;
    }
    let q = RatMat::from_matrix(g)?;
    let n = q.rows;
    let mut data = Vec::with_capacity((n + 1) * n);
    for j in 0..n {
        for i in 0..n {
            data.push(q.get(i, j).clone());
        }
    }
    data.extend((0..n).map(|_| Rat::one()));
    let mut rhs = vec![Rat::zero(); n + 1];
    rhs[n] = Rat::one();
    unique_solution(&RatMat { rows: n + 1, cols: n, data }, &rhs).map(row_matrix)
}

/// Outcome of the reference `list(linsolve(...))[0]` on a stationary or
/// limiting system.
#[derive(Clone, Debug, PartialEq)]
pub enum Linear {
    /// The unique solution, one value per unknown.
    Unique(Vec<Expr>),
    /// `linsolve` is the EmptySet: the reference raises IndexError.
    Inconsistent,
}

/// The augmented reference system of `stationary_distribution`: the
/// equations of `(P - I)^T` with the first row replaced by ones, `b = e_0`.
pub fn stationary_system(t: &Matrix) -> Vec<Vec<Expr>> {
    let n = t.rows();
    let one = Expr::one();
    (0..n)
        .map(|i| {
            let mut row: Vec<Expr> = if i == 0 {
                vec![one.clone(); n]
            } else {
                (0..n).map(|j| if i == j { ops::sub(t.get(j, i), &one) } else { t.get(j, i).clone() }).collect()
            };
            row.push(if i == 0 { one.clone() } else { Expr::zero() });
            row
        })
        .collect()
}

/// The augmented reference system of the CTMC `limiting_distribution`: the
/// entries of `Matrix([w]) * G` (coefficient of `w_k` in equation `j` is
/// `G[k, j]`) followed by `sum(w) - 1`.
pub fn limiting_system(g: &Matrix) -> Vec<Vec<Expr>> {
    let n = g.rows();
    let mut rows: Vec<Vec<Expr>> = (0..n)
        .map(|j| {
            let mut row: Vec<Expr> = (0..n).map(|k| g.get(k, j).clone()).collect();
            row.push(Expr::zero());
            row
        })
        .collect();
    let mut last = vec![Expr::one(); n];
    last.push(Expr::one());
    rows.push(last);
    rows
}

use crate::bigfloat::BigFloat;
use crate::number::Number;

/// The precision of `RR` when every entry is a number and one is a float
/// (SymPy's `construct_domain`: the largest float precision).
fn real_precision(rows: &[Vec<Expr>]) -> Option<u32> {
    let mut precision = None;
    for e in rows.iter().flatten() {
        let n = e.as_number()?;
        if let Some(p) = n.float_precision() {
            precision = Some(precision.map_or(p, |q: u32| q.max(p)));
        }
    }
    precision
}

/// `solve_lin_sys` over `RR` (SymPy 1.14): equations without coefficients
/// are dropped (inconsistent when their constant is nonzero), then the dense
/// Gauss-Jordan elimination of `ddm_irref` with partial pivoting
/// (`max(range(i, m), key=abs)`, the first maximal row), with the row
/// normalised by the reciprocal of its pivot and the elimination
/// `ak[l] -= akj * ai[l]`, all in mpmath arithmetic at the domain precision.
/// `None` when the solution is not unique (free parameters).
pub fn real_solve(rows: &[Vec<Expr>]) -> Option<Linear> {
    let prec = real_precision(rows)?;
    let width = rows.first()?.len();
    let n = width - 1;
    let convert = |e: &Expr| -> BigFloat {
        match e.as_number().expect("number") {
            Number::Float(f) => BigFloat::from_f64(*f),
            Number::MpFloat(v, _) => v.clone(),
            other => other.to_bigfloat(prec),
        }
    };
    let mut a: Vec<Vec<BigFloat>> = Vec::with_capacity(rows.len());
    for row in rows {
        let values: Vec<BigFloat> = row.iter().map(convert).collect();
        if values[..n].iter().all(BigFloat::is_zero) {
            if values[n].is_zero() {
                continue;
            }
            return Some(Linear::Inconsistent);
        }
        a.push(values);
    }
    let m = a.len();
    let mut pivots = Vec::new();
    let mut i = 0;
    for j in 0..width {
        let mut ip = i;
        for r in i + 1..m {
            if a[r][j].abs().cmp(&a[ip][j].abs()) == std::cmp::Ordering::Greater {
                ip = r;
            }
        }
        a.swap(i, ip);
        if a[i][j].is_zero() {
            continue;
        }
        let inverse = BigFloat::one().div(&a[i][j], prec)?;
        for l in j..width {
            a[i][l] = a[i][l].mul(&inverse, prec);
        }
        let pivot_row = a[i].clone();
        for (k, ak) in a.iter_mut().enumerate() {
            if k == i || ak[j].is_zero() {
                continue;
            }
            let akj = ak[j].clone();
            ak[j] = BigFloat::zero();
            for l in j + 1..width {
                ak[l] = ak[l].sub(&akj.mul(&pivot_row[l], prec), prec);
            }
        }
        pivots.push(j);
        i += 1;
        if i >= m {
            break;
        }
    }
    if pivots.last() == Some(&n) {
        return Some(Linear::Inconsistent);
    }
    if pivots.len() != n {
        return None;
    }
    // `solve_lin_sys(_raw=False)` gives an exact zero as the integer 0.
    Some(Linear::Unique(
        (0..n).map(|r| if a[r][n].is_zero() { Expr::zero() } else { Expr::num(Number::mp_float(a[r][n].clone(), prec)) }).collect(),
    ))
}

use crate::matrix::mpoly::MPoly;
use crate::matrix::ratfunc::{as_numer_denom, cancel_polys, cancel_qpolys, PolyCtx};

type Frac = (MPoly, MPoly);

fn frac_add(x: &Frac, y: &Frac) -> Frac {
    cancel_polys(&x.0.mul(&y.1).add(&y.0.mul(&x.1)), &x.1.mul(&y.1))
}

fn frac_mul(x: &Frac, y: &Frac) -> Frac {
    cancel_polys(&x.0.mul(&y.0), &x.1.mul(&y.1))
}

/// `solve_lin_sys` over the fraction field `QQ(a, b, ...)` of the symbols
/// of the coefficients: every element is kept as SymPy's `FracElement`
/// (`PolyElement.cancel`: gcd removed, positive leading coefficient of the
/// denominator in the sorted generators) and the unique reduced row echelon
/// form gives the values, converted by `numer.as_expr()/denom.as_expr()`.
/// `None` outside rational functions of symbols, for floats, and when the
/// solution is not unique.
pub fn fraction_solve<E>(rows: &[Vec<Expr>], check: &mut dyn FnMut() -> Result<(), E>) -> Result<Option<Linear>, E> {
    let Some(width) = rows.first().map(Vec::len) else { return Ok(None) };
    let n = width - 1;
    let parts: Vec<(Expr, Expr)> = rows.iter().flatten().map(as_numer_denom).collect();
    let exprs: Vec<Expr> = parts.iter().flat_map(|(p, q)| [p.clone(), q.clone()]).collect();
    let Some((ctx, polys)) = PolyCtx::from_exprs(&exprs, false) else { return Ok(None) };
    if !ctx.gens.iter().all(Expr::is_symbol) {
        return Ok(None);
    }
    let mut values = Vec::with_capacity(parts.len());
    for pair in polys.chunks(2) {
        if pair[1].p.is_zero() {
            return Ok(None);
        }
        values.push(cancel_qpolys(&pair[0], &pair[1]));
    }
    let mut a: Vec<Vec<Frac>> = values.chunks(width).map(<[Frac]>::to_vec).collect();
    let m = a.len();
    let mut pivots = Vec::new();
    let mut r = 0;
    for j in 0..width {
        check()?;
        let Some(p) = (r..m).find(|&k| !a[k][j].0.is_zero()) else { continue };
        a.swap(r, p);
        let inverse = cancel_polys(&a[r][j].1, &a[r][j].0);
        let pivot_row: Vec<Frac> = a[r].iter().map(|v| frac_mul(v, &inverse)).collect();
        a[r] = pivot_row.clone();
        for (k, ak) in a.iter_mut().enumerate() {
            if k == r || ak[j].0.is_zero() {
                continue;
            }
            let factor = ak[j].clone();
            let negated = (factor.0.neg(), factor.1.clone());
            for l in j..width {
                ak[l] = frac_add(&ak[l], &frac_mul(&negated, &pivot_row[l]));
            }
        }
        pivots.push(j);
        r += 1;
        if r >= m {
            break;
        }
    }
    if pivots.last() == Some(&n) {
        return Ok(Some(Linear::Inconsistent));
    }
    if pivots.len() != n {
        return Ok(None);
    }
    Ok(Some(Linear::Unique((0..n).map(|k| ops::div(&ctx.to_expr(&a[k][n].0), &ctx.to_expr(&a[k][n].1))).collect())))
}

/// The reference `list(linsolve(system))[0]` of a stationary or limiting
/// system with float or symbolic entries; `None` when the kernel declines.
pub fn solve_system<E>(rows: &[Vec<Expr>], check: &mut dyn FnMut() -> Result<(), E>) -> Result<Option<Linear>, E> {
    if rows.iter().flatten().any(|e| e.has_host() || e.has_matrix()) {
        return Ok(None);
    }
    if real_precision(rows).is_some() {
        return Ok(real_solve(rows));
    }
    // Floats inside symbolic entries (RR(a, b)) are declined by PolyCtx.
    fraction_solve(rows, check)
}

/// Outcome of the exact fundamental matrix.
pub enum Fundamental {
    Matrix(Matrix),
    /// The reference raises ValueError("The fundamental matrix doesn't exist.").
    Singular,
}

/// `fundamental_matrix` for exact rational transition matrices:
/// `(I - C)^-1` when there are transient states, else
/// `(I - P + W)^-1` with the rows of `W` the (unique) stationary vector.
pub fn fundamental_rational(t: &Matrix, c: &Classes) -> Option<Fundamental> {
    if !t.data().iter().all(Expr::is_rational_number) {
        return None;
    }
    let (_, _, _, q) = decompose(t, c);
    let n = t.rows();
    let (m, size) = if q.rows() > 0 {
        let q = RatMat::from_matrix(&q)?;
        let k = q.rows;
        let data = (0..k * k)
            .map(|x| {
                let (i, j) = (x / k, x % k);
                let v = q.get(i, j).neg();
                if i == j { v.add(&Rat::one()) } else { v }
            })
            .collect();
        (RatMat { rows: k, cols: k, data }, k)
    } else {
        let w = stationary_rational(t)?;
        let w = RatMat::from_matrix(&w)?;
        let p = RatMat::from_matrix(t)?;
        let data = (0..n * n)
            .map(|x| {
                let (i, j) = (x / n, x % n);
                let v = w.get(0, j).sub(p.get(i, j));
                if i == j { v.add(&Rat::one()) } else { v }
            })
            .collect();
        (RatMat { rows: n, cols: n, data }, n)
    };
    if size == 0 {
        // A 0x0 chain: the reference inverts the empty matrix.
        return Some(Fundamental::Matrix(Matrix::zeros(0, 0)));
    }
    Some(match m.inv() {
        Some(inv) => Fundamental::Matrix(inv.to_matrix()),
        None => Fundamental::Singular,
    })
}

/// `absorbing_probabilities` (`N * R`) for exact rational chains with
/// transient states; `None` otherwise (the reference product of the ergodic
/// branch is a shape error, raised by the caller's matrix product).
pub fn absorbing_rational(t: &Matrix, c: &Classes) -> Option<Fundamental> {
    let (_, _, r, q) = decompose(t, c);
    if q.rows() == 0 {
        return None;
    }
    Some(match fundamental_rational(t, c)? {
        Fundamental::Matrix(nm) => {
            let nm = RatMat::from_matrix(&nm)?;
            let r = RatMat::from_matrix(&r)?;
            Fundamental::Matrix(nm.matmul(&r).to_matrix())
        }
        Fundamental::Singular => Fundamental::Singular,
    })
}

/// `true`/`false` atoms.
pub fn truth(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Bool(b) => Some(*b),
        _ => None,
    }
}

/// One P/E query over the paths of a finite chain started in `initial`.
pub struct PathQuery<'a> {
    /// The transition kernel of each gap between consecutive observation times.
    pub kernels: &'a [Matrix],
    pub initial: usize,
    /// Value substituted for the observation of state `i` in the gate and in
    /// a probability event (the state or its code).
    pub coded: &'a [Expr],
    /// Value substituted in an expectation (the state itself).
    pub labels: &'a [Expr],
    /// Each observed random symbol with the index of its observation time.
    pub atoms: &'a [(Expr, usize)],
    pub gate: Option<&'a Expr>,
    pub event: &'a Expr,
    pub expectation: bool,
}

/// The reference enumeration: for each path `(initial,) + tail` with the
/// tails in `itertools.product` order, `weight = 1 * K0[a0, a1] * ...`,
/// `denominator += weight * [gate]`, and `numerator += weight * [gate] *
/// [event]` (probability) or `weight * [gate] * expr(path)` (expectation),
/// with the same binary operations. `Ok(None)` when an indicator is not
/// decided (the caller keeps its Piecewise form).
pub fn enumerate_paths<E>(
    q: &PathQuery<'_>,
    subs: &mut dyn FnMut(&Expr, &[(Expr, Expr)]) -> Result<Expr, E>,
    check: &mut dyn FnMut() -> Result<(), E>,
) -> Result<Option<(Expr, Expr, bool)>, E> {
    let n = q.coded.len();
    let k = q.kernels.len();
    if n == 0 || q.labels.len() != n || q.initial >= n || q.kernels.iter().any(|m| m.shape() != (n, n)) || q.atoms.iter().any(|(_, s)| *s > k) {
        return Ok(None);
    }
    let one = Expr::one();
    let mut numerator = Expr::zero();
    let mut denominator = Expr::zero();
    let mut path = vec![q.initial; k + 1];
    for slot in path.iter_mut().skip(1) {
        *slot = 0;
    }
    let mut count = 0usize;
    let mut certain = true;
    loop {
        count += 1;
        if count % 256 == 0 {
            check()?;
        }
        let mut weight = one.clone();
        for (s, kernel) in q.kernels.iter().enumerate() {
            weight = ops::mul2(&weight, kernel.get(path[s], path[s + 1]));
        }
        let coded: Vec<(Expr, Expr)> = q.atoms.iter().map(|(a, s)| (a.clone(), q.coded[path[*s]].clone())).collect();
        let indicator = match q.gate {
            None => one.clone(),
            Some(gate) => match truth(&subs(gate, &coded)?) {
                Some(b) => {
                    certain &= b;
                    Expr::int(b as i64)
                }
                None => return Ok(None),
            },
        };
        let gated = ops::mul2(&weight, &indicator);
        denominator = ops::add2(&denominator, &gated);
        let term = if q.expectation {
            let labels: Vec<(Expr, Expr)> = q.atoms.iter().map(|(a, s)| (a.clone(), q.labels[path[*s]].clone())).collect();
            subs(q.event, &labels)?
        } else {
            match truth(&subs(q.event, &coded)?) {
                Some(b) => Expr::int(b as i64),
                None => return Ok(None),
            }
        };
        numerator = ops::add2(&numerator, &ops::mul2(&gated, &term));
        // Next tail in lexicographic order (last position fastest).
        let mut pos = k;
        loop {
            if pos == 0 {
                return Ok(Some((numerator, denominator, certain)));
            }
            path[pos] += 1;
            if path[pos] < n {
                break;
            }
            path[pos] = 0;
            pos -= 1;
        }
    }
}

/// The value of a query from its path sums, as SymPy assembles it. When the
/// gate holds on every path (the condition is the initial observation, of
/// probability one by construction) SymPy returns `gprob * sum(K[gstate, s]
/// for s in states)` (or the weighted sum of an expectation) without
/// normalising or simplifying, which is the numerator, whatever the row sums
/// are (Float rows summing to 1.0000000000000002, unnormalised symbolic
/// rows). `None` when the caller must form the simplified quotient.
pub fn assemble(numerator: &Expr, certain: bool) -> Option<Expr> {
    certain.then(|| numerator.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rat(p: i64, q: i64) -> Expr {
        ops::div(&Expr::int(p), &Expr::int(q))
    }

    fn m(rows: &[&[(i64, i64)]]) -> Matrix {
        Matrix::from_rows(rows.iter().map(|r| r.iter().map(|&(p, q)| rat(p, q)).collect()).collect()).unwrap()
    }

    #[test]
    fn cpython_set_order_matches_small_tables() {
        assert_eq!(set_order(&[3, 1, 2]), vec![1, 2, 3]);
        // 9 collides with 1 in the 8-slot table and takes slot 2.
        assert_eq!(set_order(&[9, 2]), vec![9, 2]);
        assert_eq!(set_order(&[1, 9]), vec![1, 9]);
        // Five entries resize to 32 slots: sorted again.
        assert_eq!(set_order(&[9, 1, 2, 3, 4]), vec![1, 2, 3, 4, 9]);
    }

    #[test]
    fn classes_of_reference_examples() {
        let t = m(&[
            &[(0, 1), (1, 2), (1, 2), (0, 1), (0, 1)],
            &[(0, 1), (0, 1), (0, 1), (1, 1), (0, 1)],
            &[(1, 2), (0, 1), (1, 2), (0, 1), (0, 1)],
            &[(0, 1), (1, 1), (0, 1), (0, 1), (0, 1)],
            &[(0, 1), (3, 10), (0, 1), (3, 10), (2, 5)],
        ]);
        let c = communication_classes(&t, &mut || true).unwrap();
        assert_eq!(c.classes, vec![vec![1, 3], vec![0, 2], vec![4]]);
        assert_eq!(c.recurrent, vec![true, false, false]);
        assert_eq!(c.periods, vec![2, 1, 1]);
        let (states, a, b, cc) = decompose(&t, &c);
        assert_eq!(states, vec![1, 3, 0, 2, 4]);
        assert_eq!(a.shape(), (2, 2));
        assert_eq!(b.shape(), (3, 2));
        assert_eq!(cc.shape(), (3, 3));
        assert_eq!(a.data(), &[Expr::zero(), Expr::one(), Expr::one(), Expr::zero()]);
    }

    #[test]
    fn non_terminating_reference_search_gets_the_true_period() {
        let t = m(&[
            &[(0, 1), (1, 2), (1, 2), (0, 1), (0, 1)],
            &[(0, 1), (0, 1), (0, 1), (1, 1), (0, 1)],
            &[(1, 1), (0, 1), (0, 1), (0, 1), (0, 1)],
            &[(0, 1), (0, 1), (0, 1), (0, 1), (1, 1)],
            &[(1, 1), (0, 1), (0, 1), (0, 1), (0, 1)],
        ]);
        let c = communication_classes(&t, &mut || true).unwrap();
        assert_eq!(c.classes.len(), 1);
        assert_eq!(c.periods, vec![2]);
    }

    #[test]
    fn stationary_and_fundamental() {
        let t = m(&[&[(1, 5), (2, 5), (2, 5)], &[(1, 10), (1, 2), (2, 5)], &[(3, 5), (3, 10), (1, 10)]]);
        let w = stationary_rational(&t).unwrap();
        assert_eq!(w.data(), &[rat(11, 39), rat(16, 39), rat(4, 13)]);
        let t6 = m(&[
            &[(1, 1), (0, 1), (0, 1), (0, 1), (0, 1)],
            &[(1, 2), (0, 1), (1, 2), (0, 1), (0, 1)],
            &[(0, 1), (1, 2), (0, 1), (1, 2), (0, 1)],
            &[(0, 1), (0, 1), (1, 2), (0, 1), (1, 2)],
            &[(0, 1), (0, 1), (0, 1), (0, 1), (1, 1)],
        ]);
        let c = communication_classes(&t6, &mut || true).unwrap();
        let Some(Fundamental::Matrix(nm)) = fundamental_rational(&t6, &c) else { panic!() };
        assert_eq!(nm.data()[0..3], [rat(3, 2), rat(1, 1), rat(1, 2)]);
        let Some(Fundamental::Matrix(ab)) = absorbing_rational(&t6, &c) else { panic!() };
        assert_eq!(ab.data(), &[rat(3, 4), rat(1, 4), rat(1, 2), rat(1, 2), rat(1, 4), rat(3, 4)]);
        let g = m(&[&[(-1, 1), (1, 1)], &[(1, 1), (-1, 1)]]);
        let x = Expr::symbol("x");
        let gate = crate::relational::eq(&x, &Expr::int(0));
        let kernels = [g.clone()];
        let states = [Expr::int(0), Expr::int(1)];
        let q = PathQuery { kernels: &kernels, initial: 0, coded: &states, labels: &states, atoms: &[(x.clone(), 1)], gate: None, event: &gate, expectation: false };
        let mut subs = |e: &Expr, p: &[(Expr, Expr)]| -> Result<Expr, ()> { Ok(crate::subs::simultaneous(e, p)) };
        let (num, den, _) = enumerate_paths(&q, &mut subs, &mut || Ok(())).unwrap().unwrap();
        assert_eq!((num, den), (Expr::int(-1), Expr::zero()));
        assert_eq!(limiting_rational(&g).unwrap().data(), &[rat(1, 2), rat(1, 2)]);
    }

    fn floats(v: &[Expr]) -> Vec<f64> {
        v.iter().map(|e| e.as_number().map_or(f64::NAN, |n| n.to_f64())).collect()
    }

    #[test]
    fn float_systems_follow_the_partial_pivot_rref() {
        let t = Matrix::from_rows(vec![vec![Expr::float(0.3), Expr::float(0.7)], vec![Expr::float(0.9), Expr::float(0.1)]]).unwrap();
        let Ok(Some(Linear::Unique(w))) = solve_system::<()>(&stationary_system(&t), &mut || Ok(())) else { panic!() };
        assert_eq!(floats(&w), vec![0.5625, 0.4375]);
        let g = Matrix::from_rows(vec![vec![Expr::float(-0.75), Expr::float(0.75)], vec![Expr::float(0.125), Expr::float(-0.125)]]).unwrap();
        let Ok(Some(Linear::Unique(w))) = solve_system::<()>(&limiting_system(&g), &mut || Ok(())) else { panic!() };
        assert_eq!(floats(&w), vec![0.1428571428571429, 0.8571428571428571]);
        // A pivot in the right-hand side column: SymPy's EmptySet.
        let rows = vec![vec![Expr::float(1.0), Expr::float(1.0)], vec![Expr::float(0.0), Expr::float(2.0)]];
        assert_eq!(real_solve(&rows), Some(Linear::Inconsistent));
        // Non-unique: declined to the parametric fallback.
        let rows = vec![vec![Expr::float(1.0), Expr::float(1.0), Expr::float(1.0)]];
        assert_eq!(real_solve(&rows), None);
    }

    #[test]
    fn symbolic_systems_give_cancelled_quotients() {
        let (a, b) = (Expr::symbol("a"), Expr::symbol("b"));
        let g = Matrix::from_rows(vec![vec![ops::neg(&a), a.clone()], vec![b.clone(), ops::neg(&b)]]).unwrap();
        let Ok(Some(Linear::Unique(w))) = solve_system::<()>(&limiting_system(&g), &mut || Ok(())) else { panic!() };
        let s = ops::add2(&a, &b);
        assert_eq!(w, vec![ops::div(&b, &s), ops::div(&a, &s)]);
    }

    #[test]
    fn certain_gates_keep_the_unnormalised_sum() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let k = Matrix::from_rows(vec![vec![Expr::float(0.3), Expr::float(0.7)], vec![Expr::float(0.9), Expr::float(0.1)]]).unwrap();
        let kernels = [k.clone()];
        let states = [Expr::int(0), Expr::int(1)];
        let gate = crate::relational::eq(&x, &Expr::int(0));
        let event = crate::relational::eq(&y, &Expr::int(1));
        let atoms = [(x.clone(), 0), (y.clone(), 1)];
        let q = PathQuery { kernels: &kernels, initial: 0, coded: &states, labels: &states, atoms: &atoms, gate: Some(&gate), event: &event, expectation: false };
        let mut subs = |e: &Expr, p: &[(Expr, Expr)]| -> Result<Expr, ()> { Ok(crate::subs::simultaneous(e, p)) };
        let (num, _, certain) = enumerate_paths(&q, &mut subs, &mut || Ok(())).unwrap().unwrap();
        assert!(certain);
        assert_eq!(assemble(&num, certain), Some(num.clone()));
        assert_eq!(floats(&[num]), vec![0.7]);
        let gate = crate::relational::eq(&y, &Expr::int(0));
        let q = PathQuery { gate: Some(&gate), ..q };
        let (num, _, certain) = enumerate_paths(&q, &mut subs, &mut || Ok(())).unwrap().unwrap();
        assert!(!certain);
        assert_eq!(assemble(&num, certain), None);
    }
}
