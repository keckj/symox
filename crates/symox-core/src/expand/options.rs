//! Checked expansion options and a shared allocation/work budget.
use crate::expr::{Expr, Kind};
use rustc_hash::FxHashMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct ExpandOptions {
    pub deep: bool,
    pub basic: bool,
    pub mul: bool,
    pub multinomial: bool,
    pub power_base: bool,
    pub power_exp: bool,
    pub log: bool,
    pub trig: bool,
    pub func: bool,
    pub force: bool,
}
impl Default for ExpandOptions {
    fn default() -> Self {
        Self {
            deep: true,
            basic: true,
            mul: true,
            multinomial: true,
            power_base: true,
            power_exp: true,
            log: true,
            trig: false,
            func: false,
            force: false,
        }
    }
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExpandError(pub String);
impl std::fmt::Display for ExpandError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for ExpandError {}
type Result<T> = std::result::Result<T, ExpandError>;

/// Shared by scalar coefficients and matrix products. Reservations happen
/// before allocating Cartesian products or expanding large integer powers.
///
/// SymPy's `expand` has no size limit, so the default budget is unbounded
/// (only arithmetic overflow of the counts fails). A term limit and a work
/// budget apply only when requested explicitly (`ExpandBudget::new`, the
/// `max_terms` hint of the Python layer); exceeding them raises an error and
/// never returns a partially expanded result.
pub struct ExpandBudget {
    max_terms: usize,
    remaining_work: usize,
    initial_work: usize,
    scalar_cache: FxHashMap<(usize, ExpandOptions), (Expr, Expr)>,
    matrix_cache: FxHashMap<
        (usize, ExpandOptions),
        (crate::matrixexpr::MatrixExpr, crate::matrixexpr::MatrixExpr),
    >,
}
impl Default for ExpandBudget {
    fn default() -> Self {
        Self::unlimited()
    }
}
impl ExpandBudget {
    /// An explicit resource guard: at most `max_terms` terms per sum and a
    /// work budget of 5 000 000 elementary products.
    pub fn new(max_terms: usize) -> Self {
        Self::with_work(max_terms, 5_000_000)
    }
    /// SymPy's behavior: no term limit and no work budget.
    pub fn unlimited() -> Self {
        Self::with_work(usize::MAX, usize::MAX)
    }
    fn with_work(max_terms: usize, work: usize) -> Self {
        Self {
            max_terms,
            remaining_work: work,
            initial_work: work,
            scalar_cache: FxHashMap::default(),
            matrix_cache: FxHashMap::default(),
        }
    }
    pub(crate) fn matrix_result(
        &self,
        e: &crate::matrixexpr::MatrixExpr,
        options: &ExpandOptions,
    ) -> Option<crate::matrixexpr::MatrixExpr> {
        self.matrix_cache
            .get(&(e.dag_id(), *options))
            .map(|(_, v)| v.clone())
    }
    pub(crate) fn remember_matrix(
        &mut self,
        e: &crate::matrixexpr::MatrixExpr,
        options: &ExpandOptions,
        value: &crate::matrixexpr::MatrixExpr,
    ) {
        self.matrix_cache
            .insert((e.dag_id(), *options), (e.clone(), value.clone()));
    }
    pub fn claim_terms(&mut self, n: usize) -> Result<()> {
        if n > self.max_terms {
            Err(ExpandError(format!(
                "expansion term limit exceeded (maximum {})",
                self.max_terms
            )))
        } else {
            Ok(())
        }
    }
    pub fn claim_product(&mut self, left: usize, right: usize) -> Result<()> {
        let work = left
            .checked_mul(right)
            .ok_or_else(|| ExpandError("expansion work budget exceeded".into()))?;
        if work > self.remaining_work {
            return Err(ExpandError("expansion work budget exceeded".into()));
        }
        self.remaining_work -= work;
        Ok(())
    }
    fn power_terms(&mut self, terms: usize, exponent: u64, support_cap: Option<usize>) -> Result<usize> {
        if terms <= 1 || exponent == 0 {
            return Ok(1);
        }
        // C(exponent+terms-1,terms-1), stopping at the configured bound.
        // Loop over the smaller argument, never over a huge exponent.
        let n = exponent
            .checked_add((terms - 1) as u64)
            .ok_or_else(|| ExpandError("expansion term limit exceeded".into()))?;
        let k = exponent.min((terms - 1) as u64);
        let mut count = 1u128;
        for i in 1..=k {
            count = count
                .checked_mul((n - k + i) as u128)
                .ok_or_else(|| ExpandError("expansion term limit exceeded".into()))?
                / (i as u128);
            if support_cap.is_some_and(|cap| count >= cap as u128) {
                count = support_cap.unwrap() as u128;
                break;
            }
            if count > self.max_terms as u128 {
                return Err(ExpandError(format!(
                    "expansion term limit exceeded (maximum {})",
                    self.max_terms
                )));
            }
        }
        let count = count as usize;
        let exponent = usize::try_from(exponent)
            .map_err(|_| ExpandError("expansion work budget exceeded".into()))?;
        self.claim_product(
            exponent,
            terms
                .checked_mul(count)
                .ok_or_else(|| ExpandError("expansion work budget exceeded".into()))?,
        )?;
        Ok(count)
    }
}

/// A polynomial in one physical symbol has at most degree+1 monomials,
/// regardless of how many distributive products combine into each coefficient.
/// Reject other grammars conservatively; they retain the existing estimator.
fn univariate_support_cap(e: &Expr) -> Option<usize> {
    fn merge(a: &mut Option<Expr>, b: Option<Expr>) -> Option<()> {
        if let Some(b) = b {
            if let Some(a) = a { if !a.same_presentation(&b) { return None; } }
            else { *a = Some(b); }
        }
        Some(())
    }
    fn degree(e: &Expr) -> Option<(Option<Expr>,usize)> {
        match e.kind() {
            Kind::Num(_) => Some((None,0)),
            Kind::Sym(_) if crate::decompose::commutativity(e) != Some(false) => Some((Some(e.clone()),1)),
            Kind::Add(_,terms) => {
                let mut variable = None; let mut maximum = 0;
                for (term,_) in terms {
                    let (symbol,power) = degree(term)?;
                    merge(&mut variable,symbol)?; maximum = maximum.max(power);
                }
                Some((variable,maximum))
            }
            Kind::Mul(_,factors) => {
                let mut variable = None; let mut sum = 0usize;
                for (base,exponent) in factors {
                    if !exponent.is_integer_number() { return None; }
                    let exponent = usize::try_from(exponent.as_number()?.to_i64()?).ok()?;
                    let (symbol,power) = degree(base)?;
                    merge(&mut variable,symbol)?; sum = sum.checked_add(power.checked_mul(exponent)?)?;
                }
                Some((variable,sum))
            }
            Kind::Pow(base,exponent) => {
                if !exponent.is_integer_number() { return None; }
                let exponent = usize::try_from(exponent.as_number()?.to_i64()?).ok()?;
                let (symbol,power) = degree(base)?;
                Some((symbol,power.checked_mul(exponent)?))
            }
            _ => None,
        }
    }
    degree(e)?.1.checked_add(1)
}

/// Check a conservative sparse-polynomial bound before the arithmetic pass.
/// Physical identities preserve the presentation of alpha-equal binders.
fn reserve_arithmetic(
    e: &Expr,
    options: &ExpandOptions,
    budget: &mut ExpandBudget,
    memo: &mut FxHashMap<usize, (Expr, usize)>,
) -> Result<usize> {
    if let Some((_, count)) = memo.get(&e.dag_id()) {
        return Ok(*count);
    }
    let mut child = |e: &Expr| {
        if options.deep {
            reserve_arithmetic(e, options, budget, memo)
        } else {
            Ok(if e.is_add() { e.args().len() } else { 1 })
        }
    };
    let count = match e.kind() {
        Kind::Add(c, terms) => {
            let mut count = usize::from(!c.is_zero());
            for (term, _) in terms {
                count = count
                    .checked_add(child(term)?)
                    .ok_or_else(|| ExpandError("expansion term limit exceeded".into()))?;
            }
            count.max(1)
        }
        Kind::Mul(..) => {
            let support_cap = univariate_support_cap(e);
            let counts = e
                .args()
                .iter()
                .map(&mut child)
                .collect::<Result<Vec<_>>>()?;
            if options.mul {
                let mut count = 1usize;
                for next in counts {
                    budget.claim_product(count, next)?;
                    count = count.checked_mul(next).unwrap_or(usize::MAX)
                        .min(support_cap.unwrap_or(usize::MAX));
                    budget.claim_terms(count)?;
                }
                count
            } else {
                1
            }
        }
        Kind::Pow(base, exponent) => {
            let terms = child(base)?;
            child(exponent)?;
            if options.multinomial && terms > 1 {
                if let Some(n) = exponent.as_number() {
                    let integer = n.floor().to_i64();
                    if let Some(k) = integer {
                        if k != 0 {
                            let count = budget.power_terms(terms, k.unsigned_abs(), univariate_support_cap(e))?;
                            if k > 0 {
                                count
                            } else {
                                1
                            }
                        } else {
                            1
                        }
                    } else if n.is_integer() {
                        return Err(ExpandError("expansion exponent limit exceeded".into()));
                    } else {
                        1
                    }
                } else {
                    1
                }
            } else {
                1
            }
        }
        _ => {
            if options.deep {
                for argument in e.args() {
                    child(&argument)?;
                }
            }
            1
        }
    };
    budget.claim_terms(count)?;
    memo.insert(e.dag_id(), (e.clone(), count));
    Ok(count)
}

#[derive(Clone, Copy)]
enum Hint {
    Trig,
    Log,
    PowerBase,
    PowerExp,
}
fn is_matrix_metadata(e: &Expr) -> bool {
    matches!(
        e.kind(),
        Kind::Func(
            crate::Func::MatrixProxy(_)
                | crate::Func::MatrixElement(_)
                | crate::Func::MatrixTrace(_)
                | crate::Func::MatrixDeterminant(_)
                | crate::Func::MatrixDotProduct(_, _),
            _
        )
    )
}
fn matrix_metadata(e: &Expr, options: &ExpandOptions, budget: &mut ExpandBudget) -> Result<Expr> {
    if !options.deep {
        return Ok(e.clone());
    }
    let key = (e.dag_id(), *options);
    if let Some((_, value)) = budget.scalar_cache.get(&key) {
        return Ok(value.clone());
    }
    let Kind::Func(function, args) = e.kind() else {
        return Ok(e.clone());
    };
    if let crate::Func::MatrixDotProduct(a,b)=function {
        let na=a.expand_with_options(options,budget).map_err(|e|ExpandError(e.to_string()))?;
        let nb=b.expand_with_options(options,budget).map_err(|e|ExpandError(e.to_string()))?;
        let value=if na.dag_id()==a.dag_id()&&nb.dag_id()==b.dag_id(){e.clone()}else{
            crate::matrixcontraction::dot_product(&na,&nb,false).map_err(|e|ExpandError(e.to_string()))?
        };
        budget.scalar_cache.insert(key,(e.clone(),value.clone()));return Ok(value);
    }
    let (crate::Func::MatrixProxy(parent)
    | crate::Func::MatrixElement(parent)
    | crate::Func::MatrixTrace(parent)
    | crate::Func::MatrixDeterminant(parent)) = function
    else {
        return Ok(e.clone());
    };
    let mapped = parent
        .expand_with_options(options, budget)
        .map_err(|e| ExpandError(e.to_string()))?;
    let arguments = args
        .iter()
        .map(|e| expand_with_options(e, options, budget))
        .collect::<Result<Vec<_>>>()?;
    let value = if mapped.dag_id() == parent.dag_id()
        && arguments
            .iter()
            .zip(args)
            .all(|(a, b)| a.same_presentation(b))
    {
        e.clone()
    } else {
        let parent = std::sync::Arc::new(mapped);
        let function = match function {
            crate::Func::MatrixProxy(_) => crate::Func::MatrixProxy(parent),
            crate::Func::MatrixElement(_) => crate::Func::MatrixElement(parent),
            crate::Func::MatrixTrace(_) => crate::Func::MatrixTrace(parent),
            _ => crate::Func::MatrixDeterminant(parent),
        };
        Expr::from_kind(Kind::Func(function, arguments))
    };
    budget.scalar_cache.insert(key, (e.clone(), value.clone()));
    Ok(value)
}
fn hint_pass(
    e: &Expr,
    options: &ExpandOptions,
    hint: Hint,
    budget: &mut ExpandBudget,
    memo: &mut FxHashMap<usize, (Expr, Expr)>,
) -> Result<Expr> {
    if e.is_atom() {
        return Ok(e.clone());
    }
    if let Some((_, value)) = memo.get(&e.dag_id()) {
        return Ok(value.clone());
    }
    if is_matrix_metadata(e) {
        return matrix_metadata(e, options, budget);
    }
    let mut error = None;
    let mut visit = |child: &Expr| match hint_pass(child, options, hint, budget, memo) {
        Ok(value) => value,
        Err(e) => {
            error = Some(e);
            child.clone()
        }
    };
    let mapped = if !options.deep {
        e.clone()
    } else if matches!(hint, Hint::PowerExp) {
        match e.kind() {
            Kind::Add(c, terms) => {
                let transformed: Vec<_> = terms
                    .iter()
                    .map(|(term, k)| (visit(term), k.clone()))
                    .collect();
                if transformed
                    .iter()
                    .zip(terms)
                    .all(|((a, _), (b, _))| a.same_presentation(b))
                {
                    e.clone()
                } else {
                    Expr::from_unevaluated_kind(Kind::Add(c.clone(), transformed))
                }
            }
            Kind::Mul(..) => {
                let args = e.args();
                let transformed: Vec<_> = args.iter().map(&mut visit).collect();
                if args
                    .iter()
                    .zip(&transformed)
                    .all(|(a, b)| a.same_presentation(b))
                {
                    e.clone()
                } else {
                    crate::simp::expanded_product(&transformed)
                }
            }
            _ => crate::traverse::map_children(e, &mut visit),
        }
    } else if matches!(hint, Hint::PowerBase) && e.is_mul() {
        // The reference visits each factor of a product as a power node
        // (`(a*x)**n` in `a*(a*x)**n`), not its base and exponent apart.
        let args = e.args();
        let transformed: Vec<_> = args.iter().map(&mut visit).collect();
        if args.iter().zip(&transformed).all(|(a, b)| a.same_presentation(b)) {
            e.clone()
        } else {
            crate::ops::mul(&transformed)
        }
    } else {
        crate::traverse::map_children(e, &mut visit)
    };
    if let Some(error) = error {
        return Err(error);
    }
    budget.claim_product(1, mapped.args().len())?;
    if matches!(hint, Hint::Trig) {
        if let Kind::Func(
            crate::Func::Sin
            | crate::Func::Cos
            | crate::Func::Tan
            | crate::Func::Cot
            | crate::Func::Sinh
            | crate::Func::Cosh
            | crate::Func::Tanh
            | crate::Func::Coth,
            args,
        ) = mapped.kind()
        {
            if args[0].is_add() {
                let terms = 1usize
                    .checked_shl(args[0].args().len().saturating_sub(1) as u32)
                    .ok_or_else(|| {
                        ExpandError("trigonometric expansion term limit exceeded".into())
                    })?;
                budget.claim_terms(terms)?;
                budget.claim_product(terms, args[0].args().len())?;
            }
        }
    }
    let result = match hint {
        Hint::Trig => crate::simp::expand_trig_deep(&mapped, false),
        Hint::Log => crate::simp::expand_log_deep(&mapped, options.force, false),
        Hint::PowerBase => crate::simp::expand_power_base_deep(&mapped, options.force, false),
        Hint::PowerExp => crate::simp::expand_power_exp_deep(&mapped, false),
    };
    memo.insert(e.dag_id(), (e.clone(), result.clone()));
    Ok(result)
}

/// Whether a logarithm occurs, each shared node visited once.
fn has_log_node(e: &Expr, seen: &mut rustc_hash::FxHashSet<usize>) -> bool {
    if !seen.insert(e.dag_id()) {
        return false;
    }
    e.is_func(&crate::Func::Log) || e.args().iter().any(|a| has_log_node(a, seen))
}

pub fn expand_with_options(
    e: &Expr,
    options: &ExpandOptions,
    budget: &mut ExpandBudget,
) -> Result<Expr> {
    let key = (e.dag_id(), *options);
    if let Some((_, value)) = budget.scalar_cache.get(&key) {
        return Ok(value.clone());
    }
    // A fresh top-level budget has a reproducible work cost. Nested calls can
    // share matrix/scalar reservations, so they use only the per-call cache.
    let fresh = budget.scalar_cache.is_empty() && budget.matrix_cache.is_empty()
        && budget.remaining_work == budget.initial_work;
    if fresh {
        if let Some((value,work)) = super::memo::get(e,options,(budget.max_terms,budget.initial_work)) {
            budget.claim_product(1,work)?;
            budget.scalar_cache.insert(key,(e.clone(),value.clone()));
            return Ok(value);
        }
    }
    let initial_work = budget.remaining_work;
    let mut value = e.clone();
    if options.func {
        value = crate::special::expand_func_deep(&value, options.deep);
    }
    if options.log {
        value = hint_pass(&value, options, Hint::Log, budget, &mut FxHashMap::default())?;
    }
    reserve_arithmetic(&value, options, budget, &mut FxHashMap::default())?;
    let arithmetic = |value: &Expr, budget: &mut ExpandBudget| -> Result<Expr> {
        let mut error = None;
        let value = super::arithmetic_expand_with(value, options, &mut |e| match matrix_metadata(
            e, options, budget,
        ) {
            Ok(value) => Some(value),
            Err(failure) => {
                error = Some(failure);
                Some(e.clone())
            }
        });
        error.map_or(Ok(value), Err)
    };
    value = arithmetic(&value, budget)?;
    // SymPy applies the hints in sorted order (multinomial and mul before
    // power_base), then repeats multinomial and mul until nothing changes.
    if options.power_base {
        let next = hint_pass(&value, options, Hint::PowerBase, budget, &mut FxHashMap::default())?;
        if next != value {
            reserve_arithmetic(&next, options, budget, &mut FxHashMap::default())?;
            value = arithmetic(&next, budget)?;
        }
    }
    if options.trig {
        value = hint_pass(
            &value,
            options,
            Hint::Trig,
            budget,
            &mut FxHashMap::default(),
        )?;
        reserve_arithmetic(&value, options, budget, &mut FxHashMap::default())?;
        value = arithmetic(&value, budget)?;
    }
    // the reference repeats the mul, multinomial and log hints until
    // nothing changes: a logarithm produced by the first pass is expanded too
    if options.log {
        for _ in 0..16 {
            if !has_log_node(&value, &mut rustc_hash::FxHashSet::default()) {
                break;
            }
            let next = hint_pass(&value, options, Hint::Log, budget, &mut FxHashMap::default())?;
            if next == value {
                break;
            }
            reserve_arithmetic(&next, options, budget, &mut FxHashMap::default())?;
            value = arithmetic(&next, budget)?;
        }
    }
    if options.power_exp {
        value = hint_pass(
            &value,
            options,
            Hint::PowerExp,
            budget,
            &mut FxHashMap::default(),
        )?;
    }
    budget.scalar_cache.insert(key, (e.clone(), value.clone()));
    if fresh {
        super::memo::insert(e,&value,options,(budget.max_terms,budget.initial_work),initial_work-budget.remaining_work);
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn reused_expansions_replay_work_and_keep_distinct_limits() {
        let expression = parse_str("(x+1)*(x+2)*(x+3)*(x+4)*(x+5)").unwrap();
        let options = ExpandOptions::default();
        let mut first = ExpandBudget::new(6);
        let expected = expand_with_options(&expression,&options,&mut first).unwrap();
        let mut repeated = ExpandBudget::new(6);
        assert_eq!(expand_with_options(&expression,&options,&mut repeated).unwrap(),expected);
        assert_eq!(first.remaining_work,repeated.remaining_work);
        assert!(first.remaining_work < 5_000_000);
        assert!(expand_with_options(&expression,&options,&mut ExpandBudget::new(5)).is_err());
        let mut shallow = options; shallow.mul = false;
        assert_eq!(expand_with_options(&expression,&shallow,&mut ExpandBudget::default()).unwrap(),expression);
    }
    #[test]
    fn functional_hint_respects_depth_and_arithmetic_controls() {
        for (source, deep, wanted) in [
            ("gamma(x+2)", true, "x*(x+1)*gamma(x)"),
            ("gamma(x+2)", false, "x*(x+1)*gamma(x)"),
            ("sin(gamma(x+2))", true, "sin(x*(x+1)*gamma(x))"),
            ("sin(gamma(x+2))", false, "sin(gamma(x+2))"),
        ] {
            let options = ExpandOptions { func: true, deep, mul: false, ..Default::default() };
            let value = expand_with_options(&parse_str(source).unwrap(), &options,
                &mut ExpandBudget::default()).unwrap();
            assert_eq!(value, parse_str(wanted).unwrap());
        }
    }
    #[test]
    fn shallow_and_independent_arithmetic_hints() {
        for (source, options, wanted) in [
            (
                "sin((x+y)^2)",
                ExpandOptions {
                    deep: false,
                    ..Default::default()
                },
                "sin((x+y)^2)",
            ),
            (
                "(x+y)^2",
                ExpandOptions {
                    deep: false,
                    ..Default::default()
                },
                "x^2+2*x*y+y^2",
            ),
            (
                "z*(x+y)^2",
                ExpandOptions {
                    mul: false,
                    ..Default::default()
                },
                "z*(x^2+2*x*y+y^2)",
            ),
            (
                "z*(x+y)^2",
                ExpandOptions {
                    multinomial: false,
                    ..Default::default()
                },
                "z*(x+y)^2",
            ),
        ] {
            assert_eq!(
                expand_with_options(
                    &parse_str(source).unwrap(),
                    &options,
                    &mut ExpandBudget::default()
                )
                .unwrap(),
                parse_str(wanted).unwrap()
            );
        }
    }
    #[test]
    fn huge_powers_fail_before_allocation() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        for exponent in [1_000_000_000, i64::MIN, i64::MAX] {
            let value = crate::ops::pow(&crate::ops::add2(&x, &y), &Expr::int(exponent));
            let before = value.hash64();
            assert!(expand_with_options(
                &value,
                &ExpandOptions::default(),
                &mut ExpandBudget::new(1000)
            )
            .is_err());
            assert_eq!(value.hash64(), before);
        }
    }
    #[test]
    fn sparse_exponents_do_not_overflow_machine_i32() {
        let value = parse_str("(x^2147483647+y)^2").unwrap();
        let wanted = parse_str("x^4294967294+2*x^2147483647*y+y^2").unwrap();
        assert_eq!(
            expand_with_options(
                &value,
                &ExpandOptions::default(),
                &mut ExpandBudget::default()
            )
            .unwrap(),
            wanted
        );
    }
}
