//! Exact held derivative orders, separate from the existing u32 engine.
//!
//! Arguments have the public Derivative shape: body, Tuple(variable, count)...
//! A count is never truncated or inferred integral from a numerical value.
//! All-machine calls retain the original path; symbolic/noninteger orders stay
//! held when the native engine has no exact n-th-order evaluation rule.
use crate::expr::{Expr, Func, Kind};
use crate::ops::add2;

#[derive(Clone, Debug)]
pub enum Error {
    Value(String),
    Type(String),
    Unsupported(String),
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Value(s) | Self::Type(s) | Self::Unsupported(s) => f.write_str(s),
        }
    }
}
thread_local! {static ERRORS:std::cell::RefCell<Vec<Option<Error>>>=const{std::cell::RefCell::new(Vec::new())};}
struct Scope(bool);
impl Drop for Scope {
    fn drop(&mut self) {
        if self.0 {
            ERRORS.with(|e| {
                e.borrow_mut().pop();
            });
        }
    }
}
/// A bounded call-local collector. Nested API calls consume their own errors;
/// catching one never poisons its caller, and unwinding always pops the frame.
pub fn capture<R>(operation: impl FnOnce() -> R) -> (R, Option<Error>) {
    ERRORS.with(|e| e.borrow_mut().push(None));
    let mut scope = Scope(true);
    let result = operation();
    let error = ERRORS.with(|e| e.borrow_mut().pop().unwrap());
    scope.0 = false;
    (result, error)
}
pub fn record(error: Error) {
    ERRORS.with(|e| {
        if let Some(slot) = e.borrow_mut().last_mut() {
            if slot.is_none() {
                *slot = Some(error);
            }
        }
    });
}
/// NaN is an internal cancellation barrier only inside a live error frame.
/// Outside a frame the caller must retain its complete unevaluated operation.
pub(crate) fn failed(error: Error, held: impl FnOnce() -> Expr) -> Expr {
    let active = ERRORS.with(|errors| !errors.borrow().is_empty());
    if active {
        record(error);
        Expr::nan()
    } else {
        held()
    }
}
pub fn failure_pending() -> bool {
    ERRORS.with(|e| e.borrow().iter().any(Option::is_some))
}
pub fn checked<R>(operation: impl FnOnce() -> R) -> Result<R, Error> {
    let (value, error) = capture(operation);
    error.map_or(Ok(value), Err)
}

pub fn is_derivative(e: &Expr) -> bool {
    matches!(
        e.kind(),
        Kind::Derivative(..) | Kind::Func(Func::SymbolicDerivative, _)
    )
}
pub fn parts(e: &Expr) -> Option<(&Expr, Vec<(Expr, Expr)>)> {
    let Kind::Func(Func::SymbolicDerivative, args) = e.kind() else {
        return None;
    };
    if args.len() < 2 {
        return None;
    }
    let body = args.first()?;
    let vars = args[1..]
        .iter()
        .map(|pair| {
            let pair = crate::special::tuple_items(pair)?;
            (pair.len() == 2).then(|| (pair[0].clone(), pair[1].clone()))
        })
        .collect::<Option<Vec<_>>>()?;
    Some((body, vars))
}
pub fn data(e: &Expr) -> Option<(&Expr, Vec<(Expr, Expr)>)> {
    if let Kind::Derivative(body, vars) = e.kind() {
        Some((
            body,
            vars.iter()
                .map(|(v, n)| (v.clone(), Expr::int(*n as i64)))
                .collect(),
        ))
    } else {
        parts(e)
    }
}
pub fn literal(e: &Expr) -> Option<u32> {
    let n = e.as_number()?;
    if n.is_float() || !n.is_integer() {
        return None;
    }
    e.as_i64().and_then(|n| u32::try_from(n).ok())
}
pub fn zero(e: &Expr) -> bool {
    literal(e) == Some(0)
}
pub fn one(e: &Expr) -> bool {
    literal(e) == Some(1)
}
pub(crate) fn raw(body: &Expr, vars: &[(Expr, Expr)]) -> Expr {
    let mut args = vec![body.clone()];
    args.extend(
        vars.iter()
            .map(|(v, n)| crate::special::tuple(vec![v.clone(), n.clone()])),
    );
    Expr::from_kind(Kind::Func(Func::SymbolicDerivative, args))
}
pub fn validate_orders(vars: &[(Expr, Expr)], evaluate: bool) -> Result<(), Error> {
    crate::assumptions::exact_queries(|| {
        for (v, n) in vars {
            if zero(n) {
                continue;
            }
            if crate::assumptions::is_negative(n) == Some(true) {
                return Err(Error::Value(
                    "order of differentiation must be nonnegative".into(),
                ));
            }
            if evaluate && matches!(n.kind(), Kind::Bool(_)) {
                return Err(Error::Type(
                    "BooleanAtom not allowed in this context.".into(),
                ));
            }
            let valid = if parts(v).is_some() {
                valid_variable(v)?
            } else {
                super::is_valid_variable(v)
            };
            if !valid {
                return Err(Error::Value(format!(
                    "\nCan't calculate derivative wrt {}.",
                    crate::printer::to_str(v)
                )));
            }
            if evaluate {
                if n.is_constant(crate::expr::Constant::NaN) {
                    return Err(Error::Type("Invalid NaN comparison".into()));
                }
                if crate::assumptions::fact(n, "extended_real") == Some(false) {
                    return Err(Error::Type(format!(
                        "Invalid comparison of non-real {}",
                        crate::printer::to_str(n)
                    )));
                }
            }
        }
        Ok(())
    })
}
pub fn valid_variable(e: &Expr) -> Result<bool, Error> {
    let Some((body, vars)) = parts(e) else {
        return Ok(false);
    };
    if matches!(body.kind(),Kind::Func(Func::Undefined(head),_) if head.host().is_some()) {
        return Err(Error::Unsupported(
            "Host _diff_wrt protocols for symbolic-order derivative variables are not implemented"
                .into(),
        ));
    }
    if !super::applied_undef(body) {
        return Ok(false);
    }
    // Reference Derivative._diff_wrt evaluates the held derivative and checks
    // that it remains a Derivative. Errors from this check are observable.
    Ok(is_derivative(&construct(body, vars, true)?))
}
pub fn validate_tree(e: &Expr, evaluate: bool) -> Result<(), Error> {
    for (node, _) in crate::inspection::scalar_nodes(e) {
        if matches!(node.kind(), Kind::Func(Func::SymbolicDerivative, _)) {
            let (_, vars) = parts(&node)
                .ok_or_else(|| Error::Value("invalid symbolic Derivative arguments".into()))?;
            validate_orders(&vars, evaluate)?;
        }
    }
    Ok(())
}
fn count_add(a: &Expr, b: &Expr) -> Expr {
    match (a.as_number(), b.as_number()) {
        (Some(a), Some(b)) => Expr::num(a.add(b)),
        _ => add2(a, b),
    }
}
fn merge(vars: Vec<(Expr, Expr)>) -> Result<Vec<(Expr, Expr)>, Error> {
    let mut merged: Vec<(Expr, Expr)> = Vec::new();
    for (v, n) in vars {
        if zero(&n) {
            continue;
        }
        if let Some((w, m)) = merged.last_mut() {
            if *w == v {
                if matches!(m.kind(), Kind::Bool(_)) || matches!(n.kind(), Kind::Bool(_)) {
                    return Err(Error::Type(
                        "BooleanAtom not allowed in this context.".into(),
                    ));
                }
                *m = count_add(m, &n);
                if zero(m) {
                    merged.pop();
                }
                continue;
            }
        }
        merged.push((v, n));
    }
    Ok(merged)
}
fn sorted(vars: Vec<(Expr, Expr)>) -> Result<Vec<(Expr, Expr)>, Error> {
    let order = super::variable_order(&vars);
    merge(order.into_iter().map(|i| vars[i].clone()).collect())
}
pub fn construct(body: &Expr, vars: Vec<(Expr, Expr)>, evaluate: bool) -> Result<Expr, Error> {
    validate_orders(&vars, false)?;
    if (parts(body).is_some() || vars.iter().any(|(_, n)| literal(n).is_none()))
        && vars.iter().any(|(v, n)| {
            !zero(n)
                && v.as_sym().is_some_and(|symbol| {
                    v.host_handle().is_some() || symbol.scalar_identity().is_some()
                })
        })
    {
        return Err(Error::Unsupported("Symbolic derivative orders with host or specialized Symbol variables are not implemented".into()));
    }
    if evaluate {
        // Count metadata is free for substitution, but formal derivatives do
        // not differentiate it. A fresh Symbol absent from the actual body
        // therefore annihilates the derivative, even if it occurs as a count.
        if let Some((inner, held)) = parts(body) {
            if vars.iter().all(|(v, _)| v.is_symbol()) {
                let free = crate::traverse::free_symbols(inner);
                if vars.iter().any(|(v, n)| {
                    literal(n).is_some_and(|n| n > 0)
                        && !free.contains(v)
                        && !held.iter().any(|(w, _)| w == v)
                }) {
                    return Ok(Expr::zero());
                }
            }
        }
        let free = crate::traverse::free_symbols(body);
        if crate::assumptions::exact_queries(|| {
            vars.iter().any(|(v, n)| {
                v.is_symbol()
                    && !free.contains(v)
                    && crate::assumptions::is_positive(n) == Some(true)
            })
        }) {
            return Ok(Expr::zero());
        }
        validate_orders(&vars, true)?;
    }
    // The original literal-order route, including its overflow-preserving
    // adjacent groups, is deliberately unchanged.
    if parts(body).is_none() {
        if let Some(machine) = vars
            .iter()
            .map(|(v, n)| Some((v.clone(), literal(n)?)))
            .collect::<Option<Vec<_>>>()
        {
            return Ok(if evaluate {
                super::diff_evaluated(body, &machine)
            } else {
                super::derivative_unevaluated(body, machine)
            });
        }
    }
    let (base, mut all) = data(body)
        .map(|(e, v)| (e.clone(), v))
        .unwrap_or_else(|| (body.clone(), Vec::new()));
    all.extend(vars);
    if all.windows(2).any(|pair| {
        pair[0].0 == pair[1].0
            && (matches!(pair[0].1.kind(), Kind::Bool(_))
                || matches!(pair[1].1.kind(), Kind::Bool(_)))
    }) {
        return Err(Error::Type(
            "BooleanAtom not allowed in this context.".into(),
        ));
    }
    let mut all = merge(all)?;
    if evaluate {
        all = sorted(all)?;
        validate_orders(&all, true)?;
    }
    if all.is_empty() {
        return Ok(base);
    }
    if evaluate{if let Some(value)=super::integer_orders::groups(&base,&all){return Ok(value);}}
    if let Some(machine) = all
        .iter()
        .map(|(v, n)| Some((v.clone(), literal(n)?)))
        .collect::<Option<Vec<_>>>()
    {
        return Ok(if evaluate {
            super::diff_evaluated(&base, &machine)
        } else {
            super::derivative_unevaluated(&base, machine)
        });
    }
    if evaluate {
        if base.is_zero() {
            return Ok(base);
        }
        let free = crate::traverse::free_symbols(&base);
        if crate::assumptions::exact_queries(|| {
            all.iter().any(|(v, n)| {
                crate::assumptions::is_positive(n) == Some(true)
                    && v.is_symbol()
                    && !free.contains(v)
            })
        }) {
            return Ok(Expr::zero());
        }
        let prefix = all.iter().take_while(|(_, n)| literal(n).is_some()).count();
        if prefix > 0 {
            let machine = all[..prefix]
                .iter()
                .map(|(v, n)| (v.clone(), literal(n).unwrap()))
                .collect::<Vec<_>>();
            let evaluated = super::diff_evaluated(&base, &machine);
            return construct(&evaluated, all[prefix..].to_vec(), false);
        }
    }
    Ok(raw(&base, &all))
}
/// Structural construction for trusted native callers. Evaluated construction
/// and changed counters must use construct -> Result. Generic transformations
/// hold this head rather than rebuilding its metadata.
pub(crate) fn rebuild(args: Vec<Expr>) -> Expr {
    Expr::from_kind(Kind::Func(Func::SymbolicDerivative, args))
}

pub fn canonical(e: &Expr) -> Expr {
    match data(e) {
        Some((body, vars)) => sorted(vars)
            .and_then(|vars| construct(body, vars, false))
            .unwrap_or_else(|error| {
                record(error);
                e.clone()
            }),
        None => e.clone(),
    }
}
pub fn doit(e: &Expr, deep: bool) -> Expr {
    let Some((body, vars)) = parts(e) else {
        return e.clone();
    };
    let body = if deep {
        crate::functions::doit(body, true)
    } else {
        body.clone()
    };
    construct(&body, vars, true).unwrap_or_else(|error| {
        record(error);
        e.clone()
    })
}
pub fn differentiate(e: &Expr, x: &Expr) -> Expr {
    let Some((body, mut vars)) = parts(e) else {
        return e.clone();
    };
    if x.is_symbol() && !crate::traverse::free_symbols(body).contains(x) {
        return Expr::zero();
    }
    if vars.iter().any(|(v, _)| v == x) {
        vars.push((x.clone(), Expr::one()));
        return construct(body, vars, false)
            .and_then(|result| {
                if let Some((_, vars)) = parts(&result) {
                    validate_orders(&vars, true)?;
                }
                Ok(result)
            })
            .unwrap_or_else(|error| {
                record(error);
                super::derivative_unevaluated(e, vec![(x.clone(), 1)])
            });
    }
    let differentiated = super::diff(body, x);
    if differentiated.is_zero() {
        return differentiated;
    }
    if let Some((inner, more)) = data(&differentiated) {
        vars.extend(more);
        return construct(inner, vars, false)
            .and_then(|result| {
                if let Some((_, vars)) = parts(&result) {
                    validate_orders(&vars, true)?;
                }
                Ok(result)
            })
            .unwrap_or_else(|error| {
                record(error);
                super::derivative_unevaluated(e, vec![(x.clone(), 1)])
            });
    }
    construct(&differentiated, vars, true).unwrap_or_else(|error| {
        record(error);
        super::derivative_unevaluated(e, vec![(x.clone(), 1)])
    })
}
pub fn total(vars: &[(Expr, Expr)]) -> Result<Expr, Error> {
    if vars.iter().any(|(_, n)| matches!(n.kind(), Kind::Bool(_))) {
        return Err(Error::Type(
            "BooleanAtom not allowed in this context.".into(),
        ));
    }
    Ok(vars
        .iter()
        .map(|(_, n)| n.clone())
        .reduce(|a, b| count_add(&a, &b))
        .unwrap_or_else(Expr::zero))
}

/// Derivative._eval_subs, with counts substituted before differentiation
/// variables. Scalar replacements are Subs operations, not a derivative wrt1.
pub(crate) fn substitute(
    e: &Expr,
    old: &Expr,
    new: &Expr,
    transform: &mut dyn FnMut(&Expr) -> Expr,
) -> Option<Expr> {
    let (body, vars) = parts(e)?;
    let changed = vars
        .iter()
        .map(|(v, n)| (v.clone(), transform(n)))
        .collect::<Vec<_>>();
    let base = match construct(body, changed, false) {
        Ok(e) => e,
        Err(error) => {
            return Some(failed(error, || {
                crate::subs::held_subs(e, vec![old.clone()], vec![new.clone()])
            }));
        }
    };
    let Some((body, vars)) = data(&base) else {
        return Some(transform(&base));
    };
    if let Some((old_body, old_vars)) = data(old) {
        if body == old_body && canonical(&base) == canonical(old) {
            return Some(new.clone());
        }
        // Compare a sub-multiset only when its ordering is unambiguous.
        if body == old_body
            && vars
                .iter()
                .enumerate()
                .all(|(i, (v, _))| !vars[..i].iter().any(|(w, _)| w == v))
            && old_vars
                .iter()
                .enumerate()
                .all(|(i, (v, _))| !old_vars[..i].iter().any(|(w, _)| w == v))
        {
            let count = |list: &[(Expr, Expr)], v: &Expr| {
                list.iter()
                    .find(|(w, _)| w == v)
                    .map(|(_, n)| n.clone())
                    .unwrap_or_else(Expr::zero)
            };
            let subset = crate::assumptions::exact_queries(|| {
                old_vars
                    .iter()
                    .all(|(v, n)| crate::relational::le(n, &count(&vars, v)) == Expr::boolean(true))
            });
            if subset {
                let rest = vars
                    .iter()
                    .map(|(v, n)| (v.clone(), crate::ops::sub(n, &count(&old_vars, v))))
                    .collect();
                return Some(construct(new, rest, false).unwrap_or_else(|error| {
                    failed(error, || {
                        crate::subs::held_subs(e, vec![old.clone()], vec![new.clone()])
                    })
                }));
            }
        }
    }
    if vars.iter().any(|(v, _)| v == old) && !super::is_valid_variable(new) {
        if old.is_symbol() {
            return Some(crate::subs::held_subs(
                &base,
                vec![old.clone()],
                vec![new.clone()],
            ));
        }
        let dummy = Expr::dummy("xi", crate::Assumptions::NONE);
        let masked =
            crate::subs::xreplace(&base, &[(old.clone(), dummy.clone())].into_iter().collect());
        return Some(crate::subs::held_subs(
            &masked,
            vec![dummy],
            vec![new.clone()],
        ));
    }
    if crate::traverse::has(&base, old)
        && vars.iter().any(|(v, _)| {
            crate::traverse::has(new, v)
                && !crate::traverse::has(old, v)
                && crate::traverse::has(body, v)
        })
    {
        return Some(crate::subs::held_subs(
            &base,
            vec![old.clone()],
            vec![new.clone()],
        ));
    }
    let mapped = vars
        .iter()
        .map(|(v, n)| (transform(v), n.clone()))
        .collect::<Vec<_>>();
    if mapped
        .iter()
        .zip(&vars)
        .any(|((v, _), (old, _))| v != old && crate::traverse::has(body, v))
    {
        return Some(crate::subs::held_subs(
            &base,
            vec![old.clone()],
            vec![new.clone()],
        ));
    }
    let mut masked = body.clone();
    let mut variables = Vec::new();
    let mut points = Vec::new();
    let mut final_vars = Vec::new();
    for (i, ((variable, n), (before, _))) in mapped.into_iter().zip(&vars).enumerate() {
        if super::is_valid_variable(&variable) {
            final_vars.push((variable, n));
        } else {
            let dummy = Expr::dummy(&format!("xi_{i}"), crate::Assumptions::NONE);
            masked = crate::subs::xreplace(
                &masked,
                &[(before.clone(), dummy.clone())].into_iter().collect(),
            );
            final_vars.push((dummy.clone(), n));
            variables.push(dummy);
            points.push(variable);
        }
    }
    let expression = construct(&transform(&masked), final_vars, false).unwrap_or_else(|error| {
        failed(error, || {
            crate::subs::held_subs(e, vec![old.clone()], vec![new.clone()])
        })
    });
    Some(if variables.is_empty() {
        expression
    } else {
        crate::subs::held_subs(&expression, variables, points)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn literal_machine_path_and_exact_held_counts() {
        let x = Expr::symbol("x");
        let n = Expr::symbol("n");
        let body = crate::functions::sin(&x);
        assert_eq!(
            construct(&body, vec![(x.clone(), Expr::int(2))], false).unwrap(),
            super::super::derivative_unevaluated(&body, vec![(x.clone(), 2)])
        );
        let held = construct(&body, vec![(x.clone(), n.clone())], false).unwrap();
        assert_eq!(parts(&held).unwrap().1, vec![(x.clone(), n.clone())]);
        assert_eq!(crate::traverse::free_symbols(&held).len(), 2);
        let large = Expr::int(4294967296);
        let held = construct(&body, vec![(x.clone(), large.clone())], false).unwrap();
        assert_eq!(parts(&held).unwrap().1, vec![(x, large)]);
        assert_eq!(
            node_weight(&held).cmp_value(&crate::number::Number::Int(4294967298)),
            Some(std::cmp::Ordering::Equal)
        );
    }
    #[test]
    fn errors_are_call_scoped_and_cannot_cancel() {
        let x = Expr::symbol("x");
        let n = Expr::symbol("n");
        let body = crate::functions::sin(&x);
        let a = construct(&body, vec![(x.clone(), n.clone())], false).unwrap();
        let b = construct(
            &body,
            vec![(
                x.clone(),
                crate::ops::sub(&crate::ops::neg(&n), &Expr::int(2)),
            )],
            false,
        )
        .unwrap();
        let expression = crate::ops::sub(&a, &b);
        assert!(checked(|| crate::subs::subs1(&expression, &n, &Expr::minus_one())).is_err());
        assert_eq!(
            checked(|| crate::subs::subs1(&a, &n, &Expr::zero())).unwrap(),
            body
        );
        let (value, error) = capture(|| {
            assert!(checked(|| crate::subs::subs1(&a, &n, &Expr::minus_one())).is_err());
            crate::subs::subs1(&a, &n, &Expr::zero())
        });
        assert!(error.is_none());
        assert_eq!(value, body);
    }
    #[test]
    fn opaque_walkers_preserve_an_uncomputed_invalid_order_without_a_frame() {
        let a = Expr::symbol("a");
        let x = Expr::symbol("x");
        let count = crate::ops::sub(
            &crate::ops::sub(
                &crate::ops::sub(&crate::ops::add2(&a, &Expr::one()).powi(2), &a.powi(2)),
                &crate::ops::mul2(&Expr::int(2), &a),
            ),
            &Expr::int(2),
        );
        assert_eq!(crate::expand::expand(&count), Expr::minus_one());
        let held = construct(&crate::functions::sin(&x), vec![(x, count)], false).unwrap();
        assert_eq!(crate::expand::expand(&held), held);
        assert_eq!(crate::simplify::simplify(&held), held);
        assert_eq!(crate::evalf::evalf(&held), held);
        assert!(!failure_pending());
    }
    #[test]
    fn unguarded_native_substitution_keeps_the_whole_failed_operation() {
        let x = Expr::symbol("x");
        let n = Expr::symbol("n");
        let body = crate::functions::sin(&x);
        let a = construct(&body, vec![(x.clone(), n.clone())], false).unwrap();
        let b = construct(
            &body,
            vec![(x, crate::ops::sub(&crate::ops::neg(&n), &Expr::int(2)))],
            false,
        )
        .unwrap();
        let expression = crate::ops::sub(&a, &b);
        let held = crate::subs::subs1(&expression, &n, &Expr::minus_one());
        assert!(matches!(held.kind(), Kind::Func(Func::Subs, _)));
        assert!(!held.is_zero() && !held.is_constant(crate::expr::Constant::NaN));
        let held = failed(Error::Value("outside".into()), || expression.clone());
        assert_eq!(held, expression);
        assert!(!failure_pending());
    }
    #[test]
    fn scope_pops_on_unwind() {
        let result = std::panic::catch_unwind(|| {
            capture(|| {
                record(Error::Value("test".into()));
                panic!("interrupt")
            })
        });
        assert!(result.is_err());
        assert!(checked(|| Expr::one()).is_ok());
        assert!(!failure_pending());
    }
}

/// The reference sorting rank for a Derivative counts its order, not the
/// (variable,count) metadata tree. Rational/Float counts remain exact native
/// Numbers for comparison; they are never rounded to an integer rank.
pub fn node_weight(e: &Expr) -> crate::number::Number {
    use crate::number::Number;
    if let Some((body, vars)) = data(e) {
        vars.iter().fold(node_weight(body), |sum, (_, count)| {
            let weight = count
                .as_number()
                .cloned()
                .unwrap_or_else(|| match count.kind() {
                    Kind::Const(crate::expr::Constant::Infinity) => Number::float(f64::INFINITY),
                    Kind::Const(crate::expr::Constant::NaN) => Number::float(f64::NAN),
                    _ => Number::from_bigint(crate::traverse::node_count(count).into()),
                });
            sum.add(&weight)
        })
    } else {
        Number::from_bigint(crate::traverse::node_count(e).into())
    }
}

/// Exact closed form for the complete held literal-order request, without
/// applying the evaluated API's final simplify/factor_terms option.
pub fn literal_evaluation(e:&Expr)->Option<Expr>{let(body,vars)=data(e)?;super::integer_orders::groups(body,&vars)}
