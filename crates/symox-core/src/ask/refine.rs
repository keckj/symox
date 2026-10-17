//! Assumption-dependent scalar rewrites, following SymPy's refine handlers.
use super::{
    engine::{is_boolean_kind, Cx},
    AskError, AskOptions, Pred,
};
use crate::{
    expr::{Expr, Func, Kind},
    functions::func,
    ops,
};
use rustc_hash::FxHashMap;

pub fn refine(e: &Expr, assumptions: &Expr) -> Result<Expr, AskError> {
    refine_with(e, assumptions, &AskOptions::default())
}

pub fn refine_with(
    e: &Expr,
    assumptions: &Expr,
    options: &AskOptions<'_>,
) -> Result<Expr, AskError> {
    let cx = Cx::new(&options.context, options.external);
    visit(e, assumptions, &cx, &mut FxHashMap::default())
}

fn visit(
    e: &Expr,
    a: &Expr,
    cx: &Cx<'_>,
    memo: &mut FxHashMap<usize, (Expr, Expr)>,
) -> Result<Expr, AskError> {
    if let Some((_, value)) = memo.get(&e.dag_id()) {
        return Ok(value.clone());
    }
    let mut failure = None;
    let node = crate::traverse::map_children(e, &mut |child| {
        if failure.is_some() {
            return child.clone();
        }
        match visit(child, a, cx, memo) {
            Ok(value) => value,
            Err(error) => {
                failure = Some(error);
                child.clone()
            }
        }
    });
    if let Some(error) = failure {
        return Err(error);
    }
    let result = rewrite(&node, a, cx)?;
    let result = if result.same_presentation(&node) {
        result
    } else {
        visit(&result, a, cx, memo)?
    };
    memo.insert(e.dag_id(), (e.clone(), result.clone()));
    Ok(result)
}

fn rewrite(e: &Expr, a: &Expr, cx: &Cx<'_>) -> Result<Expr, AskError> {
    use Pred::*;
    if is_boolean_kind(e) {
        return Ok(cx
            .ask(e, a)?
            .map(Expr::boolean)
            .unwrap_or_else(|| e.clone()));
    }
    match e.kind() {
        Kind::Func(Func::Abs, args) => {
            let x = &args[0];
            if cx.is(Real, x, a)? && cx.askp(Negative, x, a)? == Some(false) {
                return Ok(x.clone());
            }
            if cx.is(Negative, x, a)? {
                return Ok(ops::neg(x));
            }
            if matches!(x.kind(), Kind::Mul(..)) {
                let mut outside = Vec::new();
                let mut inside = Vec::new();
                for factor in x.args() {
                    let r = rewrite(&func(Func::Abs, vec![factor.clone()]), a, cx)?;
                    match r.kind() {
                        Kind::Func(Func::Abs, args) => inside.push(args[0].clone()),
                        _ => outside.push(r),
                    }
                }
                outside.push(func(Func::Abs, vec![ops::mul(&inside)]));
                return Ok(ops::mul(&outside));
            }
        }
        Kind::Pow(b, exponent) => {
            if let Kind::Func(Func::Abs, args) = b.kind() {
                if cx.is(Real, &args[0], a)? && cx.is(Even, exponent, a)? {
                    return Ok(ops::pow(&args[0], exponent));
                }
            }
            if b.as_number().is_some_and(|n| n.is_negative()) {
                if cx.is(Even, exponent, a)? {
                    return Ok(ops::pow(&ops::neg(b), exponent));
                }
                if cx.is(Odd, exponent, a)? {
                    return Ok(ops::neg(&ops::pow(&ops::neg(b), exponent)));
                }
            }
            if let Kind::Pow(inner, power) = b.kind() {
                if exponent == &Expr::rational(1, 2)
                    && cx.is(Real, inner, a)?
                    && power.as_i64().is_some_and(|n| n >= 0 && n % 2 == 0)
                {
                    let abs = rewrite(&func(Func::Abs, vec![inner.clone()]), a, cx)?;
                    return Ok(ops::pow(&abs, &ops::mul2(power, exponent)));
                }
            }
            if b == &Expr::int(-1) && matches!(exponent.kind(), Kind::Add(..)) {
                let mut terms = Vec::new();
                let mut odd = 0;
                let mut changed = false;
                for term in exponent.args() {
                    if cx.is(Even, &term, a)? {
                        changed = true;
                    } else if cx.is(Odd, &term, a)? {
                        odd += 1;
                        changed = true;
                    } else {
                        terms.push(term);
                    }
                }
                if changed {
                    if odd % 2 == 1 {
                        terms.push(Expr::one());
                    }
                    return Ok(ops::pow(b, &ops::add(&terms)));
                }
            }
        }
        Kind::Func(f @ (Func::Re | Func::Im), args) => {
            let x = &args[0];
            if cx.is(Real, x, a)? {
                return Ok(if *f == Func::Re {
                    x.clone()
                } else {
                    Expr::zero()
                });
            }
            if cx.is(Imaginary, x, a)? {
                return Ok(if *f == Func::Re {
                    Expr::zero()
                } else {
                    ops::neg(&ops::mul2(&Expr::i(), x))
                });
            }
        }
        Kind::Func(Func::Arg, args) => {
            if cx.is(Positive, &args[0], a)? {
                return Ok(Expr::zero());
            }
            if cx.is(Negative, &args[0], a)? {
                return Ok(Expr::pi());
            }
        }
        Kind::Func(Func::Sign, args) => {
            let x = &args[0];
            if cx.is(Zero, x, a)? {
                return Ok(Expr::zero());
            }
            // Like SymPy, the real/imaginary gate uses intrinsic facts.
            let t = Expr::boolean(true);
            if cx.is(Real, x, &t)? {
                if cx.is(Positive, x, a)? {
                    return Ok(Expr::one());
                }
                if cx.is(Negative, x, a)? {
                    return Ok(Expr::int(-1));
                }
            }
            if cx.is(Imaginary, x, &t)? {
                let im = func(Func::Im, vec![x.clone()]);
                if cx.is(Positive, &im, a)? {
                    return Ok(Expr::i());
                }
                if cx.is(Negative, &im, a)? {
                    return Ok(ops::neg(&Expr::i()));
                }
            }
        }
        Kind::Func(Func::Atan2, args) => {
            let (y, x) = (&args[0], &args[1]);
            let angle = || func(Func::Atan, vec![ops::div(y, x)]);
            if cx.is(Real, y, a)? && cx.is(Positive, x, a)? {
                return Ok(angle());
            }
            if cx.is(Negative, x, a)? {
                if cx.is(Negative, y, a)? {
                    return Ok(ops::sub(&angle(), &Expr::pi()));
                }
                if cx.is(Positive, y, a)? {
                    return Ok(ops::add2(&angle(), &Expr::pi()));
                }
                if cx.is(Zero, y, a)? {
                    return Ok(Expr::pi());
                }
            }
            if cx.is(Zero, x, a)? {
                if cx.is(Positive, y, a)? {
                    return Ok(ops::div(&Expr::pi(), &Expr::int(2)));
                }
                if cx.is(Negative, y, a)? {
                    return Ok(ops::neg(&ops::div(&Expr::pi(), &Expr::int(2))));
                }
                if cx.is(Zero, y, a)? {
                    return Ok(Expr::nan());
                }
            }
        }
        _ => {}
    }
    Ok(e.clone())
}
