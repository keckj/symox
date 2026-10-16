//! Generic tree traversal helpers.

use rustc_hash::FxHashSet;

use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::ops::{add, mul, mul2, pow};
use crate::relational::{piecewise, rel};

/// Rebuild `e` with every direct child mapped through `f`, using the
/// canonical constructors. Returns `e` itself when nothing changed.
pub fn map_children(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => e.clone(),
        Kind::Add(c, terms) => {
            let mut changed = false;
            let mut v: Vec<Expr> = Vec::with_capacity(terms.len() + 1);
            let mut keep: Vec<(Expr, crate::number::Number)> = Vec::new();
            for (t, k) in terms {
                let nt = f(t);
                if nt.same_presentation(t) {
                    keep.push((t.clone(), k.clone()));
                } else {
                    changed = true;
                    v.push(mul2(&Expr::num(k.clone()), &nt));
                }
            }
            if !changed {
                return e.clone();
            }
            v.push(crate::ops::add_from_terms(c.clone(), keep));
            add(&v)
        }
        Kind::Mul(c, facs) => {
            let mut changed = false;
            let mut v: Vec<Expr> = Vec::with_capacity(facs.len() + 1);
            v.push(Expr::num(c.clone()));
            for (b, x) in facs {
                let nb = f(b);
                let nx = if x.is_number() { x.clone() } else { f(x) };
                if !(nb.same_presentation(b) && nx.same_presentation(x)) {
                    changed = true;
                }
                v.push(if x.is_one() && nx.is_one() { nb } else { pow(&nb, &nx) });
            }
            if !changed {
                return e.clone();
            }
            mul(&v)
        }
        Kind::Pow(b, x) => {
            let (nb, nx) = (f(b), f(x));
            if nb.same_presentation(b) && nx.same_presentation(x) {
                return e.clone();
            }
            pow(&nb, &nx)
        }
        Kind::Func(fun, args) => {
            if let Func::MatrixElement(parent) | Func::MatrixProxy(parent) | Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent) = fun {
                let proxy = |matrix: crate::matrixexpr::MatrixExpr| Expr::from_kind(Kind::Func(Func::MatrixProxy(std::sync::Arc::new(matrix)), vec![]));
                let mapped = if !matches!(fun, Func::MatrixProxy(_)) {
                    match f(&proxy((**parent).clone())).kind() {
                        Kind::Func(Func::MatrixProxy(new), _) => Ok((**new).clone()),
                        _ => Ok((**parent).clone()),
                    }
                } else {
                    // Both callbacks descend through f; do not map scalars in
                    // a child again after f has already transformed it.
                    let callback = std::cell::RefCell::new(&mut *f);
                    parent.map_children(&mut |child| {
                        match callback.borrow_mut()(&proxy(child.clone())).kind() {
                            Kind::Func(Func::MatrixProxy(new), _) => (**new).clone(),
                            _ => child.clone(),
                        }
                    }, &mut |scalar| callback.borrow_mut()(scalar))
                };
                if let Ok(matrix) = mapped {
                    let na: Vec<Expr> = args.iter().map(&mut *f).collect();
                    if matrix == **parent && na == *args { return e.clone(); }
                    if matches!(fun, Func::MatrixElement(_)) && na.len() == 2 {
                        return matrix.element(&na[0], &na[1]).unwrap_or_else(|_| e.clone());
                    }
                    if matches!(fun, Func::MatrixTrace(_)) {
                        return Expr::from_kind(Kind::Func(Func::MatrixTrace(std::sync::Arc::new(matrix)), vec![]));
                    }
                    if matches!(fun, Func::MatrixDeterminant(_)) {
                        return Expr::from_kind(Kind::Func(Func::MatrixDeterminant(std::sync::Arc::new(matrix)), vec![]));
                    }
                    return Expr::from_kind(Kind::Func(Func::MatrixProxy(std::sync::Arc::new(matrix)), vec![]));
                }
                return e.clone();
            }
            let na: Vec<Expr> = args.iter().map(|a| f(a)).collect();
            if na.iter().zip(args).all(|(a, b)| a.same_presentation(b)) {
                return e.clone();
            }
            func(fun.clone(), na)
        }
        Kind::Derivative(inner, vars) => {
            let ni = f(inner);
            let nv:Vec<_>=vars.iter().map(|(v,n)|(f(v),*n)).collect();
            if ni.same_presentation(inner) && nv.iter().zip(vars).all(|((a,n),(b,m))|n==m && a.same_presentation(b)) {
                return e.clone();
            }
            crate::diff::derivative_unevaluated(&ni,nv)
        }
        Kind::Rel(op, a, b) => {
            let (na, nb) = (f(a), f(b));
            if na.same_presentation(a) && nb.same_presentation(b) {
                return e.clone();
            }
            rel(*op, &na, &nb)
        }
        Kind::Piecewise(p) => {
            let np: Vec<(Expr, Expr)> = p.iter().map(|(x, c)| (f(x), f(c))).collect();
            if np.iter().zip(p).all(|((a,c),(b,d))|a.same_presentation(b) && c.same_presentation(d)) {
                return e.clone();
            }
            piecewise(np)
        }
    }
}

/// Pre-order iterator over all sub-expressions (with repetitions).
pub fn preorder(e: &Expr, f: &mut dyn FnMut(&Expr)) {
    f(e);
    for a in e.args() {
        preorder(&a, f);
    }
}

/// Set of free symbols.
pub fn free_symbols(e: &Expr) -> Vec<Expr> {
    e.cached_free_symbols(|| collect_free_symbols(e))
}

fn collect_free_symbols(e: &Expr) -> Vec<Expr> {
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut out = Vec::new();
    fn rec(e: &Expr, seen: &mut FxHashSet<Expr>, out: &mut Vec<Expr>, visited: &mut FxHashSet<Expr>) {
        if e.mask() == 0 {
            return;
        }
        if e.is_symbol() {
            if seen.insert(e.clone()) {
                out.push(e.clone());
            }
            return;
        }
        if !visited.insert(e.clone()) {
            return;
        }
        match e.kind() {
            Kind::Func(Func::Limit,args)=>{
                for symbol in free_symbols(&args[0]).into_iter().filter(|s|*s!=args[1]){
                    if seen.insert(symbol.clone()){out.push(symbol);}
                }
                rec(&args[2],seen,out,visited);
            },
            Kind::Func(Func::Subs,args)=>{
                let variables=args[1].args();
                for symbol in free_symbols(&args[0]).into_iter().filter(|s|!variables.contains(s)){
                    if seen.insert(symbol.clone()){out.push(symbol);}
                }
                rec(&args[2],seen,out,visited);
            },
            Kind::Func(Func::MatrixProxy(parent), _) => {
                for matrix in parent.matrix_symbols() {
                    let proxy = Expr::from_kind(Kind::Func(Func::MatrixProxy(std::sync::Arc::new(matrix)), vec![]));
                    if seen.insert(proxy.clone()) { out.push(proxy); }
                }
                for symbol in parent.scalar_free_symbols() { rec(&symbol, seen, out, visited); }
            }
            Kind::Func(Func::MatrixElement(parent) | Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent), args) => {
                rec(&Expr::from_kind(Kind::Func(Func::MatrixProxy(parent.clone()), vec![])), seen, out, visited);
                for index in args { rec(index, seen, out, visited); }
            }
            Kind::Func(crate::expr::Func::Set(op), args)
                if matches!(op, crate::sets::SetOp::Lambda | crate::sets::SetOp::ConditionSet)
                    && args.len() >= 2 => {
                // The binder applies to the body/condition, but not the base set.
                for symbol in free_symbols(&args[1]).into_iter().filter(|s| *s != args[0]) {
                    if seen.insert(symbol.clone()) {
                        out.push(symbol);
                    }
                }
                for base in &args[2..] {
                    rec(base, seen, out, visited);
                }
            }
            Kind::Add(_, terms) => {
                for (t, _) in terms {
                    rec(t, seen, out, visited)
                }
            }
            Kind::Mul(_, facs) => {
                for (b, x) in facs {
                    rec(b, seen, out, visited);
                    rec(x, seen, out, visited);
                }
            }
            _ => {
                for a in e.args() {
                    rec(&a, seen, out, visited)
                }
            }
        }
    }
    let mut visited = FxHashSet::default();
    rec(e, &mut seen, &mut out, &mut visited);
    out.sort_by(|a, b| a.cmp_struct(b));
    out
}

/// Does `e` contain the sub-expression `s` (structurally)?
pub fn has(e: &Expr, s: &Expr) -> bool {
    if !s.has_binding_identity() && s.mask() & !e.mask() != 0 {
        return false;
    }
    if e == s {
        return true;
    }
    match e.kind() {
        Kind::Func(Func::MatrixElement(parent) | Func::MatrixProxy(parent) | Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent), args) => {
            let matrix_match = if let Kind::Func(Func::MatrixProxy(target), _) = s.kind() {
                let mut pending = vec![parent.as_ref().clone()];
                let mut seen = FxHashSet::default();
                let mut found = false;
                while let Some(matrix) = pending.pop() {
                    if !seen.insert(matrix.dag_id()) { continue; }
                    if &matrix == target.as_ref() { found = true; break; }
                    pending.extend(matrix.children());
                }
                found
            } else { false };
            if matrix_match { return true; }
            let mut fields = Vec::new();
            parent.scalar_fields(&mut fields);
            fields.iter().any(|scalar| has(scalar, s)) || args.iter().any(|a| has(a, s))
        }
        Kind::Add(c, terms) => {
            if let Kind::Add(sc, subterms) = s.kind() {
                // `x + y` is contained in `x + y + z`
                // Structural containment does not allow the sign changes or
                // constant offsets supported by algebraic substitution.
                if (sc.is_zero() || c == sc)
                    && subterms.iter().all(|term| terms.contains(term))
                {
                    return true;
                }
            }
            terms.iter().any(|(t, _)| has(t, s))
        }
        Kind::Mul(_, facs) => facs.iter().any(|(b, x)| has(b, s) || has(x, s)),
        _ => e.args().iter().any(|a| has(a, s)),
    }
}

/// Number of nodes (SymPy's `count_ops`-like size measure).
pub fn node_count(e: &Expr) -> usize {
    1 + e.args().iter().map(node_count).sum::<usize>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbol_cache_preserves_assumptions_and_replacements() {
        let x = Expr::symbol("x");
        let y = Expr::symbol_with("x", crate::symbol::Assumptions::from_kwargs([("positive", true)]).unwrap());
        let e = add(&[x.clone(), y.clone()]);
        assert_eq!(free_symbols(&e).len(), 2);
        let mut returned = free_symbols(&e);
        returned.clear();
        assert_eq!(free_symbols(&e).len(), 2);
        let map = [(x, Expr::one())].into_iter().collect();
        let replaced = crate::subs::xreplace(&e, &map);
        assert_eq!(free_symbols(&replaced), vec![y]);
    }

    #[test]
    fn structural_sum_containment_differs_from_substitution() {
        let parse = |s| crate::parser::parse_str(s).unwrap();
        assert!(has(&parse("x+y+z"), &parse("x+y")));
        assert!(!has(&parse("x+y"), &parse("x+y+1")));
        assert!(!has(&parse("x-y"), &parse("y-x")));
        assert!(!has(&parse("2*x+2*y"), &parse("x+y")));
    }
}

/// Ordering used by SymPy's immutable `.args`, which differs from the
/// presentation ordering used by default_sort_key/as_ordered_terms.
fn args_cmp(a: &Expr, b: &Expr) -> std::cmp::Ordering {
    fn rank(e: &Expr) -> (u8, String) {
        match e.kind() {
            Kind::Num(_) => (0, String::new()),
            Kind::Const(_) => (1, e.class_name().into()),
            Kind::Sym(_) => (2, String::new()),
            Kind::Pow(..) => (3, String::new()),
            Kind::Mul(..) => (4, String::new()),
            Kind::Add(..) => (5, String::new()),
            _ => (6, e.class_name().into()),
        }
    }
    let ordering = rank(a).cmp(&rank(b));
    if ordering != std::cmp::Ordering::Equal { return ordering; }
    match (a.kind(), b.kind()) {
        (Kind::Num(x), Kind::Num(y)) => x.cmp_value(y).unwrap_or_else(|| a.cmp_struct(b)),
        (Kind::Sym(x), Kind::Sym(y)) => x.cmp(y),
        _ => {
            let aa = ordered_args(a);
            let bb = ordered_args(b);
            aa.len().cmp(&bb.len()).then_with(|| {
                for (x, y) in aa.iter().zip(&bb) {
                    let comparison = args_cmp(x, y);
                    if comparison != std::cmp::Ordering::Equal { return comparison; }
                }
                a.cmp_struct(b)
            })
        }
    }
}

pub fn ordered_args(e: &Expr) -> Vec<Expr> {
    e.cached_ordered_args(|| {
        if let Kind::Derivative(inner,vars)=e.kind(){
            let mut args=vec![inner.clone()];
            args.extend(vars.iter().map(|(v,n)|crate::functions::undefined("Tuple",vec![v.clone(),Expr::int(*n as i64)])));
            return args;
        }
        let mut args = e.args();
        let canonical = e.is_canonical() && (e.is_add() || e.is_mul());
        if canonical { args.sort_by(args_cmp); }
        args
    })
}
