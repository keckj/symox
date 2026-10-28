//! Generic tree traversal helpers.

use rustc_hash::FxHashSet;

use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::ops::{add, mul, mul2, pow};
use crate::relational::{piecewise, rel};

/// Rebuild `e` with every direct child mapped through `f`, using the
/// canonical constructors. Returns `e` itself when nothing changed.
pub fn map_children(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    if matches!(e.kind(),Kind::Func(Func::SymbolicDerivative,_)){return e.clone();}
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
            if e.is_deferred_numeric_power(){crate::ops::deferred::rebuild(&nb,&nx)}else{pow(&nb,&nx)}
        }
        Kind::Func(fun, args) => {
            if let Func::MatrixDotProduct(a,b)=fun {
                let mut mapped=|parent:&std::sync::Arc<crate::matrixexpr::MatrixExpr>|{
                    let proxy=Expr::from_kind(Kind::Func(Func::MatrixProxy(parent.clone()),vec![]));
                    if let Kind::Func(Func::MatrixProxy(new),_)=f(&proxy).kind(){new.clone()}else{parent.clone()}
                };
                let na=mapped(a);let nb=mapped(b);
                if na.dag_id()==a.dag_id()&&nb.dag_id()==b.dag_id(){return e.clone();}
                return crate::matrixcontraction::dot_product(&na,&nb,false).unwrap_or_else(|_|e.clone());
            }

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
                        return matrix.element_constructed(&na[0], &na[1]).unwrap_or_else(|_| e.clone());
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
            Kind::Func(f @ (Func::Integral|Func::Sum|Func::Product),args)=>{
                for symbol in crate::binders::free_symbols(f,args){if seen.insert(symbol.clone()){out.push(symbol);}}
            },
            Kind::Func(Func::Special(f),args) if f.is_transform() && args.len()>=3 => {
                for symbol in free_symbols(&args[0]).into_iter().chain(free_symbols(&args[2])).filter(|s|*s!=args[1]){
                    if seen.insert(symbol.clone()){out.push(symbol);}
                }
            },
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
            Kind::Func(Func::MatrixDotProduct(a,b),_)=>{
                for parent in [a,b]{rec(&Expr::from_kind(Kind::Func(Func::MatrixProxy(parent.clone()),vec![])),seen,out,visited);}
            }
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
            Kind::Func(Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction),_)=>{
                if let Some((variables,body))=crate::sets::lambda_parts(e){
                    for symbol in free_symbols(&body).into_iter().filter(|s|!variables.contains(s)){
                        if seen.insert(symbol.clone()){out.push(symbol);}
                    }
                }
            },
            Kind::Func(crate::expr::Func::Set(op), args)
                if matches!(op, crate::sets::SetOp::ConditionSet)
                    && args.len() >= 2 => {
                // The binder applies to the body/condition, but not the base set.
                let bound = if crate::special::tuple_items(&args[0]).is_some() {
                    crate::sets::lambda_checked(&args[0], &args[1]).ok()
                        .and_then(|lambda| crate::sets::lambda_parts(&lambda).map(|(variables, _)| variables))
                        .unwrap_or_else(|| vec![args[0].clone()])
                } else { vec![args[0].clone()] };
                for symbol in free_symbols(&args[1]).into_iter().filter(|s| !bound.contains(s)) {
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
            // Differentiation variables do not become free merely by being
            // present in a held derivative. Read its compact body directly;
            // Expr::args() expands each variable according to its order.
            Kind::Derivative(body,_) => rec(body,seen,out,visited),
            Kind::Func(Func::SymbolicDerivative,_)=>{
                if let Some((body,variables))=crate::diff::symbolic::parts(e){
                    rec(body,seen,out,visited);
                    for(_,count)in variables{rec(&count,seen,out,visited);}
                }
            },
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

/// Does `e` contain the sub-expression `s` (structurally)? The generator
/// of a root object is bound: algorithms that ask for dependence do not see
/// it. `has_structural` also visits it, as the reference's `Basic.has`.
pub fn has(e: &Expr, s: &Expr) -> bool {
    has_in(e, s, false)
}

/// `has` that also visits the bound polynomial argument of a root object.
pub fn has_structural(e: &Expr, s: &Expr) -> bool {
    has_in(e, s, true)
}

fn has_in(e: &Expr, s: &Expr, structural: bool) -> bool {
    if !s.has_binding_identity() && !(structural && e.has_bound_generator()) && s.mask() & !e.mask() != 0 {
        return false;
    }
    if e == s {
        return true;
    }
    if !structural && crate::poly::rootof::is_rootof(e) {
        return false;
    }
    match e.kind() {
        Kind::Derivative(body,variables)=>has_in(body,s,structural) || variables.iter().any(|(variable,count)|
            has_in(&crate::special::tuple(vec![variable.clone(),Expr::int(*count as i64)]),s,structural)),
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>[a,b].iter().any(|parent|has_in(&Expr::from_kind(Kind::Func(Func::MatrixProxy((*parent).clone()),vec![])),s,structural)),
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
            fields.iter().any(|scalar| has_in(scalar, s, structural)) || args.iter().any(|a| has_in(a, s, structural))
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
            // Visit the public structural arguments, including the additive
            // constant and each term's explicit coefficient. Internal storage
            // keeps these numbers apart from the term bases.
            e.args().iter().any(|a| has_in(a, s, structural))
        }
        // The public arguments include an explicit coefficient, but exclude
        // the implicit exponent one in the compact factor representation.
        Kind::Mul(..) => (s.is_mul() && (e.has_noncommutative_host() || s.has_noncommutative_host()) && noncommutative_mul_contains(e, s))
            || e.args().iter().any(|a| has_in(a, s, structural)),
        _ => e.args().iter().any(|a| has_in(a, s, structural)),
    }
}

/// SymPy's `AssocOp._has_matcher` for products: the commutative arguments
/// of `s` form a subset of those of `e` and its noncommutative arguments a
/// consecutive run of those of `e`.
fn noncommutative_mul_contains(e: &Expr, s: &Expr) -> bool {
    let split = |m: &Expr| -> (Vec<Expr>, Vec<Expr>) {
        m.args().into_iter().partition(|a| crate::decompose::commutativity(a) == Some(true))
    };
    let ((c, nc), (ec, enc)) = (split(s), split(e));
    if !c.iter().all(|a| ec.contains(a)) { return false; }
    nc.is_empty() || (nc.len() <= enc.len() && enc.windows(nc.len()).any(|w| w == nc.as_slice()))
}

/// Estimate cache retention from the native DAG, without constructing
/// public Add terms or Mul powers. Physical nodes and their backing vectors
/// are counted once, including large integer/rational/Float mantissas. Opaque
/// matrix graphs and Python wildcard callbacks bypass cross-call caches.
pub(crate) fn retention_weight(expressions:&[&Expr],limit:usize)->Option<usize>{
    use crate::number::Number;
    fn number(n:&Number)->usize{
        let bits=match n{
            Number::BigInt(n)=>n.bits(),
            Number::Rational(n)=>n.numer().bits().saturating_add(n.denom().bits()),
            Number::MpFloat(n,_)=>n.man.bits(),
            _=>0,
        };
        1usize.saturating_add(usize::try_from(bits.saturating_add(511)/512).unwrap_or(usize::MAX))
    }
    // Host nodes dispatch user hooks (`_eval_subs`, `_eval_evalf`...) on
    // every call and own Python instances that must stay collectable.
    if expressions.iter().any(|e|e.has_matrix()||e.has_wild_properties()||e.has_host()){return None;}
    let mut stack=expressions.to_vec();
    let mut seen=FxHashSet::default();
    let mut weight=0usize;
    while let Some(e)=stack.pop(){
        if !seen.insert(e.dag_id()){continue;}
        weight=weight.saturating_add(2);
        match e.kind(){
            Kind::Num(n)=>weight=weight.saturating_add(number(n)),
            Kind::Sym(s)=>weight=weight.saturating_add(s.name().len()/64+1),
            Kind::Add(c,terms)=>{
                weight=weight.saturating_add(number(c)).saturating_add(terms.len());
                if weight>limit{return None;}
                for(term,c)in terms{weight=weight.saturating_add(number(c));stack.push(term);}
            },
            Kind::Mul(c,factors)=>{
                weight=weight.saturating_add(number(c)).saturating_add(factors.len()/4+1);
                if weight>limit{return None;}
                for(base,exponent)in factors{stack.extend([base,exponent]);}
            },
            Kind::Pow(a,b)|Kind::Rel(_,a,b)=>stack.extend([a,b]),
            Kind::Func(head,args)=>{
                weight=weight.saturating_add(head.name().len()/64+1).saturating_add(args.len()/8+1);
                if weight>limit{return None;}
                stack.extend(args);
            },
            Kind::Derivative(body,variables)=>{
                weight=weight.saturating_add(variables.len()/4+1);
                if weight>limit{return None;}
                stack.push(body);stack.extend(variables.iter().map(|(v,_)|v));
            },
            Kind::Piecewise(pairs)=>{
                weight=weight.saturating_add(pairs.len()/4+1);
                if weight>limit{return None;}
                for(value,condition)in pairs{stack.extend([value,condition]);}
            },
            Kind::Const(_)|Kind::Bool(_)=>{},
        }
        if weight>limit{return None;}
    }
    Some(weight)
}

/// Number of nodes (SymPy's `count_ops`-like size measure).
pub fn node_count(e: &Expr) -> usize {
    fn count(e:&Expr,depth:usize,memo:&mut rustc_hash::FxHashMap<usize,(Expr,usize)>)->usize {
        if let Some((_,size))=memo.get(&e.dag_id()) {return *size;}
        // A size heuristic must remain finite on deeply nested/shared DAGs.
        if depth>=512 {return usize::MAX;}
        let size=if let Kind::Derivative(body,variables)=e.kind() {
            variables.iter().fold(count(body,depth+1,memo).saturating_add(1),|size,(variable,_)|
                size.saturating_add(count(variable,depth+1,memo)).saturating_add(2))
        } else {
            e.args().iter().fold(1usize,|size,child|size.saturating_add(count(child,depth+1,memo)))
        };
        memo.insert(e.dag_id(),(e.clone(),size));
        size
    }
    count(e,0,&mut rustc_hash::FxHashMap::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retention_accounting_is_bounded_by_the_physical_graph_and_literal_size(){
        let mut value=Expr::symbol("x");
        for _ in 0..80{value=crate::functions::undefined("pair",vec![value.clone(),value]);}
        assert!(retention_weight(&[&value],1024).is_some());
        assert!(retention_weight(&[&value],10).is_none());
        let wide=Expr::num(crate::number::Number::from_bigint(num_bigint::BigInt::from(1)<<1_000_000usize));
        assert!(retention_weight(&[&wide],128).is_none());
        assert!(retention_weight(&[&Expr::one()],128).is_some());
    }

    #[test]
    fn node_size_memoizes_shared_graphs_and_compact_derivatives() {
        let x=Expr::symbol("x");
        let body=crate::functions::undefined("f",vec![x.clone()]);
        let derivative=Expr::from_kind(Kind::Derivative(body,vec![(x.clone(),1_000_000_000)]));
        assert_eq!(node_count(&derivative),6);
        let mut expression=x;
        for _ in 0..80 {expression=crate::functions::undefined("f",vec![expression.clone(),expression]);}
        assert_eq!(node_count(&expression),usize::MAX);
    }

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

    #[test]
    fn containment_visits_explicit_arithmetic_arguments() {
        let parse = |s| crate::parser::parse_str(s).unwrap();
        for text in ["3*x", "x+3", "3*x+y", "x+3*y"] {
            assert!(has(&parse(text), &Expr::int(3)), "{text}");
        }
        for text in ["x*y", "3*x", "x+y", "x**2"] {
            assert!(!has(&parse(text), &Expr::one()), "{text}");
        }
        assert!(has(&parse("x+1"), &Expr::one()));
        assert!(has(&parse("x**3"), &Expr::int(3)));
        assert!(!has(&parse("x/3"), &Expr::int(3)));
    }
}

/// Ordering used by SymPy's immutable `.args`, which differs from the
/// presentation ordering used by default_sort_key/as_ordered_terms.
/// SymPy compares two Floats by their hashable content, the mpmath tuple
/// `(sign, man, exp, bc)` (odd mantissa): positive values precede negative
/// ones and magnitudes are not ordered by value.
fn mpf_cmp(x: &crate::number::Number, y: &crate::number::Number) -> Option<std::cmp::Ordering> {
    use crate::number::Number;
    use num_traits::{Signed, Zero};
    fn tuple(n: &Number) -> Option<(bool, num_bigint::BigInt, crate::bigfloat::Exponent, u64)> {
        let f = match n {
            Number::Float(v) if v.is_finite() && *v != 0.0 => crate::bigfloat::BigFloat::from_f64(*v),
            Number::MpFloat(v, _) if !v.is_zero() => v.clone(),
            _ => return None,
        };
        let (mut man, mut exp) = (f.man.abs(), f.exp.clone());
        if man.is_zero() { return None; }
        let shift = man.trailing_zeros().unwrap_or(0);
        man >>= shift;
        exp += shift;
        let bc = man.bits();
        Some((f.man.is_negative(), man, exp, bc))
    }
    let ordering = tuple(x)?.cmp(&tuple(y)?);
    (ordering != std::cmp::Ordering::Equal).then_some(ordering)
}

pub(crate) fn args_cmp(a: &Expr, b: &Expr) -> std::cmp::Ordering {
    if a.ptr_eq(b) { return std::cmp::Ordering::Equal; }
    fn rank(e: &Expr) -> (usize, &str) {
        // SymPy core.basic.ordering_of_classes / _cmp_name. Names absent
        // from this historical table follow every known class, then compare
        // lexically: notably ComplexInfinity, Dummy and lowercase sin.
        const CLASSES: &[&str] = &[
            "Zero", "One", "Half", "Infinity", "NaN", "NegativeOne", "NegativeInfinity",
            "Integer", "Rational", "Float", "Exp1", "Pi", "ImaginaryUnit", "Symbol", "Wild",
            "Pow", "Mul", "Add", "Derivative", "Integral", "Abs", "Sign", "Sqrt", "Floor",
            "Ceiling", "Re", "Im", "Arg", "Conjugate", "Exp", "Log", "Sin", "Cos", "Tan",
            "Cot", "ASin", "ACos", "ATan", "ACot", "Sinh", "Cosh", "Tanh", "Coth",
            "ASinh", "ACosh", "ATanh", "ACoth", "RisingFactorial", "FallingFactorial",
            "factorial", "binomial", "Gamma", "LowerGamma", "UpperGamma", "PolyGamma", "Erf",
            "Chebyshev", "Chebyshev2", "Function", "WildFunction", "Lambda", "Order",
            "Equality", "Unequality", "StrictGreaterThan", "StrictLessThan", "GreaterThan", "LessThan",
        ];
        // SymPy ranks type(self).__name__, the host class of a host node.
        let name = if e.is_host() { crate::host::class_name(e).unwrap_or(e.class_name()) } else { e.class_name() };
        (CLASSES.iter().position(|class| *class == name).unwrap_or(CLASSES.len()), name)
    }
    let ordering = rank(a).cmp(&rank(b));
    if ordering != std::cmp::Ordering::Equal { return ordering; }
    match (a.kind(), b.kind()) {
        (Kind::Num(x), Kind::Num(y)) => mpf_cmp(x, y).unwrap_or_else(|| x.cmp_value(y).unwrap_or_else(|| a.cmp_struct(b))),
        (Kind::Sym(x), Kind::Sym(y)) => x.cmp_basic(y),
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
        if canonical && e.is_mul() && e.has_noncommutative_host() {
            // Factors that may not commute keep their order after the
            // sorted commutative part, as in SymPy's Mul._from_args.
            let boundary = args.iter().position(|factor| crate::decompose::commutativity(factor) != Some(true)).unwrap_or(args.len());
            args[..boundary].sort_by(args_cmp);
        } else if canonical { args.sort_by(args_cmp); }
        if canonical && e.has_host() && args.iter().any(coefficient_host) {
            // SymPy Add.flatten and Mul.flatten hold an AccumBounds term in
            // the coefficient slot, which precedes the sorted terms.
            let (mut front, rest): (Vec<Expr>, Vec<Expr>) = args.into_iter().partition(coefficient_host);
            front.extend(rest);
            return front;
        }
        args
    })
}

/// Whether a canonical sum or product holds an AccumBounds term next to
/// a Number or to other bounds, which SymPy folds into one coefficient
/// through the host arithmetic of the bounds.
pub fn has_unfolded_bounds(e: &Expr) -> bool {
    if !e.has_host() || !e.is_canonical() || !(e.is_add() || e.is_mul()) { return false; }
    let args = ordered_args(e);
    let bounds = args.iter().take_while(|a| coefficient_host(a)).count();
    // A distributed Number times bounds, as in 2*(x + AccumBounds(1, 2)).
    if bounds == 0 { return e.is_add() && args.iter().any(|a| a.is_mul() && has_unfolded_bounds(a)); }
    bounds > 1 || args.iter().any(|a| a.is_number() || matches!(a.kind(), Kind::Const(
        crate::expr::Constant::Infinity | crate::expr::Constant::NegInfinity
        | crate::expr::Constant::ComplexInfinity | crate::expr::Constant::NaN)))
}

/// An AccumBounds node, which SymPy's Add and Mul keep as their coefficient.
fn coefficient_host(e: &Expr) -> bool {
    e.is_host() && crate::host::class_name(e) == Some("AccumulationBounds")
}
