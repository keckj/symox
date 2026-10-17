//! Structural inspection without expansion or evaluation of held arithmetic.
use rustc_hash::{FxHashMap, FxHashSet};

use crate::expr::{Constant, Expr, Func, Kind};

/// Public scalar children, with compact coefficients/powers made explicit.
/// These presentation nodes must never pass through the active constructors:
/// inspecting a held `2**1000000000` must not allocate that integer.
pub fn scalar_children(e: &Expr) -> Vec<Expr> {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => Vec::new(),
        Kind::Add(c, terms) => {
            let mut children = Vec::with_capacity(terms.len()+1);
            if !c.is_zero() { children.push(Expr::num(c.clone())); }
            children.extend(terms.iter().map(|(term, coefficient)| {
                if coefficient.is_one() { term.clone() } else {
                    Expr::from_unevaluated_kind(Kind::Mul(coefficient.clone(), vec![(term.clone(), Expr::one())]))
                }
            }));
            children
        }
        Kind::Mul(c, factors) => {
            let mut children = Vec::with_capacity(factors.len()+1);
            if !c.is_one() { children.push(Expr::num(c.clone())); }
            children.extend(factors.iter().map(|(base, exponent)| {
                if exponent.is_one() { base.clone() } else {
                    Expr::from_unevaluated_kind(Kind::Pow(base.clone(), exponent.clone()))
                }
            }));
            children
        }
        Kind::Pow(base, exponent) | Kind::Rel(_, base, exponent) => vec![base.clone(), exponent.clone()],
        Kind::Derivative(body, variables) => {
            let mut children = vec![body.clone()];
            children.extend(variables.iter().map(|(variable, count)| {
                crate::functions::undefined("Tuple", vec![variable.clone(), Expr::int(*count as i64)])
            }));
            children
        }
        Kind::Func(_, children) => children.clone(),
        Kind::Piecewise(branches) => branches.iter().flat_map(|(value, condition)| [value.clone(), condition.clone()]).collect(),
    }
}

/// Unique physical nodes and their public scalar children. Parent records keep
/// synthesized presentation nodes alive while their pointer IDs are memoized.
pub fn scalar_nodes(e: &Expr) -> Vec<(Expr, Vec<Expr>)> {
    let mut seen = FxHashSet::default();
    let mut pending = vec![e.clone()];
    let mut nodes = Vec::new();
    while let Some(node) = pending.pop() {
        if !seen.insert(node.dag_id()) { continue; }
        let children = scalar_children(&node);
        pending.extend(children.iter().cloned());
        nodes.push((node, children));
    }
    nodes
}

#[derive(Clone, Copy, Debug)]
pub enum PolynomialKind { Polynomial, Rational, Algebraic }

/// Structural classification, matching the three-valued expression queries.
/// A generic function depending on a generator remains unknown. Arithmetic
/// containers require every child to be classified positively, like SymPy.
pub fn classify(e: &Expr, generators: &[Expr], kind: PolynomialKind) -> Option<bool> {
    if generators.is_empty() {
        return Some(!matches!(kind, PolynomialKind::Rational) || !illegal(e));
    }
    let mut context = Classifier {
        generators, kind,
        dependencies: FxHashMap::default(),
        classifications: FxHashMap::default(),
        children: FxHashMap::default(),
        free_dependencies: FxHashMap::default(),
    };
    context.classify(e)
}

fn illegal(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN))
}

struct Classifier<'a> {
    generators: &'a [Expr],
    kind: PolynomialKind,
    dependencies: FxHashMap<usize, bool>,
    classifications: FxHashMap<usize, Option<bool>>,
    children: FxHashMap<usize, (Expr, Vec<Expr>)>,
    free_dependencies: FxHashMap<(usize, Vec<Expr>), bool>,
}
impl Classifier<'_> {
    fn children(&mut self, e: &Expr) -> Vec<Expr> {
        self.children.entry(e.dag_id()).or_insert_with(|| (e.clone(), scalar_children(e))).1.clone()
    }
    fn contains(&mut self, e: &Expr) -> bool {
        if let Some(answer) = self.dependencies.get(&e.dag_id()) { return *answer; }
        let mut found = self.generators.iter().any(|generator| generator == e);
        if !found {
            // SymPy's structural has() includes a subset of additive terms.
            if let Kind::Add(constant, terms) = e.kind() {
                found = self.generators.iter().any(|generator| {
                    if let Kind::Add(c, ts) = generator.kind() {
                        (c.is_zero() || c == constant) && ts.iter().all(|term| terms.contains(term))
                    } else { false }
                });
            }
        }
        if !found&&e.is_mul(){
            let children=self.children(e);
            found=self.generators.iter().any(|generator|generator.is_mul()&&scalar_children(generator).iter().all(|child|children.contains(child)));
        }
        if !found {
            found = self.children(e).iter().any(|child| self.contains(child));
        }
        self.dependencies.insert(e.dag_id(), found);
        found
    }
    fn free_dependency(&mut self, e: &Expr) -> bool {
        self.free_contains(e, &[])
    }
    fn free_contains(&mut self, e:&Expr, bound:&[Expr])->bool{
        let key=(e.dag_id(),bound.to_vec());
        if let Some(found)=self.free_dependencies.get(&key){return *found;}
        let children=if e.is_mul(){self.children(e)}else{Vec::new()};
        // Only free candidate nodes participate in a match. A partial pattern
        // y+z is not a free argument of Integral(x+y+z, (x,0,1)): that Add also
        // depends on bound x. Its individual free children remain searchable.
        let eligible=bound.is_empty()||crate::traverse::free_symbols(e).iter().all(|symbol|!bound.contains(symbol));
        let matches=eligible&&self.generators.iter().any(|generator|{
            let exact=generator==e;
            let partial=!matches!(self.kind,PolynomialKind::Rational);
            let additive=if let (Kind::Add(c,terms),Kind::Add(tc,tts))=(e.kind(),generator.kind()){
                (tc.is_zero()||c==tc)&&tts.iter().all(|term|terms.contains(term))
            }else{false};
            let multiplicative=e.is_mul()&&generator.is_mul()&&scalar_children(generator).iter().all(|child|children.contains(child));
            (exact||partial&&(additive||multiplicative))&&crate::traverse::free_symbols(generator).iter().all(|symbol|!bound.contains(symbol))
        });
        let found=if matches{true}else{match e.kind(){
            Kind::Func(function @ (Func::Integral|Func::Sum),args)=>{
                let limits=crate::binders::limits(args);
                let mut bound=bound.to_vec();
                let mut found=false;
                for limit in limits.iter().rev(){
                    if limit[1..].iter().any(|endpoint|self.free_contains(endpoint,&bound)){found=true;break;}
                    if *function==Func::Integral&&limit.len()==1{
                        if self.free_contains(&limit[0],&bound){found=true;break;}
                    }else{bound.push(limit[0].clone());}
                }
                found||self.free_contains(&args[0],&bound)
            },
            Kind::Func(Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction),_)=>{
                if let Some((variables,body))=crate::sets::lambda_parts(e){
                    let mut bound=bound.to_vec();bound.extend(variables);
                    self.free_contains(&body,&bound)
                }else{false}
            },
            Kind::Func(Func::Set(crate::sets::SetOp::ConditionSet),args) if args.len()>=2=>{
                let mut local=bound.to_vec();local.push(args[0].clone());
                self.free_contains(&args[1],&local)||args[2..].iter().any(|base|self.free_contains(base,bound))
            },
            _=>self.children(e).iter().any(|child|self.free_contains(child,bound)),
        }};
        self.free_dependencies.insert(key,found);
        found
    }
    fn power(&mut self, base: &Expr, exponent: &Expr) -> Option<bool> {
        if self.contains(exponent) { return Some(false); }
        if !self.contains(base) { return Some(true); }
        let valid_exponent = match self.kind {
            PolynomialKind::Polynomial => exponent.as_number().is_some_and(|number| number.is_integer() && !number.is_negative()),
            PolynomialKind::Rational => exponent.is_integer_number(),
            PolynomialKind::Algebraic => exponent.is_rational_number(),
        };
        if matches!(self.kind,PolynomialKind::Polynomial){
            return Some(self.classify(base)==Some(true)&&valid_exponent);
        }
        match self.classify(base) {
            Some(true) => Some(valid_exponent),
            answer => answer,
        }
    }
    fn classify(&mut self, e: &Expr) -> Option<bool> {
        if let Some(answer) = self.classifications.get(&e.dag_id()) { return *answer; }
        let answer = match e.kind() {
            Kind::Num(_) | Kind::Sym(_) => Some(true),
            Kind::Const(_) => Some(!matches!(self.kind, PolynomialKind::Rational) || !illegal(e)),
            Kind::Add(_, _) | Kind::Mul(_, _) => {
                Some(self.children(e).iter().all(|child| self.classify(child) == Some(true)))
            }
            Kind::Pow(base, exponent) => self.power(base, exponent),
            _ if self.generators.iter().any(|generator| generator == e) => Some(true),
            _ if !self.free_dependency(e) => Some(true),
            _ => None,
        };
        self.classifications.insert(e.dag_id(), answer);
        answer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::number::Number;
    #[test]
    fn polynomial_queries_preserve_held_powers_and_shared_graphs() {
        let x = Expr::symbol("x");
        let power = Expr::from_unevaluated_kind(Kind::Pow(x.clone(), Expr::int(1_000_000_000)));
        assert_eq!(classify(&power, &[x.clone()], PolynomialKind::Polynomial), Some(true));
        let coefficient = Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2), Expr::int(1_000_000_000)));
        let held = Expr::from_unevaluated_kind(Kind::Mul(Number::ONE, vec![(coefficient, Expr::one()), (x.clone(), Expr::one())]));
        assert_eq!(classify(&held, &[x.clone()], PolynomialKind::Polynomial), Some(true));
        assert!(scalar_nodes(&held).len() < 10);
        let mut shared = x.clone();
        for _ in 0..80 { shared = crate::functions::undefined("f", vec![shared.clone(), shared]); }
        assert!(scalar_nodes(&shared).len() < 100);
        assert_eq!(classify(&shared, &[x], PolynomialKind::Polynomial), None);
    }
}
