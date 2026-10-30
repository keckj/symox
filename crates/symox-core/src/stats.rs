//! Exact scalar probability kernels. Unsupported inputs retain the caller's
//! general probability model; None never denotes a zero expectation.
use crate::{Expr, Number, ops};
pub mod markov;
use crate::poly::expression_domain::dict_from_expr;
use num_bigint::BigInt;

#[derive(Clone)]
pub enum Law {
    Binomial { n: Expr, p: Expr, success: Expr, failure: Expr },
    Normal { mean: Expr, sigma: Expr },
    Poisson { rate: Expr },
    Finite(Vec<(Expr, Expr)>),
}

/// Port of the complete raw-moment formula, including arbitrary success and
/// failure values. The caller validates distribution parameters.
pub fn binomial_moment<E>(n: &Expr, p: &Expr, success: &Expr, failure: &Expr,
    order: usize, check: &mut impl FnMut() -> Result<(), E>) -> Result<Expr, E> {
    let mut stirling = vec![BigInt::from(1)];
    let mut base = Vec::new();
    for q in 0..=order {
        check()?;
        if q > 0 {
            let previous = std::mem::replace(&mut stirling, vec![BigInt::from(0); q + 1]);
            for j in 1..=q {
                stirling[j] = &previous[j - 1]
                    + previous.get(j).map(|s| s * j).unwrap_or_default();
            }
        }
        let terms = stirling.iter().enumerate().map(|(j, s)| {
            let coefficient = Expr::num(Number::from_bigint(s.clone()));
            let falling = crate::special::ff(n, &Expr::int(j as i64));
            ops::mul2(&ops::mul2(&coefficient, &falling), &ops::powi(p, j as i64))
        }).collect::<Vec<_>>();
        base.push(ops::add(&terms));
    }
    let location = ops::mul2(n, failure);
    let scale = ops::sub(success, failure);
    let mut terms = Vec::new();
    let mut choose = BigInt::from(1);
    for (q, moment) in base.into_iter().enumerate() {
        check()?;
        let coefficient = Expr::num(Number::from_bigint(choose.clone()));
        let term = ops::mul2(&coefficient, &ops::powi(&location, (order-q) as i64));
        terms.push(ops::mul2(&term, &ops::mul2(&ops::powi(&scale, q as i64), &moment)));
        if q < order { choose = choose * (order-q) / (q+1); }
    }
    Ok(ops::add(&terms))
}

/// Same term order and arithmetic grouping as the original Gaussian formula.
pub fn normal_moment<E>(mean: &Expr, sigma: &Expr, order: usize,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Expr, E> {
    let mut terms = Vec::new();
    for j in 0..=order/2 {
        check()?;
        let choose = crate::ntheory::binomial(&BigInt::from(order), &BigInt::from(2*j));
        let term = ops::mul2(&Expr::num(Number::from_bigint(choose)), &ops::powi(mean, (order-2*j) as i64));
        let term = ops::mul2(&term, &ops::powi(sigma, (2*j) as i64));
        let term = ops::mul2(&term, &Expr::num(Number::from_bigint(crate::ntheory::factorial((2*j) as u64))));
        let denominator = (BigInt::from(1) << j) * crate::ntheory::factorial(j as u64);
        terms.push(ops::div(&term, &Expr::num(Number::from_bigint(denominator))));
    }
    Ok(ops::add(&terms))
}

/// Bound only the optional native recurrence: the Python implementation
/// remains the fallback for every declined state space.
pub fn normal_raw_budget(powers:&[u32])->bool {
    powers.iter().try_fold(1usize,|states,power| {
        let width=usize::try_from(*power).ok()?.checked_add(1)?;
        states.checked_mul(width).filter(|states|*states<=100_000)
    }).is_some()
}

fn joint_parameter(value:&Expr)->bool {
    use crate::Kind;
    if value.has_host()||value.has_matrix()||value.has_nonfinite()||value.has_held_arithmetic()
        ||value.has_noncommutative_host() {return false;}
    let mut pending=vec![value.clone()];
    let mut seen=rustc_hash::FxHashSet::default();
    while let Some(value)=pending.pop() {
        if !seen.insert(value.dag_id()) {continue;}
        match value.kind() {
            Kind::Num(_)=>{},
            Kind::Sym(symbol) if symbol.host().is_none()&&symbol.scalar_identity().is_none()=>{},
            Kind::Add(_,terms)=>pending.extend(terms.iter().map(|(term,_)|term.clone())),
            Kind::Mul(_,factors)=>{
                if factors.iter().any(|(_,power)|!power.is_integer_number()) {return false;}
                pending.extend(factors.iter().map(|(base,_)|base.clone()));
            },
            Kind::Pow(base,power) if power.is_integer_number()=>pending.push(base.clone()),
            _=>return false,
        }
    }
    true
}

/// Gaussian integration by parts, with a memoized DAG of exponent vectors.
/// This preserves the original first-active-coordinate, product grouping and
/// increasing covariance-index sum order. The explicit stack accepts deep
/// moments without a Python or Rust recursive call chain.
/// The Boolean says whether the result has already received exact native
/// simplification; Float results retain the caller's existing final simplifier.
pub fn normal_raw<E>(powers:&[u32],mean:&[Expr],covariance:&[Expr],
    check:&mut impl FnMut()->Result<(),E>)->Result<Option<(Expr,bool)>,E> {
    check()?;
    if powers.iter().all(|power|*power==0) {return Ok(Some((Expr::one(),true)));}
    let dimension=powers.len();
    if mean.len()!=dimension||dimension.checked_mul(dimension)!=Some(covariance.len())
        ||!normal_raw_budget(powers)||!mean.iter().chain(covariance).all(joint_parameter) {
        return Ok(None);
    }
    let exact=mean.iter().chain(covariance).all(|value|value.float_info().max_precision.is_none());
    let root=powers.to_vec();
    let mut memo:rustc_hash::FxHashMap<Vec<u32>,Expr>=rustc_hash::FxHashMap::default();
    let mut pending=vec![(root.clone(),false)];
    while let Some((counts,visited))=pending.pop() {
        check()?;
        if memo.contains_key(&counts) {continue;}
        let Some(i)=counts.iter().position(|power|*power>0) else {
            memo.insert(counts,Expr::one());continue;
        };
        let mut remainder=counts.clone();remainder[i]-=1;
        if !visited {
            let mut dependencies=vec![remainder.clone()];
            for (j,power) in remainder.iter().enumerate() {
                if *power>0 {let mut lower=remainder.clone();lower[j]-=1;dependencies.push(lower);}
            }
            pending.push((counts,true));
            pending.extend(dependencies.into_iter().rev().map(|counts|(counts,false)));
            continue;
        }
        let mut value=ops::mul2(&mean[i],&memo[&remainder]);
        for (j,power) in remainder.iter().enumerate() {
            if *power>0 {
                let mut lower=remainder.clone();lower[j]-=1;
                let weight=ops::mul2(&Expr::int(i64::from(*power)),&covariance[i*dimension+j]);
                let term=ops::mul2(&weight,&memo[&lower]);
                value=ops::add2(&value,&term);
            }
        }
        memo.insert(counts,value);
    }
    check()?;
    let value=memo.remove(&root).unwrap();
    let value=if exact{crate::simplify::simplify(&value)}else{value};
    check()?;
    Ok(Some((value,exact)))
}

/// Touchard polynomial T_n(rate) = sum_j S(n,j) rate^j. One Stirling
/// row replaces n+1 independent evaluations of Stirling numbers. The
/// term order and arithmetic grouping remain the distribution's original ones.
pub fn poisson_moment<E>(rate:&Expr,order:usize,
    check:&mut impl FnMut()->Result<(),E>)->Result<Expr,E> {
    check()?;
    let mut row=vec![BigInt::from(0);order+1];
    row[0]=BigInt::from(1);
    for degree in 1..=order {
        check()?;
        for j in (1..=degree).rev() {
            row[j]=&row[j-1]+&row[j]*j;
        }
        row[0]=BigInt::from(0);
    }
    let mut terms=Vec::with_capacity(order+1);
    for (j,coefficient) in row.into_iter().enumerate() {
        check()?;
        let coefficient=Expr::num(Number::from_bigint(coefficient));
        terms.push(ops::mul2(&coefficient,&ops::powi(rate,j as i64)));
    }
    Ok(ops::add(&terms))
}

fn exact_scalar(e: &Expr) -> bool {
    fn rec(e:&Expr,seen:&mut rustc_hash::FxHashMap<Expr,bool>)->bool {
        use crate::expr::Kind;
        if let Some(result)=seen.get(e) {return *result;}
        if e.has_nonfinite() {return false;}
        let result=match e.kind() {
            Kind::Num(n) => n.is_rational(),
            Kind::Sym(s) => s.host().is_none() && s.scalar_identity().is_none() && crate::decompose::commutes(e),
            Kind::Add(c, terms) => c.is_rational() && terms.iter().all(|(t,c)| c.is_rational() && rec(t,seen)),
            Kind::Mul(c, terms) => c.is_rational() && terms.iter().all(|(b,p)| rec(b,seen) && p.is_integer_number()),
            Kind::Pow(b,p) => rec(b,seen) && p.is_integer_number(),
            _ => false,
        };
        seen.insert(e.clone(),result);result
    }
    rec(e,&mut rustc_hash::FxHashMap::default())
}

fn law_exact(law: &Law) -> bool {
    match law {
        Law::Binomial { n,p,success,failure } => [n,p,success,failure].into_iter().all(exact_scalar),
        Law::Normal { mean,sigma } => exact_scalar(mean) && exact_scalar(sigma),
        Law::Poisson { rate } => exact_scalar(rate),
        Law::Finite(table) => table.iter().all(|(v,p)| exact_scalar(v) && exact_scalar(p)),
    }
}

/// A conservative expansion budget, before conversion can expand powers or
/// allocate a dense multivariate grid. Declining only selects the existing
/// general expectation path; it does not limit the mathematical query.
fn expansion_budget(e: &Expr, atoms: &[Expr], seen: &mut rustc_hash::FxHashMap<Expr,Option<usize>>) -> Option<usize> {
    use crate::Kind;
    if let Some(bound)=seen.get(e) { return *bound; }
    let power = |base:&Expr, exponent:&Expr, seen:&mut rustc_hash::FxHashMap<Expr,Option<usize>>| {
        let n=exponent.as_number().filter(|n|n.is_integer())?.to_i64()?;
        if n.unsigned_abs()>256 {return None;}
        let terms=expansion_budget(base,atoms,seen)?;
        if n<0 {
            return (!atoms.iter().any(|a|crate::traverse::has(base,a))).then_some(1);
        }
        terms.checked_pow(n as u32).filter(|n|*n<=4096)
    };
    let result=match e.kind() {
        Kind::Num(_) | Kind::Sym(_)=>Some(1),
        Kind::Add(c,terms)=>{
            let mut bound=usize::from(!c.is_zero());
            for (term,_) in terms {
                bound=bound.checked_add(expansion_budget(term,atoms,seen)?)?;
                if bound>4096 {return None;}
            }
            Some(bound)
        },
        Kind::Mul(_,terms)=>{
            let mut bound=1usize;
            for (base,exponent) in terms {
                bound=bound.checked_mul(power(base,exponent,seen)?)?;
                if bound>4096 {return None;}
            }
            Some(bound)
        },
        Kind::Pow(base,exponent)=>power(base,exponent,seen),
        _=>None,
    };
    seen.insert(e.clone(),result);result
}

fn raw_moment<E>(law: &Law, order: usize, check: &mut impl FnMut() -> Result<(), E>) -> Result<Expr, E> {
    match law {
        Law::Binomial { n,p,success,failure } => binomial_moment(n,p,success,failure,order,check),
        Law::Normal { mean,sigma } => normal_moment(mean,sigma,order,check),
        Law::Poisson { rate } => poisson_moment(rate,order,check),
        Law::Finite(table) => {
            let mut terms = Vec::with_capacity(table.len());
            for (v,p) in table { check()?; terms.push(ops::mul2(&ops::powi(v,order as i64),p)); }
            Ok(ops::add(&terms))
        }
    }
}

/// Integrate every monomial over independent scalar laws. The polynomial
/// conversion, moment products, accumulation and simplification are native.
/// Floating, host-defined or nonpolynomial cases retain established routes.
pub fn polynomial_query<E>(expression: &Expr, atoms: &[Expr], laws: &[Law], variance: bool,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Option<Expr>, E> {
    if atoms.is_empty() || atoms.len()!=laws.len() || !laws.iter().all(law_exact)
        || expression.float_info().max_precision.is_some() || expression.has_nonfinite() {
        return Ok(None);
    }
    // Distribution descriptors assert independence. A rate containing one
    // of the random generators is a hierarchical model, not an independent law.
    if laws.iter().any(|law|matches!(law,Law::Poisson{rate} if atoms.iter().any(|atom|crate::traverse::has(rate,atom)))) {
        return Ok(None);
    }
    let expectation = |e: &Expr, check: &mut dyn FnMut() -> Result<(), E>| -> Result<Option<Expr>, E> {
        if expansion_budget(e,atoms,&mut rustc_hash::FxHashMap::default()).is_none() { return Ok(None); }
        let Ok(polynomial) = dict_from_expr(e, atoms) else { return Ok(None); };
        let mut terms = Vec::new();
        let mut moments: std::collections::HashMap<(usize,u32),Expr> = std::collections::HashMap::new();
        for (powers, coefficient) in polynomial {
            check()?;
            if !exact_scalar(&coefficient) { return Ok(None); }
            let mut term = coefficient;
            for (index, degree) in powers.into_iter().enumerate() {
                if degree == 0 && !matches!(laws[index],Law::Finite(_)) { continue; }
                // Optimization budget only: the existing general computation
                // remains available at every unsupported order.
                if degree > 256 { return Ok(None); }
                let key = (index,degree);
                let moment = if let Some(value)=moments.get(&key) { value.clone() } else {
                    let value = raw_moment(&laws[index],degree as usize,&mut ||check())?;
                    moments.insert(key,value.clone()); value
                };
                term = ops::mul2(&term,&moment);
            }
            terms.push(term);
        }
        Ok(Some(crate::simplify::simplify(&ops::add(&terms))))
    };
    let Some(mean)=expectation(expression,check)? else { return Ok(None); };
    if !variance { return Ok(Some(mean)); }
    let Some(second)=expectation(&ops::powi(expression,2),check)? else { return Ok(None); };
    check()?;
    Ok(Some(crate::simplify::simplify(&ops::sub(&second,&ops::powi(&mean,2)))))
}


#[cfg(test)]
mod poisson_tests {
    use super::*;
    #[test]
    fn touchard_recurrence_and_degenerate_mass_are_exact() {
        let rate=Expr::symbol("lambda");
        let mut check=||Ok::<(),()>(());
        for order in 0..=8 {
            let value=poisson_moment(&rate,order,&mut check).unwrap();
            let next=poisson_moment(&rate,order+1,&mut check).unwrap();
            let recurrence=ops::mul2(&rate,&ops::add2(&value,&crate::diff::diff(&value,&rate)));
            assert!(crate::expand::expand(&ops::sub(&next,&recurrence)).is_zero());
            assert_eq!(poisson_moment(&Expr::zero(),order,&mut check).unwrap(),if order==0{Expr::one()}else{Expr::zero()});
        }
    }
    #[test]
    fn moment_interruption_and_dependent_rates_do_not_return_values() {
        let mut calls=0;
        let result=poisson_moment(&Expr::int(2),8,&mut ||{calls+=1;if calls==3{Err("stop")}else{Ok(())}});
        assert_eq!(result,Err("stop"));
        let x=Expr::symbol("X");
        let law=Law::Poisson{rate:x.clone()};
        assert!(polynomial_query(&x,&[x.clone()],&[law],false,&mut ||Ok::<(),()>(())).unwrap().is_none());
    }
}


#[cfg(test)]
mod joint_normal_tests {
    use super::*;
    #[test]
    fn correlated_pairing_and_deep_univariate_moments_are_exact() {
        let zero=Expr::zero();let one=Expr::one();
        let result=normal_raw(&[2,2],&[zero.clone(),zero.clone()],&[Expr::int(2),one.clone(),one.clone(),Expr::int(3)],&mut ||Ok::<(),()>(())).unwrap().unwrap().0;
        assert_eq!(result,Expr::int(8)); // C00*C11 + 2*C01^2
        let moment=normal_raw(&[1200],&[zero],&[one],&mut ||Ok::<(),()>(())).unwrap().unwrap().0;
        let expected=(1..1200).step_by(2).fold(BigInt::from(1),|product,n|product*n);
        assert_eq!(moment,Expr::num(Number::from_bigint(expected)));
    }
    #[test]
    fn unsupported_states_and_interruptions_remain_distinct_from_zero() {
        assert!(!normal_raw_budget(&[1000,1000]));
        assert!(normal_raw(&[1],&[],&[],&mut ||Ok::<(),()>(())).unwrap().is_none());
        let mut checks=0;
        let result=normal_raw(&[8],&[Expr::zero()],&[Expr::one()],&mut ||{checks+=1;if checks==3{Err("stop")}else{Ok(())}});
        assert_eq!(result,Err("stop"));
        assert_eq!(normal_raw(&[0],&[],&[],&mut ||Ok::<(),()>(())).unwrap(),Some((Expr::one(),true)));
    }
}
