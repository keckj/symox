//! Exact two-phase primal simplex. All pivots and comparisons use QQ;
//! unsupported scalar domains retain the caller's general implementation.
use crate::poly::Rat;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Relation { Le, Ge, Eq }

#[derive(Clone)]
pub struct Constraint { pub coefficients: Vec<Rat>, pub relation: Relation, pub rhs: Rat }

#[derive(Debug, PartialEq, Eq)]
pub enum Failure { Infeasible, Unbounded, InvalidShape, InvalidBounds }

#[derive(Clone)]
pub struct Bounds { pub lower: Option<Rat>, pub upper: Option<Rat> }

pub type Solution = (Rat, Vec<Rat>);

/// Prepare an arbitrary-bound rational program over nonnegative variables.
/// Columns retain the original variable order, including the positive then
/// negative columns of a free variable, so Bland's deterministic solution is
/// the same as the general caller's transformation.
pub fn linprog<E>(cost: &[Rat], constraints: &[Constraint], bounds: &[Bounds],
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Result<Solution, Failure>, E> {
    check()?;
    let n=cost.len();
    if bounds.len()!=n||constraints.iter().any(|row|row.coefficients.len()!=n) {
        return Ok(Err(Failure::InvalidShape));
    }
    let mut offsets=Vec::with_capacity(n);
    let mut mapping:Vec<Vec<(usize,Rat)>>=Vec::with_capacity(n);
    let mut upper=Vec::new();
    let mut columns=0;
    for (i,bound) in bounds.iter().enumerate() {
        check()?;
        if let (Some(lower),Some(higher))=(&bound.lower,&bound.upper) {
            if higher<lower {return Ok(Err(Failure::InvalidBounds));}
        }
        if let Some(lower)=&bound.lower {
            offsets.push(lower.clone());mapping.push(vec![(columns,Rat::one())]);columns+=1;
            if let Some(higher)=&bound.upper {upper.push((i,higher.sub(lower)));}
        }else if let Some(higher)=&bound.upper {
            offsets.push(higher.clone());mapping.push(vec![(columns,Rat::one().neg())]);columns+=1;
        }else {
            offsets.push(Rat::zero());mapping.push(vec![(columns,Rat::one()),(columns+1,Rat::one().neg())]);columns+=2;
        }
    }
    let mut transformed=Vec::with_capacity(constraints.len()+upper.len());
    for row in constraints {
        check()?;
        let mut coefficients=vec![Rat::zero();columns];
        let mut rhs=row.rhs.clone();
        for (i,a) in row.coefficients.iter().enumerate() {
            for (j,factor) in &mapping[i] {coefficients[*j]=a.mul(factor);}
            rhs=rhs.sub(&a.mul(&offsets[i]));
        }
        transformed.push(Constraint{coefficients,relation:row.relation,rhs});
    }
    for (i,rhs) in upper {
        check()?;
        let mut coefficients=vec![Rat::zero();columns];
        for (j,factor) in &mapping[i] {coefficients[*j]=factor.clone();}
        transformed.push(Constraint{coefficients,relation:Relation::Le,rhs});
    }
    let mut objective=vec![Rat::zero();columns];
    for (i,a) in cost.iter().enumerate() {
        for (j,factor) in &mapping[i] {objective[*j]=a.mul(factor);}
    }
    let (_,solution)=match solve(&objective,&transformed,check)? {
        Ok(solution)=>solution,Err(error)=>return Ok(Err(error)),
    };
    let mut original=Vec::with_capacity(n);
    let mut optimum=Rat::zero();
    for ((offset,terms),a) in offsets.into_iter().zip(mapping).zip(cost) {
        check()?;
        let value=terms.into_iter().fold(offset,|value,(j,factor)|value.add(&factor.mul(&solution[j])));
        optimum=optimum.add(&a.mul(&value));original.push(value);
    }
    check()?;
    Ok(Ok((optimum,original)))
}

fn pivot<E>(rows: &mut [Vec<Rat>], basis: &mut [usize], leaving: usize, entering: usize,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<(), E> {
    let divisor = rows[leaving][entering].clone();
    for entry in &mut rows[leaving] { *entry = entry.div(&divisor); }
    for i in 0..rows.len() {
        check()?;
        if i == leaving { continue; }
        let coefficient = rows[i][entering].clone();
        if coefficient.is_zero() { continue; }
        for j in 0..rows[i].len() {
            rows[i][j] = rows[i][j].sub(&coefficient.mul(&rows[leaving][j]));
        }
    }
    basis[leaving] = entering;
    Ok(())
}

fn optimize<E>(rows: &mut [Vec<Rat>], basis: &mut [usize], cost: &[Rat],
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Result<Solution, Failure>, E> {
    loop {
        check()?;
        let entering = (0..cost.len()).find(|&j| {
            if basis.contains(&j) { return false; }
            let reduced = rows.iter().zip(basis.iter()).fold(cost[j].clone(),
                |c, (row, &b)| c.sub(&cost[b].mul(&row[j])));
            reduced.is_negative()
        });
        let Some(entering) = entering else { break; };
        let leaving = (0..rows.len()).filter(|&i| rows[i][entering] > Rat::zero())
            .min_by(|&i, &j| {
                let ri = rows[i].last().unwrap().div(&rows[i][entering]);
                let rj = rows[j].last().unwrap().div(&rows[j][entering]);
                ri.cmp(&rj).then(basis[i].cmp(&basis[j]))
            });
        let Some(leaving) = leaving else { return Ok(Err(Failure::Unbounded)); };
        pivot(rows, basis, leaving, entering, check)?;
    }
    let mut solution = vec![Rat::zero(); cost.len()];
    for (row, &b) in rows.iter().zip(basis.iter()) { solution[b] = row.last().unwrap().clone(); }
    let optimum = cost.iter().zip(&solution).fold(Rat::zero(), |v, (c, x)| v.add(&c.mul(x)));
    Ok(Ok((optimum, solution)))
}

/// Minimize `cost*x` over nonnegative variables with exact linear constraints.
/// The callback propagates interruptions separately from mathematical failures.
pub fn solve<E>(cost: &[Rat], constraints: &[Constraint],
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Result<Solution, Failure>, E> {
    let n = cost.len();
    if constraints.iter().any(|c| c.coefficients.len() != n) { return Ok(Err(Failure::InvalidShape)); }
    let mut width = n;
    let mut basis = Vec::new();
    let mut artificial = Vec::new();
    let mut records = Vec::new();
    for constraint in constraints {
        check()?;
        let mut coefficients = constraint.coefficients.clone();
        let mut rhs = constraint.rhs.clone();
        let mut relation = constraint.relation;
        if rhs.is_negative() {
            rhs = rhs.neg();
            coefficients.iter_mut().for_each(|c| *c = c.neg());
            relation = match relation { Relation::Le => Relation::Ge, Relation::Ge => Relation::Le, Relation::Eq => Relation::Eq };
        }
        let mut extras = Vec::new();
        if relation == Relation::Le {
            extras.push((width, Rat::one()));
            basis.push(width); width += 1;
        } else {
            if relation == Relation::Ge { extras.push((width, Rat::one().neg())); width += 1; }
            extras.push((width, Rat::one()));
            basis.push(width); artificial.push(width); width += 1;
        }
        records.push((coefficients, extras, rhs));
    }
    let mut rows = records.into_iter().map(|(mut coefficients, extras, rhs)| {
        coefficients.resize(width, Rat::zero());
        for (j, value) in extras { coefficients[j] = value; }
        coefficients.push(rhs); coefficients
    }).collect::<Vec<_>>();
    let phase_one = (0..width).map(|j| if artificial.contains(&j) { Rat::one() } else { Rat::zero() }).collect::<Vec<_>>();
    let (optimum, _) = match optimize(&mut rows, &mut basis, &phase_one, check)? {
        Ok(result) => result,
        Err(error) => return Ok(Err(error)),
    };
    if !optimum.is_zero() { return Ok(Err(Failure::Infeasible)); }
    let mut redundant = Vec::new();
    for i in 0..basis.len() {
        check()?;
        if !artificial.contains(&basis[i]) { continue; }
        let entering = (0..width).find(|j| !artificial.contains(j) && !basis.contains(j) && !rows[i][*j].is_zero());
        if let Some(entering) = entering { pivot(&mut rows, &mut basis, i, entering, check)?; }
        else { redundant.push(i); }
    }
    for &i in redundant.iter().rev() { rows.remove(i); basis.remove(i); }
    let columns = (0..width).filter(|j| !artificial.contains(j)).collect::<Vec<_>>();
    let mut renumber = vec![0; width];
    for (i, &j) in columns.iter().enumerate() { renumber[j] = i; }
    rows = rows.into_iter().map(|row| {
        let mut out = columns.iter().map(|&j| row[j].clone()).collect::<Vec<_>>();
        out.push(row.last().unwrap().clone()); out
    }).collect();
    for b in &mut basis { *b = renumber[*b]; }
    let fullcost = columns.iter().map(|&j| if j < n { cost[j].clone() } else { Rat::zero() }).collect::<Vec<_>>();
    Ok(optimize(&mut rows, &mut basis, &fullcost, check)?.map(|(optimum, mut solution)| {
        solution.truncate(n); (optimum, solution)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn q(n: i64) -> Rat { Rat::from(crate::poly::Int::from(n)) }
    fn constraint(values: &[i64], relation: Relation, rhs: i64) -> Constraint {
        Constraint { coefficients: values.iter().map(|&v| q(v)).collect(), relation, rhs: q(rhs) }
    }
    #[test]
    fn two_phases_redundancy_and_unboundedness() {
        let constraints = [constraint(&[1,1],Relation::Ge,3), constraint(&[1,-1],Relation::Eq,1), constraint(&[2,-2],Relation::Eq,2)];
        let solution = solve(&[q(1),q(2)], &constraints, &mut || Ok::<_, ()>(())).unwrap().unwrap();
        assert_eq!(solution, (q(4),vec![q(2),q(1)]));
        assert_eq!(solve(&[q(-1)], &[], &mut || Ok::<_, ()>(())).unwrap(), Err(Failure::Unbounded));
        assert_eq!(solve(&[q(1)], &[constraint(&[1],Relation::Le,-1)], &mut || Ok::<_, ()>(())).unwrap(), Err(Failure::Infeasible));
    }
    #[test]
    fn interrupt_is_not_a_mathematical_failure() {
        assert_eq!(solve(&[q(1)], &[], &mut || Err("interrupted")), Err("interrupted"));
    }
    #[test]
    fn variable_transformations_preserve_primal_and_objective() {
        let bounds=[Bounds{lower:None,upper:None},Bounds{lower:Some(q(1)),upper:Some(q(3))}];
        let rows=[constraint(&[1,1],Relation::Le,6)];
        assert_eq!(linprog(&[q(-1),q(-2)],&rows,&bounds,&mut ||Ok::<_,()>(())).unwrap(),Ok((q(-9),vec![q(3),q(3)])));
        let upper_only=[Bounds{lower:None,upper:Some(q(-2))}];
        assert_eq!(linprog(&[q(-3)],&[],&upper_only,&mut ||Ok::<_,()>(())).unwrap(),Ok((q(6),vec![q(-2)])));
        assert_eq!(linprog(&[q(1)],&[],&upper_only,&mut ||Ok::<_,()>(())).unwrap(),Err(Failure::Unbounded));
    }
    #[test]
    fn bound_failures_empty_programs_and_interruptions_are_distinct() {
        let bounds=[Bounds{lower:Some(q(2)),upper:Some(q(1))}];
        assert_eq!(linprog(&[q(1)],&[],&bounds,&mut ||Ok::<_,()>(())).unwrap(),Err(Failure::InvalidBounds));
        assert_eq!(linprog(&[],&[],&[],&mut ||Ok::<_,()>(())).unwrap(),Ok((q(0),vec![])));
        assert_eq!(linprog(&[],&[constraint(&[],Relation::Eq,1)],&[],&mut ||Ok::<_,()>(())).unwrap(),Err(Failure::Infeasible));
        assert_eq!(linprog(&[q(1)],&[],&bounds,&mut ||Err("stop")),Err("stop"));
    }
}
