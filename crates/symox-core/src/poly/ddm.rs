//! Exact dense domain-matrix products, independent of Python representations.
use super::Coeff;

/// Accumulate A += B*C in row-major order. The caller proves disjoint mutable
/// rows and immutable scalar inputs. Each completed entry is published in the
/// same order as the general implementation; interruptions preserve that
/// prefix rather than being converted to a successful refusal.
/// False means a nonstandard shape, before any mutation or publication.
pub fn imatmul<C:Coeff,E>(a:&mut [Vec<C>],b:&[Vec<C>],c:&[Vec<C>],
    check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(usize,usize,&C)->Result<(),E>)->Result<bool,E> {
    check()?;
    let Some(first)=a.first() else{return Ok(false)};
    let n=first.len();let k=c.len();
    if n==0||k==0||a.len()!=b.len()||a.iter().any(|row|row.len()!=n)
        ||b.iter().any(|row|row.len()!=k)||c.iter().any(|row|row.len()!=n) {
        return Ok(false);
    }
    for (i,(ai,bi)) in a.iter_mut().zip(b).enumerate() {
        for (j,aij) in ai.iter_mut().enumerate() {
            check()?;
            let mut value=aij.clone();
            for (k,bik) in bi.iter().enumerate() {
                if k%64==0 {check()?;}
                value=value.add(&bik.mul(&c[k][j]));
            }
            publish(i,j,&value)?;
            *aij=value;
        }
    }
    check()?;
    Ok(true)
}

/// Publication events preserve both outer row permutations and the order of
/// individual cell assignments of the exact Gauss-Jordan implementation.
pub enum RrefEvent<'a> { Swap(usize,usize), Set(usize,usize,&'a super::Rat) }

/// The direct QQ reduction, without promoting integer scalar representations
/// or substituting fraction-free/modular mutation sequences. Partial pivoting
/// compares absolute rational values exactly and keeps the first maximum.
pub fn irref<E>(a:&mut [Vec<super::Rat>],partial:bool,
    check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(RrefEvent<'_>)->Result<(),E>)->Result<Option<Vec<usize>>,E> {
    if !partial {
        return field_irref(a,false,check,&mut |event| publish(match event {
            FieldRrefEvent::Swap(i,j)=>RrefEvent::Swap(i,j),
            FieldRrefEvent::Set(i,j,value)=>RrefEvent::Set(i,j,value),
        }));
    }
    check()?;
    let m=a.len();
    let Some(first)=a.first() else{return Ok(Some(Vec::new()))};
    let n=first.len();
    if a.iter().any(|row|row.len()!=n) {return Ok(None);}
    let mut pivots=Vec::new();let mut i=0;
    for j in 0..n {
        check()?;
        if partial {
            let mut best=i;let mut magnitude=a[i][j].abs();
            for (ip,row) in a.iter().enumerate().skip(i+1) {
                check()?;
                let value=row[j].abs();
                if value>magnitude {best=ip;magnitude=value;}
            }
            publish(RrefEvent::Swap(i,best))?;a.swap(i,best);
        }
        if a[i][j].is_zero() {
            let mut replacement=None;
            for (ip,row) in a.iter().enumerate().skip(i+1) {
                check()?;
                if !row[j].is_zero() {replacement=Some(ip);break;}
            }
            let Some(ip)=replacement else{continue};
            publish(RrefEvent::Swap(i,ip))?;a.swap(i,ip);
        }
        let inverse=a[i][j].inv();
        for l in j..n {
            check()?;
            let value=a[i][l].mul(&inverse);
            publish(RrefEvent::Set(i,l,&value))?;a[i][l]=value;
        }
        for k in 0..m {
            check()?;
            if k==i||a[k][j].is_zero() {continue;}
            let factor=a[k][j].clone();
            let zero=a[k][j].sub(&factor);
            publish(RrefEvent::Set(k,j,&zero))?;a[k][j]=zero;
            for l in (j+1)..n {
                check()?;
                let value=a[k][l].sub(&factor.mul(&a[i][l]));
                publish(RrefEvent::Set(k,l,&value))?;a[k][l]=value;
            }
        }
        pivots.push(j);i+=1;
        if i>=m {break;}
    }
    check()?;
    Ok(Some(pivots))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{Int,Rat};
    fn q(n:i64,d:i64)->Rat {Rat::new(Int::from(n),Int::from(d))}
    #[test]
    fn exact_accumulation_preserves_nonzero_target_and_rectangular_products() {
        let mut a=vec![vec![q(1,7),q(2,3)]];
        let b=vec![vec![q(2,3),q(-1,4),q(5,2)]];
        let c=vec![vec![q(1,2),q(1,3)],vec![q(3,2),q(2,5)],vec![q(4,3),q(-1,2)]];
        let mut published=Vec::new();
        assert!(imatmul(&mut a,&b,&c,&mut ||Ok::<_,()>(()),&mut |i,j,v|{published.push((i,j,v.clone()));Ok(())}).unwrap());
        assert_eq!(a,vec![vec![q(577,168),q(-83,180)]]);
        assert_eq!(published,vec![(0,0,a[0][0].clone()),(0,1,a[0][1].clone())]);
    }
    #[test]
    fn invalid_shape_is_untouched_and_interruptions_keep_completed_prefix() {
        let mut a=vec![vec![Int::from(3),Int::from(4)]];
        let original=a.clone();
        assert!(!imatmul(&mut a,&[vec![Int::one()]],&[],&mut ||Ok::<_,()>(()),&mut |_,_,_|panic!("publication")).unwrap());
        assert_eq!(a,original);
        let mut publications=0;
        let result=imatmul(&mut a,&[vec![Int::one()]],&[vec![Int::from(2),Int::from(5)]],&mut ||Ok::<_,&str>(()),&mut |_,_,_|{publications+=1;if publications==2 {Err("stop")}else{Ok(())}});
        assert_eq!(result,Err("stop"));
        assert_eq!(a,vec![vec![Int::from(5),Int::from(4)]]);
    }
    #[test]
    fn rref_exact_pivots_swaps_and_sub_float_magnitude_comparisons() {
        let mut a=vec![vec![q(0,1),q(2,1)],vec![q(1,1),q(3,1)]];
        let mut first=None;
        let pivots=irref(&mut a,false,&mut ||Ok::<_,()>(()),&mut |event| {
            if let RrefEvent::Swap(i,j)=event {if first.is_none(){first=Some((i,j));}}
            Ok(())
        }).unwrap().unwrap();
        assert_eq!(pivots,vec![0,1]);assert_eq!(first,Some((0,1)));
        assert_eq!(a,vec![vec![q(1,1),q(0,1)],vec![q(0,1),q(1,1)]]);
        let tiny=Rat::new(Int::one(),Int::from(num_bigint::BigInt::from(10u8).pow(350)));
        for (lower,expected) in [(q(1,1).add(&tiny),(0,1)),(q(-1,1),(0,0))] {
            let mut a=vec![vec![q(1,1),q(0,1)],vec![lower,q(1,1)]];let mut first=None;
            irref(&mut a,true,&mut ||Ok::<_,()>(()),&mut |event| {
                if let RrefEvent::Swap(i,j)=event {if first.is_none(){first=Some((i,j));}}
                Ok(())
            }).unwrap().unwrap();
            assert_eq!(first,Some(expected));
        }
    }
    #[test]
    fn rref_empty_invalid_shapes_and_interrupted_publication_prefix() {
        assert_eq!(irref(&mut [],false,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap(),Some(vec![]));
        let mut invalid=vec![vec![q(1,1)],vec![q(1,1),q(2,1)]];let before=invalid.clone();
        assert_eq!(irref(&mut invalid,false,&mut ||Ok::<_,()>(()),&mut |_|panic!("publication")).unwrap(),None);
        assert_eq!(invalid,before);
        let mut a=vec![vec![q(2,1),q(4,1)],vec![q(1,1),q(3,1)]];let mut checks=0;
        let result=irref(&mut a,false,&mut ||{checks+=1;if checks==4{Err("stop")}else{Ok(())}},&mut |_|Ok(()));
        assert_eq!(result,Err("stop"));
        assert_eq!(a,vec![vec![q(1,1),q(4,1)],vec![q(1,1),q(3,1)]]);
        let mut a=vec![vec![q(0,1),q(2,1)],vec![q(1,1),q(3,1)]];let mut checks=0;
        let result=irref(&mut a,false,&mut ||{checks+=1;if checks==4{Err("stop")}else{Ok(())}},&mut |_|Ok(()));
        assert_eq!(result,Err("stop"));
        assert_eq!(a,vec![vec![q(1,1),q(3,1)],vec![q(0,1),q(2,1)]]);
    }

}

/// Observable mutations of the DDM Bareiss determinant. Lower triangular
/// entries remain unchanged; the pivot search swaps whole existing rows.
pub enum DetEvent<'a, C> { Swap(usize, usize), Set(usize, usize, &'a C) }

/// The in-place DDM Bareiss algorithm, rather than the compacting Matrix
/// determinant. Returning None means a non-square shape, before any mutation.
/// Exact quotient failures and interruptions retain the published prefix.
pub fn idet<C: Coeff, E>(a: &mut [Vec<C>],
    check: &mut impl FnMut() -> Result<(), E>,
    publish: &mut impl FnMut(DetEvent<'_, C>) -> Result<(), E>,
) -> Result<Option<C>, super::fraction_free::Failure<E, C>> {
    use super::fraction_free::Failure;
    check().map_err(Failure::Interrupted)?;
    let n = a.len();
    if a.iter().any(|row| row.len() != n) { return Ok(None); }
    if n == 0 { return Ok(Some(C::one())); }
    let mut sign = C::one();
    for k in 0..n-1 {
        check().map_err(Failure::Interrupted)?;
        if a[k][k].is_zero() {
            let mut replacement = None;
            for (i, row) in a.iter().enumerate().skip(k+1) {
                check().map_err(Failure::Interrupted)?;
                if !row[k].is_zero() { replacement = Some(i); break; }
            }
            let Some(i) = replacement else { return Ok(Some(C::zero())); };
            publish(DetEvent::Swap(k,i)).map_err(Failure::Interrupted)?;
            a.swap(k,i);
            sign = sign.neg();
        }
        let divisor = if k == 0 { C::one() } else { a[k-1][k-1].clone() };
        for i in k+1..n {
            for j in k+1..n {
                check().map_err(Failure::Interrupted)?;
                let numerator = a[i][j].mul(&a[k][k]).sub(&a[i][k].mul(&a[k][j]));
                let value = numerator.div_exact(&divisor)
                    .ok_or_else(|| Failure::NotExact(numerator, divisor.clone()))?;
                publish(DetEvent::Set(i,j,&value)).map_err(Failure::Interrupted)?;
                a[i][j] = value;
            }
        }
    }
    check().map_err(Failure::Interrupted)?;
    Ok(Some(sign.mul(&a[n-1][n-1])))
}

#[cfg(test)]
mod determinant_tests {
    use super::*;
    use crate::poly::Int;
    #[test]
    fn determinant_keeps_lower_triangle_and_swap_order() {
        let mut a = vec![vec![Int::from(0),Int::from(1),Int::from(2)],
                         vec![Int::from(2),Int::from(3),Int::from(5)],
                         vec![Int::from(4),Int::from(6),Int::from(9)]];
        let mut swaps = Vec::new(); let mut writes = Vec::new();
        let d = idet(&mut a, &mut || Ok::<_,()>(()), &mut |event| {
            match event { DetEvent::Swap(i,j) => swaps.push((i,j)),
                          DetEvent::Set(i,j,value) => writes.push((i,j,value.clone())) }
            Ok(())
        }).unwrap().unwrap();
        assert_eq!(d,Int::from(2));
        assert_eq!(swaps,vec![(0,1)]);
        assert_eq!(a[1][0],Int::from(0)); assert_eq!(a[2][0],Int::from(4));
        assert_eq!(writes.iter().map(|&(i,j,_)|(i,j)).collect::<Vec<_>>(), vec![(1,1),(1,2),(2,1),(2,2),(2,2)]);
    }
    #[test]
    fn determinant_empty_shapes_and_interruptions() {
        assert_eq!(idet::<Int,()>(&mut [], &mut ||Ok(()), &mut |_|Ok(())).unwrap(),Some(Int::one()));
        let mut a=vec![vec![Int::one(),Int::one()]];
        assert_eq!(idet(&mut a,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap(),None);
        let mut a=vec![vec![Int::one(),Int::one()],vec![Int::one(),Int::one()]];
        let result=idet(&mut a,&mut ||Ok::<_,u8>(()),&mut |_|Err(9));
        assert!(matches!(result,Err(crate::poly::fraction_free::Failure::Interrupted(9))));
        assert_eq!(a[1][1],Int::one());
    }
}

#[derive(Debug)]
pub enum InexactRrefError<E> { Interrupted(E), Arithmetic(super::convert::PolyError) }

/// The augmented Gauss-Jordan loop used by solve_lin_sys over an inexact
/// coefficient domain. DomainElement operations are the existing scalar
/// operations: first nonzero pivot, one/pivot followed by multiply, and
/// a-factor*b for every column (including the pivot and earlier columns).
/// Never replace these operations by exact rationalization, fused updates or
/// forced zero/one assignments; all change the RR/CC rounded computation.
pub fn inexact_linear_rref<E>(rows: &mut [Vec<super::domains::DomainElement>],
    domain: &super::domains::ExactDomain, width: usize,
    check: &mut impl FnMut() -> Result<(),E>) -> Result<Option<Vec<usize>>,InexactRrefError<E>> {
    use InexactRrefError::{Interrupted,Arithmetic};
    check().map_err(Interrupted)?;
    if domain.inexact_precision().is_none() && domain.inexact_ground().is_none() { return Ok(None); }
    if width==0 || rows.iter().any(|row|row.len()!=width || row.iter().any(|v|v.domain!=*domain)) { return Ok(None); }
    if rows.iter().flatten().any(|v|!super::expression_domain::native_coefficient_payload(&v.as_expr())) { return Ok(None); }
    let mut pivots=Vec::new(); let mut r=0;
    if rows.is_empty() { return Ok(Some(pivots)); }
    for col in 0..width {
        check().map_err(Interrupted)?;
        let mut pivot=None;
        for (i,row) in rows.iter().enumerate().skip(r) {
            check().map_err(Interrupted)?;
            if !row[col].is_zero() { pivot=Some(i);break; }
        }
        let Some(pivot)=pivot else { continue; };
        rows.swap(r,pivot);
        let reciprocal=domain.convert(&crate::Expr::one()).map_err(Arithmetic)?
            .binary(&rows[r][col],'/').map_err(Arithmetic)?;
        let mut normalized=Vec::with_capacity(width);
        for value in &rows[r] {
            check().map_err(Interrupted)?;
            normalized.push(value.binary(&reciprocal,'*').map_err(Arithmetic)?);
        }
        rows[r]=normalized;
        for i in 0..rows.len() {
            check().map_err(Interrupted)?;
            if i==r || rows[i][col].is_zero() { continue; }
            let factor=rows[i][col].clone();
            let mut reduced=Vec::with_capacity(width);
            for (a,b) in rows[i].iter().zip(&rows[r]) {
                check().map_err(Interrupted)?;
                let product=factor.binary(b,'*').map_err(Arithmetic)?;
                reduced.push(a.binary(&product,'-').map_err(Arithmetic)?);
            }
            rows[i]=reduced;
        }
        pivots.push(col);r+=1;
        if r==rows.len() { break; }
    }
    check().map_err(Interrupted)?;
    Ok(Some(pivots))
}

#[cfg(test)]
mod inexact_linear_tests {
    use super::*;
    use crate::poly::domains::ExactDomain;
    use crate::Expr;
    use std::sync::Arc;
    #[test]
    fn first_nonzero_gauss_jordan_and_composite_domain() {
        for domain in [ExactDomain::RR,ExactDomain::Inexact{precision:140,complex:true},
            ExactDomain::Composite{ground:Arc::new(ExactDomain::Inexact{precision:140,complex:false}),ring:false,
                                   parameters:vec![Expr::symbol("t")],order:"lex".into()}] {
            let mut rows=[vec![0,2,4],vec![1,3,5]].map(|row|row.into_iter().map(|n|domain.convert(&Expr::int(n)).unwrap()).collect::<Vec<_>>());
            let pivots=inexact_linear_rref(&mut rows,&domain,3,&mut ||Ok::<_,()>(())).unwrap().unwrap();
            assert_eq!(pivots,vec![0,1]);
            let expected=[[1,0,-1],[0,1,2]];
            for (row,expected) in rows.iter().zip(expected) {
                for (value,expected) in row.iter().zip(expected) {
                    assert_eq!(value.domain,domain);
                    assert!(value.binary(&domain.convert(&Expr::int(expected)).unwrap(),'-').unwrap().is_zero());
                }
            }
        }
    }
    #[test]
    fn refusal_and_interruptions_are_not_inconsistency() {
        let domain=ExactDomain::RR;
        let mut rows=vec![vec![domain.convert(&Expr::one()).unwrap()]];
        assert!(matches!(inexact_linear_rref(&mut rows,&domain,1,&mut ||Err(19)),Err(InexactRrefError::Interrupted(19))));
        assert!(inexact_linear_rref(&mut rows,&ExactDomain::QQ,1,&mut ||Ok::<_,()>(())).unwrap().is_none());
        assert_eq!(inexact_linear_rref(&mut [],&domain,usize::MAX,&mut ||Ok::<_,()>(())).unwrap(),Some(vec![]));
        assert_eq!(inexact_linear_rref(&mut rows,&domain,1,&mut ||Ok::<_,()>(())).unwrap(),Some(vec![0]));
    }
}


/// Exact field row reduction shared by DDM first-pivot elimination and
/// PolyMatrix's ground-field route. `full_rows` retains PolyMatrix's e/p and
/// a-c*b on every column; the DDM mode retains its published pivot-zero and
/// left-column omissions. Partial pivoting remains in the existing Rat path.
pub enum FieldRrefEvent<'a,C> { Swap(usize,usize), Set(usize,usize,&'a C) }

pub fn field_irref<C:super::int::Field,E>(a:&mut [Vec<C>],full_rows:bool,
    check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(FieldRrefEvent<'_,C>)->Result<(),E>) -> Result<Option<Vec<usize>>,E> {
    check()?;
    let m=a.len();
    let Some(first)=a.first() else{return Ok(Some(Vec::new()))};
    let n=first.len();
    if a.iter().any(|row|row.len()!=n) {return Ok(None);}
    let mut pivots=Vec::new();let mut i=0;
    for j in 0..n {
        check()?;
        if a[i][j].is_zero() {
            let mut replacement=None;
            for (k,row) in a.iter().enumerate().skip(i+1) {
                check()?;
                if !row[j].is_zero() {replacement=Some(k);break;}
            }
            let Some(k)=replacement else{continue};
            publish(FieldRrefEvent::Swap(i,k))?;a.swap(i,k);
        }
        let pivot=a[i][j].clone();
        let inverse=if full_rows {None}else{Some(pivot.inv())};
        for (l,cell) in a[i].iter_mut().enumerate().skip(if full_rows {0}else{j}) {
            check()?;
            let value=if let Some(inverse)=&inverse {cell.mul(inverse)}else{cell.div(&pivot)};
            publish(FieldRrefEvent::Set(i,l,&value))?;*cell=value;
        }
        for k in 0..m {
            check()?;
            if k==i {continue;}
            let (row,pivot_row)=if k<i {
                let (left,right)=a.split_at_mut(i);(&mut left[k],&right[0])
            }else {
                let (left,right)=a.split_at_mut(k);(&mut right[0],&left[i])
            };
            if row[j].is_zero() {continue;}
            let factor=row[j].clone();
            if !full_rows {
                let zero=row[j].sub(&factor);
                publish(FieldRrefEvent::Set(k,j,&zero))?;row[j]=zero;
            }
            for (l,(cell,pivot_value)) in row.iter_mut().zip(pivot_row.iter()).enumerate().skip(if full_rows {0}else{j+1}) {
                check()?;
                let value=cell.sub(&factor.mul(pivot_value));
                publish(FieldRrefEvent::Set(k,l,&value))?;*cell=value;
            }
        }
        pivots.push(j);i+=1;
        if i>=m {break;}
    }
    check()?;
    Ok(Some(pivots))
}

/// Canonical free-column basis of an already reduced exact field matrix.
/// Each vector has its free entry one and pivot entries -R[i,free]. An empty
/// basis therefore means proven full column rank, never a declined reduction.
pub fn field_nullspace<C:super::int::Field,E>(a:&[Vec<C>],cols:usize,pivots:&[usize],
    check:&mut impl FnMut()->Result<(),E>) -> Result<Option<Vec<Vec<C>>>,E> {
    check()?;
    if a.iter().any(|row|row.len()!=cols)||pivots.len()>a.len()
        ||pivots.iter().any(|&p|p>=cols)||pivots.windows(2).any(|p|p[0]>=p[1]) {return Ok(None);}
    if pivots.len()==cols {return Ok(Some(Vec::new()));}
    // Representation limits, not a work quota: the caller retains its full
    // fallback when a single Rust vector cannot represent the actual output.
    if !matches!(cols.checked_mul(std::mem::size_of::<C>()),Some(n) if n<=isize::MAX as usize)
        ||!matches!((cols-pivots.len()).checked_mul(std::mem::size_of::<Vec<C>>()),Some(n) if n<=isize::MAX as usize) {return Ok(None);}
    let mut basis=Vec::new();let mut next=0;
    for free in 0..cols {
        check()?;
        if pivots.get(next)==Some(&free) {next+=1;continue;}
        let mut vector=vec![C::zero();cols];vector[free]=C::one();
        for (i,&pivot) in pivots.iter().enumerate() {
            check()?;
            if !a[i][free].is_zero() {vector[pivot]=a[i][free].neg();}
        }
        basis.push(vector);
    }
    check()?;
    Ok(Some(basis))
}

#[cfg(test)]
mod ground_field_tests {
    use super::*;
    use crate::poly::{Coeff,Int,Rat,domains::{ExactDomain,ElementValue},ratfunc::{RF,Scope}};
    use crate::Expr;
    fn q(n:i64)->Rat {Rat::from(Int::from(n))}
    #[test]
    fn exact_rational_function_rref_and_nullspace_certificate() {
        let _scope=Scope::new(1);let y=Expr::symbol("y");
        let k=ExactDomain::functions(false,false,vec![y.clone()],"lex").unwrap();
        let make=|e:&Expr|match k.convert(e).unwrap().value {ElementValue::Function(r)=>r,_=>panic!("RF")};
        let original=vec![vec![make(&y),make(&Expr::one()),RF::zero()],
                          vec![RF::zero(),make(&Expr::one()),make(&y)]];
        let mut rows=original.clone();
        let pivots=field_irref(&mut rows,true,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap().unwrap();
        assert_eq!(pivots,vec![0,1]);
        let basis=field_nullspace(&rows,3,&pivots,&mut ||Ok::<_,()>(())).unwrap().unwrap();
        assert_eq!(basis.len(),1);assert!(basis[0][2].is_one());
        for row in original {
            let residual=row.iter().zip(&basis[0]).fold(RF::zero(),|v,(a,b)|v.add(&a.mul(b)));
            assert!(residual.is_zero());
        }
    }
    #[test]
    fn publication_modes_keep_their_distinct_column_sequences() {
        let original=vec![vec![q(0),q(1)],vec![q(1),q(1)]];
        let mut ddm=original.clone();let mut de=Vec::new();
        irref(&mut ddm,false,&mut ||Ok::<_,()>(()),&mut |e|{match e {
            RrefEvent::Swap(i,j)=>de.push((i,j,true)),RrefEvent::Set(i,j,_)=>de.push((i,j,false))};Ok(())}).unwrap();
        assert_eq!(de,vec![(0,1,true),(0,0,false),(0,1,false),(1,1,false),(0,1,false)]);
        let mut full=original;let mut fe=Vec::new();
        field_irref(&mut full,true,&mut ||Ok::<_,()>(()),&mut |e|{match e {
            FieldRrefEvent::Swap(i,j)=>fe.push((i,j,true)),FieldRrefEvent::Set(i,j,_)=>fe.push((i,j,false))};Ok(())}).unwrap();
        assert_eq!(ddm,full);
        assert_eq!(fe,vec![(0,1,true),(0,0,false),(0,1,false),(1,0,false),(1,1,false),(0,0,false),(0,1,false)]);
    }
    #[test]
    fn empty_rectangular_matrices_and_invalid_storage_are_distinct() {
        let basis=field_nullspace::<Rat,()>(&[],3,&[],&mut ||Ok(())).unwrap().unwrap();
        assert_eq!(basis,vec![vec![q(1),q(0),q(0)],vec![q(0),q(1),q(0)],vec![q(0),q(0),q(1)]]);
        assert!(field_nullspace(&[vec![q(1)]],2,&[],&mut ||Ok::<_,()>(())).unwrap().is_none());
        assert!(field_nullspace(&[vec![q(1)]],1,&[1],&mut ||Ok::<_,()>(())).unwrap().is_none());
        let mut rows=vec![vec![q(1)],vec![q(1),q(2)]];let original=rows.clone();
        assert!(field_irref(&mut rows,true,&mut ||Ok::<_,()>(()),&mut |_|panic!("publication")).unwrap().is_none());
        assert_eq!(rows,original);
        assert_eq!(field_nullspace(&[vec![q(1)]],1,&[0],&mut ||Ok::<_,()>(())).unwrap(),Some(vec![]));
        assert!(field_nullspace::<Rat,()>(&[],usize::MAX,&[],&mut ||Ok(())).unwrap().is_none());
    }
    #[test]
    fn interrupts_are_errors_not_successful_empty_spaces() {
        let mut rows=vec![vec![q(1),q(2)]];
        assert_eq!(field_irref(&mut rows,true,&mut ||Err("stop"),&mut |_|Ok(())),Err("stop"));
        let mut checks=0;
        assert_eq!(field_nullspace(&rows,2,&[],&mut ||{checks+=1;if checks==2{Err("stop")}else{Ok(())}}),Err("stop"));
    }
}
