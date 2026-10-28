//! Complete rank-revealing fraction-free PLDU, including skipped columns.
use super::Coeff;
pub use super::fraction_free::Failure;

pub enum Eliminate<'a,C> { Swap(usize,usize),Set(usize,usize,&'a C),Retain(usize,usize) }

/// Packed fraction-free LU. Shape refusal precedes every publication; an
/// interruption keeps the successfully published prefix. Every column is
/// considered and exact division uses the preceding actual pivot minor.
pub fn eliminate<C:Coeff,E>(a:&mut [Vec<C>],n:usize,check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(Eliminate<'_,C>)->Result<(),E>)->Result<Option<Vec<usize>>,Failure<E,C>> {
    check().map_err(Failure::Interrupted)?;
    if a.iter().any(|row|row.len()!=n){return Ok(None);}
    let m=a.len();let mut rank=0;let mut previous=None;let mut columns=Vec::new();
    for j in 0..n {
        check().map_err(Failure::Interrupted)?;
        if rank==m{break;}
        // Preserve the established full-column zero observation before pivot
        // search; only the scalar arithmetic, not Python callbacks, is here.
        if a.iter().all(|row|row[j].is_zero()){continue;}
        let mut replacement=None;
        for (i,row) in a.iter().enumerate().skip(rank){
            check().map_err(Failure::Interrupted)?;
            if !row[j].is_zero(){replacement=Some(i);break;}
        }
        let Some(pivot_row)=replacement else{continue};
        if pivot_row!=rank{
            publish(Eliminate::Swap(rank,pivot_row)).map_err(Failure::Interrupted)?;
            a.swap(rank,pivot_row);
        }
        let pivot=a[rank][j].clone();
        for i in (rank+1)..m {
            check().map_err(Failure::Interrupted)?;
            let multiplier=a[i][j].clone();
            for k in (j+1)..n {
                check().map_err(Failure::Interrupted)?;
                let numerator=pivot.mul(&a[i][k]).sub(&a[rank][k].mul(&multiplier));
                let value=match previous.as_ref(){
                    Some(divisor)=>numerator.div_exact(divisor).ok_or_else(||Failure::NotExact(numerator,divisor.clone()))?,
                    None=>numerator,
                };
                publish(Eliminate::Set(i,k,&value)).map_err(Failure::Interrupted)?;a[i][k]=value;
            }
            // Python's packed-multiplier assignment retains the same object.
            publish(Eliminate::Retain(i,j)).map_err(Failure::Interrupted)?;
        }
        columns.push(j);rank+=1;previous=Some(pivot);
    }
    check().map_err(Failure::Interrupted)?;
    Ok(Some(columns))
}

/// Write P only after the caller's actual zeros allocation. Invalid input
/// declines before mutation so the complete Python behavior remains available.
pub fn permutation<E>(indices:&[usize],m:usize,check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(usize,usize)->Result<(),E>)->Result<bool,E> {
    check()?;
    if indices.len()!=m||indices.iter().any(|&j|j>=m){return Ok(false);}
    let mut seen=vec![false;m];
    for &j in indices{if seen[j]{return Ok(false);}seen[j]=true;}
    for (i,&j) in indices.iter().enumerate(){check()?;publish(i,j)?;}
    check()?;Ok(true)
}

pub enum Lower<'a,C> {
    Copy(usize,usize,usize,usize,&'a C),
    ZeroUpper(usize,usize),
    OneLower(usize,usize),
}

/// Extract L's pivots/multipliers and clear U in their original write order.
/// Source indices accompany Copy so bindings preserve actual object aliases.
pub fn lower<C:Coeff,E>(u:&mut [Vec<C>],l:&mut [Vec<C>],n:usize,
    check:&mut impl FnMut()->Result<(),E>,publish:&mut impl FnMut(Lower<'_,C>)->Result<(),E>)->Result<bool,E> {
    check()?;let m=u.len();
    if l.len()!=m||l.iter().any(|row|row.len()!=m)||u.iter().any(|row|row.len()!=n){return Ok(false);}
    let mut i=0;
    for j in 0..n {
        check()?;if i==m{break;}
        if !u[i][j].is_zero(){
            let value=u[i][j].clone();publish(Lower::Copy(i,i,i,j,&value))?;l[i][i]=value;
            for k in (i+1)..m{
                check()?;let value=u[k][j].clone();publish(Lower::Copy(k,i,k,j,&value))?;l[k][i]=value;
                publish(Lower::ZeroUpper(k,j))?;u[k][j]=C::zero();
            }
            i+=1;
        }
    }
    for (k,row) in l.iter_mut().enumerate().skip(i){check()?;publish(Lower::OneLower(k,k))?;row[k]=C::one();}
    check()?;Ok(true)
}

/// D[0][0]=L[0][0] and the unused initial di=K.one remain in the real Python
/// prefix. Only the subsequent coefficient products are computed here.
pub fn diagonal<C:Coeff,E>(l:&[Vec<C>],d:&mut [Vec<C>],check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(usize,&C)->Result<(),E>)->Result<bool,E> {
    check()?;let m=l.len();
    if d.len()!=m||d.iter().any(|row|row.len()!=m)||l.iter().any(|row|row.len()!=m){return Ok(false);}
    for i in 1..m {check()?;let value=l[i-1][i-1].mul(&l[i][i]);publish(i,&value)?;d[i][i]=value;}
    check()?;Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{Int,Rat};
    fn ints(rows:&[&[i64]])->Vec<Vec<Int>>{rows.iter().map(|row|row.iter().copied().map(Int::from).collect()).collect()}
    fn full<C:Coeff>(source:&[Vec<C>],n:usize)->(Vec<usize>,Vec<Vec<C>>,Vec<Vec<C>>,Vec<Vec<C>>){
        let m=source.len();let mut u=source.to_vec();let mut perm=(0..m).collect::<Vec<_>>();
        eliminate(&mut u,n,&mut ||Ok::<_,()>(()),&mut |event|{if let Eliminate::Swap(i,j)=event{perm.swap(i,j);}Ok(())}).unwrap();
        let mut l=vec![vec![C::zero();m];m];let mut d=l.clone();
        assert!(lower(&mut u,&mut l,n,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap());
        if m>0{d[0][0]=l[0][0].clone();}
        assert!(diagonal(&l,&mut d,&mut ||Ok::<_,()>(()),&mut |_,_|Ok(())).unwrap());
        for i in 0..m{for j in 0..n{
            let result=(0..m).fold(Rat::zero(),|acc,k|acc.add(&l[i][k].to_rat().mul(&u[k][j].to_rat()).div(&d[k][k].to_rat())));
            assert_eq!(result,source[perm[i]][j].to_rat());
        }}
        (perm,l,d,u)
    }
    #[test]
    fn late_pivots_and_skipped_denominator_satisfy_full_exact_factorization(){
        let a=ints(&[&[0,0,1,2],&[0,0,2,3]]);let (_,l,d,u)=full(&a,4);
        assert_eq!(l,ints(&[&[1,0],&[2,-1]]));assert_eq!(d,ints(&[&[1,0],&[0,-1]]));
        assert_eq!(u,ints(&[&[0,0,1,2],&[0,0,0,-1]]));
        let a=ints(&[&[0,1,2,3],&[0,2,3,4],&[0,3,4,6]]);let (_,l,d,_)=full(&a,4);
        assert_eq!(l,ints(&[&[1,0,0],&[2,-1,0],&[3,-2,-1]]));assert_eq!(d,ints(&[&[1,0,0],&[0,-1,0],&[0,0,1]]));
    }
    #[test]
    fn rational_empty_deficient_and_row_cycle_conventions(){
        let a=ints(&[&[0,0,2],&[3,0,1],&[4,5,6]]);assert_eq!(full(&a,3).0,vec![1,2,0]);
        for a in [ints(&[&[0,0,0,0],&[0,0,2,3]]),ints(&[&[1,2],&[2,4],&[0,0]]),vec![vec![];3]]{let n=a[0].len();full(&a,n);}
        full::<Int>(&[],17);
        let a=vec![vec![Rat::zero(),Rat::zero(),Rat::new(Int::from(1),Int::from(2)),Rat::new(Int::from(2),Int::from(3))],vec![Rat::zero(),Rat::zero(),Rat::new(Int::from(2),Int::from(5)),Rat::new(Int::from(3),Int::from(7))]];full(&a,4);
    }
    #[test]
    fn shape_refusal_and_interrupt_preserve_publication_prefix(){
        let mut a=ints(&[&[1],&[2,3]]);let original=a.clone();
        assert!(eliminate(&mut a,2,&mut ||Ok::<_,()>(()),&mut |_|panic!("ragged publication")).unwrap().is_none());assert_eq!(a,original);
        let mut a=ints(&[&[2,1,3],&[4,3,5]]);let mut writes=0;
        let result=eliminate(&mut a,3,&mut ||Ok(()),&mut |e|{if let Eliminate::Set(..)=e{writes+=1;if writes==2{return Err("interrupt");}}Ok(())});
        assert!(matches!(result,Err(Failure::Interrupted("interrupt"))));assert_eq!(a,ints(&[&[2,1,3],&[4,2,5]]));
    }
    #[test]
    fn permutation_and_lower_events_preserve_source_tokens(){
        let mut writes=Vec::new();assert!(permutation(&[1,2,0],3,&mut ||Ok::<_,()>(()),&mut |i,j|{writes.push((i,j));Ok(())}).unwrap());assert_eq!(writes,vec![(0,1),(1,2),(2,0)]);
        assert!(!permutation(&[0,0],2,&mut ||Ok::<_,()>(()),&mut |_,_|panic!("bad permutation")).unwrap());
        let mut u=ints(&[&[0,1,2],&[0,3,4]]);let mut l=ints(&[&[0,0],&[0,0]]);let mut events=Vec::new();
        lower(&mut u,&mut l,3,&mut ||Ok::<_,()>(()),&mut |e|{events.push(match e{Lower::Copy(i,j,r,c,_)=>('c',i,j,r,c),Lower::ZeroUpper(i,j)=>('z',i,j,0,0),Lower::OneLower(i,j)=>('o',i,j,0,0)});Ok(())}).unwrap();
        assert_eq!(events,vec![('c',0,0,0,1),('c',1,0,1,1),('z',1,1,0,0),('c',1,1,1,2)]);
    }
}
