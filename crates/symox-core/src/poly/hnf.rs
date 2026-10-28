//! Exact ordinary column Hermite reduction, in the established operation order.
use super::Int;

pub enum Event<'a> {
    Operation(usize,usize,&'a Int,&'a Int,&'a Int,&'a Int),
    Set(usize,usize,&'a Int),
}

fn columns<E>(matrix:&mut [Vec<Int>],i:usize,j:usize,a:&Int,b:&Int,c:&Int,d:&Int,
    check:&mut impl FnMut()->Result<(),E>,publish:&mut impl FnMut(Event<'_>)->Result<(),E>)->Result<(),E> {
    publish(Event::Operation(i,j,a,b,c,d))?;
    for (k,row) in matrix.iter_mut().enumerate(){
        check()?;
        let original=row[i].clone();
        let value=a.mul(&original).add(&b.mul(&row[j]));
        publish(Event::Set(k,i,&value))?;row[i]=value;
        // The second expression sees the first assignment when i==j.
        let value=c.mul(&original).add(&d.mul(&row[j]));
        publish(Event::Set(k,j,&value))?;row[j]=value;
    }
    Ok(())
}

/// Reduce in place and return the first retained column. None is a shape
/// refusal before publication, never a claim that the lattice has rank zero.
pub fn ordinary<E>(a:&mut [Vec<Int>],n:usize,check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(Event<'_>)->Result<(),E>)->Result<Option<usize>,E> {
    check()?;
    if a.iter().any(|row|row.len()!=n){return Ok(None);}
    let mut k=n;
    let zero=Int::zero();let one=Int::one();let minus_one=one.neg();
    for i in (0..a.len()).rev(){
        check()?;
        if k==0 {break;}
        k-=1;
        for j in (0..k).rev(){
            check()?;
            if !a[i][j].is_zero(){
                let left=a[i][k].clone();let right=a[i][j].clone();
                let (g,mut u,mut v)=left.gcdex(&right);
                if !left.is_zero()&&right.div_rem_trunc(&left).1.is_zero(){
                    v=zero.clone();u=if left.is_negative(){minus_one.clone()}else{one.clone()};
                }
                let r=left.div_exact_unchecked(&g);let s=right.div_exact_unchecked(&g);
                columns(a,k,j,&u,&v,&s.neg(),&r,check,publish)?;
            }
        }
        let mut pivot=a[i][k].clone();
        if pivot.is_negative(){
            columns(a,k,k,&minus_one,&zero,&minus_one,&zero,check,publish)?;
            pivot=pivot.neg();
        }
        if pivot.is_zero(){k+=1;}
        else {
            for j in (k+1)..n {
                check()?;
                let q=a[i][j].div_mod_floor(&pivot).0;
                columns(a,j,k,&one,&q.neg(),&zero,&one,check,publish)?;
            }
        }
    }
    check()?;
    Ok(Some(k))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ints(rows:&[&[i64]])->Vec<Vec<Int>> {rows.iter().map(|row|row.iter().copied().map(Int::from).collect()).collect()}
    #[test]
    fn canonical_example_and_exact_column_lattice_transformation() {
        let original=ints(&[&[12,6,4],&[3,9,6],&[2,16,14]]);let mut a=original.clone();
        let mut transform=ints(&[&[1,0,0],&[0,1,0],&[0,0,1]]);
        let k=ordinary(&mut a,3,&mut ||Ok::<_,()>(()),&mut |e|{
            if let Event::Operation(i,j,p,q,r,s)=e{
                for row in &mut transform{
                    let old=row[i].clone();row[i]=p.mul(&old).add(&q.mul(&row[j]));row[j]=r.mul(&old).add(&s.mul(&row[j]));
                }
            }Ok(())
        }).unwrap().unwrap();
        assert_eq!(k,0);assert_eq!(a,ints(&[&[10,0,2],&[0,15,3],&[0,0,2]]));
        for i in 0..3 {for j in 0..3 {
            let sum=(0..3).fold(Int::zero(),|v,l|v.add(&original[i][l].mul(&transform[l][j])));
            assert_eq!(sum,a[i][j]);
        }}
    }
    #[test]
    fn negative_pivots_and_floor_remainders_have_reference_representatives() {
        let mut a=ints(&[&[-2,-3]]);let mut writes=Vec::new();
        let k=ordinary(&mut a,2,&mut ||Ok::<_,()>(()),&mut |e|{if let Event::Set(i,j,v)=e{writes.push((i,j,v.clone()));}Ok(())}).unwrap();
        assert_eq!(k,Some(1));assert_eq!(a,ints(&[&[0,1]]));assert!(!writes.is_empty());
        let mut b=ints(&[&[-3]]);let mut writes=Vec::new();
        assert_eq!(ordinary(&mut b,1,&mut ||Ok::<_,()>(()),&mut |e|{if let Event::Set(_,_,v)=e{writes.push(v.clone());}Ok(())}).unwrap(),Some(0));
        assert_eq!(writes,vec![Int::from(3),Int::from(3)]);
    }
    #[test]
    fn deficient_rank_tall_shapes_and_empty_rows_preserve_column_crop() {
        let mut a=ints(&[&[1,2,3],&[2,4,6]]);
        let k=ordinary(&mut a,3,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap().unwrap();
        assert_eq!(k,2);assert_eq!(a,ints(&[&[0,0,1],&[0,0,2]]));
        let mut a=ints(&[&[0,0],&[0,0],&[0,0]]);
        assert_eq!(ordinary(&mut a,2,&mut ||Ok::<_,()>(()),&mut |_|panic!("zero publication")).unwrap(),Some(2));
        assert_eq!(ordinary(&mut [],7,&mut ||Ok::<_,()>(()),&mut |_|panic!("empty publication")).unwrap(),Some(7));
    }
    #[test]
    fn malformed_refusal_and_interruption_keep_the_published_prefix() {
        let mut a=ints(&[&[1],&[2,3]]);let saved=a.clone();
        assert_eq!(ordinary(&mut a,2,&mut ||Ok::<_,()>(()),&mut |_|panic!("shape publication")).unwrap(),None);assert_eq!(a,saved);
        let mut a=ints(&[&[-3]]);let mut writes=0;
        let result=ordinary(&mut a,1,&mut ||Ok(()),&mut |e|{if matches!(e,Event::Set(..)){writes+=1;if writes==2{return Err("interrupt");}}Ok(())});
        assert_eq!(result,Err("interrupt"));assert_eq!(a,ints(&[&[3]]));assert_eq!(writes,2);
    }
}
