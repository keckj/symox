//! Exact field forward/back substitution (QQ, GF(p), QQ<a>) in the
//! established column-major order.
use super::Field;
#[cfg(test)]
use super::Rat;

#[derive(Debug,PartialEq,Eq)]
pub enum Failure<E> { Interrupted(E), NonInvertible }

pub enum Event<'a,C=super::Rat> {
    Y(usize,usize,&'a C,bool), // true copies the original RHS object
    X(usize,usize,&'a C),
}

/// Shape refusal happens before any publication. The caller has performed the
/// real public prefix, including row permutations and scratch allocation.
pub fn solve<C:Field,E>(l:&[Vec<C>],u:&[Vec<C>],b:&[Vec<C>],columns:usize,
    check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(Event<'_,C>)->Result<(),E>)->Result<bool,Failure<E>> {
    check().map_err(Failure::Interrupted)?;
    let m=u.len();
    let Some(first)=u.first() else{return Ok(true)};
    let n=first.len();
    if m<n||l.len()!=m||b.len()!=m||u.iter().any(|row|row.len()!=n)
        ||l.iter().any(|row|row.len()!=m)||b.iter().any(|row|row.len()!=columns){return Ok(false);}
    // Every read cell was published earlier in the same column: no domain
    // zero is needed (GF(p)/QQ<a> zeros depend on the element's field).
    let mut y:Vec<Vec<Option<C>>>=vec![vec![None;columns];m];
    let mut x:Vec<Vec<Option<C>>>=vec![vec![None;columns];n];
    for k in 0..columns {
        for i in 0..m {
            check().map_err(Failure::Interrupted)?;
            let mut rhs=b[i][k].clone();
            for j in 0..i {
                check().map_err(Failure::Interrupted)?;
                rhs=rhs.sub(&l[i][j].mul(y[j][k].as_ref().unwrap()));
            }
            publish(Event::Y(i,k,&rhs,i==0)).map_err(Failure::Interrupted)?;
            y[i][k]=Some(rhs);
        }
    }
    if m>n {
        for row in y.iter().skip(n) {
            for value in row {
                check().map_err(Failure::Interrupted)?;
                if !value.as_ref().unwrap().is_zero(){return Err(Failure::NonInvertible);}
            }
        }
    }
    for k in 0..columns {
        for i in (0..n).rev() {
            check().map_err(Failure::Interrupted)?;
            if u[i][i].is_zero(){return Err(Failure::NonInvertible);}
            let mut rhs=y[i][k].clone().unwrap();
            for j in (i+1)..n {
                check().map_err(Failure::Interrupted)?;
                rhs=rhs.sub(&u[i][j].mul(x[j][k].as_ref().unwrap()));
            }
            let value=rhs.mul(&u[i][i].inv());
            publish(Event::X(i,k,&value)).map_err(Failure::Interrupted)?;
            x[i][k]=Some(value);
        }
    }
    check().map_err(Failure::Interrupted)?;
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Int;
    fn q(n:i64,d:i64)->Rat {Rat::new(Int::from(n),Int::from(d))}
    fn rows(values:&[&[i64]])->Vec<Vec<Rat>> {values.iter().map(|row|row.iter().map(|&v|q(v,1)).collect()).collect()}

    #[test]
    fn exact_multiple_rhs_preserves_forward_then_descending_publication() {
        let l=rows(&[&[1,0],&[2,1]]);let u=rows(&[&[3,4],&[0,5]]);let b=rows(&[&[7,8],&[9,10]]);
        let mut events=Vec::new();
        assert_eq!(solve(&l,&u,&b,2,&mut ||Ok::<_,()>(()),&mut |event| {
            events.push(match event {Event::Y(i,k,v,copy)=>(false,i,k,v.clone(),copy),Event::X(i,k,v)=>(true,i,k,v.clone(),false)});Ok(())
        }),Ok(true));
        assert_eq!(events,vec![(false,0,0,q(7,1),true),(false,1,0,q(-5,1),false),(false,0,1,q(8,1),true),(false,1,1,q(-6,1),false),
            (true,1,0,q(-1,1),false),(true,0,0,q(11,3),false),(true,1,1,q(-6,5),false),(true,0,1,q(64,15),false)]);
    }

    #[test]
    fn tall_consistency_errors_follow_complete_forward_publication() {
        let l=rows(&[&[1,0,0],&[2,1,0],&[3,4,1]]);let u=rows(&[&[2,1],&[0,3],&[0,0]]);
        let mut output=Vec::new();
        assert_eq!(solve(&l,&u,&rows(&[&[4],&[14],&[36]]),1,&mut ||Ok::<_,()>(()),&mut |e| {if let Event::X(i,_,v)=e {output.push((i,v.clone()));}Ok(())}),Ok(true));
        assert_eq!(output,vec![(1,q(2,1)),(0,q(1,1))]);
        let mut count=0;
        assert_eq!(solve(&l,&u,&rows(&[&[4],&[14],&[37]]),1,&mut ||Ok::<_,()>(()),&mut |e| {assert!(matches!(e,Event::Y(..)));count+=1;Ok(())}),Err(Failure::NonInvertible));
        assert_eq!(count,3);
    }

    #[test]
    fn zero_pivot_keeps_completed_x_prefix_and_empty_rhs_does_not_probe_pivots() {
        let l=rows(&[&[1,0],&[0,1]]);let u=rows(&[&[0,1],&[0,2]]);let mut output=Vec::new();
        assert_eq!(solve(&l,&u,&rows(&[&[1,2],&[4,6]]),2,&mut ||Ok::<_,()>(()),&mut |e| {if let Event::X(i,k,v)=e {output.push((i,k,v.clone()));}Ok(())}),Err(Failure::NonInvertible));
        assert_eq!(output,vec![(1,0,q(2,1))]);
        assert_eq!(solve(&l,&u,&[vec![],vec![]],0,&mut ||Ok::<_,()>(()),&mut |_|panic!("empty RHS publication")),Ok(true));
    }

    #[test]
    fn shape_refusal_and_interrupted_publication_do_not_report_success() {
        let l=rows(&[&[1,0],&[0,1]]);let u=rows(&[&[2,1],&[0,3]]);
        assert_eq!(solve(&l,&u,&[vec![q(1,1)]],1,&mut ||Ok::<_,()>(()),&mut |_|panic!("shape publication")),Ok(false));
        let mut completed=0;
        let result=solve(&l,&u,&rows(&[&[1],&[2]]),1,&mut ||Ok(()),&mut |_| {if completed==2 {Err("interrupt")}else{completed+=1;Ok(())}});
        assert_eq!(result,Err(Failure::Interrupted("interrupt")));assert_eq!(completed,2);
    }
}
