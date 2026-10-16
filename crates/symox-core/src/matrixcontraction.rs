//! Scalar, bilinear contractions of typed matrix vectors.
use std::sync::Arc;
use crate::expr::{Expr,Func,Kind};
use crate::matrix::{MatrixError,MResult};
use crate::matrixexpr::MatrixExpr;
use crate::ops::{add,sub};

fn vector(matrix:&MatrixExpr)->MResult<(Expr,bool)>{
    if matrix.rows().is_one(){Ok((matrix.cols().clone(),true))}
    else if matrix.cols().is_one(){Ok((matrix.rows().clone(),false))}
    else{Err(MatrixError::Type("DotProduct arguments must be row or column vectors".into()))}
}

pub fn validate(a:&MatrixExpr,b:&MatrixExpr)->MResult<Expr>{
    let(n,_)=vector(a)?;let(m,_)=vector(b)?;
    if n!=m&&crate::matrix::ratfunc::is_zero_assume(&sub(&n,&m))!=Some(true){return Err(MatrixError::Type("DotProduct arguments must have the same length".into()));}
    Ok(n)
}

/// Construct a held, ordered DotProduct or explicitly evaluate its scalar sum.
/// Complex vectors are paired bilinearly: callers request conjugation explicitly.
pub fn dot_product(a:&MatrixExpr,b:&MatrixExpr,evaluate:bool)->MResult<Expr>{
    let n=validate(a,b)?;
    if !evaluate{return Ok(Expr::from_kind(Kind::Func(Func::MatrixDotProduct(Arc::new(a.clone()),Arc::new(b.clone())),vec![])));}
    let(_,arow)=vector(a)?;let(_,brow)=vector(b)?;
    let entry=|matrix:&MatrixExpr,row:bool,i:&Expr|if row{matrix.element(&Expr::zero(),i)}else{matrix.element(i,&Expr::zero())};
    let product=|i:&Expr|->MResult<Expr>{
        let x=entry(a,arow,i)?;let y=entry(b,brow,i)?;
        // Use the same nonfinite policy as dense and entrywise products.
        Ok(MatrixExpr::entry_product(&x,&y))
    };
    if let Some(length)=n.as_i64(){
        let length=usize::try_from(length).map_err(|_|MatrixError::Value("Vector length must be nonnegative".into()))?;
        if length>100_000{return Err(MatrixError::NotImplemented("Explicit DotProduct evaluation exceeds 100000 terms".into()));}
        return Ok(add(&(0..length).map(|i|product(&Expr::int(i as i64))).collect::<MResult<Vec<_>>>()?));
    }
    let index=Expr::dummy("_i",crate::Assumptions::INTEGER);
    Ok(crate::summation::sum_unevaluated(&product(&index)?,&index,&Expr::zero(),&sub(&n,&Expr::one())))
}

#[cfg(test)]
mod tests{
    use super::*;
    use crate::matrix::Matrix;
    #[test]
    fn bilinear_orientations_empty_and_symbolic_sum(){
        let a=MatrixExpr::explicit(Matrix::from_rows(vec![vec![Expr::i(),Expr::int(2)]]).unwrap());
        let b=MatrixExpr::explicit(Matrix::from_rows(vec![vec![Expr::i(),Expr::int(3)]]).unwrap());
        for x in [a.clone(),a.transpose(true).unwrap()]{for y in [b.clone(),b.transpose(true).unwrap()]{assert_eq!(dot_product(&x,&y,true).unwrap(),Expr::int(5));}}
        let a=MatrixExpr::zero(Expr::one(),Expr::zero()).unwrap();let b=MatrixExpr::zero(Expr::zero(),Expr::one()).unwrap();assert_eq!(dot_product(&a,&b,true).unwrap(),Expr::zero());
        let n=Expr::symbol_with("n",crate::Assumptions::INTEGER.with(crate::Assumptions::NONNEGATIVE));let a=MatrixExpr::symbol("a",Expr::one(),n.clone()).unwrap();let b=MatrixExpr::symbol("b",n,Expr::one()).unwrap();
        assert!(matches!(dot_product(&a,&b,true).unwrap().kind(),Kind::Func(Func::Sum,_)));
    }
    #[test]
    fn scalar_derivatives_metadata_and_shared_expansion(){
        let t=Expr::symbol_with("t",crate::Assumptions::REAL);
        let a=MatrixExpr::explicit(Matrix::from_rows(vec![vec![crate::ops::add2(&Expr::one(),&t),crate::ops::pow(&t,&Expr::int(2))]]).unwrap());
        let b=MatrixExpr::explicit(Matrix::from_rows(vec![vec![t.clone()],vec![crate::ops::add2(&Expr::one(),&t)]]).unwrap());
        let value=dot_product(&a,&b,false).unwrap();
        assert_eq!(crate::traverse::free_symbols(&value),vec![t.clone()]);
        let exact=dot_product(&a,&b,true).unwrap();
        for order in 1..=3{
            let actual=crate::functions::doit(&crate::diff::diff_many(&value,&[(t.clone(),order)]),true);
            let expected=crate::diff::diff_many(&exact,&[(t.clone(),order)]);
            assert!(crate::expand::expand(&sub(&actual,&expected)).is_zero());
        }
        let mut dag=t;
        for _ in 0..35{let parent=MatrixExpr::explicit(Matrix::from_rows(vec![vec![dag]]).unwrap());dag=dot_product(&parent,&parent,false).unwrap();}
        assert_eq!(crate::expand::expand(&dag),dag);
        let zero=MatrixExpr::explicit(Matrix::from_rows(vec![vec![Expr::zero()]]).unwrap());
        let infinity=MatrixExpr::explicit(Matrix::from_rows(vec![vec![crate::ops::mul2(&Expr::oo(),&Expr::symbol("x"))]]).unwrap());
        assert_eq!(dot_product(&zero,&infinity,true).unwrap(),Expr::nan());
    }
}
