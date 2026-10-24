//! Structural instance queries across the scalar and matrix DAGs.
use rustc_hash::FxHashSet;
use crate::expr::{Expr, Func, Kind};
use super::{MatrixExpr, MatrixKind, MatrixOperand};

enum Node { Matrix(MatrixExpr), Scalar(Expr) }

/// Search public expression arguments without materializing matrix entries.
/// Matrix dimensions are arguments of leaves, but not implicit arguments of
/// products; the invisible scalar unit of a MatMul is not a public argument.
fn contains(root:Node,scalar:Option<&Expr>,matrix:Option<&MatrixExpr>)->bool{
    let mut pending=vec![root];
    let mut scalars=FxHashSet::default();let mut matrices=FxHashSet::default();
    // Keep owners alive: scalar args can manufacture temporary canonical terms.
    let mut owners=Vec::new();
    while let Some(node)=pending.pop(){match node{
        Node::Scalar(value)=>{
            if !scalars.insert(value.dag_id()){continue;}
            if scalar.is_some_and(|target|&value==target){return true;}
            if let (Kind::Add(c,terms),Some(target))=(value.kind(),scalar){
                if let Kind::Add(tc,tts)=target.kind(){
                    if (tc.is_zero()||c==tc)&&tts.iter().all(|t|terms.contains(t)){return true;}
                }
            }
            if let (Kind::Mul(..),Some(target))=(value.kind(),scalar){
                if matches!(target.kind(),Kind::Mul(..)){
                    let args=value.args();let wanted=target.args();
                    let (comm,nc):(Vec<_>,Vec<_>)=args.iter().partition(|e|crate::decompose::commutativity(e)==Some(true));
                    let (wc,wn):(Vec<_>,Vec<_>)=wanted.iter().partition(|e|crate::decompose::commutativity(e)==Some(true));
                    if wc.iter().all(|e|comm.contains(e))&&(wn.is_empty()||nc.windows(wn.len()).any(|w|w==wn.as_slice())){return true;}
                }
            }
            match value.kind(){
                Kind::Func(Func::MatrixDotProduct(a,b),_)=>{
                    pending.push(Node::Matrix(a.as_ref().clone()));pending.push(Node::Matrix(b.as_ref().clone()));
                },
                Kind::Func(Func::MatrixElement(parent)|Func::MatrixProxy(parent)|Func::MatrixTrace(parent)|Func::MatrixDeterminant(parent),_)=>{
                    pending.push(Node::Matrix(parent.as_ref().clone()));
                },_=>{},
            }
            pending.extend(value.args().into_iter().map(Node::Scalar));
            owners.push(value);
        },
        Node::Matrix(value)=>{
            if !matrices.insert(value.dag_id()){continue;}
            if matrix.is_some_and(|target|&value==target){return true;}
            if let Some(target)=matrix{
                match (value.kind(),target.kind()){
                    (MatrixKind::Add(args),MatrixKind::Add(wanted))=>if wanted.iter().all(|e|args.contains(e)){return true;},
                    (MatrixKind::Mul(c,args),MatrixKind::Mul(tc,wanted))=>{
                        if (tc.is_one()||c==tc)&&!wanted.is_empty()&&args.windows(wanted.len()).any(|w|w==wanted.as_slice()){return true;}
                    },_=>{},
                }
            }
            pending.extend(value.children().into_iter().map(Node::Matrix));
            let mut fields=Vec::new();
            match value.kind(){
                MatrixKind::Permute(_,_,axis)=>fields.push(Expr::int(*axis as i64)),
                MatrixKind::Symbol(_)|MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>fields.extend([value.rows().clone(),value.cols().clone()]),
                MatrixKind::Dft(_)=>fields.push(value.rows().clone()),
                MatrixKind::Function(lambda)=>fields.extend([value.rows().clone(),value.cols().clone(),lambda.clone()]),
                MatrixKind::ApplyFunction(function,_)=>fields.push(function.clone()),
                MatrixKind::Explicit(m)=>fields.extend(m.data().iter().cloned()),
                MatrixKind::Mul(coefficient,_)=>if !coefficient.is_one(){fields.push(coefficient.clone());},
                MatrixKind::Pow(_,exponent)=>fields.push(exponent.clone()),
                MatrixKind::Inverse(_)=>fields.push(Expr::minus_one()),
                MatrixKind::Slice(_,rows,cols)=>for spec in [rows,cols]{fields.extend(spec.args().into_iter().flatten());},
                MatrixKind::HadamardPower(a,b)=>for operand in [a,b]{if let MatrixOperand::Scalar(e)=operand{fields.push(e.clone());}},
                _=>{},
            }
            pending.extend(fields.into_iter().map(Node::Scalar));
        },
    }}
    false
}

impl MatrixExpr{
    /// Structural scalar instance containment, including typed scalar parents.
    pub fn has_scalar(&self,target:&Expr)->bool{contains(Node::Matrix(self.clone()),Some(target),None)}
    /// Structural matrix instance containment, including matrices in Trace,
    /// Determinant, MatrixElement and DotProduct scalar arguments.
    pub fn has_matrix(&self,target:&Self)->bool{contains(Node::Matrix(self.clone()),None,Some(target))}
}

pub fn expression_has_scalar(value:&Expr,target:&Expr)->bool{contains(Node::Scalar(value.clone()),Some(target),None)}
pub fn expression_has_matrix(value:&Expr,target:&MatrixExpr)->bool{contains(Node::Scalar(value.clone()),None,Some(target))}

#[cfg(test)]
mod tests{
    use super::*;
    use crate::ops;
    #[test]
    fn public_arguments_and_shared_scalar_graph(){
        let x=Expr::symbol("x");let y=Expr::symbol("y");
        let a=MatrixExpr::symbol("A",Expr::int(2),Expr::int(2)).unwrap();
        let product=MatrixExpr::mul(Expr::one(),&[a.clone(),a.clone()],false).unwrap();
        assert!(product.has_matrix(&a));assert!(product.has_scalar(&Expr::int(2)));
        assert!(!product.has_scalar(&Expr::one()));
        let mut shared=x.clone();
        for _ in 0..80{shared=Expr::from_unevaluated_kind(Kind::Add(crate::number::Number::ZERO,vec![(shared.clone(),crate::number::Number::ONE),(shared.clone(),crate::number::Number::ONE)]));}
        let product=MatrixExpr::mul(shared,&[a],false).unwrap();
        assert!(product.has_scalar(&x));assert!(!product.has_scalar(&y));
        let coefficient=ops::mul(&[Expr::int(7),x]);
        let product=MatrixExpr::mul(coefficient,&[product],false).unwrap();
        assert!(product.has_scalar(&Expr::int(7)));
    }
}
