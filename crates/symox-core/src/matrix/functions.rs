//! Exact matrix functions through polynomial evaluation and Hermite interpolation.
//!
//! Interpolation modulo the characteristic polynomial handles repeated and
//! defective eigenvalues without a choice of eigenvector/Jordan basis.
use rustc_hash::{FxHashMap,FxHashSet};
use super::{Matrix,MatrixError,MResult};
use crate::expr::{Constant,Expr,Func,Kind};
use crate::number::Number;
use crate::ops::mul2;
#[path="function_budget.rs"]
pub(in crate::matrix) mod budget;
use budget::Budget;
#[path="skew_exp.rs"]
mod skew_exp;

fn finite(e:&Expr,seen:&mut FxHashSet<usize>)->bool{
    if !seen.insert(e.dag_id()){return true;}
    match e.kind(){
        Kind::Sym(symbol) if symbol.facts().get("finite")==Some(false)=>return false,
        Kind::Func(Func::Undefined(head),_) if !head.internal&&head.facts().get("finite")==Some(false)=>return false,
        Kind::Const(Constant::Infinity|Constant::NegInfinity|Constant::ComplexInfinity|Constant::NaN)=>return false,
        Kind::Num(Number::Float(value)) if !value.is_finite()=>return false,
        Kind::Add(c,terms) if matches!(c,Number::Float(v) if !v.is_finite())||terms.iter().any(|(_,c)|matches!(c,Number::Float(v) if !v.is_finite()))=>return false,
        Kind::Mul(c,_) if matches!(c,Number::Float(v) if !v.is_finite())=>return false,
        _=>{},
    }
    budget::raw_children(e).iter().all(|child|finite(child,seen))
}
fn ensure_finite(e:&Expr)->MResult<()>{
    if finite(e,&mut FxHashSet::default()){Ok(())}
    else{Err(MatrixError::Value("Matrix function is not analytic at the spectrum: a pole or nonfinite value was encountered".into()))}
}
fn check_matrix(matrix:&Matrix)->MResult<()>{
    if !matrix.is_square(){return Err(MatrixError::NonSquare("Matrix functions require a square matrix".into()));}
    let mut seen=FxHashSet::default();
    if matrix.data().iter().all(|e|finite(e,&mut seen)){Ok(())}
    else{Err(MatrixError::Value("Matrix functions require finite matrix entries".into()))}
}
fn unavailable(error:MatrixError)->MatrixError{
    MatrixError::NotImplemented(format!("Cannot construct an exact analytic matrix function: {error}"))
}

// Recognize/evaluate the polynomial AST directly. This bypasses the spectral
// solver completely, including for matrices with unsolvable characteristic
// polynomials, and bounds every multiplication during binary exponentiation.
fn polynomial(e:&Expr,x:&Expr,matrix:&Matrix,memo:&mut FxHashMap<usize,(Expr,Option<Matrix>)>,budget:&mut Budget)->MResult<Option<Matrix>>{
    if let Some((_,value))=memo.get(&e.dag_id()){return Ok(value.clone());}
    budget.expression(e)?;
    let result=if e==x{Some(matrix.clone())}
    else if !e.may_contain(x)||!crate::traverse::free_symbols(e).contains(x){let identity=budget.identity(matrix.rows())?;Some(budget.scale(&identity,e)?)}
    else{match e.kind(){
        Kind::Add(c,terms)=>{
            let identity=budget.identity(matrix.rows())?;
            let mut value=budget.scale(&identity,&Expr::num(c.clone()))?;
            for(term,coefficient)in terms{
                let Some(term)=polynomial(term,x,matrix,memo,budget)?else{return Ok(None);};
                let term=budget.scale(&term,&Expr::num(coefficient.clone()))?;
                value=budget.add(&value,&term)?;
            }
            Some(value)
        },
        Kind::Mul(c,factors)=>{
            let identity=budget.identity(matrix.rows())?;
            let mut value=budget.scale(&identity,&Expr::num(c.clone()))?;
            for(base,exponent)in factors{
                let Some(n)=budget::exponent(exponent)else{return Ok(None);};
                let Some(base)=polynomial(base,x,matrix,memo,budget)?else{return Ok(None);};
                let power=budget.power(&base,n)?;
                value=budget.mul(&value,&power)?;
            }
            Some(value)
        },
        Kind::Pow(base,exponent)=>{
            if let Some(n)=budget::exponent(exponent){
                polynomial(base,x,matrix,memo,budget)?.map(|base|budget.power(&base,n)).transpose()?
            }else{None}
        },
        _=>None,
    }};
    memo.insert(e.dag_id(),(e.clone(),result.clone()));
    Ok(result)
}

impl Matrix{
    /// Evaluate a scalar analytic function on the matrix spectrum exactly.
    /// `x` is an independent scalar symbol, not a parameter of the matrix.
    pub fn analytic_func(&self,function:&Expr,x:&Expr)->MResult<Matrix>{
        check_matrix(self)?;
        if !x.is_symbol(){return Err(MatrixError::Value("The analytic function parameter must be a symbol".into()));}
        if !crate::traverse::free_symbols(function).contains(x){return Err(MatrixError::Value("The symbol must be a parameter of the analytic function".into()));}
        if self.data().iter().any(|entry|crate::traverse::free_symbols(entry).contains(x)){
            return Err(MatrixError::Value("The analytic function symbol must not be a parameter of the matrix".into()));
        }
        if matches!(function.kind(),Kind::Bool(_)|Kind::Rel(..)|Kind::Func(Func::Logic(_)|Func::Set(_)|Func::MatrixProxy(_),_)){
            return Err(MatrixError::Value("The analytic function must be a scalar expression".into()));
        }
        ensure_finite(function)?;
        if self.rows()==0{return Ok(Matrix::zeros(0,0));}
        let mut budget=Budget::new();budget.expression(function)?;budget.matrix(self)?;
        if let Some(value)=polynomial(function,x,self,&mut FxHashMap::default(),&mut budget)?{
            for entry in value.data(){ensure_finite(entry)?;}return Ok(value);
        }
        if self.is_diagonal(){
            let mut values=Vec::with_capacity(self.rows());
            for i in 0..self.rows(){values.push(sample(function,x,self.get(i,i),&mut budget)?);}
            return budget.diagonal(&values);
        }
        let input_height=budget.matrix(self)?;budget.interpolation(self.rows(),input_height)?;
        budget.spectral_preflight(self)?;
        let eigenvalues=self.eigenvals().map_err(unavailable)?;
        if eigenvalues.iter().map(|(_,n)|*n).sum::<usize>()!=self.rows(){
            return Err(MatrixError::NotImplemented("The exact matrix spectrum is incomplete".into()));
        }
        let max_order=eigenvalues.iter().map(|(_,n)|*n).max().unwrap_or(0);
        let mut derivatives=vec![function.clone()];
        for order in 1..max_order{budget.work(1)?;let derivative=crate::diff::diff(&derivatives[order-1],x);budget.expression(&derivative)?;derivatives.push(derivative);}
        let n=self.rows();let mut equations=budget.zeros(n)?;let mut samples=Vec::with_capacity(n);let mut row=0;
        for(lambda,multiplicity)in eigenvalues{
            budget.expression(&lambda)?;ensure_finite(&lambda)?;
            for(order,derivative)in derivatives.iter().enumerate().take(multiplicity){
                samples.push(sample(derivative,x,&lambda,&mut budget)?);
                for degree in order..n{
                    let mut coefficient=Expr::one();
                    for k in 0..order{coefficient=mul2(&coefficient,&Expr::int((degree-k) as i64));}
                    let power=budget.scalar_power(&lambda,&Expr::int((degree-order) as i64))?;
                    equations.set(row,degree,mul2(&coefficient,&power));
                }
                row+=1;
            }
        }
        let values=Matrix::column(samples);
        let height=budget.matrix(&equations)?.max(budget.matrix(&values)?);budget.interpolation(n,height)?;
        let solved=equations.solve(&values).map_err(unavailable)?;
        let coefficients=budget.cancel(&solved)?;
        budget.matrix(&coefficients)?;
        let residual=budget.mul(&equations,&coefficients)?.sub(&values)?;
        if !budget.cancel(&residual)?.is_zero_matrix(){
            return Err(MatrixError::NotImplemented("Cannot certify the exact Hermite interpolation equations".into()));
        }
        let identity=budget.identity(n)?;let mut value=budget.zeros(n)?;
        for coefficient in coefficients.data().iter().rev(){let product=budget.mul(&value,self)?;let term=budget.scale(&identity,coefficient)?;let sum=budget.add(&product,&term)?;value=budget.cancel(&sum)?;}
        for entry in value.data(){ensure_finite(entry)?;}
        Ok(value)
    }

    /// Primary matrix logarithm obtained by applying the existing scalar log
    /// at every eigenvalue. Negative real eigenvalues use scalar log's +I*pi
    /// determination; no global log-product or log-exp identity is assumed.
    pub fn log(&self)->MResult<Matrix>{self.log_with_simplification(true)}

    /// Optionally cancel the final result. Necessary simplifications inside
    /// Hermite interpolation and its exact certificates are always retained.
    pub fn log_with_simplification(&self,simplify:bool)->MResult<Matrix>{
        crate::assumptions::exact_queries(||{
            check_matrix(self)?;
            let x=Expr::dummy("_matrix_log",crate::Assumptions::default());
            let result=self.analytic_func(&crate::functions::func(Func::Log,vec![x.clone()]),&x)?;
            if simplify{Budget::new().cancel(&result)}else{Ok(result)}
        })
    }

    /// Matrix exponential. Exact scalar exponentials remain symbolic; this
    /// never approximates unresolved eigenvalues or Jordan blocks numerically.
    pub fn exp(&self)->MResult<Matrix>{
        check_matrix(self)?;
        // Diagonal matrices retain their existing scalar evaluation path,
        // including precision ranges beyond the new rotation evaluator.
        if !self.is_diagonal(){
            if let Some(result)=skew_exp::evaluate(self)?{return Ok(result);}
        }
        let x=Expr::dummy("_matrix_exp",crate::Assumptions::default());
        let result=self.analytic_func(&crate::functions::func(Func::Exp,vec![x.clone()]),&x)?;
        if crate::assumptions::exact_queries(||self.data().iter().all(|e|crate::assumptions::is_real(e)==Some(true))){
            let projected=result.applyfunc(|e|crate::functions::func(Func::Re,vec![e.clone()]));
            Budget::new().cancel(&projected)
        }else{Ok(result)}
    }
}
fn sample(function:&Expr,x:&Expr,point:&Expr,budget:&mut Budget)->MResult<Expr>{
    budget.scalar_preflight(function,Some((x,point)))?;
    let value=crate::functions::doit(&crate::subs::subs(function,&[(x.clone(),point.clone())]),true);
    budget.expression(&value)?;ensure_finite(&value)?;Ok(value)
}

#[cfg(test)]
mod tests{
    use super::*;
    fn m(data:&[&[i64]])->Matrix{Matrix::from_rows(data.iter().map(|row|row.iter().map(|n|Expr::int(*n)).collect()).collect()).unwrap()}
    fn equal(a:&Matrix,b:&Matrix){assert_eq!(a.shape(),b.shape());assert!(a.sub(b).unwrap().cancel().is_zero_matrix());}
    #[test]
    fn polynomials_resolvent_and_repeated_jordan_exponential(){
        let x=Expr::symbol("x");let a=m(&[&[1,2,0],&[0,1,3],&[0,0,1]]);
        let f=crate::ops::add2(&crate::ops::pow(&x,&Expr::int(7)),&crate::ops::mul2(&Expr::int(3),&x));
        equal(&a.analytic_func(&f,&x).unwrap(),&a.powi(7).unwrap().add(&a.scale(&Expr::int(3))).unwrap());
        let resolvent=crate::ops::pow(&crate::ops::sub(&x,&Expr::int(3)),&Expr::minus_one());
        equal(&a.analytic_func(&resolvent,&x).unwrap(),&a.sub(&Matrix::eye(3).scale(&Expr::int(3))).unwrap().inv().unwrap());
        let exponential=a.exp().unwrap();let expected=m(&[&[1,2,3],&[0,1,3],&[0,0,1]]).scale(&crate::functions::func(Func::Exp,vec![Expr::one()]));
        equal(&exponential,&expected);
        assert!(matches!(a.analytic_func(&crate::ops::pow(&crate::ops::sub(&x,&Expr::one()),&Expr::minus_one()),&x),Err(MatrixError::Value(_))));
        assert_eq!(Matrix::zeros(0,0).exp().unwrap().shape(),(0,0));
    }
    #[test]
    fn logarithm_nilpotent_series_scalar_branches_and_similarity(){
        let nil=m(&[&[0,2,0],&[0,0,3],&[0,0,0]]);
        for lambda in [Expr::one(),Expr::int(2),Expr::int(-2),Expr::i(),Expr::symbol_with("p",crate::Assumptions::POSITIVE)]{
            let a=Matrix::eye(3).scale(&lambda).add(&nil).unwrap();
            let inverse=crate::ops::pow(&lambda,&Expr::minus_one());
            let expected=Matrix::eye(3).scale(&crate::functions::log(&lambda))
                .add(&nil.scale(&inverse)).unwrap()
                .sub(&nil.powi(2).unwrap().scale(&crate::ops::mul2(&Expr::rational(1,2),&crate::ops::pow(&inverse,&Expr::int(2))))).unwrap();
            let actual=a.log().unwrap();equal(&actual,&expected);
            equal(&a.matmul(&actual).unwrap(),&actual.matmul(&a).unwrap());
            equal(&actual.exp().unwrap(),&a);
        }
        let p=m(&[&[1,2,-1],&[0,1,3],&[0,0,1]]);
        let a=Matrix::eye(3).add(&nil).unwrap();
        let transformed=p.matmul(&a).unwrap().matmul(&p.inv().unwrap()).unwrap();
        equal(&transformed.log().unwrap(),&p.matmul(&a.log().unwrap()).unwrap().matmul(&p.inv().unwrap()).unwrap());
    }

    #[test]
    fn logarithm_domain_empty_and_spectral_preflight(){
        assert_eq!(Matrix::zeros(0,0).log().unwrap().shape(),(0,0));
        assert!(matches!(Matrix::zeros(0,2).log(),Err(MatrixError::NonSquare(_))));
        assert!(matches!(Matrix::zeros(2,2).log(),Err(MatrixError::Value(_))));
        assert!(matches!(m(&[&[0,1],&[0,0]]).log(),Err(MatrixError::Value(_))));
        for value in [Expr::oo(),Expr::zoo(),Expr::nan()]{assert!(matches!(Matrix::new(1,1,vec![value]).unwrap().log(),Err(MatrixError::Value(_))));}
        let entries=(0..36).map(|i|crate::ops::add(&(0..6).map(|j|Expr::symbol(&format!("a{i}_{j}"))).collect::<Vec<_>>())).collect();
        assert!(matches!(Matrix::new(6,6,entries).unwrap().log(),Err(MatrixError::NotImplemented(_))));
        let held=Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2),Expr::int(1_000_000_000)));
        assert!(matches!(Matrix::new(1,1,vec![held]).unwrap().log(),Err(MatrixError::NotImplemented(_))));
    }

    #[test]
    fn nonfinite_declared_symbols_are_not_finite_matrix_parameters(){
        let x=Expr::symbol("x");
        for flags in [vec![("finite",false)],vec![("infinite",true)],vec![("extended_real",true),("finite",false)]] {
            let facts=crate::assumption_facts::FactSet::from_kwargs(flags).unwrap();
            let symbol=crate::Symbol::with_facts("z",facts).unwrap();
            let z=Expr::from_kind(Kind::Sym(symbol));
            let matrix=Matrix::new(1,1,vec![z.clone()]).unwrap();
            assert!(matches!(matrix.exp(),Err(MatrixError::Value(_))));
            assert!(matches!(matrix.analytic_func(&x,&x),Err(MatrixError::Value(_))));
            let f=Expr::from_unevaluated_kind(Kind::Add(Number::int(0),vec![(x.clone(),Number::int(1)),(z.clone(),Number::int(1))]));
            assert!(matches!(Matrix::eye(1).analytic_func(&f,&x),Err(MatrixError::Value(_))));
            assert!(matches!(sample(&x,&x,&z,&mut Budget::new()),Err(MatrixError::Value(_))));
        }
        for flags in [vec![],vec![("extended_real",true)]] {
            let facts=crate::assumption_facts::FactSet::from_kwargs(flags).unwrap();
            let symbol=crate::Symbol::with_facts("z",facts).unwrap();
            let z=Expr::from_kind(Kind::Sym(symbol));
            assert!(Matrix::new(1,1,vec![z]).unwrap().exp().is_ok());
        }
    }
    #[test]
    fn functional_calculus_resource_guards_precede_evaluation(){
        let x=Expr::symbol("x");
        let held=Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2),Expr::int(1_000_000_000)));
        let raw=Expr::from_unevaluated_kind(Kind::Mul(Number::int(1),vec![(Expr::int(2),Expr::int(1_000_000_000))]));
        assert!(matches!(Matrix::new(1,1,vec![raw]).unwrap().exp(),Err(MatrixError::NotImplemented(_))));
        let parent=crate::matrixexpr::MatrixExpr::explicit(Matrix::new(1,1,vec![held.clone()]).unwrap());
        let trace=Expr::from_kind(Kind::Func(Func::MatrixTrace(std::sync::Arc::new(parent)),vec![]));
        assert!(matches!(Matrix::new(1,1,vec![trace]).unwrap().exp(),Err(MatrixError::NotImplemented(_))));
        let a=Matrix::new(1,1,vec![held.clone()]).unwrap();
        assert!(matches!(a.exp(),Err(MatrixError::NotImplemented(_))));
        let f=Expr::from_unevaluated_kind(Kind::Mul(Number::int(1),vec![(x.clone(),Expr::one()),(held,Expr::one())]));
        let a=m(&[&[2,1],&[0,2]]);
        assert!(matches!(a.analytic_func(&f,&x),Err(MatrixError::NotImplemented(_))));
        let power=Expr::from_unevaluated_kind(Kind::Pow(x.clone(),Expr::int(1_000_000_000)));
        assert!(matches!(a.analytic_func(&power,&x),Err(MatrixError::NotImplemented(_))));
        let nested=crate::functions::func(Func::Sin,vec![power]);
        assert!(matches!(m(&[&[2]]).analytic_func(&nested,&x),Err(MatrixError::NotImplemented(_))));
        let gamma=crate::functions::func(Func::Gamma,vec![x.clone()]);
        assert!(matches!(m(&[&[1_000_000_000]]).analytic_func(&gamma,&x),Err(MatrixError::NotImplemented(_))));
        let huge=Expr::num(Number::from_bigint(num_bigint::BigInt::from(1u32)<<100));
        let f=Expr::from_unevaluated_kind(Kind::Pow(x.clone(),huge));
        assert_eq!(Matrix::eye(2).analytic_func(&f,&x).unwrap(),Matrix::eye(2));
        assert!(m(&[&[0,1],&[0,0]]).analytic_func(&f,&x).unwrap().is_zero_matrix());
        assert_eq!(m(&[&[0,1],&[1,0]]).analytic_func(&f,&x).unwrap(),Matrix::eye(2));
        let mut dag=x.clone();for _ in 0..50{dag=Expr::from_unevaluated_kind(Kind::Mul(Number::int(1),vec![(dag.clone(),Expr::one()),(dag,Expr::one())]));}
        assert!(matches!(a.analytic_func(&dag,&x),Err(MatrixError::NotImplemented(_))));
    }

}
