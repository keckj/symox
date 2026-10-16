//! Typed, noncommutative symbolic matrix expressions over scalar Expr dimensions.
use std::cmp::Ordering;
use std::sync::Arc;
use rustc_hash::{FxHashMap, FxHashSet};
use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::matrix::{Matrix, MatrixError, MResult};
use crate::matrix::ratfunc::is_zero_assume;
use crate::ops::{add, mul2, neg, pow, sub};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct MatrixExpr(Arc<MatrixNode>);
#[derive(Debug, PartialEq, Eq, Hash)]
struct MatrixNode { kind: MatrixKind, rows: Expr, cols: Expr }
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MatrixKind {
    Symbol(String), Explicit(Matrix), Identity, Zero, One,
    Add(Vec<MatrixExpr>), Mul(Expr, Vec<MatrixExpr>), Pow(MatrixExpr, Expr),
    Transpose(MatrixExpr), Adjoint(MatrixExpr), Inverse(MatrixExpr), Conjugate(MatrixExpr),
    HadamardProduct(Vec<MatrixExpr>), HadamardPower(MatrixOperand,MatrixOperand),
    Kronecker(Vec<MatrixExpr>),
}
#[derive(Clone,Debug,PartialEq,Eq,Hash,PartialOrd,Ord)]
pub enum MatrixOperand {Scalar(Expr),Matrix(MatrixExpr)}
impl MatrixOperand {
    pub fn matrix(&self)->Option<&MatrixExpr>{if let Self::Matrix(m)=self{Some(m)}else{None}}
    fn scalar(&self)->Option<&Expr>{if let Self::Scalar(e)=self{Some(e)}else{None}}
    fn entry(&self,i:&Expr,j:&Expr)->MResult<Expr>{match self{Self::Scalar(e)=>Ok(e.clone()),Self::Matrix(m)=>m.element(i,j)}}
    fn map(&self,matrices:&mut dyn FnMut(&MatrixExpr)->MatrixExpr,scalars:&mut dyn FnMut(&Expr)->Expr)->Self{match self{Self::Scalar(e)=>Self::Scalar(scalars(e)),Self::Matrix(m)=>Self::Matrix(matrices(m))}}
}
impl Ord for MatrixExpr {
    fn cmp(&self, other: &Self) -> Ordering {
        if Arc::ptr_eq(&self.0,&other.0){return Ordering::Equal;}
        let kind=self.kind_name().cmp(other.kind_name());
        if kind!=Ordering::Equal{return kind;}
        match (self.kind(),other.kind()){
            (MatrixKind::Symbol(a),MatrixKind::Symbol(b))=>a.cmp(b).then_with(||self.rows().cmp(other.rows())).then_with(||self.cols().cmp(other.cols())),
            (MatrixKind::Explicit(a),MatrixKind::Explicit(b))=>a.shape().cmp(&b.shape()).then_with(||a.data().cmp(b.data())),
            (MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One,MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One)=>self.rows().cmp(other.rows()).then_with(||self.cols().cmp(other.cols())),
            (MatrixKind::Add(a),MatrixKind::Add(b))|(MatrixKind::HadamardProduct(a),MatrixKind::HadamardProduct(b))|(MatrixKind::Kronecker(a),MatrixKind::Kronecker(b))=>a.cmp(b),
            (MatrixKind::HadamardPower(a,aa),MatrixKind::HadamardPower(b,bb))=>a.cmp(b).then_with(||aa.cmp(bb)),
            (MatrixKind::Mul(a,aa),MatrixKind::Mul(b,bb))=>a.cmp(b).then_with(||aa.cmp(bb)),
            (MatrixKind::Pow(a,aa),MatrixKind::Pow(b,bb))=>a.cmp(b).then_with(||aa.cmp(bb)),
            (MatrixKind::Transpose(a),MatrixKind::Transpose(b))|(MatrixKind::Adjoint(a),MatrixKind::Adjoint(b))|(MatrixKind::Inverse(a),MatrixKind::Inverse(b))|(MatrixKind::Conjugate(a),MatrixKind::Conjugate(b))=>a.cmp(b),
            _=>unreachable!("distinct matrix kinds have distinct names"),
        }
    }
}
impl PartialOrd for MatrixExpr { fn partial_cmp(&self, other:&Self)->Option<Ordering>{Some(self.cmp(other))} }

fn dimension(e:&Expr)->MResult<()> {
    if crate::assumptions::is_integer(e) == Some(false) || crate::assumptions::is_negative(e) == Some(true) {
        return Err(MatrixError::Value("Matrix dimensions must be nonnegative integers".into()));
    }
    Ok(())
}
fn compatible(a:&Expr,b:&Expr)->MResult<()> {
    if a != b && is_zero_assume(&sub(a,b)) == Some(false) {
        return Err(MatrixError::Shape(format!("Incompatible matrix dimensions: {a} and {b}")));
    }
    Ok(())
}
fn finite_dim(e:&Expr)->MResult<usize> {
    e.as_i64().and_then(|n| usize::try_from(n).ok()).ok_or_else(||MatrixError::NotImplemented("Explicit matrices require finite integer dimensions".into()))
}
impl MatrixExpr {
    fn new(kind:MatrixKind,rows:Expr,cols:Expr)->Self {Self(Arc::new(MatrixNode{kind,rows,cols}))}
    pub fn kind(&self)->&MatrixKind {&self.0.kind}
    pub fn rows(&self)->&Expr {&self.0.rows}
    pub fn cols(&self)->&Expr {&self.0.cols}
    pub fn shape(&self)->(Expr,Expr) {(self.rows().clone(),self.cols().clone())}
    pub fn symbol(name:&str,rows:Expr,cols:Expr)->MResult<Self> {
        dimension(&rows)?;dimension(&cols)?;
        Ok(Self::new(MatrixKind::Symbol(name.to_string()),rows,cols))
    }
    pub fn identity(n:Expr)->MResult<Self> { dimension(&n)?;Ok(Self::new(MatrixKind::Identity,n.clone(),n)) }
    pub fn zero(rows:Expr,cols:Expr)->MResult<Self> {dimension(&rows)?;dimension(&cols)?;Ok(Self::new(MatrixKind::Zero,rows,cols))}
    pub fn one(rows:Expr,cols:Expr)->MResult<Self> {dimension(&rows)?;dimension(&cols)?;Ok(Self::new(MatrixKind::One,rows,cols))}
    pub fn explicit(m:Matrix)->Self {Self::new(MatrixKind::Explicit(m.clone()),Expr::int(m.rows() as i64),Expr::int(m.cols() as i64))}
    fn scalar_nonfinite(e:&Expr)->bool{
        match e.kind(){Kind::Const(crate::expr::Constant::Infinity|crate::expr::Constant::NegInfinity|crate::expr::Constant::ComplexInfinity|crate::expr::Constant::NaN)=>true,Kind::Add(..)|Kind::Mul(..)=>e.args().iter().any(Self::scalar_nonfinite),_=>false}
    }
    fn has_nonfinite(&self)->bool{
        let mut scalars=Vec::new();self.scalar_fields(&mut scalars);
        scalars.iter().any(Self::scalar_nonfinite)
    }
    fn entry_product(a:&Expr,b:&Expr)->Expr{
        if (a.is_zero()&&Self::scalar_nonfinite(b))||(b.is_zero()&&Self::scalar_nonfinite(a)){Expr::nan()}else{mul2(a,b)}
    }
    pub fn hadamard_product(factors:&[Self],evaluate:bool)->MResult<Self>{
        let first=factors.first().ok_or_else(||MatrixError::Value("HadamardProduct needs at least one argument".into()))?;
        let (rows,cols)=first.shape();
        for factor in factors{compatible(&rows,factor.rows())?;compatible(&cols,factor.cols())?;}
        if !evaluate{return Ok(Self::new(MatrixKind::HadamardProduct(factors.to_vec()),rows,cols));}
        let mut flat=Vec::new();
        for factor in factors{if let MatrixKind::HadamardProduct(args)=factor.kind(){flat.extend(args.iter().cloned());}else{flat.push(factor.clone());}}
        if flat.iter().any(|f|matches!(f.kind(),MatrixKind::Zero))&&!flat.iter().any(Self::has_nonfinite){return Self::zero(rows,cols);}
        let mut explicit:Option<Matrix>=None;let mut rest=Vec::new();
        for factor in flat{match factor.kind(){
            MatrixKind::One=>{},
            MatrixKind::Explicit(matrix)=>{explicit=Some(if let Some(previous)=explicit{Matrix::new(matrix.rows(),matrix.cols(),previous.data().iter().zip(matrix.data()).map(|(a,b)|Self::entry_product(a,b)).collect())?}else{matrix.clone()});},
            _=>rest.push(factor),
        }}
        if let Some(matrix)=explicit{rest.push(Self::explicit(matrix));}
        rest.sort();
        let mut grouped=Vec::new();let mut i=0;
        while i<rest.len(){let mut end=i+1;while end<rest.len()&&rest[end]==rest[i]{end+=1;}
            grouped.push(if end-i==1{rest[i].clone()}else{Self::hadamard_power(MatrixOperand::Matrix(rest[i].clone()),MatrixOperand::Scalar(Expr::int((end-i) as i64)),true)?});i=end;
        }
        grouped.sort();
        Ok(match grouped.len(){0=>Self::one(rows,cols)?,1=>grouped.remove(0),_=>Self::new(MatrixKind::HadamardProduct(grouped),rows,cols)})
    }
    pub fn hadamard_power(base:MatrixOperand,exponent:MatrixOperand,evaluate:bool)->MResult<Self>{
        let matrix=base.matrix().or_else(||exponent.matrix()).ok_or_else(||MatrixError::Type("HadamardPower requires a matrix operand".into()))?;
        let (rows,cols)=matrix.shape();
        if let (Some(a),Some(b))=(base.matrix(),exponent.matrix()){compatible(a.rows(),b.rows())?;compatible(a.cols(),b.cols())?;}
        if evaluate{
            if exponent.scalar().is_some_and(Expr::is_zero){return Self::one(rows,cols);}
            if exponent.scalar().is_some_and(Expr::is_one){return Ok(base.matrix().unwrap().clone());}
            let explicit=|operand:&MatrixOperand|match operand{MatrixOperand::Scalar(_)=>true,MatrixOperand::Matrix(m)=>matches!(m.kind(),MatrixKind::Explicit(_))};
            if explicit(&base)&&explicit(&exponent){
                let mut entries=Vec::new();for i in 0..finite_dim(&rows)?{for j in 0..finite_dim(&cols)?{let(i,j)=(Expr::int(i as i64),Expr::int(j as i64));entries.push(pow(&base.entry(&i,&j)?,&exponent.entry(&i,&j)?));}}
                return Ok(Self::explicit(Matrix::new(finite_dim(&rows)?,finite_dim(&cols)?,entries)?));
            }
        }
        Ok(Self::new(MatrixKind::HadamardPower(base,exponent),rows,cols))
    }
    pub fn kronecker(factors:&[Self],evaluate:bool)->MResult<Self>{
        if factors.is_empty(){return Ok(Self::explicit(Matrix::eye(1)));}
        let rows=crate::ops::mul(&factors.iter().map(|m|m.rows().clone()).collect::<Vec<_>>());
        let cols=crate::ops::mul(&factors.iter().map(|m|m.cols().clone()).collect::<Vec<_>>());
        if factors.iter().all(|m|matches!(m.kind(),MatrixKind::Identity)){return Self::identity(rows);}
        if !evaluate{return Ok(Self::new(MatrixKind::Kronecker(factors.to_vec()),rows,cols));}
        let mut flat=Vec::new();let mut coefficient=Expr::one();
        for factor in factors{match factor.kind(){
            MatrixKind::Kronecker(args)=>flat.extend(args.iter().cloned()),
            MatrixKind::Mul(c,args)=>{coefficient=Self::entry_product(&coefficient,c);flat.push(Self::mul(Expr::one(),args,true)?);},
            _=>flat.push(factor.clone()),
        }}
        if flat.iter().any(|m|matches!(m.kind(),MatrixKind::Zero))&&!flat.iter().any(Self::has_nonfinite){return Self::zero(rows,cols);}
        if flat.len()==1{return Self::mul(coefficient,&flat,true);}
        if flat.iter().all(|m|matches!(m.kind(),MatrixKind::Explicit(_))){
            let node=Self::new(MatrixKind::Kronecker(flat),rows,cols);
            return Ok(Self::explicit(node.as_explicit()?.scale(&coefficient)));
        }
        let node=Self::new(MatrixKind::Kronecker(flat),rows,cols);
        if coefficient.is_one(){Ok(node)}else{Self::mul(coefficient,&[node],true)}
    }
    pub fn add(terms:&[Self],evaluate:bool)->MResult<Self> {
        let first=terms.first().ok_or_else(||MatrixError::Type("Matrix addition requires a matrix".into()))?;
        let (rows,cols)=first.shape();
        for term in terms {compatible(&rows,term.rows())?;compatible(&cols,term.cols())?;}
        if !evaluate {return Ok(Self::new(MatrixKind::Add(terms.to_vec()),rows,cols));}
        let mut flat=Vec::new();
        for term in terms {if let MatrixKind::Add(args)=term.kind(){flat.extend(args.iter().cloned());}else{flat.push(term.clone());}}
        let mut groups:Vec<(Self,Expr)>=Vec::new();
        let mut explicit:Option<Matrix>=None;
        for term in flat {
            if matches!(term.kind(),MatrixKind::Zero){continue;}
            if let MatrixKind::Explicit(m)=term.kind(){explicit=Some(if let Some(prev)=explicit{prev.add(m)?}else{m.clone()});continue;}
            let (base,coefficient)=match term.kind(){
                MatrixKind::Mul(c,factors)=>(Self::mul(Expr::one(),factors,true)?,c.clone()),
                _=>(term.clone(),Expr::one()),
            };
            if let Some((_,c))=groups.iter_mut().find(|(b,_)|b==&base){*c=crate::ops::add2(c,&coefficient);}else{groups.push((base,coefficient));}
        }
        let mut result=Vec::new();
        if let Some(m)=explicit {if !m.is_zero_matrix(){result.push(Self::explicit(m));}}
        for (base,c) in groups {if !c.is_zero(){result.push(Self::mul(c,&[base],true)?);}}
        result.sort();
        Ok(match result.len(){0=>Self::zero(rows,cols)?,1=>result.remove(0),_=>Self::new(MatrixKind::Add(result),rows,cols)})
    }
    pub fn mul(coefficient:Expr,factors:&[Self],evaluate:bool)->MResult<Self> {
        let first=factors.first().ok_or_else(||MatrixError::Type("Matrix multiplication requires a matrix".into()))?;
        let rows=first.rows().clone(); let cols=factors.last().unwrap().cols().clone();
        for pair in factors.windows(2){compatible(pair[0].cols(),pair[1].rows())?;}
        if !evaluate{return Ok(Self::new(MatrixKind::Mul(coefficient,factors.to_vec()),rows,cols));}
        let mut scalar=coefficient;
        let mut flat=Vec::new();
        for factor in factors {if let MatrixKind::Mul(c,fs)=factor.kind(){scalar=mul2(&scalar,c);flat.extend(fs.iter().cloned());}else{flat.push(factor.clone());}}
        if scalar.is_zero() || flat.iter().any(|f|matches!(f.kind(),MatrixKind::Zero)) {return Self::zero(rows,cols);}
        let mut result:Vec<Self>=Vec::new();
        for factor in flat {
            if matches!(factor.kind(),MatrixKind::Identity){continue;}
            if let Some(prev)=result.last().cloned(){
                if let (MatrixKind::Explicit(a),MatrixKind::Explicit(b))=(prev.kind(),factor.kind()){
                    result.pop();result.push(Self::explicit(a.matmul(b)?));continue;
                }
                let (a,ae)=prev.as_base_exp();let(b,be)=factor.as_base_exp();
                if a==b {
                    result.pop();let merged=a.pow(crate::ops::add2(&ae,&be),true)?;
                    if !matches!(merged.kind(),MatrixKind::Identity){result.push(merged);}continue;
                }
            }
            result.push(factor);
        }
        if result.is_empty(){return Ok(if scalar.is_one(){Self::identity(rows)?}else{Self::new(MatrixKind::Mul(scalar,vec![Self::identity(rows.clone())?]),rows,cols)});}
        if result.len()==1 {
            if scalar.is_one(){return Ok(result.remove(0));}
            if let MatrixKind::Explicit(m)=result[0].kind(){return Ok(Self::explicit(m.scale(&scalar)));}
            if let MatrixKind::Add(args)=result[0].kind(){return Self::add(&args.iter().map(|a|Self::mul(scalar.clone(),&[a.clone()],true)).collect::<MResult<Vec<_>>>()?,true);}
        }
        Ok(Self::new(MatrixKind::Mul(scalar,result),rows,cols))
    }
    pub fn as_base_exp(&self)->(Self,Expr){match self.kind(){MatrixKind::Pow(b,e)=>(b.clone(),e.clone()),MatrixKind::Inverse(b)=>(b.clone(),Expr::minus_one()),_=>(self.clone(),Expr::one())}}
    pub fn pow(&self, exponent:Expr, evaluate:bool)->MResult<Self>{
        compatible(self.rows(),self.cols())?;
        if evaluate {
            if exponent.is_zero(){return Self::identity(self.rows().clone());}
            if exponent.is_one(){return Ok(self.clone());}
            if exponent==Expr::minus_one(){return self.inverse(true);}
            if matches!(self.kind(),MatrixKind::Identity){return Ok(self.clone());}
            if let Some(n)=exponent.as_i64(){
                if let MatrixKind::Explicit(m)=self.kind(){return Ok(Self::explicit(m.powi(n)?));}
                if matches!(self.kind(),MatrixKind::Zero){if n<0{return Err(MatrixError::NonInvertible("Zero matrix is not invertible".into()));}return Ok(self.clone());}
                if let MatrixKind::Pow(b,e)=self.kind(){if e.is_integer_number(){return b.pow(mul2(e,&exponent),true);}}
                if let MatrixKind::Inverse(b)=self.kind(){return b.pow(neg(&exponent),true);}
            }
        }
        Ok(Self::new(MatrixKind::Pow(self.clone(),exponent),self.rows().clone(),self.cols().clone()))
    }
    pub fn inverse(&self,evaluate:bool)->MResult<Self>{
        compatible(self.rows(),self.cols())?;
        if evaluate {match self.kind(){
            MatrixKind::Identity=>return Ok(self.clone()),MatrixKind::Zero=>return Err(MatrixError::NonInvertible("Zero matrix is not invertible".into())),
            MatrixKind::Inverse(b)=>return Ok(b.clone()),
            MatrixKind::Explicit(m)=>return Ok(Self::explicit(m.inv()?)),
            MatrixKind::One=>{if let Some(n)=self.rows().as_i64(){return if n==0{Self::identity(Expr::zero())}else if n==1{Ok(self.clone())}else{Err(MatrixError::NonInvertible("An all-ones matrix of size greater than one is singular".into()))};}},
            MatrixKind::Kronecker(factors) if factors.iter().all(|m|m.rows()==m.cols())=>return Self::kronecker(&factors.iter().map(|m|m.inverse(true)).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Mul(c,fs) if fs.iter().all(|f|f.rows()==f.cols())=>return Self::mul(pow(c,&Expr::minus_one()),&fs.iter().rev().map(|f|f.inverse(true)).collect::<MResult<Vec<_>>>()?,true),
            _=>{},
        }}
        Ok(Self::new(MatrixKind::Inverse(self.clone()),self.rows().clone(),self.cols().clone()))
    }
    pub fn transpose(&self,evaluate:bool)->MResult<Self>{self.unary(0,evaluate)}
    pub fn adjoint(&self,evaluate:bool)->MResult<Self>{self.unary(1,evaluate)}
    pub fn conjugate(&self,evaluate:bool)->MResult<Self>{self.unary(2,evaluate)}
    fn require_square(&self)->MResult<()>{
        compatible(self.rows(),self.cols()).map_err(|_|MatrixError::NonSquare("Trace and determinant require a square matrix".into()))
    }
    pub fn trace_unevaluated(&self)->MResult<Expr>{
        self.require_square()?;
        Ok(Expr::from_kind(Kind::Func(Func::MatrixTrace(Arc::new(self.clone())),vec![])))
    }
    pub fn determinant_unevaluated(&self)->MResult<Expr>{
        self.require_square()?;
        Ok(Expr::from_kind(Kind::Func(Func::MatrixDeterminant(Arc::new(self.clone())),vec![])))
    }
    /// Evaluate certified trace identities, retaining a typed scalar trace
    /// when the matrix is abstract. No implicit expansion into scalar entries.
    pub fn trace(&self)->MResult<Expr>{
        self.require_square()?;
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),MatrixKind::Identity|MatrixKind::One=>Ok(self.rows().clone()),
            MatrixKind::Kronecker(factors) if factors.iter().all(|m|m.rows()==m.cols())=>Ok(crate::ops::mul(&factors.iter().map(Self::trace).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Explicit(matrix)=>matrix.trace(),
            MatrixKind::Add(terms)=>Ok(add(&terms.iter().map(Self::trace).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Mul(coefficient,factors)=>{
                if factors.len()==1{return Ok(mul2(coefficient,&factors[0].trace()?));}
                // The trace is invariant under cyclic rotations, including
                // rectangular compatible factors forming a closed chain.
                let transposed=factors.iter().rev().map(|factor|factor.transpose(true)).collect::<MResult<Vec<_>>>()?;
                let mut canonical=factors.clone();
                for orientation in [factors,&transposed]{
                    for i in 0..orientation.len(){
                        let rotation=orientation[i..].iter().chain(&orientation[..i]).cloned().collect::<Vec<_>>();
                        if rotation<canonical{canonical=rotation;}
                    }
                }
                let parent=Self::mul(Expr::one(),&canonical,false)?;
                Ok(mul2(coefficient,&parent.trace_unevaluated()?))
            },
            MatrixKind::Transpose(base)=>base.trace(),
            MatrixKind::Adjoint(base)|MatrixKind::Conjugate(base)=>Ok(func(Func::Conjugate,vec![base.trace()?])),
            _=>self.trace_unevaluated(),
        }
    }
    pub fn determinant(&self)->MResult<Expr>{
        self.require_square()?;
        if self.rows().is_zero(){return Ok(Expr::one());}
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),MatrixKind::Identity=>Ok(Expr::one()),
            MatrixKind::One=>match self.rows().as_i64(){Some(1)=>Ok(Expr::one()),Some(n) if n>1=>Ok(Expr::zero()),_=>self.determinant_unevaluated()},
            MatrixKind::Kronecker(factors) if factors.iter().all(|m|m.rows()==m.cols())=>{
                let mut terms=Vec::new();for (index,factor) in factors.iter().enumerate(){let exponent=crate::ops::mul(&factors.iter().enumerate().filter(|(i,_)|*i!=index).map(|(_,m)|m.rows().clone()).collect::<Vec<_>>());terms.push(pow(&factor.determinant()?,&exponent));}Ok(crate::ops::mul(&terms))
            },
            MatrixKind::Explicit(matrix)=>matrix.det(),
            MatrixKind::Mul(coefficient,factors) if factors.iter().all(|factor|factor.rows()==factor.cols())=>{
                let mut values=vec![pow(coefficient,self.rows())];
                values.extend(factors.iter().map(Self::determinant).collect::<MResult<Vec<_>>>()?);
                Ok(crate::ops::mul(&values))
            },
            MatrixKind::Transpose(base)=>base.determinant(),
            MatrixKind::Adjoint(base)|MatrixKind::Conjugate(base)=>Ok(func(Func::Conjugate,vec![base.determinant()?])),
            MatrixKind::Inverse(base)=>{
                let determinant=base.determinant()?;
                if determinant.is_zero(){return Err(MatrixError::NonInvertible("A singular matrix has no inverse".into()));}
                Ok(pow(&determinant,&Expr::minus_one()))
            },
            MatrixKind::Pow(base,exponent) if crate::assumptions::is_integer(exponent)==Some(true)=>{
                let determinant=base.determinant()?;
                if determinant.is_zero()&&crate::assumptions::is_negative(exponent)==Some(true){return Err(MatrixError::NonInvertible("A singular matrix has no negative power".into()));}
                Ok(pow(&determinant,exponent))
            },
            _=>self.determinant_unevaluated(),
        }
    }
    pub fn trace_explicit(&self)->MResult<Expr>{
        self.require_square()?;
        let n=finite_dim(self.rows())?;
        Ok(add(&(0..n).map(|i|self.element(&Expr::int(i as i64),&Expr::int(i as i64))).collect::<MResult<Vec<_>>>()?))
    }
    pub fn determinant_explicit(&self)->MResult<Expr>{self.require_square()?;self.as_explicit()?.det()}
    /// Jacobi's adjugate formula is valid even at singular matrices. In
    /// contrast, det(A)*trace(A^-1*A') loses those cases through division.
    pub fn determinant_diff(&self,variable:&Expr)->MResult<Expr>{
        self.require_square()?;
        let derivative=self.diff(variable)?;
        if matches!(derivative.kind(),MatrixKind::Zero){return Ok(Expr::zero());}
        self.as_explicit()?.adjugate()?.matmul(&derivative.as_explicit()?)?.trace()
    }
    /// Derivative with respect to a scalar variable. Matrix factors retain
    /// their order; dimensions describe the result and are not scalar entries.
    pub fn diff(&self,variable:&Expr)->MResult<Self>{
        if !matches!(variable.kind(),Kind::Sym(_)|Kind::Func(Func::Undefined(_),_)){
            return Err(MatrixError::NotImplemented("MatrixExpr differentiation requires a scalar symbol or scalar function; matrix-variable derivatives require tensor expressions".into()));
        }
        let mut cache=FxHashMap::default();
        self.diff_cached(variable,&mut cache)
    }
    pub fn diff_many(&self,variables:&[(Expr,u32)])->MResult<Self>{
        let mut result=self.clone();
        for (variable,order) in variables {for _ in 0..*order {result=result.diff(variable)?;if matches!(result.kind(),MatrixKind::Zero){break;}}}
        Ok(result)
    }
    fn scalar_diff(e:&Expr,variable:&Expr)->MResult<Expr>{
        if crate::assumptions::is_real(variable)!=Some(true){
            let mut nonholomorphic=false;
            crate::traverse::preorder(e,&mut |node|{
                if let Kind::Func(Func::Conjugate|Func::Re|Func::Im,args)=node.kind(){
                    if args.iter().any(|a|!crate::diff::diff(a,variable).is_zero()){nonholomorphic=true;}
                }
            });
            if nonholomorphic{return Err(MatrixError::NotImplemented("Differentiating conjugate, real or imaginary matrix coefficients requires a real differentiation variable".into()));}
        }
        Ok(crate::diff::diff(e,variable))
    }
    fn diff_cached(&self,variable:&Expr,cache:&mut FxHashMap<Self,Self>)->MResult<Self>{
        if let Some(result)=cache.get(self){return Ok(result.clone());}
        let zero=||Self::zero(self.rows().clone(),self.cols().clone());
        let result=match self.kind(){
            MatrixKind::Symbol(_)|MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>zero()?,
            MatrixKind::Explicit(matrix)=>{
                let derived=matrix.try_applyfunc::<MatrixError>(|e|Self::scalar_diff(e,variable))?;
                if derived.is_zero_matrix(){zero()?}else{Self::explicit(derived)}
            },
            MatrixKind::Add(terms)=>Self::add(&terms.iter().map(|a|a.diff_cached(variable,cache)).collect::<MResult<Vec<_>>>()?,true)?,
            MatrixKind::HadamardProduct(factors)|MatrixKind::Kronecker(factors)=>{
                let mut terms=Vec::new();
                for (i,factor) in factors.iter().enumerate(){let derivative=factor.diff_cached(variable,cache)?;if matches!(derivative.kind(),MatrixKind::Zero){continue;}let mut product=factors.clone();product[i]=derivative;
                    terms.push(if matches!(self.kind(),MatrixKind::Kronecker(_)){Self::kronecker(&product,true)?}else{Self::hadamard_product(&product,true)?});}
                if terms.is_empty(){zero()?}else{Self::add(&terms,true)?}
            },
            MatrixKind::HadamardPower(base,exponent)=>{
                if exponent.scalar().is_some_and(Expr::is_zero){return zero();}
                let mut derive=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>Self::scalar_diff(e,variable).map(MatrixOperand::Scalar),MatrixOperand::Matrix(m)=>m.diff_cached(variable,cache).map(MatrixOperand::Matrix)};
                let db=derive(base)?;let de=derive(exponent)?;
                let is_zero=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>e.is_zero(),MatrixOperand::Matrix(m)=>matches!(m.kind(),MatrixKind::Zero)};
                if is_zero(&db)&&is_zero(&de){zero()?}
                else if is_zero(&de){
                    let previous=match exponent{MatrixOperand::Scalar(e)=>MatrixOperand::Scalar(sub(e,&Expr::one())),MatrixOperand::Matrix(m)=>MatrixOperand::Matrix(Self::add(&[m.clone(),Self::mul(Expr::minus_one(),&[Self::one(m.rows().clone(),m.cols().clone())?],true)?],true)?)};
                    let mut factors=vec![Self::hadamard_power(base.clone(),previous,true)?];let mut coefficient=Expr::one();
                    for operand in [exponent,&db]{match operand{MatrixOperand::Scalar(e)=>coefficient=mul2(&coefficient,e),MatrixOperand::Matrix(m)=>factors.push(m.clone())}}
                    Self::mul(coefficient,&[Self::hadamard_product(&factors,true)?],true)?
                }else if let MatrixOperand::Scalar(a)=base{
                    let MatrixOperand::Matrix(exponent)=exponent else{unreachable!()};let MatrixOperand::Matrix(de)=de else{unreachable!()};let MatrixOperand::Scalar(da)=db else{unreachable!()};
                    let mut terms=vec![Self::mul(func(Func::Log,vec![a.clone()]),&[Self::hadamard_product(&[self.clone(),de],true)?],true)?];
                    if !da.is_zero(){terms.push(Self::mul(crate::ops::div(&da,a),&[Self::hadamard_product(&[self.clone(),exponent.clone()],true)?],true)?);}
                    Self::add(&terms,true)?
                }else{
                    let explicit=self.as_explicit().map_err(|_|MatrixError::NotImplemented("Variable Hadamard exponents with a matrix base require finite dimensions for exact entrywise differentiation".into()))?;
                    Self::explicit(explicit.try_applyfunc::<MatrixError>(|e|Self::scalar_diff(e,variable))?)
                }
            },
            MatrixKind::Mul(coefficient,factors)=>{
                if coefficient.is_zero(){return zero();}
                let mut terms=Vec::new();
                let dc=Self::scalar_diff(coefficient,variable)?;
                if !dc.is_zero(){terms.push(Self::mul(dc,factors,true)?);}
                for (index,factor) in factors.iter().enumerate(){
                    let derivative=factor.diff_cached(variable,cache)?;
                    if matches!(derivative.kind(),MatrixKind::Zero){continue;}
                    let mut product=factors.clone();product[index]=derivative;
                    terms.push(Self::mul(coefficient.clone(),&product,true)?);
                }
                if terms.is_empty(){zero()?}else{Self::add(&terms,true)?}
            },
            MatrixKind::Pow(base,exponent)=>{
                if exponent.is_zero(){return zero();}
                let derivative=base.diff_cached(variable,cache)?;
                let de=Self::scalar_diff(exponent,variable)?;
                if matches!(derivative.kind(),MatrixKind::Zero)&&de.is_zero(){zero()?}
                else{
                    let n=exponent.as_i64().ok_or_else(||MatrixError::NotImplemented("Derivatives of matrix powers require a constant integer exponent".into()))?;
                    if n<0{
                        let inverse=base.inverse(true)?;
                        let dinverse=Self::mul(Expr::minus_one(),&[inverse.clone(),derivative,inverse.clone()],true)?;
                        Self::diff_power(&inverse,&dinverse,n.unsigned_abs(),&mut FxHashMap::default())?
                    }else{Self::diff_power(base,&derivative,n as u64,&mut FxHashMap::default())?}
                }
            },
            MatrixKind::Inverse(base)=>{
                let derivative=base.diff_cached(variable,cache)?;
                if matches!(derivative.kind(),MatrixKind::Zero){zero()?}
                else{let inverse=base.inverse(true)?;Self::mul(Expr::minus_one(),&[inverse.clone(),derivative,inverse],true)?}
            },
            MatrixKind::Transpose(base)=>base.diff_cached(variable,cache)?.transpose(true)?,
            MatrixKind::Adjoint(base)|MatrixKind::Conjugate(base)=>{
                let derivative=base.diff_cached(variable,cache)?;
                if matches!(derivative.kind(),MatrixKind::Zero){zero()?}
                else{
                    if crate::assumptions::is_real(variable)!=Some(true){return Err(MatrixError::NotImplemented("Differentiating a dependent adjoint or conjugate requires a real scalar variable".into()));}
                    if matches!(self.kind(),MatrixKind::Adjoint(_)){derivative.adjoint(true)?}else{derivative.conjugate(true)?}
                }
            },
        };
        cache.insert(self.clone(),result.clone());Ok(result)
    }
    /// A balanced product rule keeps very large integer powers compact.
    fn diff_power(base:&Self,derivative:&Self,n:u64,cache:&mut FxHashMap<u64,Self>)->MResult<Self>{
        if n==0{return Self::zero(base.rows().clone(),base.cols().clone());}
        if n==1{return Ok(derivative.clone());}
        if let Some(result)=cache.get(&n){return Ok(result.clone());}
        let left=n/2;let right=n-left;
        let exponent=|n:u64|Expr::num(crate::number::Number::from_bigint(num_bigint::BigInt::from(n)));
        let dl=Self::diff_power(base,derivative,left,cache)?;
        let dr=Self::diff_power(base,derivative,right,cache)?;
        let result=Self::add(&[
            Self::mul(Expr::one(),&[dl,base.pow(exponent(right),true)?],true)?,
            Self::mul(Expr::one(),&[base.pow(exponent(left),true)?,dr],true)?,
        ],true)?;
        cache.insert(n,result.clone());Ok(result)
    }
    fn unary(&self,op:u8,evaluate:bool)->MResult<Self>{
        let apply=|m:&Self|m.unary(op,true);
        if evaluate {match self.kind(){
            MatrixKind::Identity=>return Ok(self.clone()),MatrixKind::One=>return if op==2{Ok(self.clone())}else{Self::one(self.cols().clone(),self.rows().clone())},MatrixKind::Zero=>return if op==2{Ok(self.clone())}else{Self::zero(self.cols().clone(),self.rows().clone())},
            MatrixKind::HadamardProduct(factors)=>return Self::hadamard_product(&factors.iter().map(apply).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Kronecker(factors)=>return Self::kronecker(&factors.iter().map(apply).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::HadamardPower(base,exponent) if op==0 || exponent.scalar().is_some_and(|e|crate::assumptions::is_integer(e)==Some(true))=>{
                let transform=|o:&MatrixOperand|->MResult<MatrixOperand>{Ok(match o{MatrixOperand::Matrix(m)=>MatrixOperand::Matrix(apply(m)?),MatrixOperand::Scalar(e)=>MatrixOperand::Scalar(if op==0{e.clone()}else{func(Func::Conjugate,vec![e.clone()])})})};
                return Self::hadamard_power(transform(base)?,transform(exponent)?,true);
            },
            MatrixKind::Explicit(m)=>return Ok(Self::explicit(match op{0=>m.transpose(),1=>m.adjoint(),_=>m.conjugate()})),
            MatrixKind::Transpose(b) if op==0=>return Ok(b.clone()),MatrixKind::Adjoint(b) if op==1=>return Ok(b.clone()),MatrixKind::Conjugate(b) if op==2=>return Ok(b.clone()),
            MatrixKind::Transpose(b) if op==1=>return b.conjugate(true),MatrixKind::Adjoint(b) if op==0=>return b.conjugate(true),
            MatrixKind::Conjugate(b) if op==0=>return b.adjoint(true),MatrixKind::Conjugate(b) if op==1=>return b.transpose(true),
            MatrixKind::Add(args)=>return Self::add(&args.iter().map(apply).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Mul(c,fs)=>{let mut factors=fs.clone();if op!=2{factors.reverse();}let c=if op==0{c.clone()}else{func(Func::Conjugate,vec![c.clone()])};return Self::mul(c,&factors.iter().map(apply).collect::<MResult<Vec<_>>>()?,true);},
            MatrixKind::Pow(b,e)=>return apply(b)?.pow(if op==1{func(Func::Conjugate,vec![e.clone()])}else{e.clone()},true),
            MatrixKind::Inverse(b)=>return apply(b)?.inverse(true),_=>{},
        }}
        let (r,c)=if op==2{self.shape()}else{(self.cols().clone(),self.rows().clone())};
        Ok(Self::new(match op{0=>MatrixKind::Transpose(self.clone()),1=>MatrixKind::Adjoint(self.clone()),_=>MatrixKind::Conjugate(self.clone())},r,c))
    }
    pub fn children(&self)->Vec<Self>{match self.kind(){MatrixKind::Add(a)|MatrixKind::Mul(_,a)|MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)=>a.clone(),MatrixKind::HadamardPower(a,b)=>[a,b].into_iter().filter_map(|o|o.matrix().cloned()).collect(),MatrixKind::Pow(b,_)|MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Inverse(b)|MatrixKind::Conjugate(b)=>vec![b.clone()],_=>Vec::new()}}
    pub fn map_children(&self, matrices:&mut dyn FnMut(&Self)->Self, scalars:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        match self.kind(){
            MatrixKind::Symbol(name)=>{let r=scalars(self.rows());let c=scalars(self.cols());if r==*self.rows()&&c==*self.cols(){Ok(self.clone())}else{Self::symbol(name,r,c)}},
            MatrixKind::Identity=>{let n=scalars(self.rows());if n==*self.rows(){Ok(self.clone())}else{Self::identity(n)}},
            MatrixKind::Zero|MatrixKind::One=>{let r=scalars(self.rows());let c=scalars(self.cols());if r==*self.rows()&&c==*self.cols(){Ok(self.clone())}else if matches!(self.kind(),MatrixKind::One){Self::one(r,c)}else{Self::zero(r,c)}},
            MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)=>{let b:Vec<_>=a.iter().map(matrices).collect();if &b==a{Ok(self.clone())}else if matches!(self.kind(),MatrixKind::Kronecker(_)){Self::kronecker(&b,true)}else{Self::hadamard_product(&b,true)}},
            MatrixKind::HadamardPower(a,b)=>{let aa=a.map(matrices,scalars);let bb=b.map(matrices,scalars);if aa==*a&&bb==*b{Ok(self.clone())}else{Self::hadamard_power(aa,bb,true)}},
            MatrixKind::Explicit(m)=>{let entries:Vec<_>=m.data().iter().map(scalars).collect();if entries==m.data(){Ok(self.clone())}else{Ok(Self::explicit(Matrix::new(m.rows(),m.cols(),entries)?))}},
            MatrixKind::Add(a)=>{let b:Vec<_>=a.iter().map(matrices).collect();if &b==a{Ok(self.clone())}else{Self::add(&b,true)}},
            MatrixKind::Mul(c,a)=>{let coefficient=scalars(c);let b:Vec<_>=a.iter().map(matrices).collect();if &b==a&&coefficient==*c{Ok(self.clone())}else{Self::mul(coefficient,&b,true)}},
            MatrixKind::Pow(b,e)=>{let base=matrices(b);let exponent=scalars(e);if base==*b&&exponent==*e{Ok(self.clone())}else{base.pow(exponent,true)}},
            MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Inverse(b)|MatrixKind::Conjugate(b)=>{let base=matrices(b);if base==*b{return Ok(self.clone());}match self.kind(){MatrixKind::Transpose(_)=>base.transpose(true),MatrixKind::Adjoint(_)=>base.adjoint(true),MatrixKind::Inverse(_)=>base.inverse(true),_=>base.conjugate(true)}},
        }
    }
    pub fn matrix_symbols(&self)->Vec<Self>{let mut out=FxHashSet::default();self.collect_matrices(&mut out);let mut out:Vec<_>=out.into_iter().collect();out.sort();out}
    fn collect_matrices(&self,out:&mut FxHashSet<Self>){if matches!(self.kind(),MatrixKind::Symbol(_)){out.insert(self.clone());}else{for c in self.children(){c.collect_matrices(out);}}}
    pub fn scalar_free_symbols(&self)->Vec<Expr>{
        let mut out=FxHashSet::default();self.collect_scalars(&mut out);let mut out:Vec<_>=out.into_iter().collect();out.sort();out
    }
    fn collect_scalars(&self,out:&mut FxHashSet<Expr>){
        let mut push=|e:&Expr|out.extend(crate::traverse::free_symbols(e));
        match self.kind(){
            MatrixKind::Symbol(_)=>{}, MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>{push(self.rows());push(self.cols());},
            MatrixKind::Explicit(m)=>{for e in m.data(){push(e);}},MatrixKind::Mul(c,_)=>push(c),MatrixKind::Pow(_,e)=>push(e),MatrixKind::HadamardPower(a,b)=>{for operand in [a,b]{if let MatrixOperand::Scalar(e)=operand{push(e);}}},_=>{},
        }
        for c in self.children(){c.collect_scalars(out);}
    }
    pub fn matrix_map(&self,f:&mut dyn FnMut(&Self)->Self)->MResult<Self>{
        let replacement=f(self);if replacement!=*self{return Ok(replacement);}
        let old=self.children();let children=old.iter().map(|c|c.matrix_map(f)).collect::<MResult<Vec<_>>>()?;
        if children==old{return Ok(self.clone());}
        self.map_children(&mut |child|old.iter().position(|c|c==child).map(|i|children[i].clone()).unwrap_or_else(||child.clone()),&mut |e|e.clone())
    }

    fn scalar_fields(&self,out:&mut Vec<Expr>){
        match self.kind(){MatrixKind::Symbol(_)|MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>{out.push(self.rows().clone());out.push(self.cols().clone());},MatrixKind::Explicit(m)=>out.extend(m.data().iter().cloned()),MatrixKind::Mul(c,_)=>out.push(c.clone()),MatrixKind::Pow(_,e)=>out.push(e.clone()),MatrixKind::HadamardPower(a,b)=>{for operand in [a,b]{if let MatrixOperand::Scalar(e)=operand{out.push(e.clone());}}},_=>{}}
        for c in self.children(){c.scalar_fields(out);}
    }
    pub fn scalar_map(&self,f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        let mut fields=Vec::new();self.scalar_fields(&mut fields);
        let values:FxHashMap<_,_>=fields.into_iter().map(|e|{let mapped=f(&e);(e,mapped)}).collect();
        if values.iter().all(|(a,b)|a==b){return Ok(self.clone());}
        self.remap_internal(&[],&mut |e|values.get(e).cloned().unwrap_or_else(||e.clone()))
    }
    pub fn map_all(&self,replacements:&[(Self,Self)],f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        self.remap_internal(replacements,f)
    }
    fn remap_internal(&self, replacements:&[(Self,Self)], f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        if let Some((_,value))=replacements.iter().find(|(key,_)|key==self){return Ok(value.clone());}
        let children=self.children();
        let mapped=children.iter().map(|c|c.remap_internal(replacements,f)).collect::<MResult<Vec<_>>>()?;
        self.map_children(&mut |child| children.iter().position(|c|c==child).map(|i|mapped[i].clone()).unwrap_or_else(||child.clone()), f)
    }

    pub fn xreplace(&self,matrices:&[(Self,Self)],scalars:&FxHashMap<Expr,Expr>)->MResult<Self>{
        let mut rules=scalars.clone();
        for (old,new) in matrices{
            let proxy=|m:&Self|Expr::from_kind(Kind::Func(Func::MatrixProxy(Arc::new(m.clone())),vec![]));
            rules.insert(proxy(old),proxy(new));
        }
        self.map_all(matrices,&mut |e|crate::subs::xreplace(e,&rules))
    }
    pub fn doit(&self,deep:bool)->MResult<Self>{
        let children=if deep{self.children().iter().map(|c|c.doit(true)).collect::<MResult<Vec<_>>>()?}else{self.children()};
        match self.kind(){MatrixKind::HadamardProduct(_)=>Self::hadamard_product(&children,true),MatrixKind::Kronecker(_)=>Self::kronecker(&children,true),MatrixKind::HadamardPower(a,b)=>{let mut index=0;let mut map=|_:&Self|{let child=children[index].clone();index+=1;child};Self::hadamard_power(a.map(&mut map,&mut |e|e.clone()),b.map(&mut map,&mut |e|e.clone()),true)},MatrixKind::Add(_)=>Self::add(&children,true),MatrixKind::Mul(c,_)=>Self::mul(c.clone(),&children,true),MatrixKind::Pow(_,e)=>children[0].pow(e.clone(),true),MatrixKind::Inverse(_)=>children[0].inverse(true),MatrixKind::Transpose(_)=>children[0].transpose(true),MatrixKind::Adjoint(_)=>children[0].adjoint(true),MatrixKind::Conjugate(_)=>children[0].conjugate(true),_=>Ok(self.clone())}
    }
    pub fn expand(&self)->MResult<Self>{
        let node=self.doit(true)?;
        match node.kind(){
            MatrixKind::Add(args)=>Self::add(&args.iter().map(|a|a.expand()).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Mul(c,fs)=>{let mut terms=vec![Vec::new()];for factor in fs{let expanded=factor.expand()?;let choices=if let MatrixKind::Add(a)=expanded.kind(){a.clone()}else{vec![expanded]};let mut next=Vec::new();for prefix in terms{for term in &choices{let mut chain=prefix.clone();chain.push(term.clone());next.push(chain);}}terms=next;}Self::add(&terms.iter().map(|chain|Self::mul(c.clone(),chain,true)).collect::<MResult<Vec<_>>>()?,true)},_=>Ok(node)
        }
    }
    fn validate_index(index:&Expr,dim:&Expr)->MResult<()> {
        if crate::assumptions::is_integer(index)==Some(false){return Err(MatrixError::Index("Matrix indices must be integers".into()));}
        if let(Some(i),Some(n))=(index.as_i64(),dim.as_i64()){if i < -n || i>=n{return Err(MatrixError::Index("Matrix index out of range".into()));}}
        Ok(())
    }
    fn kronecker_element(&self,factors:&[Self],i:&Expr,j:&Expr)->MResult<Expr>{
        if factors.iter().any(|m|crate::assumptions::is_positive(m.rows())!=Some(true)||crate::assumptions::is_positive(m.cols())!=Some(true)){return Ok(self.opaque_element(i,j));}
        let split=|index:&Expr,size:&Expr|{
            if let (Some(i),Some(n))=(index.as_i64(),size.as_i64()){(Expr::int(i.div_euclid(n)),Expr::int(i.rem_euclid(n)))}
            else{let quotient=func(Func::Floor,vec![crate::ops::div(index,size)]);let remainder=sub(index,&mul2(size,&quotient));(quotient,remainder)}
        };
        let mut i=i.clone();let mut j=j.clone();let mut value=Expr::one();
        for factor in factors.iter().rev(){let(quotient_i,row)=split(&i,factor.rows());let(quotient_j,col)=split(&j,factor.cols());value=Self::entry_product(&value,&factor.element(&row,&col)?);i=quotient_i;j=quotient_j;}
        Ok(value)
    }
    pub fn element(&self,i:&Expr,j:&Expr)->MResult<Expr>{
        Self::validate_index(i,self.rows())?;Self::validate_index(j,self.cols())?;
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),
            MatrixKind::One=>Ok(Expr::one()),
            MatrixKind::HadamardProduct(factors)=>{let mut value=Expr::one();for factor in factors{value=Self::entry_product(&value,&factor.element(i,j)?);}Ok(value)},
            MatrixKind::HadamardPower(base,exponent)=>Ok(pow(&base.entry(i,j)?,&exponent.entry(i,j)?)),
            MatrixKind::Kronecker(factors)=>self.kronecker_element(factors,i,j),
            MatrixKind::Identity if i==j=>Ok(Expr::one()),
            MatrixKind::Identity if i.as_i64().is_some()&&j.as_i64().is_some()=>Ok(Expr::zero()),
            MatrixKind::Explicit(m)=>{if let(Some(i),Some(j))=(i.as_i64(),j.as_i64()){return Ok(m.get(i.rem_euclid(m.rows() as i64) as usize,j.rem_euclid(m.cols() as i64) as usize).clone());}Ok(self.opaque_element(i,j))},
            MatrixKind::Add(a)=>Ok(add(&a.iter().map(|e|e.element(i,j)).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Mul(c,fs)=>{if fs.len()==1{return Ok(mul2(c,&fs[0].element(i,j)?));}let first=&fs[0];if let Ok(n)=finite_dim(first.cols()){let rest=Self::mul(Expr::one(),&fs[1..],true)?;let mut terms=Vec::with_capacity(n);for k in 0..n{let k=Expr::int(k as i64);terms.push(mul2(&first.element(i,&k)?,&rest.element(&k,j)?));}Ok(mul2(c,&add(&terms)))}else{Ok(self.opaque_element(i,j))}},
            MatrixKind::Transpose(b)=>b.element(j,i),MatrixKind::Adjoint(b)=>Ok(func(Func::Conjugate,vec![b.element(j,i)?])),MatrixKind::Conjugate(b)=>Ok(func(Func::Conjugate,vec![b.element(i,j)?])),
            MatrixKind::Pow(_,e) if e.as_i64().is_none()=>Ok(self.opaque_element(i,j)),
            MatrixKind::Pow(..)|MatrixKind::Inverse(_)=>{if finite_dim(self.rows()).is_ok(){let m=self.as_explicit()?;if let(Some(i),Some(j))=(i.as_i64(),j.as_i64()){return Ok(m.get(i.rem_euclid(m.rows() as i64) as usize,j.rem_euclid(m.cols() as i64) as usize).clone());}}Ok(self.opaque_element(i,j))},
            _=>Ok(self.opaque_element(i,j)),
        }
    }
    fn opaque_element(&self,i:&Expr,j:&Expr)->Expr{Expr::from_kind(crate::expr::Kind::Func(Func::MatrixElement(Arc::new(self.clone())),vec![i.clone(),j.clone()]))}
    pub fn as_explicit(&self)->MResult<Matrix>{
        let (rows,cols)=(finite_dim(self.rows())?,finite_dim(self.cols())?);
        match self.kind(){
            MatrixKind::Explicit(m)=>Ok(m.clone()),MatrixKind::Identity=>Ok(Matrix::eye(rows)),MatrixKind::Zero=>Ok(Matrix::zeros(rows,cols)),
            MatrixKind::One=>Ok(Matrix::from_fn(rows,cols,|_,_|Expr::one())),
            MatrixKind::HadamardProduct(_)|MatrixKind::HadamardPower(..)|MatrixKind::Kronecker(_)=>{let mut entries=Vec::with_capacity(rows*cols);for i in 0..rows{for j in 0..cols{entries.push(self.element(&Expr::int(i as i64),&Expr::int(j as i64))?);}}Matrix::new(rows,cols,entries)},
            MatrixKind::Add(a)=>{let mut m=Matrix::zeros(rows,cols);for e in a{m=m.add(&e.as_explicit()?)?;}Ok(m)},
            MatrixKind::Mul(c,fs)=>{let mut m=fs[0].as_explicit()?;for e in &fs[1..]{m=m.matmul(&e.as_explicit()?)?;}Ok(m.scale(c))},
            MatrixKind::Pow(b,e)=>{let n=e.as_i64().ok_or_else(||MatrixError::NotImplemented("Explicit matrix powers require an integer exponent".into()))?;b.as_explicit()?.powi(n)},
            MatrixKind::Inverse(b)=>b.as_explicit()?.inv(),MatrixKind::Transpose(b)=>Ok(b.as_explicit()?.transpose()),MatrixKind::Adjoint(b)=>Ok(b.as_explicit()?.adjoint()),MatrixKind::Conjugate(b)=>Ok(b.as_explicit()?.conjugate()),
            MatrixKind::Symbol(_)=>{let mut entries=Vec::with_capacity(rows*cols);for i in 0..rows{for j in 0..cols{entries.push(self.opaque_element(&Expr::int(i as i64),&Expr::int(j as i64)));}}Matrix::new(rows,cols,entries)},
        }
    }
    pub fn kind_name(&self)->&'static str{match self.kind(){MatrixKind::Symbol(_)=>"MatrixSymbol",MatrixKind::Explicit(_)=>"ImmutableMatrix",MatrixKind::Identity=>"Identity",MatrixKind::Zero=>"ZeroMatrix",MatrixKind::One=>"OneMatrix",MatrixKind::HadamardProduct(_)=>"HadamardProduct",MatrixKind::HadamardPower(..)=>"HadamardPower",MatrixKind::Kronecker(_)=>"KroneckerProduct",MatrixKind::Add(_)=>"MatAdd",MatrixKind::Mul(..)=>"MatMul",MatrixKind::Pow(..)=>"MatPow",MatrixKind::Transpose(_)=>"Transpose",MatrixKind::Adjoint(_)=>"Adjoint",MatrixKind::Inverse(_)=>"Inverse",MatrixKind::Conjugate(_)=>"MatrixConjugate"}}
    pub fn to_repr(&self)->String{
        let s=crate::printer::srepr;let me=|e:&Self|e.to_repr();
        let args=match self.kind(){
            MatrixKind::Symbol(name)=>format!("{name:?}, {}, {}",s(self.rows()),s(self.cols())),MatrixKind::Identity=>s(self.rows()),MatrixKind::Zero|MatrixKind::One=>format!("{}, {}",s(self.rows()),s(self.cols())),
            MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)=>a.iter().map(me).collect::<Vec<_>>().join(", "),
            MatrixKind::HadamardPower(a,b)=>{let operand=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>s(e),MatrixOperand::Matrix(m)=>m.to_repr()};format!("{}, {}",operand(a),operand(b))},
            MatrixKind::Explicit(m)=>format!("{}, {}, [{}]",m.rows(),m.cols(),m.data().iter().map(s).collect::<Vec<_>>().join(", ")),
            MatrixKind::Add(a)=>a.iter().map(me).collect::<Vec<_>>().join(", "),MatrixKind::Mul(c,a)=>{let mut a=a.iter().map(me).collect::<Vec<_>>();if !c.is_one(){a.insert(0,s(c));}a.join(", ")},
            MatrixKind::Pow(b,e)=>format!("{}, {}",b.to_repr(),s(e)),_=>self.children()[0].to_repr(),
        };format!("{}({args})",self.kind_name())
    }
    pub fn to_str(&self)->String{
        let child=|e:&Self|if matches!(e.kind(),MatrixKind::Add(_)){format!("({})",e.to_str())}else{e.to_str()};
        match self.kind(){MatrixKind::Symbol(n)=>n.clone(),MatrixKind::Explicit(m)=>m.to_str(),MatrixKind::Identity=>"I".into(),MatrixKind::Zero=>"0".into(),MatrixKind::One=>"1".into(),
            MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)=>format!("{}({})",self.kind_name(),a.iter().map(Self::to_str).collect::<Vec<_>>().join(", ")),
            MatrixKind::HadamardPower(a,b)=>{let operand=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>e.to_string(),MatrixOperand::Matrix(m)=>m.to_str()};format!("HadamardPower({}, {})",operand(a),operand(b))},
            MatrixKind::Add(a)=>{let mut out=String::new();for e in a{let v=e.to_str();if !out.is_empty(){if let Some(positive)=v.strip_prefix('-'){out.push_str(" - ");out.push_str(positive);continue;}out.push_str(" + ");}out.push_str(&v);}out},
            MatrixKind::Mul(c,a)=>{let mut terms:Vec<_>=a.iter().map(child).collect();if !c.is_one(){let scalar=if c.is_add(){format!("({c})")}else{c.to_string()};if c==&Expr::minus_one(){return format!("-{}",terms.join("*"));}terms.insert(0,scalar);}terms.join("*")},
            MatrixKind::Pow(b,e)=>format!("{}**{}",child(b),if e.is_atom()&&!e.as_number().is_some_and(|n|n.is_negative()){e.to_string()}else{format!("({e})")}),MatrixKind::Inverse(b)=>format!("{}**(-1)",child(b)),MatrixKind::Transpose(b)=>format!("{}.T",child(b)),MatrixKind::Adjoint(b)=>format!("Adjoint({})",b.to_str()),MatrixKind::Conjugate(b)=>format!("Adjoint({}.T)",child(b)),
        }
    }
}

impl std::fmt::Display for MatrixExpr { fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{f.write_str(&self.to_str())} }

#[cfg(test)]
mod tests {
    use super::*;
    fn symbol(name:&str,rows:i64,cols:i64)->MatrixExpr {MatrixExpr::symbol(name,Expr::int(rows),Expr::int(cols)).unwrap()}
    #[test]
    fn noncommutative_canonicalization() {
        let (a,b)=(symbol("A",2,2),symbol("B",2,2));
        assert_ne!(MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap(),MatrixExpr::mul(Expr::one(),&[b.clone(),a.clone()],true).unwrap());
        assert_eq!(MatrixExpr::add(&[a.clone(),a.clone()],true).unwrap(),MatrixExpr::mul(Expr::int(2),&[a.clone()],true).unwrap());
        assert_eq!(MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone(),b.inverse(true).unwrap()],true).unwrap(),a);
        let p=MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap();
        assert_eq!(p.transpose(true).unwrap(),MatrixExpr::mul(Expr::one(),&[b.transpose(true).unwrap(),a.transpose(true).unwrap()],true).unwrap());
        assert!(MatrixExpr::mul(Expr::one(),&[a.clone(),symbol("C",3,2)],true).is_err());
    }
    #[test]
    fn finite_entries_and_explicit_matrices() {
        let a=symbol("A",2,3);let b=symbol("B",3,2);
        let product=MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap();
        let value=product.element(&Expr::zero(),&Expr::one()).unwrap();
        let expected=add(&(0..3).map(|k|mul2(&a.element(&Expr::zero(),&Expr::int(k)).unwrap(),&b.element(&Expr::int(k),&Expr::one()).unwrap())).collect::<Vec<_>>());
        assert_eq!(value,expected);
        assert_eq!(product.as_explicit().unwrap(),a.as_explicit().unwrap().matmul(&b.as_explicit().unwrap()).unwrap());
        let vars=crate::traverse::free_symbols(&value);assert_eq!(vars.len(),2);
        assert!(vars.iter().all(|e|matches!(e.kind(),crate::expr::Kind::Func(Func::MatrixProxy(_),_))));
    }
    #[test]
    fn dimension_substitutions_and_single_pass_mapping() {
        let n=Expr::symbol("n");let m=Expr::symbol("m");let p=Expr::symbol("p");
        let a=MatrixExpr::symbol("A",n.clone(),n.clone()).unwrap();
        assert!(a.scalar_free_symbols().is_empty());assert_eq!(a.matrix_symbols(),vec![a.clone()]);
        let map=FxHashMap::from_iter([(n.clone(),m.clone()),(m.clone(),p)]);
        let actual=a.xreplace(&[],&map).unwrap();assert_eq!(actual.rows(),&m);
        assert!(a.as_explicit().is_err());
        let raw=MatrixExpr::mul(Expr::one(),&[a.clone(),MatrixExpr::identity(n).unwrap()],false).unwrap();
        assert_eq!(raw.scalar_map(&mut |e|e.clone()).unwrap(),raw);
        assert_eq!(raw.map_children(&mut |e|e.clone(),&mut |e|e.clone()).unwrap(),raw);
    }
    #[test]
    fn ordered_scalar_derivatives_and_inverse() {
        let x=Expr::symbol("x");let a=symbol("A",2,2);let b=symbol("B",2,2);
        let xa=MatrixExpr::mul(x.clone(),&[a.clone()],true).unwrap();
        let base=MatrixExpr::add(&[xa.clone(),b.clone()],true).unwrap();
        let derivative=base.pow(Expr::int(3),true).unwrap().diff(&x).unwrap().expand().unwrap();
        let expected=MatrixExpr::add(&[
            MatrixExpr::mul(Expr::one(),&[a.clone(),base.pow(Expr::int(2),true).unwrap()],true).unwrap(),
            MatrixExpr::mul(Expr::one(),&[base.clone(),a.clone(),base.clone()],true).unwrap(),
            MatrixExpr::mul(Expr::one(),&[base.pow(Expr::int(2),true).unwrap(),a.clone()],true).unwrap(),
        ],true).unwrap().expand().unwrap();
        let actual=derivative.as_explicit().unwrap();let wanted=expected.as_explicit().unwrap();
        for (a,b) in actual.data().iter().zip(wanted.data()){assert!(crate::expand::expand(&sub(a,b)).is_zero());}
        let inverse=base.inverse(true).unwrap();
        assert_eq!(inverse.diff(&x).unwrap(),MatrixExpr::mul(Expr::minus_one(),&[inverse.clone(),a.clone(),inverse],true).unwrap());
        assert!(matches!(a.diff(&x).unwrap().kind(),MatrixKind::Zero));
        assert!(base.pow(x.clone(),true).unwrap().diff(&x).is_err());
        assert!(base.adjoint(false).unwrap().diff(&x).is_err());
        assert!(xa.adjoint(true).unwrap().diff(&x).is_err());
    }
    #[test]
    fn abstract_trace_determinant_and_explicit_certificates(){
        let a=symbol("A",2,3);let b=symbol("B",3,2);
        let ab=MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap();
        let ba=MatrixExpr::mul(Expr::one(),&[b,a],true).unwrap();
        assert_eq!(ab.trace().unwrap(),ba.trace().unwrap());
        assert!(matches!(ab.determinant().unwrap().kind(),Kind::Func(Func::MatrixDeterminant(_),_)));
        assert_eq!(MatrixExpr::identity(Expr::int(0)).unwrap().determinant().unwrap(),Expr::one());
        assert_eq!(MatrixExpr::zero(Expr::int(0),Expr::int(0)).unwrap().determinant().unwrap(),Expr::one());
        let a=symbol("C",2,2);let b=symbol("D",2,2);
        let product=MatrixExpr::mul(Expr::int(3),&[a.clone(),b.clone()],true).unwrap();
        assert_eq!(product.determinant().unwrap(),crate::ops::mul(&[Expr::int(9),a.determinant().unwrap(),b.determinant().unwrap()]));
        assert_eq!(a.trace_explicit().unwrap(),a.as_explicit().unwrap().trace().unwrap());
        assert_eq!(a.determinant_explicit().unwrap(),a.as_explicit().unwrap().det().unwrap());
    }
    #[test]
    fn large_integer_power_derivatives_share_subexpressions(){
        let x=Expr::symbol("x");let a=symbol("A",2,2);let b=symbol("B",2,2);
        let base=MatrixExpr::add(&[MatrixExpr::mul(x.clone(),&[a],true).unwrap(),b],true).unwrap();
        let result=base.pow(Expr::int(1_000_000_000),true).unwrap().diff(&x).unwrap();
        assert_eq!(result.shape(),(Expr::int(2),Expr::int(2)));
    }
    #[test]
    fn hadamard_operands_canonicalization_and_entries(){
        let a=symbol("A",2,3);let b=symbol("B",2,3);
        let product=MatrixExpr::hadamard_product(&[b.clone(),a.clone(),a.clone()],true).unwrap();
        let reversed=MatrixExpr::hadamard_product(&[a.clone(),b.clone(),a.clone()],true).unwrap();assert_eq!(product,reversed);
        assert_eq!(product.element(&Expr::zero(),&Expr::one()).unwrap(),mul2(&a.element(&Expr::zero(),&Expr::one()).unwrap().powi(2),&b.element(&Expr::zero(),&Expr::one()).unwrap()));
        assert_eq!(MatrixExpr::hadamard_product(&[a.clone(),MatrixExpr::one(Expr::int(2),Expr::int(3)).unwrap()],true).unwrap(),a);
        let power=MatrixExpr::hadamard_power(MatrixOperand::Scalar(Expr::int(2)),MatrixOperand::Matrix(b),false).unwrap();assert_eq!(power.shape(),(Expr::int(2),Expr::int(3)));
        assert_eq!(power.scalar_map(&mut |e|e.clone()).unwrap(),power);
        assert!(MatrixExpr::hadamard_product(&[a,symbol("C",3,2)],false).is_err());
    }
    #[test]
    fn rectangular_kronecker_and_scalar_derivative(){
        let a=symbol("A",2,3);let b=symbol("B",3,2);let x=Expr::symbol("x");
        let product=MatrixExpr::kronecker(&[a.clone(),b.clone()],true).unwrap();assert_eq!(product.shape(),(Expr::int(6),Expr::int(6)));
        assert_eq!(product.element(&Expr::zero(),&Expr::int(3)).unwrap(),mul2(&a.element(&Expr::zero(),&Expr::one()).unwrap(),&b.element(&Expr::zero(),&Expr::one()).unwrap()));
        let scaled=MatrixExpr::kronecker(&[MatrixExpr::mul(x.clone(),&[a],true).unwrap(),b],false).unwrap();assert_eq!(scaled.diff(&x).unwrap(),product);
        assert!(matches!(product.determinant().unwrap().kind(),Kind::Func(Func::MatrixDeterminant(_),_)));
        assert!(product.as_explicit().is_ok());
    }
    #[test]
    fn elementwise_nonfinite_products_keep_explicit_zero_terms(){
        let zero=MatrixExpr::zero(Expr::one(),Expr::one()).unwrap();let infinite=MatrixExpr::explicit(Matrix::new(1,1,vec![Expr::oo()]).unwrap());
        for factors in [[zero.clone(),infinite.clone()],[infinite,zero]]{
            assert_eq!(MatrixExpr::hadamard_product(&factors,true).unwrap().as_explicit().unwrap().get(0,0),&Expr::nan());
            assert_eq!(MatrixExpr::kronecker(&factors,true).unwrap().as_explicit().unwrap().get(0,0),&Expr::nan());
        }
    }
}
