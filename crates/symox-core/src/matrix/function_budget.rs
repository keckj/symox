//! Local allocation/work bounds for matrix functional calculus.
use rustc_hash::{FxHashMap,FxHashSet};
use num_bigint::{BigInt,BigUint};
use num_traits::{One,Zero,ToPrimitive};
use crate::{Expr,Number};
use crate::expr::{Kind,Func};
use super::{Matrix,MatrixError,MResult};

const MAX_BITS:u64=1_000_000;
const MAX_STORAGE:u64=64_000_000;
const MAX_WORK:u64=5_000_000;
const MAX_NODES:usize=100_000;
fn limit()->MatrixError{MatrixError::NotImplemented("Exact matrix function resource limit exceeded".into())}
fn number_bits(n:&Number)->u64{match n{
    Number::Int(n)=>u64::from(64-n.unsigned_abs().leading_zeros()),Number::BigInt(n)=>n.bits(),
    Number::Rational(n)=>n.numer().bits().saturating_add(n.denom().bits()),
    Number::Float(_)=>64,Number::MpFloat(_,p)=>u64::from(*p),
}}
fn rational_profile(m:&Matrix)->Option<(u64,u64)>{
    let mut numerator=0;let mut denominator=0u64;let mut seen:FxHashSet<&BigInt>=FxHashSet::default();
    for e in m.data(){match e.as_number()?{
        Number::Int(n)=>numerator=numerator.max(u64::from(64-n.unsigned_abs().leading_zeros())),
        Number::BigInt(n)=>numerator=numerator.max(n.bits()),
        Number::Rational(n)=>{numerator=numerator.max(n.numer().bits());if seen.insert(n.denom()){denominator=denominator.saturating_add(n.denom().bits());}},
        _=>return None,
    }}Some((numerator,denominator))
}

pub(in crate::matrix) struct Budget{
    work:u64,storage:u64,nodes:FxHashMap<usize,(Expr,u64)>,bits:u64,
}
impl Budget{
    pub(in crate::matrix) fn new()->Self{Self{work:MAX_WORK,storage:MAX_STORAGE,nodes:FxHashMap::default(),bits:MAX_BITS}}
    pub(in crate::matrix) fn work(&mut self,n:u64)->MResult<()>{self.work=self.work.checked_sub(n).ok_or_else(limit)?;Ok(())}
    fn space(&self,n:u64)->MResult<()>{if n>self.storage{Err(limit())}else{Ok(())}}
    pub(in crate::matrix) fn allocate(&mut self,cells:usize)->MResult<()>{let bits=(cells as u64).checked_mul(512).ok_or_else(limit)?;self.space(bits)?;self.storage-=bits;Ok(())}
    pub(in crate::matrix) fn expression(&mut self,e:&Expr)->MResult<u64>{
        if let Some((_,bits))=self.nodes.get(&e.dag_id()){return Ok(*bits);}
        self.scalar_preflight(e,None)?;self.expression_inner(e)
    }
    fn expression_inner(&mut self,e:&Expr)->MResult<u64>{
        if let Some((_,bits))=self.nodes.get(&e.dag_id()){return Ok(*bits);}
        if self.nodes.len()>=MAX_NODES{return Err(limit());}self.work(1)?;
        let own=match e.kind(){
            Kind::Num(n)=>number_bits(n),Kind::Add(c,terms)=>terms.iter().map(|(_,n)|number_bits(n)).chain([number_bits(c)]).max().unwrap_or(0),
            Kind::Mul(c,_)=>number_bits(c),_=>0,
        };
        if own>self.bits{return Err(limit());}
        let payload=match e.kind(){Kind::Add(c,terms)=>terms.iter().fold(number_bits(c),|sum,(_,n)|sum.saturating_add(number_bits(n))),_=>own};
        self.space(payload.saturating_add(256))?;self.storage-=payload+256;
        let mut height=own;
        for child in raw_children(e){height=height.max(self.expression_inner(&child)?);}
        self.nodes.insert(e.dag_id(),(e.clone(),height));Ok(height)
    }
    /// Inspect raw DAGs before substitution can force held numeric powers.
    pub(in crate::matrix) fn scalar_preflight(&mut self,e:&Expr,replacement:Option<(&Expr,&Expr)>)->MResult<()> {
        let mut memo=FxHashMap::default();inspect(e,replacement.map(|(a,b)|Replacement::One(a,b)),self.bits,&mut memo)?;
        self.work(memo.len() as u64)
    }
    /// Bound simultaneous substitutions before forcing numeric arithmetic.
    pub(in crate::matrix) fn substitutions_preflight(&mut self,e:&Expr,replacements:&FxHashMap<Expr,Expr>)->MResult<()> {
        let mut memo=FxHashMap::default();inspect(e,Some(Replacement::Many(replacements)),self.bits,&mut memo)?;
        self.work(memo.len() as u64)
    }
    pub(in crate::matrix) fn matrix(&mut self,m:&Matrix)->MResult<u64>{let mut height=0;for e in m.data(){height=height.max(self.expression(e)?);}Ok(height)}
    pub(in crate::matrix) fn zeros(&mut self,n:usize)->MResult<Matrix>{self.allocate(n.checked_mul(n).ok_or_else(limit)?)?;let m=Matrix::zeros(n,n);self.matrix(&m)?;Ok(m)}
    pub(in crate::matrix) fn identity(&mut self,n:usize)->MResult<Matrix>{self.allocate(n.checked_mul(n).ok_or_else(limit)?)?;let m=Matrix::eye(n);self.matrix(&m)?;Ok(m)}
    pub(in crate::matrix) fn diagonal(&mut self,values:&[Expr])->MResult<Matrix>{self.allocate(values.len().checked_mul(values.len()).ok_or_else(limit)?)?;let m=Matrix::diag(values);self.matrix(&m)?;Ok(m)}
    fn preflight(&mut self,a:&Matrix,b:&Matrix,inner:usize,cells:usize)->MResult<()>{
        self.work((cells as u64).checked_mul((inner as u64).max(1)).and_then(|n|n.checked_mul(2)).ok_or_else(limit)?)?;
        let log=if inner<=1{0}else{u64::from(usize::BITS-(inner-1).leading_zeros())};
        let height=if let(Some((an,ad)),Some((bn,bd)))=(rational_profile(a),rational_profile(b)){
            an.saturating_add(bn).saturating_add(ad.saturating_add(bd).saturating_mul(2)).saturating_add(log+2)
        }else{
            self.matrix(a)?.saturating_add(self.matrix(b)?).saturating_add(log+2).saturating_mul(inner.max(1) as u64)
        };
        if height>self.bits{return Err(limit());}
        self.space(height.saturating_add(512).saturating_mul(cells as u64))?;Ok(())
    }
    pub(in crate::matrix) fn mul(&mut self,a:&Matrix,b:&Matrix)->MResult<Matrix>{
        if a.cols()!=b.rows(){return Err(MatrixError::Shape("Matrix function product shape mismatch".into()));}
        if identity(a){return Ok(b.clone());}
        if identity(b){return Ok(a.clone());}
        if a.is_zero_matrix()||b.is_zero_matrix(){self.allocate(a.rows().checked_mul(b.cols()).ok_or_else(limit)?)?;let m=Matrix::zeros(a.rows(),b.cols());self.matrix(&m)?;return Ok(m);}
        let cells=a.rows().checked_mul(b.cols()).ok_or_else(limit)?;self.preflight(a,b,a.cols(),cells)?;self.allocate(cells)?;
        let m=a.matmul(b)?;self.matrix(&m)?;Ok(m)
    }
    pub(in crate::matrix) fn add(&mut self,a:&Matrix,b:&Matrix)->MResult<Matrix>{
        if a.is_zero_matrix(){return Ok(b.clone());}
        if b.is_zero_matrix(){return Ok(a.clone());}
        self.preflight(a,b,1,a.data().len())?;self.allocate(a.data().len())?;let m=a.add(b)?;self.matrix(&m)?;Ok(m)
    }
    pub(in crate::matrix) fn scale(&mut self,a:&Matrix,c:&Expr)->MResult<Matrix>{
        self.expression(c)?;if c.is_one(){return Ok(a.clone());}
        let scalar=Matrix::new(1,1,vec![c.clone()])?;self.preflight(a,&scalar,1,a.data().len())?;self.allocate(a.data().len())?;
        let m=a.scale(c);self.matrix(&m)?;Ok(m)
    }
    pub(in crate::matrix) fn power(&mut self,m:&Matrix,mut exponent:BigUint)->MResult<Matrix>{
        if exponent.is_zero(){return self.identity(m.rows());}
        if exponent.is_one(){return Ok(m.clone());}
        if identity(m)||m.is_zero_matrix(){return Ok(m.clone());}
        let mut base=m.clone();let mut answer=self.identity(m.rows())?;
        while !exponent.is_zero(){
            self.work(1)?;
            if exponent.bit(0){answer=self.mul(&answer,&base)?;}
            exponent>>=1;
            if !exponent.is_zero(){base=self.mul(&base,&base)?;
                if base.is_zero_matrix(){return self.zeros(m.rows());}
                if identity(&base){return Ok(answer);}
            }
        }Ok(answer)
    }
    pub(in crate::matrix) fn scalar_power(&mut self,base:&Expr,exponent:&Expr)->MResult<Expr>{
        let held=Expr::from_unevaluated_kind(Kind::Pow(base.clone(),exponent.clone()));
        self.scalar_preflight(&held,None)?;
        let value=crate::ops::pow(base,exponent);self.expression(&value)?;Ok(value)
    }
    pub(in crate::matrix) fn cancel_preflight(&mut self,m:&Matrix)->MResult<()>{
        self.matrix(m)?;let mut memo=FxHashMap::default();
        for e in m.data(){
            let terms=expansion_terms(e,&mut memo)?;self.work(terms)?;
            // A polynomial monomial stores an exponent for every variable.
            // All cached scalar nodes bound the number of possible generators,
            // including opaque function atoms; width alone does not bound it.
            self.space(terms.saturating_mul(self.nodes.len() as u64).saturating_mul(32))?;
        }Ok(())
    }
    pub(in crate::matrix) fn cancel(&mut self,m:&Matrix)->MResult<Matrix>{
        self.cancel_preflight(m)?;let result=m.cancel();self.matrix(&result)?;Ok(result)
    }
    /// The caller will expand/cancel this product. Bound its joint width
    /// before scalar multiplication or rational polynomial conversion.
    pub(in crate::matrix) fn product_preflight(&mut self,a:&Matrix,b:&Matrix)->MResult<()>{
        if a.cols()!=b.rows(){return Err(MatrixError::Shape("Matrix function product shape mismatch".into()));}
        self.matrix(a)?;self.matrix(b)?;
        let mut memo=FxHashMap::default();
        let mut left=0;let mut right=0;
        for e in a.data(){if !e.is_zero(){left=left.max(expansion_terms(e,&mut memo)?);}}
        for e in b.data(){if !e.is_zero(){right=right.max(expansion_terms(e,&mut memo)?);}}
        let terms=left.saturating_mul(right).saturating_mul(a.cols() as u64);
        if terms>100_000{return Err(limit());}
        let cells=(a.rows() as u64).saturating_mul(b.cols() as u64);
        self.work(terms.saturating_mul(cells))?;
        self.space(terms.saturating_mul(self.nodes.len() as u64).saturating_mul(32))
    }
    /// Eigenvalue routines form characteristic polynomials per connected
    /// block. Bound even its intermediate products by all ordered words of
    /// length n over the scalar entry monomials (plus the identity).
    /// Numeric blocks have coefficient growth but no polynomial monomial
    /// expansion, already bounded by interpolation().
    pub(in crate::matrix) fn spectral_preflight(&mut self,m:&Matrix)->MResult<()>{
        self.cancel_preflight(m)?;
        let mut memo=FxHashMap::default();
        for block in super::super::linalg::scc(m){
            let mut alphabet=1u64;let mut numeric=true;
            for &i in &block{for &j in &block{let e=m.get(i,j);if !e.is_zero(){alphabet=alphabet.saturating_add(expansion_terms(e,&mut memo)?);}numeric &= e.as_number().is_some();}}
            let n=block.len() as u64;
            self.work(n.saturating_pow(4))?;
            if numeric{continue;}
            let mut terms=1u64;
            for _ in 0..n{terms=terms.saturating_mul(alphabet);if terms>100_000{return Err(limit());}}
            self.work(terms.saturating_mul(n).saturating_mul(n))?;
            self.space(terms.saturating_mul(self.nodes.len() as u64).saturating_mul(32))?;
        }Ok(())
    }
    pub(in crate::matrix) fn interpolation(&mut self,n:usize,height:u64)->MResult<()>{
        let n=n as u64;self.work(n.checked_mul(n).and_then(|v|v.checked_mul(n)).ok_or_else(limit)?)?;
        let log=u64::from(64-n.max(1).leading_zeros());
        let predicted=n.saturating_mul(height.saturating_add(log+2)).saturating_mul(n+1);
        if predicted>self.bits{return Err(limit());}
        self.space(predicted.saturating_add(512).saturating_mul(n.saturating_mul(n)))
    }
}
fn identity(m:&Matrix)->bool{m.is_square()&&(0..m.rows()).all(|i|(0..m.cols()).all(|j|if i==j{m.get(i,j).is_one()}else{m.get(i,j).is_zero()}))}
pub(in crate::matrix) fn exponent(e:&Expr)->Option<BigUint>{match e.as_number()?{Number::Int(n) if *n>=0=>Some(BigUint::from(*n as u64)),Number::BigInt(n)=>n.to_biguint(),_=>None}}

/// Read stored children without map_children: that reconstruction helper can
/// canonicalize held numeric factors even when its callback returns them.
pub(in crate::matrix) fn raw_children(e:&Expr)->Vec<Expr>{
    let mut out=Vec::new();
    match e.kind(){
        Kind::Add(_,terms)=>out.extend(terms.iter().map(|(e,_)|e.clone())),
        Kind::Mul(_,factors)=>{for(b,p)in factors{out.push(b.clone());out.push(p.clone());}},
        Kind::Pow(b,p)|Kind::Rel(_,b,p)=>{out.push(b.clone());out.push(p.clone());},
        Kind::Func(f,args)=>{
            out.extend(args.iter().cloned());
            match f{
                Func::MatrixElement(p)|Func::MatrixProxy(p)|Func::MatrixTrace(p)|Func::MatrixDeterminant(p)=>p.scalar_fields(&mut out),
                Func::MatrixDotProduct(a,b)=>{a.scalar_fields(&mut out);b.scalar_fields(&mut out);},_=>{},
            }
        },
        Kind::Derivative(body,vars)=>{out.push(body.clone());out.extend(vars.iter().map(|(e,_)|e.clone()));},
        Kind::Piecewise(pieces)=>{for(a,b)in pieces{out.push(a.clone());out.push(b.clone());}},_=>{},
    }out
}

// cancel() converts algebraic expressions to expanded polynomial fractions.
// Bound that conversion before entering it, including held symbolic powers.
fn expansion_terms(e:&Expr,memo:&mut FxHashMap<usize,(Expr,u64)>)->MResult<u64>{
    if let Some((_,count))=memo.get(&e.dag_id()){return Ok(*count);}
    let count=match e.kind(){
        Kind::Add(c,terms)=>{let mut n=u64::from(!c.is_zero());for(e,_)in terms{n=n.saturating_add(expansion_terms(e,memo)?);}n},
        Kind::Mul(_,factors)=>{let mut n=1u64;for(b,p)in factors{n=n.saturating_mul(expanded_power(b,p,memo)?);if n>100_000{return Err(limit());}}n},
        Kind::Pow(b,p)=>expanded_power(b,p,memo)?,_=>1,
    };
    if count>100_000{return Err(limit());}memo.insert(e.dag_id(),(e.clone(),count));Ok(count)
}
fn expanded_power(b:&Expr,p:&Expr,memo:&mut FxHashMap<usize,(Expr,u64)>)->MResult<u64>{
    let count=expansion_terms(b,memo)?;
    if count<=1{return Ok(1);}
    match p.as_number().and_then(Number::to_bigint){
        Some(n)=>{let n=n.magnitude().to_u32().ok_or_else(limit)?;let value=count.checked_pow(n).ok_or_else(limit)?;if value>100_000{Err(limit())}else{Ok(value)}},
        None=>Ok(1),
    }
}

#[derive(Clone,Copy)]
struct Scalar { numeric:bool,bits:u64,small:Option<i64> }
impl Scalar {
    fn number(n:&Number)->Self{Self{numeric:true,bits:if n.is_zero()||n.is_one()||n.is_minus_one(){0}else{number_bits(n)},small:if n.is_zero(){Some(0)}else if n.is_one(){Some(1)}else if n.is_minus_one(){Some(-1)}else{n.to_i64()}}}
    fn unknown()->Self{Self{numeric:false,bits:0,small:None}}
}
#[derive(Clone,Copy)]
enum Replacement<'a>{One(&'a Expr,&'a Expr),Many(&'a FxHashMap<Expr,Expr>)}
impl<'a> Replacement<'a>{
    fn get(self,e:&Expr)->Option<&'a Expr>{match self{
        Self::One(key,value)=>(e==key).then_some(value),Self::Many(values)=>values.get(e),
    }}
}
fn inspect(e:&Expr,replacement:Option<Replacement<'_>>,max:u64,memo:&mut FxHashMap<usize,(Expr,Scalar)>)->MResult<Scalar>{
    if let Some((_,s))=memo.get(&e.dag_id()){return Ok(*s);}
    if memo.len()>=MAX_NODES{return Err(limit());}
    if let Some(point)=replacement.and_then(|rule|rule.get(e)){return inspect(point,None,max,&mut FxHashMap::default());}
    let mut child=|e:&Expr|inspect(e,replacement,max,memo);
    let value=match e.kind(){
        Kind::Num(n)=>Scalar::number(n),
        Kind::Const(_)=>Scalar{numeric:true,bits:1,small:None},
        Kind::Sym(_)=>Scalar::unknown(),
        Kind::Add(c,terms)=>{
            let mut v=Scalar::number(c);
            for(e,c)in terms{let a=child(e)?;let b=Scalar::number(c);
                v.numeric &= a.numeric;
                v.bits=v.bits.saturating_add(a.bits).saturating_add(b.bits).saturating_add(1);
                v.small=v.small.zip(a.small.zip(b.small)).and_then(|(v,(a,b))|a.checked_mul(b).and_then(|p|v.checked_add(p)));
            }v
        },
        Kind::Mul(c,factors)=>{
            let mut v=Scalar::number(c);
            for(b,p)in factors{let b=child(b)?;let exponent=child(p)?;let a=inspect_power(b,exponent,p,max)?;
                v.numeric &= a.numeric;v.bits=v.bits.saturating_add(a.bits);
                v.small=v.small.zip(a.small).and_then(|(a,b)|a.checked_mul(b));
            }v
        },
        Kind::Pow(b,p)=>{let b=child(b)?;let exponent=child(p)?;inspect_power(b,exponent,p,max)?},
        Kind::Func(f,args)=>{
            let mut v=Scalar{numeric:true,bits:1,small:None};
            for a in raw_children(e){let a=child(&a)?;v.numeric &= a.numeric;v.bits=v.bits.max(a.bits);}
            // Exact discrete special functions can allocate factorial-sized
            // integers when substitution turns an argument into a number.
            if matches!(f,Func::Gamma|Func::Special(_)){
                for a in args{let s=child(a)?;if s.numeric && s.small.is_none_or(|n|n.unsigned_abs()>10_000){
                    let number=replacement.and_then(|rule|rule.get(a)).unwrap_or(a).as_number();
                    if number.is_none_or(|n|!n.to_f64().is_finite()||n.to_f64().abs()>10_000.0){return Err(limit());}
                }}
            }v
        },
        _=>{for c in raw_children(e){child(&c)?;}Scalar::unknown()}

    };
    if value.bits>max{return Err(limit());}
    memo.insert(e.dag_id(),(e.clone(),value));Ok(value)
}
fn inspect_power(base:Scalar,exponent:Scalar,raw_exponent:&Expr,max:u64)->MResult<Scalar>{
    let integer=match raw_exponent.as_number(){Some(Number::Int(n))=>Some(BigUint::from(n.unsigned_abs())),Some(Number::BigInt(n))=>Some(n.magnitude().clone()),_=>exponent.small.map(|n|BigUint::from(n.unsigned_abs()))};
    let numeric=base.numeric&&exponent.numeric;
    let unit=matches!(base.small,Some(-1..=1));
    let mut bits=base.bits.saturating_add(exponent.bits);
    if !unit {
        if let Some(n)=integer.as_ref(){bits=if base.bits==0{0}else{base.bits.saturating_mul(n.to_u64().unwrap_or(u64::MAX))};if bits>max{return Err(limit());}}
        else if base.numeric && exponent.numeric && exponent.bits>64{return Err(limit());}
    }
    let small=match(base.small,exponent.small){
        (Some(0),Some(p))if p>0=>Some(0),
        (Some(1),_)=>Some(1),
        (Some(-1),Some(p))=>Some(if p%2==0{1}else{-1}),
        (_,Some(0))=>Some(1),
        (Some(b),Some(p))if (0..=63).contains(&p)=>b.checked_pow(p as u32),_=>None,
    };
    Ok(Scalar{numeric,bits:if unit{0}else{bits},small})
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn preflight_power_and_dag_growth_with_tiny_test_limits(){
        let mut budget=Budget::new();budget.bits=64;
        let a=Matrix::from_rows(vec![vec![Expr::int(2),Expr::one()],vec![Expr::zero(),Expr::int(2)]]).unwrap();
        assert!(matches!(budget.power(&a,BigUint::one()<<50),Err(MatrixError::NotImplemented(_))));
        let mut budget=Budget::new();budget.bits=64;let mut a=a;
        let mut rejected=false;for _ in 0..20{match budget.mul(&a,&a){Ok(next)=>a=next,Err(MatrixError::NotImplemented(_))=>{rejected=true;break;},Err(e)=>panic!("{e}")}}
        assert!(rejected);
        let nil=Matrix::from_rows(vec![vec![Expr::zero(),Expr::one()],vec![Expr::zero(),Expr::zero()]]).unwrap();
        assert!(Budget::new().power(&nil,BigUint::one()<<100).unwrap().is_zero_matrix());
        assert_eq!(Budget::new().power(&Matrix::eye(2),BigUint::one()<<100).unwrap(),Matrix::eye(2));
        // Raw constructors preserve the latent power: rejecting it must not
        // call pow(), cancel(), substitution or any numerical evaluator.
        let held=Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2),Expr::int(1_000_000_000)));
        assert!(Budget::new().expression(&held).is_err());
        let x=Expr::symbol("x");
        let mut mixed=Expr::from_unevaluated_kind(Kind::Mul(Number::int(2),vec![(x.clone(),Expr::one())]));
        for _ in 0..80{mixed=Expr::from_unevaluated_kind(Kind::Mul(Number::int(1),vec![(mixed.clone(),Expr::one()),(mixed,Expr::one())]));}
        assert!(Budget::new().expression(&mixed).is_err());
        let mixed=Expr::from_unevaluated_kind(Kind::Mul(Number::int(2),vec![(x.clone(),Expr::one())]));
        let mixed=Expr::from_unevaluated_kind(Kind::Pow(mixed,Expr::int(1_000_000_000)));
        assert!(Budget::new().expression(&mixed).is_err());
        let f=Expr::from_unevaluated_kind(Kind::Pow(x.clone(),Expr::int(1_000_000_000)));
        assert!(Budget::new().scalar_preflight(&f,Some((&x,&Expr::int(2)))).is_err());
        for n in [-1,0,1]{assert!(Budget::new().scalar_preflight(&f,Some((&x,&Expr::int(n)))).is_ok());}
        let y=Expr::symbol("y");
        let joint=Expr::from_unevaluated_kind(Kind::Pow(crate::ops::mul2(&x,&y),Expr::int(8)));
        let mut budget=Budget::new();budget.bits=64;
        assert!(budget.scalar_preflight(&joint,Some((&x,&Expr::int(16)))).is_ok());
        assert!(budget.scalar_preflight(&joint,Some((&y,&Expr::int(16)))).is_ok());
        let replacements=FxHashMap::from_iter([(x.clone(),Expr::int(16)),(y,Expr::int(16))]);
        assert!(budget.substitutions_preflight(&joint,&replacements).is_err());
        let shifted=crate::ops::add2(&x,&Expr::one());
        let big=Expr::from_unevaluated_kind(Kind::Pow(shifted,Expr::int(1_000_000_000)));
        assert!(Budget::new().cancel(&Matrix::new(1,1,vec![big]).unwrap()).is_err());
        let symbols:(Vec<_>,Vec<_>)=(0..32).map(|i|{let x=Expr::symbol(&format!("a{i}"));let y=Expr::symbol(&format!("b{i}"));(x,y)}).unzip();
        let a=Matrix::new(1,1,vec![crate::ops::add(&symbols.0)]).unwrap();
        let b=Matrix::new(1,1,vec![crate::ops::add(&symbols.1)]).unwrap();
        let mut budget=Budget::new();budget.work=500;
        assert!(budget.product_preflight(&a,&b).is_err());
        let width=crate::ops::add(&symbols.0);
        let wide=Matrix::new(3,3,vec![width;9]).unwrap();
        assert!(Budget::new().spectral_preflight(&wide).is_err());
        let mut budget=Budget::new();budget.storage=20_000;
        assert!(budget.cancel_preflight(&a).is_err());
        let mut budget=Budget::new();budget.work=10;assert!(budget.interpolation(3,1).is_err());
        let mut budget=Budget::new();budget.storage=10;assert!(budget.identity(2).is_err());
    }
}
