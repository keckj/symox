//! Native continuation of the hypergeometric, Mathieu and Bessel families.
use super::*;
use crate::functions::{cos, sin, exp, gamma};
use crate::ops::{add, sub, mul, mul2, div, neg, pow, sqrt};

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        MathieuC | MathieuS | MathieuCPrime | MathieuSPrime => {
            let (aa,q,z)=(&a[0],&a[1],&a[2]);
            if q.is_zero() {
                let root=sqrt(aa);let v=mul2(&root,z);
                return Some(match s {
                    MathieuC=>cos(&v),MathieuS=>sin(&v),
                    MathieuCPrime=>neg(&mul2(&root,&sin(&v))),
                    _=>mul2(&root,&cos(&v)),
                });
            }
            if could_extract_minus_sign(z) {
                let value=sf(s,vec![aa.clone(),q.clone(),neg(z)]);
                return Some(if matches!(s,MathieuS|MathieuCPrime) {neg(&value)} else {value});
            }
            None
        }
        MarcumQ => {
            let (m,x,y)=(&a[0],&a[1],&a[2]);
            if is_zero(x)==Some(true) {
                if is_zero(m)==Some(true) && is_zero(y)==Some(true) {return Some(Expr::zero());}
                return Some(div(&super::uppergamma(m,&div(&y.powi(2),&int(2))),&gamma(m)));
            }
            if is_zero(m)==Some(true) && is_zero(y)==Some(true) {
                return Some(sub(&Expr::one(),&exp(&neg(&div(&x.powi(2),&int(2))))));
            }
            if x==y && matches!(as_small(m),Some(1)|Some(2)) {
                let exponential=exp(&neg(&x.powi(2)));
                let mut value=add(&[rat(1,2),mul(&[rat(1,2),exponential.clone(),super::besseli(&Expr::zero(),&x.powi(2))])]);
                if as_small(m)==Some(2) {value=add(&[value,mul2(&exponential,&super::besseli(&Expr::one(),&x.powi(2)))]);}
                return Some(value);
            }
            None
        }
        AppellF1 => {
            if a[4].is_zero() && a[5].is_zero() {return Some(Expr::one());}
            None
        }
        _ => None,
    }
}

pub(super) fn partial(s:SpecialFn,a:&[Expr],i:usize)->Option<Expr> {
    use SpecialFn::*;
    let one=Expr::one();let two=int(2);
    let p1=|x:&Expr|add(&[x.clone(),one.clone()]);
    let m1=|x:&Expr|sub(x,&one);
    match (s,i) {
        (MathieuC,2)=>Some(sf(MathieuCPrime,a.to_vec())),
        (MathieuS,2)=>Some(sf(MathieuSPrime,a.to_vec())),
        (MathieuCPrime|MathieuSPrime,2)=>Some(mul2(&sub(&mul(&[two.clone(),a[1].clone(),cos(&mul2(&two,&a[2]))]),&a[0]),&sf(if s==MathieuCPrime{MathieuC}else{MathieuS},a.to_vec()))),
        (Hn1|Hn2,1)=>Some(sub(&sf(s,vec![m1(&a[0]),a[1].clone()]),&mul2(&div(&p1(&a[0]),&a[1]),&sf(s,a.to_vec())))),
        (MarcumQ,1)=>Some(mul2(&a[1],&sub(&sf(s,vec![p1(&a[0]),a[1].clone(),a[2].clone()]),&sf(s,a.to_vec())))),
        (MarcumQ,2)=>Some(neg(&mul(&[pow(&a[1],&sub(&one,&a[0])),pow(&a[2],&a[0]),exp(&neg(&div(&add(&[a[1].powi(2),a[2].powi(2)]),&two))),super::besseli(&m1(&a[0]),&mul2(&a[1],&a[2]))]))),
        (AppellF1,4|5)=>{
            let b=if i==4{1}else{2};let mut shifted=a.to_vec();
            for j in [0,b,3] {shifted[j]=p1(&a[j]);}
            Some(mul2(&div(&mul2(&a[0],&a[b]),&a[3]),&sf(s,shifted)))
        }
        (MeijerG,2)=>{
            let aa=tuple_items(&a[0])?;let bb=tuple_items(&a[1])?;
            let an=tuple_items(&aa[0])?;let bm=tuple_items(&bb[0])?;
            let (factor,shifted)=if !an.is_empty(){
                let mut v=an.to_vec();v[0]=m1(&v[0]);
                (mul2(&m1(&an[0]),&sf(s,a.to_vec())),sf(s,vec![tuple(vec![tuple(v),aa[1].clone()]),a[1].clone(),a[2].clone()]))
            }else if !bm.is_empty(){
                let mut v=bm.to_vec();v[0]=p1(&v[0]);
                (mul2(&bm[0],&sf(s,a.to_vec())),neg(&sf(s,vec![a[0].clone(),tuple(vec![tuple(v),bb[1].clone()]),a[2].clone()])))
            }else{return Some(Expr::zero())};
            Some(div(&add(&[factor,shifted]),&a[2]))
        }
        _=>None,
    }
}

/// A conservative Slater expansion. Equal or integer-separated pole families
/// need logarithmic residues and deliberately remain held.
fn slater(a:&[Expr])->Option<Expr>{
    let aa=tuple_items(&a[0])?;let bb=tuple_items(&a[1])?;
    let an=tuple_items(&aa[0])?;let ao=tuple_items(&aa[1])?;
    let bm=tuple_items(&bb[0])?;let bo=tuple_items(&bb[1])?;
    let ap:Vec<Expr>=an.iter().chain(ao).cloned().collect();
    let bq:Vec<Expr>=bm.iter().chain(bo).cloned().collect();
    if bq.len()>16 || ap.len()>bq.len() {return None;}
    if bm.is_empty(){return Some(Expr::zero());}
    if bm.iter().enumerate().any(|(i,b)|bm.iter().enumerate().any(|(j,c)|i!=j && is_integer(&sub(b,c))!=Some(false))){return None;}
    let mut terms=Vec::new();
    let sign=if (ap.len() as isize-bm.len() as isize-an.len() as isize)%2==0 {Expr::one()}else{Expr::minus_one()};
    for (k,b) in bm.iter().enumerate(){
        let mut numerator=Vec::new();let mut denominator=Vec::new();
        for (j,c) in bm.iter().enumerate(){if j!=k{numerator.push(gamma(&sub(c,b)));}}
        for c in an{numerator.push(gamma(&add(&[Expr::one(),sub(b,c)])));}
        for c in ao{denominator.push(gamma(&sub(c,b)));}
        for c in bo{denominator.push(gamma(&add(&[Expr::one(),sub(b,c)])));}
        let params_a:Vec<Expr>=ap.iter().map(|c|add(&[Expr::one(),sub(b,c)])).collect();
        let params_b:Vec<Expr>=bq.iter().enumerate().filter(|(j,_)|*j!=k).map(|(_,c)|add(&[Expr::one(),sub(b,c)])).collect();
        terms.push(mul(&[div(&mul(&numerator),&mul(&denominator)),pow(&a[2],b),super::hyper(&params_a,&params_b,&mul2(&sign,&a[2]))]));
    }
    Some(super::expand::hyperexpand(&add(&terms)))
}

pub(crate) fn meijer_expand(a:&[Expr])->Option<Expr>{
    let aa=tuple_items(&a[0])?;let bb=tuple_items(&a[1])?;
    let an=tuple_items(&aa[0])?;let ao=tuple_items(&aa[1])?;
    let bm=tuple_items(&bb[0])?;let bo=tuple_items(&bb[1])?;
    let p=an.len()+ao.len();let q=bm.len()+bo.len();
    if p<q{return slater(a);}
    let complement=|values:&[Expr]|tuple(values.iter().map(|v|sub(&Expr::one(),v)).collect());
    let dual=vec![tuple(vec![complement(bm),complement(bo)]),tuple(vec![complement(an),complement(ao)]),a[2].powi(-1)];
    if p>q{return slater(&dual);}
    if p==0{return None;}
    let small=slater(a)?;let big=slater(&dual)?;
    Some(crate::relational::piecewise(vec![
        (small,crate::relational::lt(&crate::functions::abs(&a[2]),&Expr::one())),
        (big,crate::relational::lt(&crate::functions::abs(&dual[2]),&Expr::one())),
        (raw(SpecialFn::MeijerG,a.to_vec()),Expr::boolean(true)),
    ]))
}

pub(crate) fn rewrite(s:SpecialFn,a:&[Expr],rule:&str)->Option<Expr>{
    use SpecialFn::*;
    if matches!(s,Hn1|Hn2){
        let n=&a[0];let z=&a[1];let sign=int(if s==Hn1{1}else{-1});
        let prefactor=sqrt(&div(&Expr::pi(),&mul2(&int(2),z)));
        let order=add(&[n.clone(),rat(1,2)]);
        return match rule{
            "jn"|"yn"=>Some(add(&[super::jn(n,z),mul(&[sign,Expr::i(),super::yn(n,z)])])),
            "hankel1" if s==Hn1=>Some(mul2(&prefactor,&super::hankel1(&order,z))),
            "hankel2" if s==Hn2=>Some(mul2(&prefactor,&super::hankel2(&order,z))),
            "besselj" if is_integer(n)==Some(true)=>Some(mul2(&prefactor,&add(&[super::besselj(&order,z),mul(&[sign,Expr::i(),pow(&Expr::minus_one(),&add(&[n.clone(),Expr::one()])),super::besselj(&neg(&order),z)])]))),
            "bessely" if is_integer(n)==Some(true)=>Some(mul2(&prefactor,&add(&[mul2(&pow(&Expr::minus_one(),n),&super::bessely(&neg(&order),z)),mul(&[sign,Expr::i(),super::bessely(&order,z)])]))),
            _=>None,
        };
    }
    None
}
