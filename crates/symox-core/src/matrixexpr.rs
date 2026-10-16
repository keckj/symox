//! Typed, noncommutative symbolic matrix expressions over scalar Expr dimensions.
use std::cmp::Ordering;
use std::sync::Arc;
use rustc_hash::{FxHashMap, FxHashSet};
use crate::expr::{Expr, Func};
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
    Symbol(String), Explicit(Matrix), Identity, Zero,
    Add(Vec<MatrixExpr>), Mul(Expr, Vec<MatrixExpr>), Pow(MatrixExpr, Expr),
    Transpose(MatrixExpr), Adjoint(MatrixExpr), Inverse(MatrixExpr), Conjugate(MatrixExpr),
}
impl Ord for MatrixExpr {
    fn cmp(&self, other: &Self) -> Ordering {
        if self == other { Ordering::Equal } else { self.to_repr().cmp(&other.to_repr()) }
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
    pub fn explicit(m:Matrix)->Self {Self::new(MatrixKind::Explicit(m.clone()),Expr::int(m.rows() as i64),Expr::int(m.cols() as i64))}
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
            MatrixKind::Mul(c,fs) if fs.iter().all(|f|f.rows()==f.cols())=>return Self::mul(pow(c,&Expr::minus_one()),&fs.iter().rev().map(|f|f.inverse(true)).collect::<MResult<Vec<_>>>()?,true),
            _=>{},
        }}
        Ok(Self::new(MatrixKind::Inverse(self.clone()),self.rows().clone(),self.cols().clone()))
    }
    pub fn transpose(&self,evaluate:bool)->MResult<Self>{self.unary(0,evaluate)}
    pub fn adjoint(&self,evaluate:bool)->MResult<Self>{self.unary(1,evaluate)}
    pub fn conjugate(&self,evaluate:bool)->MResult<Self>{self.unary(2,evaluate)}
    fn unary(&self,op:u8,evaluate:bool)->MResult<Self>{
        let apply=|m:&Self|m.unary(op,true);
        if evaluate {match self.kind(){
            MatrixKind::Identity=>return Ok(self.clone()),MatrixKind::Zero=>return if op==2{Ok(self.clone())}else{Self::zero(self.cols().clone(),self.rows().clone())},
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
    pub fn children(&self)->Vec<Self>{match self.kind(){MatrixKind::Add(a)|MatrixKind::Mul(_,a)=>a.clone(),MatrixKind::Pow(b,_)|MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Inverse(b)|MatrixKind::Conjugate(b)=>vec![b.clone()],_=>Vec::new()}}
    pub fn map_children(&self, matrices:&mut dyn FnMut(&Self)->Self, scalars:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        match self.kind(){
            MatrixKind::Symbol(name)=>{let r=scalars(self.rows());let c=scalars(self.cols());if r==*self.rows()&&c==*self.cols(){Ok(self.clone())}else{Self::symbol(name,r,c)}},
            MatrixKind::Identity=>{let n=scalars(self.rows());if n==*self.rows(){Ok(self.clone())}else{Self::identity(n)}},
            MatrixKind::Zero=>{let r=scalars(self.rows());let c=scalars(self.cols());if r==*self.rows()&&c==*self.cols(){Ok(self.clone())}else{Self::zero(r,c)}},
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
            MatrixKind::Symbol(_)=>{}, MatrixKind::Identity|MatrixKind::Zero=>{push(self.rows());push(self.cols());},
            MatrixKind::Explicit(m)=>{for e in m.data(){push(e);}},MatrixKind::Mul(c,_)=>push(c),MatrixKind::Pow(_,e)=>push(e),_=>{},
        }
        for c in self.children(){c.collect_scalars(out);}
    }
    pub fn matrix_map(&self,f:&mut dyn FnMut(&Self)->Self)->MResult<Self>{
        let replacement=f(self);if replacement!=*self{return Ok(replacement);}
        let old=self.children();let children=old.iter().map(|c|c.matrix_map(f)).collect::<MResult<Vec<_>>>()?;
        if children==old{return Ok(self.clone());}
        match self.kind(){MatrixKind::Add(_)=>Self::add(&children,true),MatrixKind::Mul(c,_)=>Self::mul(c.clone(),&children,true),MatrixKind::Pow(_,e)=>children[0].pow(e.clone(),true),MatrixKind::Transpose(_)=>children[0].transpose(true),MatrixKind::Adjoint(_)=>children[0].adjoint(true),MatrixKind::Inverse(_)=>children[0].inverse(true),MatrixKind::Conjugate(_)=>children[0].conjugate(true),_=>Ok(self.clone())}
    }
    fn scalar_fields(&self,out:&mut Vec<Expr>){
        match self.kind(){MatrixKind::Symbol(_)|MatrixKind::Identity|MatrixKind::Zero=>{out.push(self.rows().clone());out.push(self.cols().clone());},MatrixKind::Explicit(m)=>out.extend(m.data().iter().cloned()),MatrixKind::Mul(c,_)=>out.push(c.clone()),MatrixKind::Pow(_,e)=>out.push(e.clone()),_=>{}}
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

    pub fn xreplace(&self,matrices:&[(Self,Self)],scalars:&FxHashMap<Expr,Expr>)->MResult<Self>{self.map_all(matrices,&mut |e|crate::subs::xreplace(e,scalars))}
    pub fn doit(&self,deep:bool)->MResult<Self>{
        let children=if deep{self.children().iter().map(|c|c.doit(true)).collect::<MResult<Vec<_>>>()?}else{self.children()};
        match self.kind(){MatrixKind::Add(_)=>Self::add(&children,true),MatrixKind::Mul(c,_)=>Self::mul(c.clone(),&children,true),MatrixKind::Pow(_,e)=>children[0].pow(e.clone(),true),MatrixKind::Inverse(_)=>children[0].inverse(true),MatrixKind::Transpose(_)=>children[0].transpose(true),MatrixKind::Adjoint(_)=>children[0].adjoint(true),MatrixKind::Conjugate(_)=>children[0].conjugate(true),_=>Ok(self.clone())}
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
    pub fn element(&self,i:&Expr,j:&Expr)->MResult<Expr>{
        Self::validate_index(i,self.rows())?;Self::validate_index(j,self.cols())?;
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),
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
            MatrixKind::Add(a)=>{let mut m=Matrix::zeros(rows,cols);for e in a{m=m.add(&e.as_explicit()?)?;}Ok(m)},
            MatrixKind::Mul(c,fs)=>{let mut m=fs[0].as_explicit()?;for e in &fs[1..]{m=m.matmul(&e.as_explicit()?)?;}Ok(m.scale(c))},
            MatrixKind::Pow(b,e)=>{let n=e.as_i64().ok_or_else(||MatrixError::NotImplemented("Explicit matrix powers require an integer exponent".into()))?;b.as_explicit()?.powi(n)},
            MatrixKind::Inverse(b)=>b.as_explicit()?.inv(),MatrixKind::Transpose(b)=>Ok(b.as_explicit()?.transpose()),MatrixKind::Adjoint(b)=>Ok(b.as_explicit()?.adjoint()),MatrixKind::Conjugate(b)=>Ok(b.as_explicit()?.conjugate()),
            MatrixKind::Symbol(_)=>{let mut entries=Vec::with_capacity(rows*cols);for i in 0..rows{for j in 0..cols{entries.push(self.opaque_element(&Expr::int(i as i64),&Expr::int(j as i64)));}}Matrix::new(rows,cols,entries)},
        }
    }
    pub fn kind_name(&self)->&'static str{match self.kind(){MatrixKind::Symbol(_)=>"MatrixSymbol",MatrixKind::Explicit(_)=>"ImmutableMatrix",MatrixKind::Identity=>"Identity",MatrixKind::Zero=>"ZeroMatrix",MatrixKind::Add(_)=>"MatAdd",MatrixKind::Mul(..)=>"MatMul",MatrixKind::Pow(..)=>"MatPow",MatrixKind::Transpose(_)=>"Transpose",MatrixKind::Adjoint(_)=>"Adjoint",MatrixKind::Inverse(_)=>"Inverse",MatrixKind::Conjugate(_)=>"MatrixConjugate"}}
    pub fn to_repr(&self)->String{
        let s=crate::printer::srepr;let me=|e:&Self|e.to_repr();
        let args=match self.kind(){
            MatrixKind::Symbol(name)=>format!("{name:?}, {}, {}",s(self.rows()),s(self.cols())),MatrixKind::Identity=>s(self.rows()),MatrixKind::Zero=>format!("{}, {}",s(self.rows()),s(self.cols())),
            MatrixKind::Explicit(m)=>format!("{}, {}, [{}]",m.rows(),m.cols(),m.data().iter().map(s).collect::<Vec<_>>().join(", ")),
            MatrixKind::Add(a)=>a.iter().map(me).collect::<Vec<_>>().join(", "),MatrixKind::Mul(c,a)=>{let mut a=a.iter().map(me).collect::<Vec<_>>();if !c.is_one(){a.insert(0,s(c));}a.join(", ")},
            MatrixKind::Pow(b,e)=>format!("{}, {}",b.to_repr(),s(e)),_=>self.children()[0].to_repr(),
        };format!("{}({args})",self.kind_name())
    }
    pub fn to_str(&self)->String{
        let child=|e:&Self|if matches!(e.kind(),MatrixKind::Add(_)){format!("({})",e.to_str())}else{e.to_str()};
        match self.kind(){MatrixKind::Symbol(n)=>n.clone(),MatrixKind::Explicit(m)=>m.to_str(),MatrixKind::Identity=>"I".into(),MatrixKind::Zero=>"0".into(),
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
}
