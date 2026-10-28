//! Sparse execution of the recursive DMP product schedule. Coefficient
//! policies retain their own arithmetic and rounding; no flattened product.
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use crate::poly::convert::PResult;

pub(crate) trait Arithmetic {
    type C: Clone + PartialEq;
    fn zero(&self) -> Self::C;
    fn one(&self) -> Self::C;
    fn two(&self) -> Self::C;
    fn is_zero(&self, c: &Self::C) -> bool;
    fn implicit_zero(&self, c: &Self::C) -> bool;
    fn add(&self, a: &Self::C, b: &Self::C) -> PResult<Self::C>;
    fn sub(&self, a: &Self::C, b: &Self::C) -> PResult<Self::C>;
    fn mul(&self, a: &Self::C, b: &Self::C) -> PResult<Self::C>;
    fn square(&self, a: &Self::C) -> PResult<Self::C>;
    fn karatsuba(&self) -> bool { false }
    /// Only policies with pure coefficient arithmetic may share a result of
    /// the exact same ordered operands within a Karatsuba call.
    fn reuse_identical_products(&self) -> bool { false }
    fn coefficient_reusable(&self, _c:&Self::C)->bool {false}
}

#[derive(Clone, Debug, PartialEq)]
enum Tree<C> {
    Leaf(BTreeMap<u64, C>),
    Node(BTreeMap<u64, Tree<C>>),
}
impl<C> Tree<C> {
    fn empty(&self) -> bool {match self {Self::Leaf(v)=>v.is_empty(),Self::Node(v)=>v.is_empty()}}
}
fn normalized<K:Arithmetic>(k:&K, mut tree:Tree<K::C>)->Tree<K::C> {
    match &mut tree {
        Tree::Leaf(v)=>{
            while v.last_key_value().is_some_and(|(_,c)|k.is_zero(c)){v.pop_last();}
            v.retain(|_,c|!k.implicit_zero(c));
        }
        Tree::Node(v)=>v.retain(|_,c|!c.empty()),
    }
    tree
}
fn make<K:Arithmetic>(k:&K, terms:Vec<(Vec<u32>,K::C)>, level:usize, width:usize)->Tree<K::C> {
    if level+1==width {
        return normalized(k,Tree::Leaf(terms.into_iter().map(|(m,c)|(u64::from(m[level]),c)).collect()));
    }
    let mut groups:BTreeMap<u64,Vec<_>>=BTreeMap::new();
    for (m,c) in terms {groups.entry(u64::from(m[level])).or_default().push((m,c));}
    normalized(k,Tree::Node(groups.into_iter().map(|(e,t)|(e,make(k,t,level+1,width))).collect()))
}
fn one<K:Arithmetic>(k:&K,depth:usize)->Tree<K::C> {
    if depth==0 {Tree::Leaf(BTreeMap::from([(0,k.one())]))}
    else {Tree::Node(BTreeMap::from([(0,one(k,depth-1))]))}
}
fn add<K:Arithmetic>(k:&K,a:&Tree<K::C>,b:&Tree<K::C>,subtract:bool)->PResult<Tree<K::C>> {
    if b.empty(){return Ok(a.clone());}
    if a.empty()&&!subtract{return Ok(b.clone());}
    let out=match (a,b) {
        (Tree::Leaf(a),Tree::Leaf(b))=>{
            let keys:BTreeSet<_>=a.keys().chain(b.keys()).copied().collect();
            let mut out=BTreeMap::new();
            // Dense arrays visit the highest degree first.
            for e in keys.into_iter().rev(){
                let c=match(a.get(&e),b.get(&e)) {
                    (Some(x),Some(y))=>if subtract{k.sub(x,y)?}else{k.add(x,y)?},
                    (Some(x),None)=>x.clone(),
                    (None,Some(y))=>if subtract{k.sub(&k.zero(),y)?}else{y.clone()},
                    _=>unreachable!(),
                };
                out.insert(e,c);
            }
            Tree::Leaf(out)
        }
        (Tree::Node(a),Tree::Node(b))=>{
            let keys:BTreeSet<_>=a.keys().chain(b.keys()).copied().collect();
            let mut out=BTreeMap::new();
            for e in keys.into_iter().rev(){
                let c=match(a.get(&e),b.get(&e)) {
                    (Some(x),Some(y))=>add(k,x,y,subtract)?,
                    (Some(x),None)=>x.clone(),
                    (None,Some(y))=>if subtract{add(k,&zero_like(y),y,true)?}else{y.clone()},
                    _=>unreachable!(),
                };
                out.insert(e,c);
            }
            Tree::Node(out)
        }
        _=>unreachable!("sparse product depth mismatch"),
    };
    Ok(normalized(k,out))
}
fn zero_like<C>(tree:&Tree<C>)->Tree<C> {match tree{Tree::Leaf(_)=>Tree::Leaf(BTreeMap::new()),Tree::Node(_)=>Tree::Node(BTreeMap::new())}}
fn scale<K:Arithmetic>(k:&K,a:&Tree<K::C>,c:&K::C)->PResult<Tree<K::C>> {
    let out=match a {
        Tree::Leaf(v)=>Tree::Leaf(v.iter().rev().map(|(e,x)|Ok((*e,k.mul(x,c)?))).collect::<PResult<_>>()?),
        Tree::Node(v)=>Tree::Node(v.iter().rev().map(|(e,x)|Ok((*e,scale(k,x,c)?))).collect::<PResult<_>>()?),
    };
    Ok(normalized(k,out))
}
fn exponent_add(a:u64,b:u64)->PResult<u64>{a.checked_add(b).ok_or_else(crate::poly::mpoly::exponent_limit_error)}
// Merge product rows in the dense schedule without retaining a Cartesian
// table of pairs. At most one pending pair per left monomial is stored.
struct PairSums {
    left:Vec<u64>,right:Vec<u64>,heap:BinaryHeap<(u64,u64,usize,usize)>,
}
impl PairSums {
    fn new(left:Vec<u64>,right:Vec<u64>,square:bool)->PResult<Self> {
        if let(Some(a),Some(b))=(left.first(),right.first()){exponent_add(*a,*b)?;}
        let mut heap=BinaryHeap::new();
        for (ia,ea) in left.iter().enumerate(){
            let ib=if square{ia}else{0};
            if let Some(eb)=right.get(ib){heap.push((ea+eb,*ea,ia,ib));}
        }
        Ok(Self{left,right,heap})
    }
}
impl Iterator for PairSums {
    type Item=(u64,u64,u64);
    fn next(&mut self)->Option<Self::Item>{
        let (sum,ea,ia,ib)=self.heap.pop()?;
        let eb=self.right[ib];
        if let Some(next)=self.right.get(ib+1){
            // All sums were bounded by the two maximal exponents in new().
            self.heap.push((self.left[ia]+next,ea,ia,ib+1));
        }
        Some((sum,ea,eb))
    }
}
struct PairGroups(std::iter::Peekable<PairSums>);
impl PairGroups {
    fn new(left:Vec<u64>,right:Vec<u64>,square:bool)->PResult<Self>{Ok(Self(PairSums::new(left,right,square)?.peekable()))}
}
impl Iterator for PairGroups {
    type Item=(u64,Vec<(u64,u64)>);
    fn next(&mut self)->Option<Self::Item>{
        let (sum,a,b)=self.0.next()?;let mut pairs=vec![(a,b)];
        while self.0.peek().is_some_and(|(e,_,_)|*e==sum){
            let (_,a,b)=self.0.next().expect("peeked pair");pairs.push((a,b));
        }
        Some((sum,pairs))
    }
}
fn leaf_slice<C:Clone>(f:&BTreeMap<u64,C>,lo:u64,hi:u64,shift:u64)->Tree<C> {
    Tree::Leaf(f.range(lo..hi).map(|(e,c)|(e-shift,c.clone())).collect())
}
fn leaf_shift<C:Clone>(f:&Tree<C>,amount:u64)->PResult<Tree<C>> {
    let Tree::Leaf(f)=f else{unreachable!()};
    Ok(Tree::Leaf(f.iter().map(|(e,c)|Ok((exponent_add(*e,amount)?,c.clone()))).collect::<PResult<_>>()?))
}
fn reusable<K:Arithmetic>(k:&K,f:&Tree<K::C>)->bool {
    match f {Tree::Leaf(v)=>v.values().all(|c|k.coefficient_reusable(c)),Tree::Node(v)=>v.values().all(|c|reusable(k,c))}
}
fn mul<K:Arithmetic>(k:&K,f:&Tree<K::C>,g:&Tree<K::C>)->PResult<Tree<K::C>> {
    if f==g{return square(k,f);}
    if f.empty()||g.empty(){return Ok(zero_like(f));}
    let result=match(f,g){
        (Tree::Leaf(a),Tree::Leaf(b))=>{
            let n=exponent_add(*a.last_key_value().unwrap().0.max(b.last_key_value().unwrap().0),1)?;
            if k.karatsuba()&&n>=100 {
                let half=n/2;
                // The exact dup_mul split/recombine order, not a new strategy.
                let (lo_a,lo_b)=(leaf_slice(a,0,half,0),leaf_slice(b,0,half,0));
                let (hi_a,hi_b)=(leaf_slice(a,half,n,half),leaf_slice(b,half,n,half));
                let lo=mul(k,&lo_a,&lo_b)?;let hi=mul(k,&hi_a,&hi_b)?;
                let sum_a=add(k,&lo_a,&hi_a,false)?;let sum_b=add(k,&lo_b,&hi_b,false)?;
                let share=k.reuse_identical_products() && reusable(k,&sum_a) && reusable(k,&sum_b);
                let mid=if share && sum_a==lo_a && sum_b==lo_b {lo.clone()}
                    else if share && sum_a==hi_a && sum_b==hi_b {hi.clone()}
                    else {mul(k,&sum_a,&sum_b)?};
                let mid=add(k,&mid,&add(k,&lo,&hi,false)?,true)?;
                return add(k,&add(k,&lo,&leaf_shift(&mid,half)?,false)?,&leaf_shift(&hi,2*half)?,false);
            }
            let groups=PairGroups::new(a.keys().rev().copied().collect(),b.keys().rev().copied().collect(),false)?;
            let mut out=BTreeMap::new();
            for (e,pairs) in groups{
                let mut c=k.zero();
                for (ea,eb) in pairs{c=k.add(&c,&k.mul(&a[&ea],&b[&eb])?)?;}
                out.insert(e,c);
            }
            Tree::Leaf(out)
        }
        (Tree::Node(a),Tree::Node(b))=>{
            let groups=PairGroups::new(a.keys().rev().copied().collect(),b.keys().rev().copied().collect(),false)?;
            let mut out=BTreeMap::new();
            for (e,pairs) in groups{
                let mut c=zero_like(a.first_key_value().unwrap().1);
                for (ea,eb) in pairs{c=add(k,&c,&mul(k,&a[&ea],&b[&eb])?,false)?;}
                out.insert(e,c);
            }
            Tree::Node(out)
        }
        _=>unreachable!("sparse product depth mismatch"),
    };
    Ok(normalized(k,result))
}
fn square<K:Arithmetic>(k:&K,f:&Tree<K::C>)->PResult<Tree<K::C>> {
    if f.empty(){return Ok(f.clone());}
    let out=match f {
        Tree::Leaf(a)=>{
            let keys:Vec<_>=a.keys().rev().copied().collect();
            let groups=PairGroups::new(keys.clone(),keys,true)?;
            let mut out=BTreeMap::new();
            for (e,pairs) in groups{
                let mut c=k.zero();let mut diagonal=None;
                for (ea,eb) in pairs{if ea==eb{diagonal=Some(ea);}else{c=k.add(&c,&k.mul(&a[&ea],&a[&eb])?)?;}}
                c=k.add(&c,&c)?;
                if let Some(ea)=diagonal{c=k.add(&c,&k.square(&a[&ea])?)?;}
                out.insert(e,c);
            }
            Tree::Leaf(out)
        }
        Tree::Node(a)=>{
            let keys:Vec<_>=a.keys().rev().copied().collect();
            let groups=PairGroups::new(keys.clone(),keys,true)?;
            let mut out=BTreeMap::new();
            for (e,pairs) in groups{
                let mut c=zero_like(a.first_key_value().unwrap().1);let mut diagonal=None;
                for (ea,eb) in pairs{if ea==eb{diagonal=Some(ea);}else{c=add(k,&c,&mul(k,&a[&ea],&a[&eb])?,false)?;}}
                c=scale(k,&c,&k.two())?;
                if let Some(ea)=diagonal{c=add(k,&c,&square(k,&a[&ea])?,false)?;}
                out.insert(e,c);
            }
            Tree::Node(out)
        }
    };
    Ok(normalized(k,out))
}
fn flatten<K:Arithmetic>(k:&K,f:&Tree<K::C>,prefix:&mut Vec<u32>,out:&mut Vec<(Vec<u32>,K::C)>)->PResult<()> {
    match f {
        Tree::Leaf(v)=>for(e,c)in v.iter().rev(){
            let e=u32::try_from(*e).map_err(|_|crate::poly::mpoly::exponent_limit_error())?;
            let mut m=prefix.clone();m.push(e);out.push((m,c.clone()));
        },
        Tree::Node(v)=>for(e,c)in v.iter().rev(){
            prefix.push(u32::try_from(*e).map_err(|_|crate::poly::mpoly::exponent_limit_error())?);
            flatten(k,c,prefix,out)?;prefix.pop();
        },
    }
    Ok(())
}
/// Execute a DMP multiplication or binary power using sparse recursive nodes.
pub(crate) fn product<K:Arithmetic>(k:&K,width:usize,a:Vec<(Vec<u32>,K::C)>,b:Option<Vec<(Vec<u32>,K::C)>>,power:Option<u32>)->PResult<Vec<(Vec<u32>,K::C)>> {
    let mut f=make(k,a,0,width);
    let result=if let Some(mut n)=power {
        let mut g=one(k,width-1);
        if n==0{g}else if n==1||f.empty()||f==g{f}else{
            loop {
                let m=n;n/=2;
                if m%2==1{g=mul(k,&g,&f)?;if n==0{break;}}
                f=square(k,&f)?;
            }
            g
        }
    }else{let g=make(k,b.expect("product rhs"),0,width);mul(k,&f,&g)?};
    let mut out=Vec::new();flatten(k,&result,&mut Vec::new(),&mut out)?;Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Mod4;
    impl Arithmetic for Mod4 {
        type C=u8;
        fn zero(&self)->u8{0} fn one(&self)->u8{1} fn two(&self)->u8{2}
        fn is_zero(&self,c:&u8)->bool{*c==0} fn implicit_zero(&self,c:&u8)->bool{*c==0}
        fn add(&self,a:&u8,b:&u8)->PResult<u8>{Ok((a+b)%4)}
        fn sub(&self,a:&u8,b:&u8)->PResult<u8>{Ok((a+4-b)%4)}
        fn mul(&self,a:&u8,b:&u8)->PResult<u8>{Ok((a*b)%4)}
        fn square(&self,a:&u8)->PResult<u8>{Ok((a*a)%4)}
    }
    #[test]
    fn sparse_recursive_product_wide_cancellation_precedes_coordinate_bound() {
        let n=1u32<<31;
        let terms=vec![(vec![n],2),(vec![0],1)];
        assert_eq!(product(&Mod4,1,terms.clone(),Some(terms.clone()),None).unwrap(),vec![(vec![0],1)]);
        assert_eq!(product(&Mod4,1,terms,None,Some(2)).unwrap(),vec![(vec![0],1)]);
        assert!(product(&Mod4,1,vec![(vec![n],1)],None,Some(2)).is_err());
    }
}

#[cfg(test)]
mod callback_tests {
    use super::*;
    use std::cell::RefCell;
    #[derive(Clone,PartialEq,Debug)]
    struct Probe(i64,String);
    struct Callbacks(RefCell<Vec<String>>);
    impl Arithmetic for Callbacks {
        type C=Probe;
        fn zero(&self)->Probe{Probe(0,"0".into())} fn one(&self)->Probe{Probe(1,"1".into())} fn two(&self)->Probe{Probe(2,"2".into())}
        fn is_zero(&self,c:&Probe)->bool{c.0==0} fn implicit_zero(&self,c:&Probe)->bool{c.0==0}
        fn add(&self,a:&Probe,b:&Probe)->PResult<Probe>{
            if a.0==0{return Ok(b.clone());}if b.0==0{return Ok(a.clone());}
            Ok(Probe(a.0+b.0,format!("({}+{})",a.1,b.1)))
        }
        fn sub(&self,a:&Probe,b:&Probe)->PResult<Probe>{Ok(Probe(a.0-b.0,format!("({}-{})",a.1,b.1)))}
        fn mul(&self,a:&Probe,b:&Probe)->PResult<Probe>{
            if a.0==0||b.0==0{return Ok(self.zero());}
            self.0.borrow_mut().push(format!("{}*{}",a.1,b.1));
            Ok(Probe(a.0*b.0,format!("{}*{}",a.1,b.1)))
        }
        fn square(&self,a:&Probe)->PResult<Probe>{self.0.borrow_mut().push(format!("sq({})",a.1));Ok(Probe(a.0*a.0,format!("sq({})",a.1)))}
    }
    #[test]
    fn sparse_recursive_product_callbacks_follow_nested_dense_schedule() {
        let k=Callbacks(RefCell::new(Vec::new()));
        let p=|v,s:&str|Probe(v,s.into());
        let f=vec![(vec![1,1],p(1,"a")),(vec![1,0],p(2,"b")),(vec![0,0],p(3,"c"))];
        let g=vec![(vec![1,1],p(4,"d")),(vec![1,0],p(5,"e")),(vec![0,0],p(6,"f"))];
        let out=product(&k,2,f.clone(),Some(g),None).unwrap();
        assert_eq!(&*k.0.borrow(),&["a*d","a*e","b*d","b*e","a*f","b*f","c*d","c*e","c*f"]);
        assert_eq!(out.iter().find(|(m,_)|m==&[1,1]).unwrap().1.1,"(a*f+c*d)");
        k.0.borrow_mut().clear();
        let out=product(&k,2,f.clone(),Some(f),None).unwrap();
        assert_eq!(&*k.0.borrow(),&["sq(a)","a*b","sq(b)","a*c","b*c","a*c*2","b*c*2","sq(c)"]);
        assert_eq!(out.iter().find(|(m,_)|m==&[2,1]).unwrap().1.1,"(a*b+a*b)");
        assert_eq!(out.iter().find(|(m,_)|m==&[1,1]).unwrap().1.1,"a*c*2");
    }
}

#[cfg(test)]
mod interruption_tests {
    use super::*;
    struct Interrupted;
    struct Host(bool);
    impl Arithmetic for Host {
        type C=u8;
        fn zero(&self)->u8{0} fn one(&self)->u8{1} fn two(&self)->u8{2}
        fn is_zero(&self,c:&u8)->bool{*c==0} fn implicit_zero(&self,c:&u8)->bool{*c==0}
        fn add(&self,a:&u8,b:&u8)->PResult<u8>{Ok(a+b)}
        fn sub(&self,a:&u8,b:&u8)->PResult<u8>{Ok(a.saturating_sub(*b))}
        fn mul(&self,_:&u8,_:&u8)->PResult<u8>{
            if self.0 {std::panic::resume_unwind(Box::new(Interrupted));}
            Err(crate::poly::convert::PolyError::Reference{kind:"RuntimeError",message:"host callback".into()})
        }
        fn square(&self,a:&u8)->PResult<u8>{self.mul(a,a)}
    }
    #[test]
    fn sparse_recursive_product_never_turns_host_errors_or_interrupts_into_zero() {
        let run=|k:&Host|product(k,1,vec![(vec![1],2)],Some(vec![(vec![1],3)]),None);
        assert!(matches!(run(&Host(false)),Err(crate::poly::convert::PolyError::Reference{kind:"RuntimeError",message}) if message=="host callback"));
        let payload=std::panic::catch_unwind(||run(&Host(true))).unwrap_err();
        assert!(payload.is::<Interrupted>());
        assert!(!Host(false).reuse_identical_products());
    }
}

#[cfg(test)]
mod row_storage_tests {
    use super::*;
    #[test]
    fn sparse_recursive_product_pair_storage_is_linear_and_ordered() {
        let keys:Vec<u64>=(0..4096).rev().collect();
        for square in [false,true]{
            let mut rows=PairSums::new(keys.clone(),keys.clone(),square).unwrap();
            assert_eq!(rows.heap.len(),4096);
            let mut previous=(u64::MAX,u64::MAX);
            for _ in 0..10000 {
                let (sum,left,_)=rows.next().unwrap();
                assert!((sum,left)<=previous);previous=(sum,left);
                assert!(rows.heap.len()<=4096);
            }
        }
    }
}
