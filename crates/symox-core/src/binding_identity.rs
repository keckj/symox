//! Conservative alpha identity for held substitutions. Public trees stay intact.
use std::sync::Arc;
use std::hash::{Hash,Hasher};
use std::cmp::Ordering;
use rustc_hash::{FxHasher,FxHashSet};
use rustc_hash::FxHashMap;
use crate::expr::{Expr,Func,Kind};

#[derive(Debug)]
pub(crate) struct Key {hash:u64,data:Data}
#[derive(Debug,PartialEq,Eq,Hash,PartialOrd,Ord)]
enum Data {Free(Expr),Bound(u32),Node(u8,Vec<Arc<Key>>),Function(Func,Vec<Arc<Key>>)}
impl Key {
    fn make(data:Data)->Arc<Self>{let mut h=FxHasher::default();data.hash(&mut h);Arc::new(Self{hash:h.finish(),data})}
    fn free(e:Expr)->Arc<Self>{Self::make(Data::Free(e))}
    fn bound(i:u32)->Arc<Self>{Self::make(Data::Bound(i))}
    fn node(tag:u8,args:Vec<Arc<Key>>)->Arc<Self>{Self::make(Data::Node(tag,args))}
    fn function(f:Func,args:Vec<Arc<Key>>)->Arc<Self>{Self::make(Data::Function(f,args))}
}
impl Hash for Key {fn hash<H:Hasher>(&self,h:&mut H){h.write_u64(self.hash);}}
impl PartialEq for Key {
 fn eq(&self,other:&Self)->bool{
  fn same(a:&Key,b:&Key,seen:&mut FxHashSet<(usize,usize)>)->bool{
   if std::ptr::eq(a,b){return true;}if a.hash!=b.hash{return false;}
   if !seen.insert((a as *const _ as usize,b as *const _ as usize)){return true;}
   match (&a.data,&b.data){
    (Data::Free(a),Data::Free(b))=>a==b,(Data::Bound(a),Data::Bound(b))=>a==b,
    (Data::Node(ta,a),Data::Node(tb,b))=>ta==tb && a.len()==b.len() && a.iter().zip(b).all(|(a,b)|same(a,b,seen)),
    (Data::Function(fa,a),Data::Function(fb,b))=>fa==fb && a.len()==b.len() && a.iter().zip(b).all(|(a,b)|same(a,b,seen)),_=>false,
   }
  }
  same(self,other,&mut FxHashSet::default())
 }
}
impl Eq for Key{}
impl Ord for Key {fn cmp(&self,other:&Self)->Ordering{self.hash.cmp(&other.hash).then_with(||if self==other{Ordering::Equal}else{self.data.cmp(&other.data)})}}
impl PartialOrd for Key {fn partial_cmp(&self,other:&Self)->Option<Ordering>{Some(self.cmp(other))}}

struct Scope { vars:Vec<(Expr,u32)>, id:u32, mask:u64 }
struct Builder { cache:FxHashMap<(usize,u32),Arc<Key>>, next_scope:u32 }
fn safe(args:&[Expr])->bool {
    if args.len()!=3{return false;}
    let variables=args[1].args();
    variables.iter().all(Expr::is_symbol) && variables.iter().all(|v|!crate::traverse::has(&args[2],v))
}
pub(crate) fn build(args:&[Expr])->Option<Arc<Key>> {
    if !safe(args) || args[0].has_matrix(){return None;}
    fn unsupported(e:&Expr,seen:&mut FxHashMap<usize,Expr>)->bool {
        if seen.insert(e.dag_id(),e.clone()).is_some(){return false;}
        if matches!(e.kind(),Kind::Func(Func::Limit|Func::Integral|Func::Sum|Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::ConditionSet),_)){return true;}
        if matches!(e.kind(),Kind::Func(Func::Undefined(f),_) if f.is_internal("Integral") || f.is_internal("Sum")){return true;}
        e.args().iter().any(|arg|unsupported(arg,seen))
    }
    if unsupported(&args[0],&mut FxHashMap::default()){return None;}
    let mut builder=Builder{cache:FxHashMap::default(),next_scope:1};
    let scope=Scope{vars:Vec::new(),id:0,mask:0};
    Some(builder.substitution(args,&scope))
}
impl Builder {
    fn scope(&mut self,vars:Vec<(Expr,u32)>)->Scope {
        let id=self.next_scope;self.next_scope+=1;
        let mask=vars.iter().fold(0,|m,(v,_)|m|v.mask());
        Scope{vars,id,mask}
    }
    fn substitution(&mut self,args:&[Expr],outer:&Scope)->Arc<Key> {
        let variables=args[1].args();let points=args[2].args();
        let mut kept:Vec<_>=outer.vars.iter().filter(|(v,_)|!variables.contains(v)).cloned().collect();
        let point_keys=points.iter().map(|point|self.walk(point,outer)).collect();
        if !safe(args) || args[0].has_matrix(){
            let scope=self.scope(kept);
            let body=self.walk(&args[0],&scope);
            return Key::node(20,vec![body,Key::free(args[1].clone()),Key::node(21,point_keys)]);
        }
        let next=outer.vars.iter().map(|(_,i)|i+1).max().unwrap_or(0);
        kept.extend(variables.iter().enumerate().map(|(i,v)|(v.clone(),next+i as u32)));
        let scope=self.scope(kept);
        let body=self.walk(&args[0],&scope);
        // Unused variables are observable in SymPy's Subs identity.
        let mut unused:Vec<_>=variables.iter().zip(&points).filter(|(v,_)|!crate::traverse::has(&args[0],v)).map(|(v,p)|Key::node(22,vec![Key::free(v.clone()),self.walk(p,outer)])).collect();
        unused.sort();
        Key::node(19,vec![body,Key::node(21,point_keys),Key::node(23,unused)])
    }
    fn walk(&mut self,e:&Expr,scope:&Scope)->Arc<Key> {
        if let Some((_,index))=scope.vars.iter().rev().find(|(v,_)|v==e){return Key::bound(*index);}
        if e.is_atom() || e.mask() & scope.mask==0 {return Key::free(e.clone());}
        let cache_key=(e.dag_id(),scope.id);
        if let Some(key)=self.cache.get(&cache_key){return key.clone();}
        let key=match e.kind(){
            Kind::Add(c,terms)=>{
                let mut terms:Vec<_>=terms.iter().map(|(term,c)|Key::node(10,vec![self.walk(term,scope),Key::free(Expr::num(c.clone()))])).collect();terms.sort();
                terms.push(Key::free(Expr::num(c.clone())));Key::node(0,terms)
            },
            Kind::Mul(c,factors)=>{
                let mut factors:Vec<_>=factors.iter().map(|(base,power)|Key::node(10,vec![self.walk(base,scope),self.walk(power,scope)])).collect();factors.sort();
                factors.push(Key::free(Expr::num(c.clone())));Key::node(1,factors)
            },
            Kind::Pow(base,power)=>Key::node(2,vec![self.walk(base,scope),self.walk(power,scope)]),
            Kind::Func(Func::Subs,args)=>self.substitution(args,scope),
            Kind::Func(f,args)=>Key::function(f.clone(),args.iter().map(|a|self.walk(a,scope)).collect()),
            Kind::Derivative(inner,vars)=>{
                let mut children=vec![self.walk(inner,scope)];
                children.extend(vars.iter().map(|(v,n)|Key::node(10,vec![self.walk(v,scope),Key::free(Expr::int(*n as i64))])));Key::node(3,children)
            },
            Kind::Rel(op,left,right)=>Key::node(30+*op as u8,vec![self.walk(left,scope),self.walk(right,scope)]),
            Kind::Piecewise(pieces)=>Key::node(4,pieces.iter().flat_map(|(e,c)|[self.walk(e,scope),self.walk(c,scope)]).collect()),
            _=>Key::free(e.clone()),
        };
        self.cache.insert(cache_key,key.clone());key
    }
}
