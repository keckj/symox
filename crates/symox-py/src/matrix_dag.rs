//! Flat, versioned pickle records for the mixed scalar/matrix DAG.
//! No textual printer or evaluating public constructor participates in the codec.
use std::sync::Arc;
use pyo3::prelude::*;
use pyo3::BoundObject;
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyList,PyTuple};
use rustc_hash::FxHashMap;
use symox::expr::{Expr,Kind,Func,Constant,RelOp,UndefinedFunction};
use symox::number::Number;
use symox::symbol::{Symbol,Assumptions};
use symox::matrix::{Matrix,MResult};
use symox::matrixexpr::{MatrixExpr,MatrixKind,MatrixOperand,BlockGrid,SliceSpec};

const VERSION:u32=1;
type Reduction<'py>=(Bound<'py,PyAny>,Bound<'py,PyTuple>);
fn invalid(message:impl Into<String>)->PyErr{PyValueError::new_err(format!("Invalid matrix DAG pickle: {}",message.into()))}
fn obj<'py,T:IntoPyObject<'py>>(py:Python<'py>,value:T)->PyResult<PyObject> where T::Error:Into<PyErr>{Ok(value.into_pyobject(py).map_err(Into::into)?.into_bound().into_any().unbind())}
fn tuple(py:Python<'_>,items:Vec<PyObject>)->PyResult<PyObject>{Ok(PyTuple::new(py,items)?.into_any().unbind())}
fn number(py:Python<'_>,value:&Number)->PyResult<PyObject>{
    let fields=match value{
        Number::Int(n)=>vec![obj(py,"int")?,obj(py,*n)?],
        Number::BigInt(n)=>vec![obj(py,"bigint")?,obj(py,n.clone())?],
        Number::Rational(n)=>vec![obj(py,"rational")?,obj(py,n.numer().clone())?,obj(py,n.denom().clone())?],
        Number::Float(n)=>vec![obj(py,"float")?,obj(py,n.to_bits())?],
        Number::MpFloat(n,p)=>vec![obj(py,"mpfloat")?,obj(py,n.man.clone())?,obj(py,n.exp)?,obj(py,*p)?],
    };tuple(py,fields)
}
fn read_number(value:&Bound<'_,PyAny>)->PyResult<Number>{
    let p=value.downcast::<PyTuple>().map_err(|_|invalid("number record must be a tuple"))?;
    if p.is_empty(){return Err(invalid("empty number record"));}
    let tag:String=p.get_item(0)?.extract()?;
    Ok(match (tag.as_str(),p.len()){
        ("int",2)=>Number::Int(p.get_item(1)?.extract()?),
        ("bigint",2)=>Number::from_bigint(p.get_item(1)?.extract()?),
        ("rational",3)=>Number::rational(p.get_item(1)?.extract()?,p.get_item(2)?.extract()?).ok_or_else(||invalid("zero rational denominator"))?,
        ("float",2)=>Number::Float(f64::from_bits(p.get_item(1)?.extract()?)),
        ("mpfloat",4)=>{let precision:u32=p.get_item(3)?.extract()?;if precision<2{return Err(invalid("invalid floating point precision"));}Number::MpFloat(symox::bigfloat::BigFloat{man:p.get_item(1)?.extract()?,exp:p.get_item(2)?.extract()?},precision)},
        _=>return Err(invalid("unknown number record")),
    })
}
#[derive(Clone)]
enum Node{S(Expr),M(MatrixExpr)}
impl Node{
    fn key(&self)->(u8,usize){match self{Self::S(e)=>(0,e.dag_id()),Self::M(e)=>(1,e.dag_id())}}
    fn dependencies(&self)->Vec<Self>{
        match self{
            Self::M(m)=>{
                let mut out:Vec<_>=m.children().into_iter().map(Self::M).collect();
                match m.kind(){
                    MatrixKind::Symbol(_)|MatrixKind::Zero|MatrixKind::One=>{out.push(Self::S(m.rows().clone()));out.push(Self::S(m.cols().clone()));},
                    MatrixKind::Identity=>out.push(Self::S(m.rows().clone())),
                    MatrixKind::ApplyFunction(function,_)=>out.push(Self::S(function.clone())),
                    MatrixKind::Function(lambda)=>{out.push(Self::S(m.rows().clone()));out.push(Self::S(m.cols().clone()));out.push(Self::S(lambda.clone()));},
                    MatrixKind::Explicit(a)=>out.extend(a.data().iter().cloned().map(Self::S)),
                    MatrixKind::Mul(c,_)|MatrixKind::Pow(_,c)=>out.push(Self::S(c.clone())),
                    MatrixKind::HadamardPower(a,b)=>for operand in [a,b]{if let MatrixOperand::Scalar(e)=operand{out.push(Self::S(e.clone()));}},
                    MatrixKind::Slice(_,rows,cols)=>out.extend(rows.fields().into_iter().chain(cols.fields()).map(Self::S)),
                    _=>{},
                }out
            },
            Self::S(e)=>{
                let mut out=Vec::new();match e.kind(){
                    Kind::Sym(s)=>out.extend(s.exclude().iter().cloned().map(Self::S)),
                    Kind::Add(_,terms)=>out.extend(terms.iter().map(|(e,_)|Self::S(e.clone()))),
                    Kind::Mul(_,terms)=>for(a,b)in terms{out.push(Self::S(a.clone()));out.push(Self::S(b.clone()));},
                    Kind::Pow(a,b)|Kind::Rel(_,a,b)=>{out.push(Self::S(a.clone()));out.push(Self::S(b.clone()));},
                    Kind::Func(f,args)=>{if let Func::MatrixDotProduct(a,b)=f{out.push(Self::M((**a).clone()));out.push(Self::M((**b).clone()));}if let Func::MatrixElement(m)|Func::MatrixProxy(m)|Func::MatrixTrace(m)|Func::MatrixDeterminant(m)=f{out.push(Self::M((**m).clone()));}out.extend(args.iter().cloned().map(Self::S));},
                    Kind::Derivative(e,vars)=>{out.push(Self::S(e.clone()));out.extend(vars.iter().map(|(e,_)|Self::S(e.clone())));},
                    Kind::Piecewise(pairs)=>for(a,b)in pairs{out.push(Self::S(a.clone()));out.push(Self::S(b.clone()));},
                    _=>{},
                }out
            }
        }
    }
}
struct Encoder<'py>{py:Python<'py>,ids:FxHashMap<(u8,usize),usize>,records:Vec<PyObject>}
impl<'py> Encoder<'py>{
    fn new(py:Python<'py>)->Self{Self{py,ids:FxHashMap::default(),records:Vec::new()}}
    fn scalar(&self,e:&Expr)->usize{self.ids[&(0,e.dag_id())]}
    fn matrix(&self,e:&MatrixExpr)->usize{self.ids[&(1,e.dag_id())]}
    fn scalars(&self,es:&[Expr])->Vec<usize>{es.iter().map(|e|self.scalar(e)).collect()}
    fn matrices(&self,es:&[MatrixExpr])->Vec<usize>{es.iter().map(|e|self.matrix(e)).collect()}
    fn push(&mut self,tag:&str,canonical:bool,fields:Vec<PyObject>)->PyResult<usize>{
        let id=self.records.len();self.records.push(tuple(self.py,vec![obj(self.py,tag)?,obj(self.py,canonical)?,tuple(self.py,fields)?])?);Ok(id)
    }
    fn visit(&mut self,root:Node)->PyResult<usize>{
        let key=root.key();let mut pending=vec![(root,false)];
        while let Some((node,ready))=pending.pop(){if self.ids.contains_key(&node.key()){continue;}
            if !ready{pending.push((node.clone(),true));for child in node.dependencies().into_iter().rev(){if !self.ids.contains_key(&child.key()){pending.push((child,false));}}continue;}
            let id=self.emit(&node)?;self.ids.insert(node.key(),id);
        }Ok(self.ids[&key])
    }
    fn emit(&mut self,node:&Node)->PyResult<usize>{
        let py=self.py;
        let(tag,canonical,fields)=match node{
            Node::M(m)=>{let mut tag=m.kind_name();let fields=match m.kind(){
                MatrixKind::Permutation(p)=>vec![obj(py,p.clone())?],
                MatrixKind::Permute(parent,p,axis)=>vec![obj(py,self.matrix(parent))?,obj(py,p.clone())?,obj(py,*axis)?],
                MatrixKind::Symbol(name)=>vec![obj(py,name.as_str())?,obj(py,self.scalar(m.rows()))?,obj(py,self.scalar(m.cols()))?],
                MatrixKind::ApplyFunction(function,parent)=>vec![obj(py,self.scalar(function))?,obj(py,self.matrix(parent))?],
                MatrixKind::Function(lambda)=>vec![obj(py,self.scalar(m.rows()))?,obj(py,self.scalar(m.cols()))?,obj(py,self.scalar(lambda))?],
                MatrixKind::Block(grid)=>vec![obj(py,grid.shape().0)?,obj(py,grid.shape().1)?,obj(py,self.matrices(grid.data()))?],
                MatrixKind::Identity=>vec![obj(py,self.scalar(m.rows()))?],
                MatrixKind::Zero|MatrixKind::One=>vec![obj(py,self.scalar(m.rows()))?,obj(py,self.scalar(m.cols()))?],
                MatrixKind::Explicit(a)=>{tag="m.explicit";vec![obj(py,a.rows())?,obj(py,a.cols())?,obj(py,self.scalars(a.data()))?]},
                MatrixKind::Add(a)|MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)|MatrixKind::BlockDiag(a)=>vec![obj(py,self.matrices(a))?],
                MatrixKind::Mul(c,a)=>vec![obj(py,self.scalar(c))?,obj(py,self.matrices(a))?],
                MatrixKind::Pow(a,b)=>vec![obj(py,self.matrix(a))?,obj(py,self.scalar(b))?],
                MatrixKind::Transpose(a)|MatrixKind::Adjoint(a)|MatrixKind::Inverse(a)|MatrixKind::Conjugate(a)|MatrixKind::Diag(a)|MatrixKind::DiagonalOf(a)|MatrixKind::Diagonal(a)=>vec![obj(py,self.matrix(a))?],
                MatrixKind::HadamardPower(a,b)=>{let operand=|a:&MatrixOperand|match a{MatrixOperand::Matrix(m)=>(true,self.matrix(m)),MatrixOperand::Scalar(s)=>(false,self.scalar(s))};vec![obj(py,operand(a))?,obj(py,operand(b))?]},
                MatrixKind::Slice(parent,rows,cols)=>{tag="m.slice";vec![obj(py,self.matrix(parent))?,obj(py,self.scalars(&rows.fields()))?,obj(py,rows.stop_is_omitted)?,obj(py,self.scalars(&cols.fields()))?,obj(py,cols.stop_is_omitted)?]},
            };(tag,false,fields)},
            Node::S(e)=>{let(tag,fields)=match e.kind(){
                Kind::Num(n)=>("s.number",vec![number(py,n)?]),
                Kind::Sym(s)=>match s.dummy_index(){
                    Some(index)=>("s.dummy.facts",vec![obj(py,s.name())?,obj(py,s.given_facts().bits())?,obj(py,index.clone())?]),
                    None=>("s.symbol.facts",vec![obj(py,s.name())?,obj(py,s.given_facts().bits())?,obj(py,s.is_wild())?,obj(py,self.scalars(s.exclude()))?]),
                },
                Kind::Const(c)=>("s.constant",vec![obj(py,c.name())?]),
                Kind::Bool(b)=>("s.bool",vec![obj(py,*b)?]),
                Kind::Add(c,terms)=>{let terms=terms.iter().map(|(t,c)|tuple(py,vec![obj(py,self.scalar(t))?,number(py,c)?])).collect::<PyResult<Vec<_>>>()?;("s.add",vec![number(py,c)?,obj(py,terms)?])},
                Kind::Mul(c,terms)=>("s.mul",vec![number(py,c)?,obj(py,terms.iter().map(|(a,b)|(self.scalar(a),self.scalar(b))).collect::<Vec<_>>())?]),
                Kind::Pow(a,b)=>("s.pow",vec![obj(py,self.scalar(a))?,obj(py,self.scalar(b))?]),
                Kind::Derivative(a,vars)=>("s.derivative",vec![obj(py,self.scalar(a))?,obj(py,vars.iter().map(|(v,n)|(self.scalar(v),*n)).collect::<Vec<_>>())?]),
                Kind::Rel(op,a,b)=>("s.rel",vec![obj(py,op.symbol())?,obj(py,self.scalar(a))?,obj(py,self.scalar(b))?]),
                Kind::Piecewise(pairs)=>("s.piecewise",vec![obj(py,pairs.iter().map(|(a,b)|(self.scalar(a),self.scalar(b))).collect::<Vec<_>>())?]),
                Kind::Func(f,args)=>{
                    let refs=obj(py,self.scalars(args))?;
                    match f{
                        Func::MatrixDotProduct(a,b)=>("s.dot",vec![obj(py,self.matrix(a))?,obj(py,self.matrix(b))?]),
                        Func::MatrixElement(m)|Func::MatrixProxy(m)|Func::MatrixTrace(m)|Func::MatrixDeterminant(m)=>{
                            let tag=match f{Func::MatrixElement(_)=>"s.element",Func::MatrixProxy(_)=>"s.proxy",Func::MatrixTrace(_)=>"s.trace",_=>"s.determinant"};(tag,vec![obj(py,self.matrix(m))?,refs])
                        },
                        Func::Predicate(p)=>("s.predicate",vec![obj(py,p.name())?,obj(py,p.is_undefined())?]),
                        Func::Undefined(head)=>("s.undefined.facts",vec![obj(py,head.name.as_ref())?,obj(py,head.head_facts().bits())?,obj(py,head.nargs.clone())?,obj(py,head.internal)?,refs]),
                        _=>("s.function",vec![obj(py,f.name())?,refs]),
                    }
                },
            };(tag,e.is_canonical(),fields)},
        };self.push(tag,canonical,fields)
    }
    fn finish(self,root:usize)->PyResult<Reduction<'py>>{Ok((self.py.import("symox._core")?.getattr("_restore_matrix_dag")?,(VERSION,self.records,root).into_pyobject(self.py)?))}
}
pub(crate) fn reduce_expr<'py>(py:Python<'py>,e:&Expr)->PyResult<Reduction<'py>>{let mut codec=Encoder::new(py);let root=codec.visit(Node::S(e.clone()))?;codec.finish(root)}
pub(crate) fn reduce_matrix<'py>(py:Python<'py>,m:&MatrixExpr)->PyResult<Reduction<'py>>{let mut codec=Encoder::new(py);let root=codec.visit(Node::M(m.clone()))?;codec.finish(root)}
pub(crate) fn reduce_grid<'py>(py:Python<'py>,grid:&BlockGrid)->PyResult<Reduction<'py>>{let mut codec=Encoder::new(py);let mut ids=Vec::new();for m in grid.data(){ids.push(codec.visit(Node::M(m.clone()))?);}let root=codec.push("grid",false,vec![obj(py,grid.shape().0)?,obj(py,grid.shape().1)?,obj(py,ids)?])?;codec.finish(root)}
pub(crate) fn reduce_dense<'py>(py:Python<'py>,matrix:&Matrix,immutable:bool)->PyResult<Reduction<'py>>{let mut codec=Encoder::new(py);let mut ids=Vec::new();for e in matrix.data(){ids.push(codec.visit(Node::S(e.clone()))?);}let root=codec.push("dense",false,vec![obj(py,matrix.rows())?,obj(py,matrix.cols())?,obj(py,ids)?,obj(py,immutable)?])?;codec.finish(root)}

enum Decoded{S(Expr),M(MatrixExpr),G(BlockGrid),D(Matrix,bool)}
struct Reader{nodes:Vec<Decoded>}
impl Reader{
    fn scalar(&self,id:usize)->PyResult<Expr>{match self.nodes.get(id){Some(Decoded::S(e))=>Ok(e.clone()),_=>Err(invalid("scalar reference must target an earlier scalar node"))}}
    fn matrix(&self,id:usize)->PyResult<MatrixExpr>{match self.nodes.get(id){Some(Decoded::M(m))=>Ok(m.clone()),_=>Err(invalid("matrix reference must target an earlier matrix node"))}}
    fn scalars(&self,ids:Vec<usize>)->PyResult<Vec<Expr>>{ids.into_iter().map(|id|self.scalar(id)).collect()}
    fn matrices(&self,ids:Vec<usize>)->PyResult<Vec<MatrixExpr>>{ids.into_iter().map(|id|self.matrix(id)).collect()}
}
fn assumptions(bits:u32)->PyResult<Assumptions>{if bits&!Assumptions::NAMES.iter().fold(0u32,|mask,(_,flag)|mask|flag.0)!=0{return Err(invalid("unknown assumption bits"));}let a=Assumptions(bits);a.closure().map_err(invalid)?;Ok(a)}
fn constant(name:&str)->PyResult<Constant>{Ok(match name{"pi"=>Constant::Pi,"E"=>Constant::E,"I"=>Constant::I,"EulerGamma"=>Constant::EulerGamma,"Catalan"=>Constant::Catalan,"GoldenRatio"=>Constant::GoldenRatio,"TribonacciConstant"=>Constant::TribonacciConstant,"oo"=>Constant::Infinity,"-oo"=>Constant::NegInfinity,"zoo"=>Constant::ComplexInfinity,"nan"=>Constant::NaN,_=>return Err(invalid("unknown constant"))})}
fn validate_function(function:&Func,args:&[Expr])->PyResult<()>{
    if function.arity().is_some_and(|n|n!=args.len()){return Err(invalid("function arity mismatch"));}
    if matches!(function,Func::Integral|Func::Sum|Func::Product){symox::binders::validate(function,args).map_err(invalid)?;}
    if let Func::Set(op @ (symox::sets::SetOp::Lambda|symox::sets::SetOp::IdentityFunction))=function{symox::sets::validate_lambda(*op,args).map_err(invalid)?;}
    if matches!(function,Func::AppliedPredicate){symox::ask::applied_from_args(args.to_vec()).map_err(|e|invalid(e.to_string()))?;}
    if matches!(function,Func::Subs){
        if args.len()!=3{return Err(invalid("Subs requires three arguments"));}
        let tuple=|e:&Expr|matches!(e.kind(),Kind::Func(Func::Undefined(head),_) if head.is_internal("Tuple"));
        if !tuple(&args[1])||!tuple(&args[2]){return Err(invalid("Subs variables and points must be tuples"));}
        let variables=args[1].args();let points=args[2].args();
        if variables.is_empty()||variables.len()!=points.len()||variables.iter().enumerate().any(|(i,v)|!v.is_symbol()||variables[..i].contains(v)){return Err(invalid("invalid Subs variables or points"));}
    }
    if function.name()=="Limit" {
        if args.len()!=4||!args[1].is_symbol()||!matches!(args[3].kind(),Kind::Sym(s) if matches!(s.name(),"+"|"-"|"+-")){return Err(invalid("invalid held Limit arguments"));}
    }
    Ok(())
}
fn matrix_result(result:MResult<MatrixExpr>)->PyResult<Decoded>{Ok(Decoded::M(result.map_err(|e|invalid(e.to_string()))?))}
#[pyfunction]
fn _restore_matrix_dag(py:Python<'_>,version:u32,records:&Bound<'_,PyList>,root:usize)->PyResult<PyObject>{
    if version!=VERSION{return Err(invalid("unsupported format version"));}
    let mut reader=Reader{nodes:Vec::with_capacity(records.len())};
    for record in records.iter(){
        let record=record.downcast::<PyTuple>().map_err(|_|invalid("node record must be a tuple"))?;
        if record.len()!=3{return Err(invalid("node record must contain tag, canonical flag and fields"));}
        let tag:String=record.get_item(0)?.extract()?;let canonical:bool=record.get_item(1)?.extract()?;
        let f=record.get_item(2)?.downcast_into::<PyTuple>().map_err(|_|invalid("node fields must be a tuple"))?;
        let arity=|n:usize|->PyResult<()>{if f.len()==n{Ok(())}else{Err(invalid(format!("{tag} expects {n} fields")))}};
        let s=|i|reader.scalar(f.get_item(i)?.extract()?);let m=|i|reader.matrix(f.get_item(i)?.extract()?);
        let ss=|i|reader.scalars(f.get_item(i)?.extract()?);let mm=|i|reader.matrices(f.get_item(i)?.extract()?);
        let node=if tag.starts_with("s."){
            let kind=match tag.as_str(){
                "s.number"=>{arity(1)?;Kind::Num(read_number(&f.get_item(0)?)?)},
                "s.symbol"=>{arity(4)?;let name:String=f.get_item(0)?.extract()?;let given=assumptions(f.get_item(1)?.extract()?)?;let wild:bool=f.get_item(2)?.extract()?;let exclude=ss(3)?;if !wild&&!exclude.is_empty(){return Err(invalid("ordinary symbol cannot have wildcard exclusions"));}Kind::Sym(if wild{Symbol::wild_with_facts(&name,exclude,Symbol::with_assumptions(&name,given).given_facts()).map_err(invalid)?}else{Symbol::with_assumptions(&name,given)})},
                "s.dummy"=>{arity(3)?;Kind::Sym(Symbol::dummy_with_index(&f.get_item(0)?.extract::<String>()?,assumptions(f.get_item(1)?.extract()?)?,f.get_item(2)?.extract()?))},
                "s.symbol.facts"=>{arity(4)?;let name:String=f.get_item(0)?.extract()?;let(known,truth):(u32,u32)=f.get_item(1)?.extract()?;let given=symox::assumption_facts::FactSet::from_bits(known,truth).map_err(invalid)?;let wild:bool=f.get_item(2)?.extract()?;let exclude=ss(3)?;if !wild&&!exclude.is_empty(){return Err(invalid("ordinary symbol cannot have wildcard exclusions"));}Kind::Sym(if wild{Symbol::wild_with_facts(&name,exclude,given).map_err(invalid)?}else{Symbol::with_facts(&name,given).map_err(invalid)?})},
                "s.dummy.facts"=>{arity(3)?;let(known,truth):(u32,u32)=f.get_item(1)?.extract()?;let given=symox::assumption_facts::FactSet::from_bits(known,truth).map_err(invalid)?;Kind::Sym(Symbol::dummy_with_facts(&f.get_item(0)?.extract::<String>()?,given,Some(f.get_item(2)?.extract()?)).map_err(invalid)?)},
                "s.constant"=>{arity(1)?;Kind::Const(constant(&f.get_item(0)?.extract::<String>()?)?)},
                "s.bool"=>{arity(1)?;Kind::Bool(f.get_item(0)?.extract()?)},
                "s.add"=>{arity(2)?;let mut terms=Vec::new();for pair in f.get_item(1)?.try_iter()?{let pair=pair?;if pair.len()?!=2{return Err(invalid("additive term requires reference and coefficient"));}terms.push((reader.scalar(pair.get_item(0)?.extract()?)?,read_number(&pair.get_item(1)?)?));}Kind::Add(read_number(&f.get_item(0)?)?,terms)},
                "s.mul"=>{arity(2)?;let pairs:Vec<(usize,usize)>=f.get_item(1)?.extract()?;Kind::Mul(read_number(&f.get_item(0)?)?,pairs.into_iter().map(|(a,b)|Ok((reader.scalar(a)?,reader.scalar(b)?))).collect::<PyResult<_>>()?)},
                "s.pow"=>{arity(2)?;Kind::Pow(s(0)?,s(1)?)},
                "s.derivative"=>{arity(2)?;let vars:Vec<(usize,u32)>=f.get_item(1)?.extract()?;if vars.iter().any(|(_,n)|*n==0){return Err(invalid("derivative order must be positive"));}Kind::Derivative(s(0)?,vars.into_iter().map(|(id,n)|Ok((reader.scalar(id)?,n))).collect::<PyResult<_>>()?)},
                "s.rel"=>{arity(3)?;let op=match f.get_item(0)?.extract::<String>()?.as_str(){"=="=>RelOp::Eq,"!="=>RelOp::Ne,"<"=>RelOp::Lt,"<="=>RelOp::Le,">"=>RelOp::Gt,">="=>RelOp::Ge,_=>return Err(invalid("unknown relation"))};Kind::Rel(op,s(1)?,s(2)?)},
                "s.piecewise"=>{arity(1)?;let pairs:Vec<(usize,usize)>=f.get_item(0)?.extract()?;Kind::Piecewise(pairs.into_iter().map(|(a,b)|Ok((reader.scalar(a)?,reader.scalar(b)?))).collect::<PyResult<_>>()?)},
                "s.undefined"|"s.undefined.facts"=>{arity(5)?;
                    let name:String=f.get_item(0)?.extract()?;
                    let given=if tag=="s.undefined"{symox::assumption_facts::FactSet::from_legacy(assumptions(f.get_item(1)?.extract()?)?)}else{
                        let(known,truth):(u32,u32)=f.get_item(1)?.extract()?;
                        symox::assumption_facts::FactSet::from_bits(known,truth).map_err(invalid)?
                    };
                    let mut head=UndefinedFunction::with_facts(&name,given,f.get_item(2)?.extract()?).map_err(invalid)?;
                    head.internal=f.get_item(3)?.extract()?;
                    if head.internal && (!matches!(name.as_str(),"Tuple"|"UnevaluatedExpr"|"Integral"|"Sum") || !head.head_facts().is_empty() || head.nargs.is_some()){
                        return Err(invalid("invalid internal function head"));
                    }
                    let args=ss(4)?;if !head.accepts(args.len()) || (head.is_internal("UnevaluatedExpr")&&args.len()!=1){return Err(invalid("undefined function arity mismatch"));}
                    if head.internal&&matches!(head.name.as_ref(),"Integral"|"Sum"){
                        if args.len()<2||head.assumptions!=Assumptions::NONE||head.nargs.is_some(){return Err(invalid("invalid legacy calculus placeholder"));}
                        let function=if head.name.as_ref()=="Integral"{Func::Integral}else{Func::Sum};
                        let limits=args[1..].iter().map(|limit|{
                            if function==Func::Integral&&limit.is_symbol(){Ok(vec![limit.clone()])}
                            else if matches!(limit.kind(),Kind::Func(Func::Undefined(head),_) if head.is_internal("Tuple")){Ok(limit.args())}
                            else{Err(invalid("invalid legacy calculus limit"))}
                        }).collect::<PyResult<Vec<_>>>()?;
                        symox::binders::held(function,&args[0],limits).map_err(invalid)?.kind().clone()
                    }else{Kind::Func(Func::Undefined(head),args)}
                },
                "s.predicate"=>{arity(2)?;let name:String=f.get_item(0)?.extract()?;let undefined:bool=f.get_item(1)?.extract()?;
                    let p=if undefined{symox::ask::Pred::undefined(&name)}else{symox::ask::Pred::from_name(&name).ok_or_else(||invalid("unknown predicate"))?};
                    Kind::Func(Func::Predicate(p),vec![])
                },
                "s.function"=>{arity(2)?;let name:String=f.get_item(0)?.extract()?;let function=Func::from_name(&name).ok_or_else(||invalid("unknown function"))?;let args=ss(1)?;
                    if matches!(function,Func::Set(symox::sets::SetOp::Lambda)){
                        if args.len()!=2{return Err(invalid("Lambda requires a signature and a body"));}
                        // Read both old [Symbol, body] and normalized [Tuple, body]
                        // records. User functions named Lambda use s.undefined.
                        symox::sets::lambda_checked(&args[0],&args[1]).map_err(invalid)?.kind().clone()
                    }else{validate_function(&function,&args)?;Kind::Func(function,args)}
                },
                "s.dot"=>{arity(2)?;symox::matrixcontraction::dot_product(&m(0)?,&m(1)?,false).map_err(|e|invalid(e.to_string()))?.kind().clone()},
                "s.element"|"s.proxy"|"s.trace"|"s.determinant"=>{
                    arity(2)?;let matrix=m(0)?;let args=ss(1)?;let n=if tag=="s.element"{2}else{0};if args.len()!=n{return Err(invalid("matrix scalar arity mismatch"));}
                    if tag=="s.element"{matrix.element_unevaluated(&args[0],&args[1]).map_err(|e|invalid(e.to_string()))?;}
                    if tag=="s.trace"||tag=="s.determinant"{matrix.trace_unevaluated().map_err(|e|invalid(e.to_string()))?;}
                    let parent=Arc::new(matrix);let function=match tag.as_str(){"s.element"=>Func::MatrixElement(parent),"s.proxy"=>Func::MatrixProxy(parent),"s.trace"=>Func::MatrixTrace(parent),_=>Func::MatrixDeterminant(parent)};Kind::Func(function,args)
                },
                _=>return Err(invalid("unknown scalar node tag")),
            };Decoded::S(if canonical{Expr::from_kind(kind)}else{Expr::from_unevaluated_kind(kind)})
        }else{
            match tag.as_str(){
                "PermutationMatrix"=>{arity(1)?;matrix_result(MatrixExpr::permutation(f.get_item(0)?.extract()?,false))?},
                "MatrixPermute"=>{arity(3)?;let parent=m(0)?;let images:Vec<usize>=f.get_item(1)?.extract()?;let axis:u8=f.get_item(2)?.extract()?;if axis>1||(if axis==0{parent.rows()}else{parent.cols()}).as_i64()!=Some(images.len() as i64){return Err(invalid("permutation axis or size mismatch"));}matrix_result(parent.permute(images,axis,false))?},
                "MatrixSymbol"=>{arity(3)?;matrix_result(MatrixExpr::symbol(&f.get_item(0)?.extract::<String>()?,s(1)?,s(2)?))?},
                "Identity"=>{arity(1)?;matrix_result(MatrixExpr::identity(s(0)?))?},
                "ElementwiseApplyFunction"=>{arity(2)?;matrix_result(MatrixExpr::apply_function(s(0)?,&m(1)?,false))?},
                "FunctionMatrix"=>{arity(3)?;matrix_result(MatrixExpr::function(s(0)?,s(1)?,s(2)?))?},
                "ZeroMatrix"=>{arity(2)?;matrix_result(MatrixExpr::zero(s(0)?,s(1)?))?},
                "OneMatrix"=>{arity(2)?;matrix_result(MatrixExpr::one(s(0)?,s(1)?))?},
                "m.explicit"|"dense"=>{arity(if tag=="dense"{4}else{3})?;let rows:usize=f.get_item(0)?.extract()?;let cols:usize=f.get_item(1)?.extract()?;let entries=ss(2)?;if rows.checked_mul(cols)!=Some(entries.len()){return Err(invalid("dense shape does not match its entries"));}let matrix=Matrix::new(rows,cols,entries).map_err(|e|invalid(e.to_string()))?;if tag=="dense"{Decoded::D(matrix,f.get_item(3)?.extract()?)}else{Decoded::M(MatrixExpr::explicit(matrix))}},
                "MatAdd"|"HadamardProduct"|"KroneckerProduct"|"BlockDiagMatrix"=>{arity(1)?;let children=mm(0)?;matrix_result(match tag.as_str(){"MatAdd"=>MatrixExpr::add(&children,false),"HadamardProduct"=>MatrixExpr::hadamard_product(&children,false),"KroneckerProduct"=>MatrixExpr::kronecker(&children,false),_=>MatrixExpr::block_diag(&children)})?},
                "BlockMatrix"=>{arity(3)?;matrix_result(MatrixExpr::block(BlockGrid::new(f.get_item(0)?.extract()?,f.get_item(1)?.extract()?,mm(2)?).map_err(|e|invalid(e.to_string()))?))?},
                "MatMul"=>{arity(2)?;matrix_result(MatrixExpr::mul(s(0)?,&mm(1)?,false))?},
                "MatPow"=>{arity(2)?;matrix_result(m(0)?.pow(s(1)?,false))?},
                "Transpose"|"Adjoint"|"Inverse"|"MatrixConjugate"=>{arity(1)?;let matrix=m(0)?;matrix_result(match tag.as_str(){"Transpose"=>matrix.transpose(false),"Adjoint"=>matrix.adjoint(false),"Inverse"=>matrix.inverse(false),_=>matrix.conjugate(false)})?},
                "DiagMatrix"|"DiagonalOf"|"DiagonalMatrix"=>{arity(1)?;let parent=m(0)?;matrix_result(match tag.as_str(){"DiagMatrix"=>MatrixExpr::diag(&parent),"DiagonalOf"=>MatrixExpr::diagonal_of(&parent),_=>MatrixExpr::diagonal(&parent)})?},
                "HadamardPower"=>{arity(2)?;let operand=|i|->PyResult<MatrixOperand>{let(is_matrix,id):(bool,usize)=f.get_item(i)?.extract()?;Ok(if is_matrix{MatrixOperand::Matrix(reader.matrix(id)?)}else{MatrixOperand::Scalar(reader.scalar(id)?)})};matrix_result(MatrixExpr::hadamard_power(operand(0)?,operand(1)?,false))?},
                "MatrixSlice"=>{arity(3)?;let rows=ss(1)?;let cols=ss(2)?;if rows.len()!=3||cols.len()!=3{return Err(invalid("slice specifications need three values"));}matrix_result(m(0)?.slice([Some(rows[0].clone()),Some(rows[1].clone()),Some(rows[2].clone())],[Some(cols[0].clone()),Some(cols[1].clone()),Some(cols[2].clone())]))?},
                "m.slice"=>{arity(5)?;let spec=|i,j|->PyResult<SliceSpec>{let values=ss(i)?;let [start,stop,step,length]:[Expr;4]=values.try_into().map_err(|_|invalid("normalized slice needs four scalar fields"))?;SliceSpec::normalized(start,stop,step,length,f.get_item(j)?.extract()?).map_err(|e|invalid(e.to_string()))};matrix_result(m(0)?.slice_normalized(spec(1,2)?,spec(3,4)?))?},
                "grid"=>{arity(3)?;Decoded::G(BlockGrid::new(f.get_item(0)?.extract()?,f.get_item(1)?.extract()?,mm(2)?).map_err(|e|invalid(e.to_string()))?)},
                _=>return Err(invalid("unknown matrix node tag")),
            }
        };reader.nodes.push(node);
    }
    match reader.nodes.get(root).ok_or_else(||invalid("root reference out of range"))?{
        Decoded::S(e)=>{
            if matches!(e.kind(),Kind::Func(Func::Set(symox::sets::SetOp::IdentityFunction),_)){
                Ok(py.import("symox")?.getattr("IdentityFunction")?.call0()?.unbind())
            }else{crate::expr_object(py,e.clone())}
        },
        Decoded::M(m)=>Ok(crate::matrixexpr::wrap_matrixexpr(m.clone()).into_pyobject(py)?.into_any().unbind()),
        Decoded::G(grid)=>Ok(crate::matrixexpr::PyBlockGrid{grid:grid.clone()}.into_pyobject(py)?.into_any().unbind()),
        Decoded::D(matrix,immutable)=>{let value=crate::matrices::wrapm(matrix.clone()).into_pyobject(py)?;if *immutable{Ok(value.call_method0("as_immutable")?.unbind())}else{Ok(value.into_any().unbind())}},
    }
}
pub(crate) fn register(module:&Bound<'_,PyModule>)->PyResult<()>{module.add_function(wrap_pyfunction!(_restore_matrix_dag,module)?)?;Ok(())}
