//! Python interface to the native typed MatrixExpr DAG.
use pyo3::prelude::*;
use pyo3::types::{PyDict,PyList,PySet,PyTuple};
use pyo3::exceptions::{PyTypeError,PyNotImplementedError};
use rustc_hash::FxHashMap;
use symox::expr::Expr;
use symox::matrixexpr::{MatrixExpr,MatrixKind};
use crate::{wrap,to_expr,PyExpr};

#[pyclass(name="_MatrixExpr",module="symox._core",frozen)]
#[derive(Clone)]
pub(crate) struct PyMatrixExpr {pub m:MatrixExpr}
pub(crate) fn wrap_matrixexpr(m:MatrixExpr)->PyMatrixExpr{PyMatrixExpr{m}}
fn error(py:Python<'_>,e:symox::matrix::MatrixError)->PyErr{crate::matrices::err(py,e)}
fn matrix_object<'py>(py:Python<'py>,m:&MatrixExpr)->PyResult<Bound<'py,PyAny>>{
    if let MatrixKind::Explicit(value)=m.kind(){return Ok(crate::matrices::wrapm(value.clone()).into_pyobject(py)?.call_method0("as_immutable")?);}
    Ok(wrap_matrixexpr(m.clone()).into_pyobject(py)?.into_any())
}
fn is_matrix(value:&Bound<'_,PyAny>)->bool{
    value.downcast::<PyMatrixExpr>().is_ok() || value.downcast::<crate::matrices::PyMatrix>().is_ok()
        || value.getattr("is_Matrix").is_ok_and(|v|v.is_truthy().unwrap_or(false))
        || value.getattr("_symox_sparse").is_ok_and(|v|v.is_truthy().unwrap_or(false))
}
pub(crate) fn to_matrixexpr(value:&Bound<'_,PyAny>)->PyResult<MatrixExpr>{
    if let Ok(v)=value.downcast::<PyMatrixExpr>(){return Ok(v.get().m.clone());}
    if value.downcast::<crate::matrices::PyMatrix>().is_ok() || value.getattr("_symox_sparse").is_ok_and(|v|v.is_truthy().unwrap_or(false)){
        return Ok(MatrixExpr::explicit(crate::matrices::to_matrix(value)?));
    }
    if value.getattr("is_MatrixExpr").is_ok_and(|v|v.is_truthy().unwrap_or(false)){
        let head:String=value.getattr("func")?.getattr("__name__")?.extract()?;
        let args=value.getattr("args")?;
        let args=args.downcast::<PyTuple>()?;
        if head=="MatrixSymbol"{
            let name:String=value.getattr("name")?.extract()?;
            return MatrixExpr::symbol(&name,to_expr(&value.getattr("rows")?)?,to_expr(&value.getattr("cols")?)?).map_err(|e|error(value.py(),e));
        }
        return construct(&head,args,false).and_then(|result|{
            result.bind(value.py()).downcast::<PyMatrixExpr>().map(|v|v.get().m.clone())
                .map_err(|_|PyTypeError::new_err("Expected a matrix expression"))
        });
    }
    if is_matrix(value){return Ok(MatrixExpr::explicit(crate::matrices::to_matrix(value)?));}
    Err(PyTypeError::new_err("Expected a matrix expression or explicit matrix"))
}
pub(crate) fn matrixexpr_to_sympy<'py>(py:Python<'py>,m:&MatrixExpr)->PyResult<Bound<'py,PyAny>>{
    let sympy=py.import("sympy")?;
    if let MatrixKind::Explicit(value)=m.kind(){return matrix_object(py,&MatrixExpr::explicit(value.clone()))?.call_method0("to_sympy");}
    if let MatrixKind::Conjugate(parent)=m.kind(){return matrixexpr_to_sympy(py,parent)?.call_method0("conjugate");}
    let args=matrix_args(py,m,true)?;
    let kwargs=PyDict::new(py);
    if matches!(m.kind(),MatrixKind::Add(_)|MatrixKind::Mul(..)|MatrixKind::Pow(..)){kwargs.set_item("evaluate",false)?;}
    sympy.getattr(m.kind_name())?.call(args,Some(&kwargs))
}
fn matrix_args<'py>(py:Python<'py>,m:&MatrixExpr,sympy:bool)->PyResult<Bound<'py,PyTuple>>{
    let scalar=|e:&Expr|->PyResult<Bound<'py,PyAny>>{if sympy{crate::expr_to_sympy(py,e)}else{Ok(wrap(e.clone()).into_pyobject(py)?.into_any())}};
    let matrix=|e:&MatrixExpr|if sympy{matrixexpr_to_sympy(py,e)}else{matrix_object(py,e)};
    let mut args=Vec::new();
    match m.kind(){
        MatrixKind::Symbol(name)=>{args.push(name.into_pyobject(py)?.into_any());args.push(scalar(m.rows())?);args.push(scalar(m.cols())?);},
        MatrixKind::Identity=>args.push(scalar(m.rows())?),MatrixKind::Zero=>{args.push(scalar(m.rows())?);args.push(scalar(m.cols())?);},
        MatrixKind::Add(a)=>for e in a{args.push(matrix(e)?);},
        MatrixKind::Mul(c,a)=>{if !c.is_one(){args.push(scalar(c)?);}for e in a{args.push(matrix(e)?);}},
        MatrixKind::Pow(b,e)=>{args.push(matrix(b)?);args.push(scalar(e)?);},
        MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Inverse(b)|MatrixKind::Conjugate(b)=>args.push(matrix(b)?),
        MatrixKind::Explicit(value)=>{args.push(value.rows().into_pyobject(py)?.into_any());args.push(value.cols().into_pyobject(py)?.into_any());args.push(PyList::new(py,value.data().iter().map(|e|wrap(e.clone())))?.into_any());},
    }
    PyTuple::new(py,args)
}
fn construct(kind:&str,args:&Bound<'_,PyTuple>,evaluate:bool)->PyResult<PyObject>{
    let py=args.py();let fail=||PyTypeError::new_err(format!("Invalid {kind} arguments"));
    let matrix=|i|to_matrixexpr(&args.get_item(i)?);
    let scalar=|i|to_expr(&args.get_item(i)?);
    let result=match kind{
        "MatrixSymbol" if args.len()==3=>{let name:String=args.get_item(0)?.extract()?;MatrixExpr::symbol(&name,scalar(1)?,scalar(2)?)},
        "Identity" if args.len()==1=>MatrixExpr::identity(scalar(0)?),"ZeroMatrix" if args.len()==2=>MatrixExpr::zero(scalar(0)?,scalar(1)?),
        "MatAdd"=>{if args.is_empty(){return Ok(wrap(Expr::zero()).into_pyobject(py)?.into_any().unbind());}let terms=args.iter().map(|v|to_matrixexpr(&v)).collect::<PyResult<Vec<_>>>()?;MatrixExpr::add(&terms,evaluate)},
        "MatMul"=>{let mut coefficient=Expr::one();let mut factors=Vec::new();for arg in args.iter(){if is_matrix(&arg){factors.push(to_matrixexpr(&arg)?);}else{coefficient=symox::ops::mul2(&coefficient,&to_expr(&arg)?);}}if factors.is_empty(){return Ok(wrap(coefficient).into_pyobject(py)?.into_any().unbind());}MatrixExpr::mul(coefficient,&factors,evaluate)},
        "MatPow" if args.len()==2=>matrix(0)?.pow(scalar(1)?,evaluate),"Transpose" if args.len()==1=>matrix(0)?.transpose(evaluate),"Adjoint" if args.len()==1=>matrix(0)?.adjoint(evaluate),"Inverse" if args.len()==1 || (args.len()==2 && scalar(1)?==Expr::minus_one())=>matrix(0)?.inverse(evaluate),"MatrixConjugate" if args.len()==1=>matrix(0)?.conjugate(evaluate),
        _=>return Err(fail()),
    }.map_err(|e|error(py,e))?;
    Ok(wrap_matrixexpr(result).into_pyobject(py)?.into_any().unbind())
}
#[pyfunction]
#[pyo3(signature=(kind,args,evaluate=false))]
fn _matrixexpr_node(kind:&str,args:&Bound<'_,PyTuple>,evaluate:bool)->PyResult<PyObject>{construct(kind,args,evaluate)}
#[pyfunction]
fn _matrix_element(parent:&Bound<'_,PyAny>,i:&Bound<'_,PyAny>,j:&Bound<'_,PyAny>)->PyResult<PyExpr>{
    let m=to_matrixexpr(parent)?;Ok(wrap(m.element(&to_expr(i)?,&to_expr(j)?).map_err(|e|error(parent.py(),e))?))
}
#[pyfunction]
fn _as_matrixexpr(value:&Bound<'_,PyAny>)->PyResult<PyMatrixExpr>{Ok(wrap_matrixexpr(to_matrixexpr(value)?))}
#[pymethods]
impl PyMatrixExpr{
    #[getter] fn _symox_matrixexpr(&self)->bool{true}
    #[getter] fn is_Matrix(&self)->bool{true}
    #[getter] fn is_MatrixExpr(&self)->bool{true}
    #[getter] fn is_commutative(&self)->bool{false}
    #[getter] fn is_Symbol(&self)->bool{false}
    #[getter] fn is_Atom(&self)->bool{false}
    #[getter] fn is_Number(&self)->bool{false}
    #[getter] fn is_MatrixSymbol(&self)->bool{matches!(self.m.kind(),MatrixKind::Symbol(_))}
    #[getter] fn is_ZeroMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::Zero)}
    #[getter] fn is_Identity(&self)->bool{matches!(self.m.kind(),MatrixKind::Identity)}
    #[getter] fn is_MatAdd(&self)->bool{matches!(self.m.kind(),MatrixKind::Add(_))}
    #[getter] fn is_MatMul(&self)->bool{matches!(self.m.kind(),MatrixKind::Mul(..))}
    #[getter] fn is_MatPow(&self)->bool{matches!(self.m.kind(),MatrixKind::Pow(..)|MatrixKind::Inverse(_))}
    #[getter] fn kind(&self)->&'static str{self.m.kind_name()}
    #[getter] fn rows(&self)->PyExpr{wrap(self.m.rows().clone())}
    #[getter] fn cols(&self)->PyExpr{wrap(self.m.cols().clone())}
    #[getter] fn shape(&self)->(PyExpr,PyExpr){(self.rows(),self.cols())}
    #[getter] fn is_square(&self)->Option<bool>{if self.m.rows()==self.m.cols(){Some(true)}else{symox::matrix::ratfunc::is_zero_assume(&symox::ops::sub(self.m.rows(),self.m.cols()))}}
    #[getter] fn name(&self)->PyResult<String>{if let MatrixKind::Symbol(name)=self.m.kind(){Ok(name.clone())}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only MatrixSymbol has a name"))}}
    #[getter] fn base(&self)->PyResult<PyMatrixExpr>{match self.m.kind(){MatrixKind::Pow(b,_)|MatrixKind::Inverse(b)=>Ok(wrap_matrixexpr(b.clone())),_=>Err(pyo3::exceptions::PyAttributeError::new_err("Only matrix powers have a base"))}}
    #[getter] fn exp(&self)->PyResult<PyExpr>{match self.m.kind(){MatrixKind::Pow(_,e)=>Ok(wrap(e.clone())),MatrixKind::Inverse(_)=>Ok(wrap(Expr::minus_one())),_=>Err(pyo3::exceptions::PyAttributeError::new_err("Only matrix powers have an exponent"))}}
    #[getter] fn arg(&self)->PyResult<PyMatrixExpr>{let c=self.m.children();if c.len()==1{Ok(wrap_matrixexpr(c[0].clone()))}else{Err(pyo3::exceptions::PyAttributeError::new_err("Matrix expression does not have one matrix argument"))}}
    #[getter] fn args<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{matrix_args(py,&self.m,false)}
    #[getter] fn func<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{py.import("symox")?.getattr(self.m.kind_name())}
    #[getter] fn free_symbols<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PySet>>{
        let result=PySet::empty(py)?;for m in self.m.matrix_symbols(){result.add(wrap_matrixexpr(m))?;}
        for e in self.m.scalar_free_symbols(){result.add(crate::expr_object(py,e)?)?;}Ok(result)
    }
    fn __str__(&self)->String{self.m.to_str()}
    fn __repr__(&self)->String{self.m.to_str()}
    fn srepr(&self)->String{self.m.to_repr()}
    fn __hash__(&self)->isize{
        if let MatrixKind::Explicit(m)=self.m.kind(){
            let entries=m.data().iter().enumerate().filter(|(_,e)|!e.is_zero()).map(|(i,e)|((i/m.cols(),i%m.cols()),e.clone())).collect();
            return crate::matrices::matrix_value_hash(m.rows(),m.cols(),entries);
        }
        use std::hash::{Hash,Hasher};let mut h=std::collections::hash_map::DefaultHasher::new();self.m.hash(&mut h);let v=h.finish() as isize;if v== -1{-2}else{v}
    }
    fn __richcmp__(&self,other:&Bound<'_,PyAny>,op:pyo3::basic::CompareOp)->PyObject{
        use pyo3::basic::CompareOp::*;let py=other.py();let same=to_matrixexpr(other).is_ok_and(|m|m==self.m);
        match op{Eq=>same.into_pyobject(py).unwrap().to_owned().into_any().unbind(),Ne=>(!same).into_pyobject(py).unwrap().to_owned().into_any().unbind(),_=>py.NotImplemented()}
    }
    fn __add__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,false,false)}
    fn __radd__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,false,true)}
    fn __sub__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,true,false)}
    fn __rsub__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,true,true)}
    fn __mul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_mul(other,false)}
    fn __rmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_mul(other,true)}
    fn __matmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if !is_matrix(other){return Err(PyTypeError::new_err("Matrix @ requires matrix operands"));}self.binary_mul(other,false)}
    fn __rmatmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if !is_matrix(other){return Err(PyTypeError::new_err("Matrix @ requires matrix operands"));}self.binary_mul(other,true)}
    fn __truediv__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if is_matrix(other){return Err(PyTypeError::new_err("Matrix division requires a scalar divisor"));}let c=symox::ops::pow(&to_expr(other)?,&Expr::minus_one());Ok(wrap_matrixexpr(MatrixExpr::mul(c,&[self.m.clone()],true).map_err(|e|error(other.py(),e))?))}
    fn __pow__(&self,exponent:&Bound<'_,PyAny>,modulo:Option<&Bound<'_,PyAny>>)->PyResult<Self>{if modulo.is_some_and(|m|!m.is_none()){return Err(PyTypeError::new_err("Modular matrix powers are unsupported"));}if is_matrix(exponent){return Err(PyTypeError::new_err("Matrix exponent must be scalar"));}Ok(wrap_matrixexpr(self.m.pow(to_expr(exponent)?,true).map_err(|e|error(exponent.py(),e))?))}
    fn __neg__(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(MatrixExpr::mul(Expr::minus_one(),&[self.m.clone()],true).map_err(|e|error(py,e))?))}
    fn __pos__(&self)->Self{self.clone()}
    fn transpose(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.transpose(true).map_err(|e|error(py,e))?))}
    fn adjoint(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.adjoint(true).map_err(|e|error(py,e))?))}
    fn conjugate(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.conjugate(true).map_err(|e|error(py,e))?))}
    fn inverse(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.inverse(true).map_err(|e|error(py,e))?))}
    fn inv(&self,py:Python<'_>)->PyResult<Self>{self.inverse(py)}
    #[getter] fn T(&self,py:Python<'_>)->PyResult<Self>{self.transpose(py)}
    #[getter] fn H(&self,py:Python<'_>)->PyResult<Self>{self.adjoint(py)}
    #[getter] fn C(&self,py:Python<'_>)->PyResult<Self>{self.conjugate(py)}
    fn as_explicit(&self,py:Python<'_>)->PyResult<crate::matrices::PyMatrix>{Ok(crate::matrices::wrapm(self.m.as_explicit().map_err(|e|error(py,e))?).into_pyobject(py)?.call_method0("as_immutable")?.extract()?)}
    fn as_mutable(&self,py:Python<'_>)->PyResult<crate::matrices::PyMatrix>{Ok(crate::matrices::wrapm(self.m.as_explicit().map_err(|e|error(py,e))?))}
    #[pyo3(signature=(deep=true,**kwargs))]
    fn doit(&self,py:Python<'_>,deep:bool,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Self>{if kwargs.is_some_and(|d|!d.is_empty()){return Err(PyNotImplementedError::new_err("MatrixExpr doit options are unsupported"));}Ok(wrap_matrixexpr(self.m.doit(deep).map_err(|e|error(py,e))?))}
    fn expand(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.expand().map_err(|e|error(py,e))?))}
    fn __getitem__(&self,key:&Bound<'_,PyAny>)->PyResult<PyExpr>{let pair=key.downcast::<PyTuple>().map_err(|_|PyTypeError::new_err("MatrixExpr indexing requires (row, col)"))?;if pair.len()!=2{return Err(PyTypeError::new_err("MatrixExpr indexing requires two indices"));}Ok(wrap(self.m.element(&to_expr(&pair.get_item(0)?)?,&to_expr(&pair.get_item(1)?)?).map_err(|e|error(key.py(),e))?))}
    fn xreplace(&self,rule:&Bound<'_,PyDict>)->PyResult<Self>{let mut matrices=Vec::new();let mut scalars=FxHashMap::default();for (a,b) in rule.iter(){if is_matrix(&a){matrices.push((to_matrixexpr(&a)?,to_matrixexpr(&b)?));}else{scalars.insert(to_expr(&a)?,to_expr(&b)?);}}Ok(wrap_matrixexpr(self.m.xreplace(&matrices,&scalars).map_err(|e|error(rule.py(),e))?))}
    #[pyo3(signature=(*args,**kwargs))]
    fn subs(&self,args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Self>{
        if kwargs.is_some_and(|d|!d.is_empty()){return Err(PyNotImplementedError::new_err("MatrixExpr subs options are unsupported"));}
        let mut pairs=Vec::new();if args.len()==2{pairs.push((args.get_item(0)?,args.get_item(1)?));}else if args.len()==1{let rule=args.get_item(0)?;if let Ok(d)=rule.downcast::<PyDict>(){pairs.extend(d.iter());}else{for pair in rule.try_iter()?{let pair=pair?;let a=pair.get_item(0)?;let b=pair.get_item(1)?;pairs.push((a,b));}}}else{return Err(PyTypeError::new_err("subs expects a mapping, pairs or (old,new)"));}
        let mut result=self.m.clone();for (a,b) in pairs{result=if is_matrix(&a){result.map_all(&[(to_matrixexpr(&a)?,to_matrixexpr(&b)?)],&mut |e|e.clone())}else{let old=to_expr(&a)?;let new=to_expr(&b)?;result.scalar_map(&mut |e|symox::subs::subs(e,&[(old.clone(),new.clone())]))}.map_err(|e|error(args.py(),e))?;}Ok(wrap_matrixexpr(result))
    }
    fn to_sympy<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{matrixexpr_to_sympy(py,&self.m)}
    fn _sympy_<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{self.to_sympy(py)}
    fn __reduce__<'py>(&self,py:Python<'py>)->PyResult<(Bound<'py,PyAny>,Bound<'py,PyTuple>)>{Ok((py.import("symox")?.getattr(self.m.kind_name())?,matrix_args(py,&self.m,false)?))}
}
impl PyMatrixExpr{
    fn binary_add(&self,other:&Bound<'_,PyAny>,subtract:bool,reverse:bool)->PyResult<Self>{let py=other.py();let mut a=self.m.clone();let mut b=to_matrixexpr(other)?;if reverse{std::mem::swap(&mut a,&mut b);}if subtract{b=MatrixExpr::mul(Expr::minus_one(),&[b],true).map_err(|e|error(py,e))?;}Ok(wrap_matrixexpr(MatrixExpr::add(&[a,b],true).map_err(|e|error(py,e))?))}
    fn binary_mul(&self,other:&Bound<'_,PyAny>,reverse:bool)->PyResult<Self>{let py=other.py();let result=if is_matrix(other){let peer=to_matrixexpr(other)?;let terms=if reverse{vec![peer,self.m.clone()]}else{vec![self.m.clone(),peer]};MatrixExpr::mul(Expr::one(),&terms,true)}else{MatrixExpr::mul(to_expr(other)?,&[self.m.clone()],true)};Ok(wrap_matrixexpr(result.map_err(|e|error(py,e))?))}
}
pub(crate) fn register(m:&Bound<'_,PyModule>)->PyResult<()>{m.add_class::<PyMatrixExpr>()?;m.add_function(wrap_pyfunction!(_matrixexpr_node,m)?)?;m.add_function(wrap_pyfunction!(_matrix_element,m)?)?;m.add_function(wrap_pyfunction!(_as_matrixexpr,m)?)?;Ok(())}
