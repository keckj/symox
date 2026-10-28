//! QQ splitfactor publication after the actual kinv/One prefix.
use pyo3::prelude::*;
use pyo3::types::{PyInt,PyList};
use symox::expr::{Expr,Kind};
use symox::poly::{dense,polyobj::Domain};
use crate::poly::PyPoly;
fn plain(p:&symox::poly::fraction::Polynomial)->bool{
    p.field.is_none()&&p.algebraic.is_none()&&p.expression.is_none()
        &&p.base.modulus.is_none()&&!p.base.has_domain_gens()&&p.base.ngens==1
        &&p.base.gens.len()==1&&matches!(p.base.domain,Domain::ZZ|Domain::QQ)
        &&p.base.rep.degree(0)<=64&&p.base.rep.iter().all(|(_,v)|v.n.bits()<=4096&&v.d.bits()<=4096)
}
fn symbol(value:&Bound<'_,PyAny>)->PyResult<Option<Expr>>{
    if !value.is_instance_of::<crate::PyExprNode>(){return Ok(None);}
    let e=crate::to_expr(value)?;let Kind::Sym(s)=e.kind()else{return Ok(None)};
    if s.host().is_some()||s.scalar_identity().is_some()||s.is_wild()
        ||symox::decompose::commutativity(&e)!=Some(true)
        ||!crate::node_type(value.py(),&e).is_some_and(|c|value.get_type().is(&c)){return Ok(None);}
    Ok(Some(e))
}
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn _risch_splitfactor(p:&Bound<'_,PyAny>,d:&Bound<'_,PyAny>,t:&Bound<'_,PyAny>,
    variables:&Bound<'_,PyAny>,derivations:&Bound<'_,PyAny>,level:&Bound<'_,PyAny>,
    one:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<Option<(PyObject,PyObject)>>{
    use pyo3::impl_::pyclass::PyClassImpl;
    py.check_signals()?;
    let expected=<PyPoly as PyClassImpl>::items_iter().flat_map(|i|i.slots)
        .find(|slot|slot.slot==pyo3::ffi::Py_tp_new).map(|slot|slot.pfunc);
    let current=unsafe{pyo3::ffi::PyType_GetSlot(py.get_type::<PyPoly>().as_type_ptr(),pyo3::ffi::Py_tp_new)};
    if !expected.is_some_and(|p|!p.is_null()&&p==current)
        ||![p,d,one].into_iter().all(|v|v.is_exact_instance_of::<PyPoly>())
        ||!variables.is_exact_instance_of::<PyList>()||!derivations.is_exact_instance_of::<PyList>()
        ||!level.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let (variables,derivations)=(variables.downcast::<PyList>()?,derivations.downcast::<PyList>()?);
    let Ok(level)=level.extract::<isize>()else{return Ok(None)};
    if variables.is_empty()||variables.len()!=derivations.len()||level>=0||level< -(variables.len()as isize){return Ok(None);}
    let active=(variables.len()as isize+level+1)as usize;
    if !t.is(&variables.get_item(active-1)?)||!d.is(&derivations.get_item(active-1)?){return Ok(None);}
    let Some(t)=symbol(t)?else{return Ok(None)};
    let(a,b,unit)=(p.downcast::<PyPoly>()?.get(),d.downcast::<PyPoly>()?.get(),one.downcast::<PyPoly>()?.get());
    if !plain(&a.p)||a.p.base.domain!=Domain::QQ||a.p.base.gens!=[t.clone()]
        ||!plain(&b.p)||b.p.base.gens!=[t.clone()]||b.p.base.rep.is_zero()
        ||!plain(&unit.p)||unit.p.base.domain!=Domain::QQ||unit.p.base.gens!=[t.clone()]
        ||!unit.p.base.rep.is_one(){return Ok(None);}
    let mut seen=rustc_hash::FxHashSet::default();
    for i in 0..active{
        py.check_signals()?;
        let Some(v)=symbol(&variables.get_item(i)?)?else{return Ok(None)};
        if !seen.insert(v.clone()){return Ok(None);}
        let derivative=derivations.get_item(i)?;
        if !derivative.is_exact_instance_of::<PyPoly>(){return Ok(None);}
        let derivative=derivative.downcast::<PyPoly>()?.get();
        if !plain(&derivative.p)||derivative.p.base.gens!=[v]{return Ok(None);}
    }
    let input=dense::from_mpoly(&a.p.base.rep,0);let tau=dense::from_mpoly(&b.p.base.rep,0);
    let Some(result)=symox::integrate::risch_splitfactor::split(&input,&tau,active>1,&mut||py.check_signals())?else{return Ok(None)};
    let normal=if result.reuse_input{p.clone().unbind()}else{
        let mut out=a.p.base.clone();out.rep=dense::to_mpoly(&result.normal,0,1);out.domain=Domain::QQ;
        PyPoly{p:out.into()}.into_pyobject(py)?.into_any().unbind()
    };
    let special=if result.reuse_one{one.clone().unbind()}else{
        let mut out=a.p.base.clone();out.rep=dense::to_mpoly(&result.special,0,1);
        out.domain=if result.special_integer{Domain::ZZ}else{Domain::QQ};
        PyPoly{p:out.into()}.into_pyobject(py)?.into_any().unbind()
    };
    py.check_signals()?;Ok(Some((normal,special)))
}
pub(crate)fn register(m:&Bound<'_,PyModule>)->PyResult<()>{m.add_function(wrap_pyfunction!(_risch_splitfactor,m)?)}
