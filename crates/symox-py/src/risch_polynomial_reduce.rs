//! Exact QQ polynomial_reduce after the real public zero-Poly prefix.
use pyo3::prelude::*;
use pyo3::types::{PyInt,PyList};
use symox::expr::{Expr,Kind};
use symox::poly::{dense,polyobj::Domain};
use crate::poly::PyPoly;

fn plain(p:&symox::poly::fraction::Polynomial)->bool {
    p.field.is_none()&&p.algebraic.is_none()&&p.expression.is_none()
        &&p.base.modulus.is_none()&&!p.base.has_domain_gens()
        &&matches!(p.base.domain,Domain::ZZ|Domain::QQ)&&p.base.gens.len()==1
        &&p.base.ngens==1&&p.base.rep.degree(0)<=256
        &&p.base.rep.iter().all(|(_,c)|c.n.bits()<=4096&&c.d.bits()<=4096)
}
fn symbol(value:&Bound<'_,PyAny>)->PyResult<Option<Expr>> {
    if !value.is_instance_of::<crate::PyExprNode>(){return Ok(None);}
    let e=crate::to_expr(value)?;
    let Kind::Sym(variable)=e.kind()else{return Ok(None)};
    if !crate::node_type(value.py(),&e).is_some_and(|c|value.get_type().is(&c))
        ||variable.host().is_some()||variable.scalar_identity().is_some()
        ||variable.is_wild()||symox::decompose::commutativity(&e)!=Some(true){return Ok(None);}
    Ok(Some(e))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn _risch_polynomial_reduce(p:&Bound<'_,PyAny>,d:&Bound<'_,PyAny>,t:&Bound<'_,PyAny>,
    variables:&Bound<'_,PyAny>,derivations:&Bound<'_,PyAny>,level:&Bound<'_,PyAny>,
    initial:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<Option<(PyObject,PyObject)>> {
    py.check_signals()?;
    // Restoring a Python __new__ dictionary entry need not restore tp_new.
    // Read the live slot and compare it with the immutable generated table;
    // a prefix callback cannot hide a constructor error by restoring aliases.
    use pyo3::impl_::pyclass::PyClassImpl;
    let stock_new=<PyPoly as PyClassImpl>::items_iter()
        .flat_map(|items|items.slots)
        .find(|slot|slot.slot==pyo3::ffi::Py_tp_new).map(|slot|slot.pfunc);
    let current_new=unsafe{pyo3::ffi::PyType_GetSlot(py.get_type::<PyPoly>().as_type_ptr(),pyo3::ffi::Py_tp_new)};
    if !stock_new.is_some_and(|slot|!slot.is_null()&&slot==current_new){return Ok(None);}
    if ![p,d,initial].into_iter().all(|v|v.is_exact_instance_of::<PyPoly>())
        ||!variables.is_exact_instance_of::<PyList>()||!derivations.is_exact_instance_of::<PyList>()
        ||!level.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let (variables,derivations)=(variables.downcast::<PyList>()?,derivations.downcast::<PyList>()?);
    let Ok(level)=level.extract::<isize>()else{return Ok(None)};
    if variables.len()!=derivations.len()||variables.is_empty()||level>=0
        ||level< -(variables.len() as isize){return Ok(None);}
    let active=(variables.len() as isize+level+1)as usize;
    if !t.is(&variables.get_item(active-1)?)||!d.is(&derivations.get_item(active-1)?){return Ok(None);}
    let Some(t)=symbol(t)? else{return Ok(None)};
    let (p_value,d_value,initial_value)=(p.downcast::<PyPoly>()?.get(),d.downcast::<PyPoly>()?.get(),initial.downcast::<PyPoly>()?.get());
    if !plain(&p_value.p)||!plain(&d_value.p)||!plain(&initial_value.p)
        ||p_value.p.base.gens!=[t.clone()]||d_value.p.base.gens!=[t.clone()]
        ||initial_value.p.base.gens!=[t.clone()]||initial_value.p.base.domain!=Domain::ZZ
        ||!initial_value.p.base.rep.is_zero(){return Ok(None);}
    let mut seen=rustc_hash::FxHashSet::default();
    for i in 0..active {
        py.check_signals()?;
        let Some(variable)=symbol(&variables.get_item(i)?)? else{return Ok(None)};
        if !seen.insert(variable.clone()){return Ok(None);}
        let derivation=derivations.get_item(i)?;
        if !derivation.is_exact_instance_of::<PyPoly>(){return Ok(None);}
        let derivation=derivation.downcast::<PyPoly>()?.get();
        if !plain(&derivation.p)||derivation.p.base.gens!=[variable]{return Ok(None);}
    }
    let p_coeff=dense::from_mpoly(&p_value.p.base.rep,0);
    let d_coeff=dense::from_mpoly(&d_value.p.base.rep,0);
    let Some(result)=symox::integrate::risch_polynomial_reduce::reduce(&p_coeff,&d_coeff,
        p_value.p.base.domain==Domain::QQ,&mut ||py.check_signals())? else{return Ok(None)};
    if result.steps==0{return Ok(Some((initial.clone().unbind(),p.clone().unbind())));}
    let mut quotient=initial_value.p.base.clone();quotient.rep=dense::to_mpoly(&result.quotient,0,1);
    quotient.domain=if result.quotient_qq{Domain::QQ}else{Domain::ZZ};
    let mut remainder=p_value.p.base.clone();remainder.rep=dense::to_mpoly(&result.remainder,0,1);
    remainder.domain=if result.remainder_qq{Domain::QQ}else{Domain::ZZ};
    let q=PyPoly{p:quotient.into()}.into_pyobject(py)?.into_any().unbind();
    let r=PyPoly{p:remainder.into()}.into_pyobject(py)?.into_any().unbind();
    py.check_signals()?;Ok(Some((q,r)))
}

/// Read-only admission hint. True only refuses native optimization; the
/// Python caller always executes its entire established body afterwards.
#[pyfunction]
fn _risch_decline_no_step(p:&Bound<'_,PyAny>,de:&Bound<'_,PyAny>,slot:&Bound<'_,PyAny>,
    py:Python<'_>)->PyResult<bool> {
    let stored=|value:&symox::poly::fraction::Polynomial| {
        value.field.is_none()&&value.algebraic.is_none()&&value.expression.is_none()
            &&value.base.modulus.is_none()&&value.base.ngens==1
            &&matches!(value.base.domain,Domain::ZZ|Domain::QQ)
    };
    if !p.is_exact_instance_of::<PyPoly>(){return Ok(false);}
    let p=p.downcast::<PyPoly>()?.get();
    if !stored(&p.p){return Ok(false);}
    if p.p.base.rep.is_zero(){return Ok(true);}
    // The exact immutable CPython member-descriptor type supplies its C
    // descr_get slot. No DE.__getattribute__, property, __getattr__, mutable
    // Python method lookup or binding of descriptor.__get__ is performed.
    if slot.get_type().as_type_ptr()!=std::ptr::addr_of_mut!(pyo3::ffi::PyMemberDescr_Type){return Ok(false);}
    let pointer=unsafe{pyo3::ffi::PyType_GetSlot(slot.get_type().as_type_ptr(),pyo3::ffi::Py_tp_descr_get)};
    let read=if pointer.is_null() {
        // Python3.9 limits GetSlot to heap types. Clear only that unsupported
        // lookup, then invoke this exact immutable C descriptor's own getter.
        if PyErr::occurred(py) {
            let error=PyErr::fetch(py);
            if !error.is_instance_of::<pyo3::exceptions::PySystemError>(py){return Err(error);}
        }
        slot.call_method1("__get__",(de,))
    }else{
        let get: pyo3::ffi::descrgetfunc=unsafe{std::mem::transmute(pointer)};
        let value=unsafe{get(slot.as_ptr(),de.as_ptr(),std::ptr::null_mut())};
        unsafe{Bound::from_owned_ptr_or_err(py,value)}
    };
    let value=match read {
        Ok(value)=>value,
        Err(error) if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            ||error.is_instance_of::<pyo3::exceptions::PyTypeError>(py)=>return Ok(false),
        Err(error)=>return Err(error),
    };
    if !value.is_exact_instance_of::<PyPoly>(){return Ok(false);}
    let d=value.downcast::<PyPoly>()?.get();
    if !stored(&d.p){return Ok(false);}
    // Degree in the first variable is stored in the leading monomial. The
    // comparison is O(1), independent of coefficients and without algebra.
    Ok(p.p.base.rep.degree(0)<d.p.base.rep.degree(0))
}

pub(crate) fn register(m:&Bound<'_,PyModule>)->PyResult<()> {
    m.add_function(wrap_pyfunction!(_risch_decline_no_step,m)?)?;
    m.add_function(wrap_pyfunction!(_risch_polynomial_reduce,m)?)
}
