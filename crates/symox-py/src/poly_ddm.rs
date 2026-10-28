//! Representation/alias checks and in-place publication for exact dense rows.
#[path="poly_function_metadata.rs"]
mod function_metadata;
#[path="poly_fflu.rs"]
mod fflu;
#[path="poly_hnf_modular.rs"]
mod hnf_modular;
#[path="poly_hnf.rs"]
mod hnf;
#[path="poly_lll.rs"]
mod lll;
#[path="poly_smith.rs"]
mod smith;
#[path="poly_lu.rs"]
mod lu;
#[path="poly_lu_solve.rs"]
mod lu_solve;
use pyo3::prelude::*;
use pyo3::types::{PyInt,PyList,PyTuple,PyDict};
use symox::poly::{Int,Rat,domains::{DomainElement,ElementValue,ExactDomain}};
use crate::poly::domains::PyDomainElement;

#[derive(Clone,Copy,PartialEq,Eq)]
pub(super) enum Scalar { Integer, ZZ, QQ }

fn sequence<'py>(value:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,mutable:bool)->Option<Vec<Bound<'py,PyAny>>> {
    if value.is_exact_instance_of::<PyList>()
        ||stock.is_some_and(|stock|value.get_type().is(stock)) {
        Some(value.downcast::<PyList>().ok()?.iter().collect())
    }else if !mutable&&value.is_exact_instance_of::<PyTuple>() {
        Some(value.downcast::<PyTuple>().ok()?.iter().collect())
    }else {None}
}

fn entries<'py>(rows:&[Bound<'py,PyAny>],mutable:bool)->Option<Vec<Vec<Bound<'py,PyAny>>>> {
    rows.iter().map(|row|sequence(row,None,mutable)).collect()
}

pub(super) fn scalar(value:&Bound<'_,PyAny>)->Option<Scalar> {
    if value.is_exact_instance_of::<PyInt>() {return Some(Scalar::Integer);}
    if !value.is_exact_instance_of::<PyDomainElement>() {return None;}
    let element=value.downcast::<PyDomainElement>().ok()?.get();
    match (&element.value.domain,&element.value.value) {
        (ExactDomain::ZZ,ElementValue::Rational(value)) if value.d.is_one()=>Some(Scalar::ZZ),
        (ExactDomain::QQ,ElementValue::Rational(_))=>Some(Scalar::QQ),
        _=>None,
    }
}

pub(super) fn integers(values:&[Vec<Bound<'_,PyAny>>],kind:Scalar,py:Python<'_>)->PyResult<Option<Vec<Vec<Int>>>> {
    let mut result=Vec::with_capacity(values.len());
    for row in values {
        py.check_signals()?;
        let mut output=Vec::with_capacity(row.len());
        for (j,value) in row.iter().enumerate() {
            if j%64==0 {py.check_signals()?;}
            if scalar(value)!=Some(kind) {return Ok(None);}
            output.push(if kind==Scalar::Integer {
                if let Ok(value)=value.extract::<i64>() {Int::from(value)}else{Int::from(value.extract::<num_bigint::BigInt>()?)}
            }else {
                let value=value.downcast::<PyDomainElement>()?.get();
                let ElementValue::Rational(value)=&value.value.value else{unreachable!()};value.n.clone()
            });
        }
        result.push(output);
    }
    Ok(Some(result))
}

pub(super) fn rationals(values:&[Vec<Bound<'_,PyAny>>],py:Python<'_>)->PyResult<Option<Vec<Vec<Rat>>>> {
    let mut result=Vec::with_capacity(values.len());
    for row in values {
        py.check_signals()?;
        let mut output=Vec::with_capacity(row.len());
        for (j,value) in row.iter().enumerate() {
            if j%64==0 {py.check_signals()?;}
            if scalar(value)!=Some(Scalar::QQ) {return Ok(None);}
            let value=value.downcast::<PyDomainElement>()?.get();
            let ElementValue::Rational(value)=&value.value.value else{unreachable!()};output.push(value.clone());
        }
        result.push(output);
    }
    Ok(Some(result))
}

pub(super) fn integer_object(py:Python<'_>,value:&Int)->PyResult<PyObject> {
    if let Some(value)=value.as_i64() {Ok(value.into_pyobject(py)?.into_any().unbind())}
    else {Ok(value.to_big().into_pyobject(py)?.into_any().unbind())}
}
pub(super) fn element_object(py:Python<'_>,value:Rat,domain:ExactDomain)->PyResult<PyObject> {
    Ok(PyDomainElement{value:DomainElement{domain,value:ElementValue::Rational(value)}}.into_pyobject(py)?.into_any().unbind())
}

/// Exact field element payloads beyond QQ whose native Python arithmetic is
/// `DomainElement::binary` on exactly these `Field` operations: GF(p) and the
/// algebraic fields QQ<a> and QQ_I. The domain is shared by every cell.
pub(super) trait FieldPayload: symox::poly::Field {
    fn payload(value:&ElementValue)->Option<&Self>;
    fn value(self)->ElementValue;
}
impl FieldPayload for symox::poly::f5b::Fp {
    fn payload(value:&ElementValue)->Option<&Self>{if let ElementValue::Finite(v)=value{Some(v)}else{None}}
    fn value(self)->ElementValue{ElementValue::Finite(self)}
}
impl FieldPayload for symox::poly::algebraic::AN {
    fn payload(value:&ElementValue)->Option<&Self>{if let ElementValue::Algebraic(v)=value{Some(v)}else{None}}
    fn value(self)->ElementValue{ElementValue::Algebraic(self)}
}

#[derive(Clone,Copy,PartialEq,Eq)]
pub(super) enum FieldKind { Finite, Algebraic }

/// The admitted non-QQ field of a stock element: GF(p) with p < 2**32 (the
/// native prime scope), an algebraic field QQ<a>, or QQ_I (only when the
/// caller's Python operations do not route through the Gaussian wrappers).
pub(super) fn field_domain(value:&Bound<'_,PyAny>,gaussian:bool)->Option<(ExactDomain,FieldKind)> {
    if !value.is_exact_instance_of::<PyDomainElement>() {return None;}
    let element=&value.downcast::<PyDomainElement>().ok()?.get().value;
    match (&element.domain,&element.value) {
        (ExactDomain::GF(p),ElementValue::Finite(v)) if v.characteristic()==*p&&symox::poly::f5b::PrimeScope::new(*p).is_ok()=>
            Some((element.domain.clone(),FieldKind::Finite)),
        (ExactDomain::Algebraic(field),ElementValue::Algebraic(v)) if v.field==*field=>
            Some((element.domain.clone(),FieldKind::Algebraic)),
        (ExactDomain::Gaussian{integral:false},ElementValue::Algebraic(v)) if gaussian&&symox::poly::gaussian::is_gaussian_field(&v.field)=>
            Some((element.domain.clone(),FieldKind::Algebraic)),
        _=>None,
    }
}

/// Every cell must be a stock element of exactly `domain` (no unification,
/// no conversion): a mixed matrix keeps the complete Python algorithm.
pub(super) fn field_values<C:FieldPayload>(values:&[Vec<Bound<'_,PyAny>>],domain:&ExactDomain,py:Python<'_>)->PyResult<Option<Vec<Vec<C>>>> {
    let mut result=Vec::with_capacity(values.len());
    for row in values {
        py.check_signals()?;
        let mut output=Vec::with_capacity(row.len());
        for (j,value) in row.iter().enumerate() {
            if j%64==0 {py.check_signals()?;}
            let Some(value)=field_value::<C>(value,domain) else{return Ok(None)};
            output.push(value);
        }
        result.push(output);
    }
    Ok(Some(result))
}

pub(super) fn field_value<C:FieldPayload>(value:&Bound<'_,PyAny>,domain:&ExactDomain)->Option<C> {
    if !value.is_exact_instance_of::<PyDomainElement>() {return None;}
    let element=&value.downcast::<PyDomainElement>().ok()?.get().value;
    if element.domain!=*domain {return None;}
    let payload=C::payload(&element.value)?;
    match (domain,&element.value) {
        (ExactDomain::GF(p),ElementValue::Finite(v)) if v.characteristic()!=*p=>return None,
        (ExactDomain::Algebraic(field),ElementValue::Algebraic(v)) if v.field!=*field=>return None,
        (ExactDomain::Gaussian{..},ElementValue::Algebraic(v)) if !symox::poly::gaussian::is_gaussian_field(&v.field)=>return None,
        _=>{}
    }
    Some(payload.clone())
}

pub(super) fn field_object<C:FieldPayload>(py:Python<'_>,value:C,domain:&ExactDomain)->PyResult<PyObject> {
    Ok(PyDomainElement{value:DomainElement{domain:domain.clone(),value:value.value()}}.into_pyobject(py)?.into_any().unbind())
}

/// Run `body` with the native coefficient scope of `domain` installed, as
/// `DomainElement::binary` does around the same field operations.
pub(super) fn with_field_scope<T>(domain:&ExactDomain,body:impl FnOnce()->PyResult<T>)->PyResult<Option<T>> {
    match domain {
        ExactDomain::GF(p)=>{
            let Ok(_scope)=symox::poly::f5b::PrimeScope::new(*p) else{return Ok(None)};
            body().map(Some)
        },
        ExactDomain::Algebraic(field)=>{let _scope=symox::poly::algebraic::AlgebraicScope::new(field.clone());body().map(Some)},
        ExactDomain::Gaussian{integral:false}=>{let _scope=symox::poly::algebraic::AlgebraicScope::new(symox::poly::gaussian::field());body().map(Some)},
        _=>Ok(None),
    }
}

#[pyfunction]
#[pyo3(signature=(a,b,c,stock=None))]
fn _ddm_imatmul<'py>(a:&Bound<'py,PyAny>,b:&Bound<'py,PyAny>,c:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,py:Python<'py>)->PyResult<bool> {
    py.check_signals()?;
    // The optional stock class is registered by DDM itself and admitted by
    // the public adapter only while its inherited list protocols are intact.
    let Some(arows)=sequence(a,stock,true) else{return Ok(false)};
    let Some(brows)=sequence(b,stock,false) else{return Ok(false)};
    let Some(crows)=sequence(c,stock,false) else{return Ok(false)};
    let mut destinations=rustc_hash::FxHashSet::default();
    for row in &arows {
        if !row.is_exact_instance_of::<PyList>()||!destinations.insert(row.as_ptr() as usize) {return Ok(false);}
    }
    if brows.iter().chain(&crows).any(|row|destinations.contains(&(row.as_ptr() as usize))) {return Ok(false);}
    let Some(a)=entries(&arows,true) else{return Ok(false)};
    let Some(b)=entries(&brows,false) else{return Ok(false)};
    let Some(c)=entries(&crows,false) else{return Ok(false)};
    let Some(first)=a.first() else{return Ok(false)};
    let n=first.len();let k=c.len();
    if n==0||k==0||a.len()!=b.len()||a.iter().any(|row|row.len()!=n)
        ||b.iter().any(|row|row.len()!=k)||c.iter().any(|row|row.len()!=n) {return Ok(false);}
    let Some(kind)=scalar(&a[0][0]) else{return Ok(false)};
    let output_rows=arows.iter().map(|row|row.downcast::<PyList>().unwrap()).collect::<Vec<_>>();
    if kind==Scalar::QQ {
        let Some(mut a)=rationals(&a,py)? else{return Ok(false)};
        let Some(b)=rationals(&b,py)? else{return Ok(false)};
        let Some(c)=rationals(&c,py)? else{return Ok(false)};
        symox::poly::ddm::imatmul(&mut a,&b,&c,&mut ||py.check_signals(),&mut |i,j,value| {
            output_rows[i].set_item(j,element_object(py,value.clone(),ExactDomain::QQ)?)
        })
    }else {
        let Some(mut a)=integers(&a,kind,py)? else{return Ok(false)};
        let Some(b)=integers(&b,kind,py)? else{return Ok(false)};
        let Some(c)=integers(&c,kind,py)? else{return Ok(false)};
        symox::poly::ddm::imatmul(&mut a,&b,&c,&mut ||py.check_signals(),&mut |i,j,value| {
            let value=if kind==Scalar::Integer {integer_object(py,value)?}else{element_object(py,Rat::from(value.clone()),ExactDomain::ZZ)?};
            output_rows[i].set_item(j,value)
        })
    }
}

#[pyfunction]
#[pyo3(signature=(a,partial,stock=None))]
fn _ddm_irref<'py>(a:&Bound<'py,PyAny>,partial:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,py:Python<'py>)->PyResult<Option<Vec<usize>>> {
    py.check_signals()?;
    if !partial.is_exact_instance_of::<pyo3::types::PyBool>() {return Ok(None);}
    let partial=partial.extract::<bool>()?;
    let Some(arows)=sequence(a,stock,true) else{return Ok(None)};
    let mut destinations=rustc_hash::FxHashSet::default();
    for row in &arows {
        if !row.is_exact_instance_of::<PyList>()||!destinations.insert(row.as_ptr() as usize) {return Ok(None);}
    }
    let Some(entries)=entries(&arows,true) else{return Ok(None)};
    let Some(first)=entries.first() else{return Ok(Some(Vec::new()))};
    let n=first.len();
    if entries.iter().any(|row|row.len()!=n) {return Ok(None);}
    if n==0 {return Ok(Some(Vec::new()));}
    if scalar(&entries[0][0])!=Some(Scalar::QQ) {
        // GF(p), QQ<a> and QQ_I: the Python body consumes only the native
        // pow/mul/sub/bool operations checked by the caller. Partial pivoting
        // is the inexact RR/CC path and stays in Python for these domains.
        if partial {return Ok(None);}
        let Some((domain,kind))=field_domain(&entries[0][0],true) else{return Ok(None)};
        let outer=a.downcast::<PyList>()?;
        let rows=arows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
        return match kind {
            FieldKind::Finite=>field_irref_dense::<symox::poly::f5b::Fp>(&entries,outer,rows,&domain,py),
            FieldKind::Algebraic=>field_irref_dense::<symox::poly::algebraic::AN>(&entries,outer,rows,&domain,py),
        };
    }
    let Some(mut values)=rationals(&entries,py)? else{return Ok(None)};
    let outer=a.downcast::<PyList>()?;
    let mut output_rows=arows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
    symox::poly::ddm::irref(&mut values,partial,&mut ||py.check_signals(),&mut |event| {
        match event {
            symox::poly::ddm::RrefEvent::Swap(i,j)=>{
                outer.set_item(i,&output_rows[j])?;
                outer.set_item(j,&output_rows[i])?;
                output_rows.swap(i,j);
                Ok(())
            },
            symox::poly::ddm::RrefEvent::Set(i,j,value)=>output_rows[i].set_item(j,element_object(py,value.clone(),ExactDomain::QQ)?),
        }
    })
}

fn field_irref_dense<C:FieldPayload>(entries:&[Vec<Bound<'_,PyAny>>],outer:&Bound<'_,PyList>,
    mut rows:Vec<Bound<'_,PyList>>,domain:&ExactDomain,py:Python<'_>)->PyResult<Option<Vec<usize>>> {
    let Some(mut values)=field_values::<C>(entries,domain,py)? else{return Ok(None)};
    with_field_scope(domain,||symox::poly::ddm::field_irref(&mut values,false,&mut ||py.check_signals(),&mut |event| {
        match event {
            symox::poly::ddm::FieldRrefEvent::Swap(i,j)=>{
                outer.set_item(i,&rows[j])?;outer.set_item(j,&rows[i])?;rows.swap(i,j);Ok(())
            },
            symox::poly::ddm::FieldRrefEvent::Set(i,j,value)=>rows[i].set_item(j,field_object(py,value.clone(),domain)?),
        }
    })).map(Option::flatten)
}

fn reduce_fraction_free<C:symox::poly::Coeff>(mut values:Vec<Vec<C>>,outer:&Bound<'_,PyList>,
    mut rows:Vec<Bound<'_,PyList>>,domain:ExactDomain,py:Python<'_>)->PyResult<Option<(PyObject,Vec<usize>)>> {
    let mut shared:Option<(usize,PyObject)>=None;
    let result=symox::poly::fraction_free::dense(&mut values,&mut ||py.check_signals(),&mut |event| {
        use symox::poly::fraction_free::DenseEvent;
        match event {
            DenseEvent::Swap(i,j)=>{outer.set_item(i,&rows[j])?;outer.set_item(j,&rows[i])?;rows.swap(i,j);Ok(())},
            DenseEvent::Set(i,j,value,group)=>{
                let value=if let Some(group)=group {
                    if !shared.as_ref().is_some_and(|(id,_)|*id==group){shared=Some((group,element_object(py,value.to_rat(),domain.clone())?));}
                    shared.as_ref().unwrap().1.clone_ref(py)
                }else{element_object(py,value.to_rat(),domain.clone())?};
                rows[i].set_item(j,value)
            },
        }
    }).map_err(super::fraction_free::failure::<C>)?;
    let Some((denominator,pivots))=result else{return Ok(None)};
    let den=if let Some(&j)=pivots.first(){rows[0].get_item(j)?.unbind()}
        else{element_object(py,denominator.to_rat(),domain)?};
    Ok(Some((den,pivots)))
}

#[pyfunction]
#[pyo3(signature=(a,k,stock=None))]
fn _ddm_irref_den<'py>(a:&Bound<'py,PyAny>,k:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,py:Python<'py>)->PyResult<Option<(PyObject,Vec<usize>)>> {
    py.check_signals()?;
    let typed=super::fraction_free::domain(k);
    let generic=if typed.is_none(){super::ff_domains::domain(k)}else{None};
    if typed.is_none()&&generic.is_none(){return Ok(None);}
    let Some(arows)=sequence(a,stock,true) else{return Ok(None)};
    let mut seen=rustc_hash::FxHashSet::default();
    for row in &arows{if !row.is_exact_instance_of::<PyList>()||!seen.insert(row.as_ptr() as usize){return Ok(None);}}
    let Some(values)=entries(&arows,true) else{return Ok(None)};
    let n=values.first().map(Vec::len).unwrap_or(0);
    if values.iter().any(|row|row.len()!=n){return Ok(None);}
    let outer=a.downcast::<PyList>()?;
    let rows=arows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect();
    if let Some(domain)=generic {return super::ff_domains::irref_den(outer,rows,&values,domain,py);}
    let domain=typed.unwrap();
    if domain==ExactDomain::ZZ {
        let Some(values)=integers(&values,Scalar::ZZ,py)? else{return Ok(None)};
        reduce_fraction_free(values,outer,rows,domain,py)
    }else {
        let Some(values)=rationals(&values,py)? else{return Ok(None)};
        reduce_fraction_free(values,outer,rows,domain,py)
    }
}

fn determinant<C: symox::poly::Coeff>(mut values: Vec<Vec<C>>, outer: &Bound<'_, PyList>,
    mut rows: Vec<Bound<'_, PyList>>, domain: ExactDomain, py: Python<'_>) -> PyResult<Option<PyObject>> {
    let result = symox::poly::ddm::idet(&mut values, &mut || py.check_signals(), &mut |event| {
        use symox::poly::ddm::DetEvent;
        match event {
            DetEvent::Swap(i,j) => {
                outer.set_item(i,&rows[j])?; outer.set_item(j,&rows[i])?; rows.swap(i,j); Ok(())
            },
            DetEvent::Set(i,j,value) => rows[i].set_item(j,element_object(py,value.to_rat(),domain.clone())?),
        }
    }).map_err(super::fraction_free::failure::<C>)?;
    result.map(|value| element_object(py,value.to_rat(),domain)).transpose()
}

#[pyfunction]
#[pyo3(signature=(a,k,stock=None))]
fn _ddm_idet<'py>(a: &Bound<'py,PyAny>, k: &Bound<'py,PyAny>,
    stock: Option<&Bound<'py,PyAny>>, py: Python<'py>) -> PyResult<Option<PyObject>> {
    py.check_signals()?;
    let typed = super::fraction_free::domain(k);
    let generic = if typed.is_none() { super::ff_domains::domain(k) } else { None };
    if typed.is_none() && generic.is_none() { return Ok(None); }
    let Some(arows) = sequence(a,stock,true) else { return Ok(None); };
    let mut seen = rustc_hash::FxHashSet::default();
    for row in &arows {
        if !row.is_exact_instance_of::<PyList>() || !seen.insert(row.as_ptr() as usize) { return Ok(None); }
    }
    let Some(values) = entries(&arows,true) else { return Ok(None); };
    if values.iter().any(|row|row.len()!=values.len()) { return Ok(None); }
    let outer = a.downcast::<PyList>()?;
    let rows = arows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect();
    if let Some(domain) = generic { return super::ff_domains::det(outer,rows,&values,domain,py); }
    let domain = typed.unwrap();
    if domain == ExactDomain::ZZ {
        let Some(values) = integers(&values,Scalar::ZZ,py)? else { return Ok(None); };
        determinant(values,outer,rows,domain,py)
    } else {
        let Some(values) = rationals(&values,py)? else { return Ok(None); };
        determinant(values,outer,rows,domain,py)
    }
}

#[pyfunction]
#[pyo3(signature=(a,k,stock=None))]
fn _ddm_berk<'py>(a: &Bound<'py,PyAny>, k: &Bound<'py,PyAny>,
    stock: Option<&Bound<'py,PyAny>>, py: Python<'py>) -> PyResult<Option<PyObject>> {
    use symox::matrix::{domain::RatMat,int::{Int as MatrixInt,Rat as MatrixRat}};
    fn matrix_integer(value: &Int) -> MatrixInt {
        value.as_i64().map(MatrixInt::from_i64).unwrap_or_else(|| MatrixInt::from_big(value.to_big()))
    }
    py.check_signals()?;
    let typed = super::fraction_free::domain(k);
    let generic = if typed.is_none() { super::ff_domains::domain(k) } else { None };
    if typed.is_none() && generic.is_none() { return Ok(None); }
    let Some(arows) = sequence(a,stock,true) else { return Ok(None); };
    if arows.iter().any(|row| !row.is_exact_instance_of::<PyList>()) { return Ok(None); }
    let Some(entries) = entries(&arows,true) else { return Ok(None); };
    let n = entries.len();
    if entries.iter().any(|row|row.len()!=n) { return Ok(None); }
    if let Some(domain) = generic { return super::ff_domains::berk(&entries,domain,py); }
    let domain = typed.unwrap();
    let values = if domain == ExactDomain::ZZ {
        let Some(values) = integers(&entries,Scalar::ZZ,py)? else { return Ok(None); };
        values.into_iter().map(|row|row.into_iter().map(Rat::from).collect()).collect()
    } else {
        let Some(values) = rationals(&entries,py)? else { return Ok(None); }; values
    };
    let mut data = Vec::new();
    for row in values {
        py.check_signals()?;
        for (j,value) in row.into_iter().enumerate() {
            if j%64==0 { py.check_signals()?; }
            // Canonical rational inputs stay exact; no expression or Float
            // conversion enters the existing rational matrix kernel.
            data.push(MatrixRat{num:matrix_integer(&value.n),den:matrix_integer(&value.d)});
        }
    }
    let coefficients = RatMat{rows:n,cols:n,data}.charpoly_checked(&mut ||py.check_signals())?;
    let output = PyList::empty(py);
    for value in coefficients {
        py.check_signals()?;
        let value = Rat::new(Int::from(value.num.to_big()),Int::from(value.den.to_big()));
        output.append(PyList::new(py,[element_object(py,value,domain.clone())?])?)?;
    }
    Ok(Some(output.into_any().unbind()))
}

#[pyfunction]
fn _inexact_linear_rref(rows: &Bound<'_,PyAny>, domain: &Bound<'_,PyAny>, width: usize,
    py: Python<'_>) -> PyResult<Option<(PyObject,Vec<usize>)>> {
    py.check_signals()?;
    if !domain.is_exact_instance_of::<super::domains::PyDomain>() || !rows.is_exact_instance_of::<PyList>() { return Ok(None); }
    let domain=domain.downcast::<super::domains::PyDomain>()?.get().inner.clone();
    let mut values=Vec::new();
    for row in rows.downcast::<PyList>()?.iter() {
        py.check_signals()?;
        if !row.is_exact_instance_of::<PyList>() { return Ok(None); }
        let row=row.downcast::<PyList>()?;
        if row.len()!=width { return Ok(None); }
        let mut entries=Vec::with_capacity(width);
        for (index,value) in row.iter().enumerate() {
            if index%64==0 { py.check_signals()?; }
            if !value.is_exact_instance_of::<PyDomainElement>() { return Ok(None); }
            entries.push(value.downcast::<PyDomainElement>()?.get().value.clone());
        }
        values.push(entries);
    }
    let result=symox::poly::ddm::inexact_linear_rref(&mut values,&domain,width,&mut ||py.check_signals())
        .map_err(|error| match error {
            symox::poly::ddm::InexactRrefError::Interrupted(error)=>error,
            symox::poly::ddm::InexactRrefError::Arithmetic(error)=>super::domains::element_error(error),
        })?;
    let Some(pivots)=result else { return Ok(None); };
    let out=PyList::empty(py);
    for row in values {
        py.check_signals()?;
        let entries=PyList::empty(py);
        for value in row { entries.append(PyDomainElement{value}.into_pyobject(py)?)?; }
        out.append(entries)?;
    }
    Ok(Some((out.into_any().unbind(),pivots)))
}

fn polymatrix_reduce<C:symox::poly::int::Field>(mut values:Vec<Vec<C>>,cols:usize,basis:bool,
    py:Python<'_>,expression:impl Fn(&C)->symox::Expr)->PyResult<Option<(PyObject,Vec<usize>)>> {
    let result=symox::poly::mpoly::checked_operation(||Ok((||->PyResult<_>{
        let Some(pivots)=symox::poly::ddm::field_irref(&mut values,true,&mut ||py.check_signals(),&mut |_|Ok(()))?
            else{return Ok(None)};
        let out=PyList::empty(py);
        if basis {
            let Some(vectors)=symox::poly::ddm::field_nullspace(&values,cols,&pivots,&mut ||py.check_signals())?
                else{return Ok(None)};
            for vector in vectors {
                py.check_signals()?;
                let row=PyList::empty(py);
                for value in vector {row.append(crate::wrap(expression(&value)))?;}
                out.append(row)?;
            }
        }else {
            for row in values {
                py.check_signals()?;
                for value in row {out.append(crate::wrap(expression(&value)))?;}
            }
        }
        Ok(Some((out.into_any().unbind(),pivots)))
    })()));
    match result {
        Ok(result)=>result,
        Err(symox::poly::convert::PolyError::Unsupported(message)) if message.contains("exponent storage limit")=>Ok(None),
        Err(error)=>Err(super::err(error)),
    }
}

/// Stock PolyMatrix ground entries: exact rationals or formal rational
/// functions over QQ with independent commutative symbol parameters. None
/// always requests the complete original route, never an empty nullspace.
#[pyfunction]
fn _polymatrix_ground(items:&Bound<'_,PyAny>,rows:&Bound<'_,PyAny>,cols:&Bound<'_,PyAny>,domain:&Bound<'_,PyAny>,
    basis:bool,stock:&Bound<'_,PyTuple>,py:Python<'_>)->PyResult<Option<(PyObject,Vec<usize>)>> {
    py.check_signals()?;
    if !items.is_exact_instance_of::<PyList>()||!domain.is_exact_instance_of::<super::domains::PyDomain>() {return Ok(None);}
    let domain=domain.downcast::<super::domains::PyDomain>()?.get().inner.clone();
    if !rows.is_exact_instance_of::<PyInt>()||!cols.is_exact_instance_of::<PyInt>() {return Ok(None);}
    let (Ok(rows),Ok(cols))=(rows.extract::<usize>(),cols.extract::<usize>()) else{return Ok(None)};
    let parameters=domain.flat_parameters();
    let rational=domain==ExactDomain::QQ;
    if !rational && !(domain.is_field()
        &&matches!(domain.composite_base(),Some(symox::poly::domains::CompositeBase::Rational{..}))
        &&!parameters.is_empty()
        &&parameters.iter().all(|p|p.is_symbol()&&!p.has_host()&&!p.has_wild_properties()
            &&symox::assumptions::fact(p,"commutative")==Some(true))) {return Ok(None);}
    let items=items.downcast::<PyList>()?;
    if rows.checked_mul(cols)!=Some(items.len()) {return Ok(None);}
    let mut entries=Vec::with_capacity(items.len());
    for (i,value) in items.iter().enumerate() {
        if i%64==0 {py.check_signals()?;}
        if !stock.iter().any(|t|value.get_type().is(&t)) {return Ok(None);}
        let Ok(value)=value.extract::<PyRef<'_,crate::PyExprNode>>() else{return Ok(None)};
        if value.e.has_wild_properties()||!symox::poly::expression_domain::native_coefficient_payload(&value.e) {return Ok(None);}
        let mut floating=false;
        symox::traverse::preorder(&value.e,&mut |e| {
            if e.as_number().is_some_and(symox::Number::is_float) {floating=true;}
        });
        if floating {return Ok(None);}
        entries.push(value.e.clone());
    }
    if rational {
        let Some(values)=entries.iter().map(|e|e.as_number().and_then(Rat::from_number)).collect::<Option<Vec<_>>>()
            else{return Ok(None)};
        let values=if cols==0 {vec![vec![];rows]}else{values.chunks(cols).map(<[_]>::to_vec).collect()};
        polymatrix_reduce(values,cols,basis,py,|r|symox::Expr::num(r.to_number()))
    }else {
        let _scope=symox::poly::ratfunc::Scope::new(parameters.len());
        let mut values=Vec::with_capacity(entries.len());
        for e in entries {
            py.check_signals()?;
            let Ok(DomainElement{value:ElementValue::Function(r),..})=domain.convert(&e) else{return Ok(None)};
            values.push(r);
        }
        let values=if cols==0 {vec![vec![];rows]}else{values.chunks(cols).map(<[_]>::to_vec).collect()};
        polymatrix_reduce(values,cols,basis,py,|r|DomainElement{domain:domain.clone(),value:ElementValue::Function(r.clone())}.as_expr())
    }
}

/// Read only actual dictionary storage. Subclasses are admitted by the
/// stock SDM adapter only after its inherited protocols have been checked.
fn sparse_entries<C:symox::poly::Coeff>(matrix:&Bound<'_,PyDict>,n:usize,py:Python<'_>,
    mut parse:impl FnMut(&Bound<'_,PyAny>)->PyResult<Option<C>>)
    ->PyResult<Option<symox::poly::sdm::SparseRows<C>>> {
    let mut rows=std::collections::BTreeMap::new();
    // No live PyDict iterator may cross a Python signal handler: that
    // handler could mutate the original dictionary. Copy references first.
    let items:Vec<_>=matrix.iter().collect();
    for (step,(i,row)) in items.into_iter().enumerate() {
        if step%64==0{py.check_signals()?;}
        if !i.is_exact_instance_of::<PyInt>() || !row.is_exact_instance_of::<PyDict>() {return Ok(None);}
        let Ok(i)=i.extract::<usize>() else{return Ok(None)};
        if i>=n{return Ok(None);}
        let mut entries=std::collections::BTreeMap::new();
        let entries_snapshot:Vec<_>=row.downcast::<PyDict>()?.iter().collect();
        for (step,(j,value)) in entries_snapshot.into_iter().enumerate() {
            if step%64==0{py.check_signals()?;}
            if !j.is_exact_instance_of::<PyInt>(){return Ok(None);}
            let Ok(j)=j.extract::<usize>() else{return Ok(None)};
            if j>=n{return Ok(None);}
            let Some(value)=parse(&value)? else{return Ok(None)};
            if !value.is_zero(){entries.insert(j,value);}
        }
        if !entries.is_empty(){rows.insert(i,entries);}
    }
    Ok(Some(rows))
}

#[pyfunction]
#[pyo3(signature=(matrix,n,k,stock=None))]
fn _sdm_berk(matrix:&Bound<'_,PyAny>,n:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>,
    stock:Option<&Bound<'_,PyAny>>,py:Python<'_>)->PyResult<Option<PyObject>> {
    py.check_signals()?;
    let typed=super::fraction_free::domain(k);
    let generic=if typed.is_none(){super::ff_domains::domain(k)}else{None};
    if typed.is_none()&&generic.is_none(){return Ok(None);}
    if !n.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let Ok(n)=n.extract::<usize>() else{return Ok(None)};
    if !matrix.is_exact_instance_of::<PyDict>() && !stock.is_some_and(|stock|matrix.get_type().is(stock)){return Ok(None);}
    let Ok(matrix)=matrix.downcast::<PyDict>() else{return Ok(None)};
    if let Some(domain)=generic {
        let Some(scope)=super::ff_domains::scope(&domain) else{return Ok(None)};
        if n==0 {
            let out=PyDict::new(py);
            out.set_item(0,super::ff_domains::object(py,&<symox::poly::domain_coeff::DomainCoeff as symox::poly::Coeff>::one())?)?;
            return Ok(Some(out.into_any().unbind()));
        }
        let Some(rows)=sparse_entries(matrix,n,py,|value|Ok(super::ff_domains::element(value,&domain)))? else{return Ok(None)};
        return super::ff_domains::sparse_berk(&rows,n,&scope,py);
    }
    let domain=typed.unwrap();
    let out=PyDict::new(py);
    if n==0{
        out.set_item(0,element_object(py,Rat::from(Int::from(1)),domain)?)?;
        return Ok(Some(out.into_any().unbind()));
    }
    if domain==ExactDomain::ZZ {
        let Some(rows)=sparse_entries(matrix,n,py,|value|{
            if scalar(value)!=Some(Scalar::ZZ){return Ok(None);}
            let value=value.downcast::<PyDomainElement>()?.get();
            let ElementValue::Rational(value)=&value.value.value else{unreachable!()};
            Ok(Some(value.n.clone()))
        })? else{return Ok(None)};
        let Some(coefficients)=symox::poly::sdm::berkowitz_checked(&rows,n,&mut ||py.check_signals())? else{return Ok(None)};
        for (i,value) in coefficients{py.check_signals()?;out.set_item(i,element_object(py,Rat::from(value),domain.clone())?)?;}
    }else{
        let Some(rows)=sparse_entries(matrix,n,py,|value|{
            if scalar(value)!=Some(Scalar::QQ){return Ok(None);}
            let value=value.downcast::<PyDomainElement>()?.get();
            let ElementValue::Rational(value)=&value.value.value else{unreachable!()};
            Ok(Some(value.clone()))
        })? else{return Ok(None)};
        let Some(coefficients)=symox::poly::sdm::berkowitz_checked(&rows,n,&mut ||py.check_signals())? else{return Ok(None)};
        for (i,value) in coefficients{py.check_signals()?;out.set_item(i,element_object(py,value,domain.clone())?)?;}
    }
    Ok(Some(out.into_any().unbind()))
}

fn sparse_convolution<C:symox::poly::Coeff>(q:&Bound<'_,PyDict>,t:&Bound<'_,PyList>,degree:usize,
    py:Python<'_>,mut parse:impl FnMut(&Bound<'_,PyAny>)->PyResult<Option<C>>,
    publish:impl Fn(&C)->PyResult<PyObject>,failed:impl Fn()->bool)
    ->PyResult<Option<PyObject>> {
    // The legacy caller takes min/max(q): an empty dictionary must retain
    // that explicit error instead of becoming a successfully empty result.
    if q.is_empty(){return Ok(None);}
    let mut q_native=std::collections::BTreeMap::new();
    let items:Vec<_>=q.iter().collect();
    for (step,(i,value)) in items.into_iter().enumerate(){
        if step%64==0{py.check_signals()?;}
        if !i.is_exact_instance_of::<PyInt>(){return Ok(None);}
        let Ok(i)=i.extract::<usize>() else{return Ok(None)};
        if i>degree{return Ok(None);}
        let Some(value)=parse(&value)? else{return Ok(None)};
        if !value.is_zero(){q_native.insert(i,value);}
    }
    let mut t_native=std::collections::BTreeMap::new();
    for (i,value) in t.iter().enumerate(){
        if i%64==0{py.check_signals()?;}
        let Some(value)=parse(&value)? else{return Ok(None)};
        if !value.is_zero(){t_native.insert(i,value);}
    }
    let result=symox::poly::sdm::toeplitz_product_checked(&q_native,&t_native,degree,&mut ||py.check_signals())?;
    if failed(){return Ok(None);}
    let out=PyDict::new(py);
    for (i,value) in result{py.check_signals()?;out.set_item(i,publish(&value)?)?;}
    Ok(Some(out.into_any().unbind()))
}

/// Finish the actual legacy prefix/recursive results, without replaying any
/// arithmetic or conversion that may have invoked a user callback already.
#[pyfunction]
fn _sdm_berk_combine(q:&Bound<'_,PyAny>,t:&Bound<'_,PyAny>,degree:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>,
    py:Python<'_>)->PyResult<Option<PyObject>> {
    py.check_signals()?;
    let typed=super::fraction_free::domain(k);
    let generic=if typed.is_none(){super::ff_domains::domain(k)}else{None};
    if typed.is_none()&&generic.is_none(){return Ok(None);}
    if !degree.is_exact_instance_of::<PyInt>() || !q.is_exact_instance_of::<PyDict>()
        || !t.is_exact_instance_of::<PyList>(){return Ok(None);}
    let Ok(degree)=degree.extract::<usize>() else{return Ok(None)};
    let q=q.downcast::<PyDict>()?;let t=t.downcast::<PyList>()?;
    if let Some(domain)=generic {
        let Some(scope)=super::ff_domains::scope(&domain) else{return Ok(None)};
        return sparse_convolution(q,t,degree,py,|value|Ok(super::ff_domains::element(value,&domain)),
            |value|super::ff_domains::object(py,value),||scope.failed());
    }
    let domain=typed.unwrap();
    let publish=|value:&Rat|element_object(py,value.clone(),domain.clone());
    if domain==ExactDomain::ZZ{
        sparse_convolution(q,t,degree,py,|value|{
            if scalar(value)!=Some(Scalar::ZZ){return Ok(None);}
            let value=value.downcast::<PyDomainElement>()?.get();
            let ElementValue::Rational(value)=&value.value.value else{unreachable!()};
            Ok(Some(value.n.clone()))
        },|value:&Int|publish(&Rat::from(value.clone())),||false)
    }else{
        sparse_convolution(q,t,degree,py,|value|{
            if scalar(value)!=Some(Scalar::QQ){return Ok(None);}
            let value=value.downcast::<PyDomainElement>()?.get();
            let ElementValue::Rational(value)=&value.value.value else{unreachable!()};
            Ok(Some(value.clone()))
        },publish,||false)
    }
}

pub fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    function_metadata::register(module)?;
    fflu::register(module)?;
    hnf::register(module)?;
    hnf_modular::register(module)?;
    lll::register(module)?;
    smith::register(module)?;
    lu::register(module)?;
    lu_solve::register(module)?;
    module.add_function(wrap_pyfunction!(_ddm_imatmul,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_irref,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_irref_den,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_idet,module)?)?;
    module.add_function(wrap_pyfunction!(_inexact_linear_rref,module)?)?;
    module.add_function(wrap_pyfunction!(_polymatrix_ground,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_berk,module)?)?;
    module.add_function(wrap_pyfunction!(_sdm_berk,module)?)?;
    module.add_function(wrap_pyfunction!(_sdm_berk_combine,module)?)?;
    Ok(())
}
