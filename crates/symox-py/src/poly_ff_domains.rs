//! Fraction-free kernels over the exact native domains beyond ZZ/QQ.
//!
//! GF(p), ZZ_I/QQ_I, QQ<a>, ZZ[x]/QQ[x]/QQ(x) (multivariate included) and
//! the composite rings over Gaussian/algebraic grounds use `DomainCoeff`,
//! whose arithmetic is the very element arithmetic that the Python bodies
//! reach through `DomainElement` and `Domain.exquo`. Inputs must be stock
//! elements of exactly the requested domain. Any arithmetic error or non-exact
//! quotient declines before the first publication (writes are buffered), so
//! the complete Python body then runs unchanged and raises its own error.
use pyo3::prelude::*;
use pyo3::types::{PyDict,PyList};
use symox::poly::{Coeff,domain_coeff::{self,DomainCoeff,DomainScope},domains::ExactDomain,fraction_free};
use crate::poly::domains::PyDomainElement;

pub(crate) enum Stop { Py(PyErr), Decline }

/// The requested domain when it is an exact `Domain` handled here.
pub(crate) fn domain(k:&Bound<'_,PyAny>)->Option<ExactDomain> {
    if !k.is_exact_instance_of::<crate::poly::domains::PyDomain>(){return None;}
    let k=k.downcast::<crate::poly::domains::PyDomain>().ok()?.get();
    domain_coeff::admitted(&k.inner).then(||k.inner.clone())
}

pub(crate) fn scope(domain:&ExactDomain)->Option<DomainScope> {DomainScope::new(domain).ok()}

pub(crate) fn element(value:&Bound<'_,PyAny>,domain:&ExactDomain)->Option<DomainCoeff> {
    if !value.is_exact_instance_of::<PyDomainElement>(){return None;}
    let element=&value.downcast::<PyDomainElement>().ok()?.get().value;
    domain_coeff::native_member(domain,element).then(||DomainCoeff(element.clone()))
}

pub(crate) fn values(cells:&[Vec<Bound<'_,PyAny>>],domain:&ExactDomain,py:Python<'_>)->PyResult<Option<Vec<Vec<DomainCoeff>>>> {
    let mut out=Vec::with_capacity(cells.len());
    for row in cells {
        py.check_signals()?;
        let mut values=Vec::with_capacity(row.len());
        for (j,value) in row.iter().enumerate() {
            if j%64==0{py.check_signals()?;}
            let Some(value)=element(value,domain) else{return Ok(None)};
            values.push(value);
        }
        out.push(values);
    }
    Ok(Some(out))
}

pub(crate) fn object(py:Python<'_>,value:&DomainCoeff)->PyResult<PyObject> {
    Ok(PyDomainElement{value:value.0.clone()}.into_pyobject(py)?.into_any().unbind())
}

/// Signal checks, plus an early stop once an element operation has failed.
pub(crate) fn checker<'a>(py:Python<'a>,scope:&'a DomainScope)->impl FnMut()->Result<(),Stop>+'a {
    move ||{py.check_signals().map_err(Stop::Py)?;if scope.failed(){Err(Stop::Decline)}else{Ok(())}}
}

/// Kernel outcome: Ok(Some) to publish, Ok(None) to decline, Err to raise.
pub(crate) fn settle<T>(result:Result<Option<T>,fraction_free::Failure<Stop,DomainCoeff>>,scope:&DomainScope)->PyResult<Option<T>> {
    match result {
        Err(fraction_free::Failure::Interrupted(Stop::Py(error)))=>Err(error),
        Err(_)=>Ok(None),
        Ok(value)=>Ok(if scope.failed(){None}else{value}),
    }
}

enum Write { Swap(usize,usize), Set(usize,usize,DomainCoeff,Option<usize>), Retain(usize,usize) }

/// Replay the buffered writes in their original order on the actual rows.
fn replay(py:Python<'_>,outer:&Bound<'_,PyList>,rows:&mut [Bound<'_,PyList>],writes:Vec<Write>,
    perm:Option<&Bound<'_,PyList>>)->PyResult<()> {
    let mut shared:Option<(usize,PyObject)>=None;
    for (step,write) in writes.into_iter().enumerate() {
        if step%64==0{py.check_signals()?;}
        match write {
            Write::Swap(i,j)=>{
                outer.set_item(i,&rows[j])?;outer.set_item(j,&rows[i])?;rows.swap(i,j);
                if let Some(perm)=perm {
                    let left=perm.get_item(i)?;let right=perm.get_item(j)?;perm.set_item(i,right)?;perm.set_item(j,left)?;
                }
            },
            Write::Set(i,j,value,group)=>{
                let value=if let Some(group)=group {
                    if !shared.as_ref().is_some_and(|(id,_)|*id==group){shared=Some((group,object(py,&value)?));}
                    shared.as_ref().unwrap().1.clone_ref(py)
                }else{object(py,&value)?};
                rows[i].set_item(j,value)?;
            },
            Write::Retain(i,j)=>{let value=rows[i].get_item(j)?;rows[i].set_item(j,value)?;},
        }
    }
    Ok(())
}

/// `ddm_irref_den`: dense fraction-free Gauss-Jordan.
pub(crate) fn irref_den(outer:&Bound<'_,PyList>,mut rows:Vec<Bound<'_,PyList>>,cells:&[Vec<Bound<'_,PyAny>>],
    domain:ExactDomain,py:Python<'_>)->PyResult<Option<(PyObject,Vec<usize>)>> {
    let Some(scope)=scope(&domain) else{return Ok(None)};
    let Some(mut values)=values(cells,&domain,py)? else{return Ok(None)};
    let mut writes=Vec::new();
    let result=fraction_free::dense(&mut values,&mut checker(py,&scope),&mut |event|{
        use fraction_free::DenseEvent;
        writes.push(match event {DenseEvent::Swap(i,j)=>Write::Swap(i,j),DenseEvent::Set(i,j,value,group)=>Write::Set(i,j,value.clone(),group)});
        Ok(())
    });
    let Some((denominator,pivots))=settle(result,&scope)? else{return Ok(None)};
    replay(py,outer,&mut rows,writes,None)?;
    let den=if let Some(&j)=pivots.first(){rows[0].get_item(j)?.unbind()}else{object(py,&denominator)?};
    Ok(Some((den,pivots)))
}

/// `ddm_idet`: in-place Bareiss determinant.
pub(crate) fn det(outer:&Bound<'_,PyList>,mut rows:Vec<Bound<'_,PyList>>,cells:&[Vec<Bound<'_,PyAny>>],
    domain:ExactDomain,py:Python<'_>)->PyResult<Option<PyObject>> {
    let Some(scope)=scope(&domain) else{return Ok(None)};
    let Some(mut values)=values(cells,&domain,py)? else{return Ok(None)};
    let mut writes=Vec::new();
    let result=symox::poly::ddm::idet(&mut values,&mut checker(py,&scope),&mut |event|{
        use symox::poly::ddm::DetEvent;
        writes.push(match event {DetEvent::Swap(i,j)=>Write::Swap(i,j),DetEvent::Set(i,j,value)=>Write::Set(i,j,value.clone(),None)});
        Ok(())
    });
    let Some(value)=settle(result,&scope)? else{return Ok(None)};
    replay(py,outer,&mut rows,writes,None)?;
    Ok(Some(object(py,&value)?))
}

/// `ddm_berk`: the characteristic polynomial is unique, and the division-free
/// sparse Berkowitz kernel computes it over any commutative ring.
pub(crate) fn berk(cells:&[Vec<Bound<'_,PyAny>>],domain:ExactDomain,py:Python<'_>)->PyResult<Option<PyObject>> {
    let Some(scope)=scope(&domain) else{return Ok(None)};
    let Some(values)=values(cells,&domain,py)? else{return Ok(None)};
    let n=values.len();
    let mut rows=symox::poly::sdm::SparseRows::new();
    for (i,row) in values.into_iter().enumerate() {
        let entries:symox::poly::sdm::SparseVector<DomainCoeff>=row.into_iter().enumerate().filter(|(_,v)|!v.is_zero()).collect();
        if !entries.is_empty(){rows.insert(i,entries);}
    }
    let result=symox::poly::sdm::berkowitz_checked(&rows,n,&mut checker(py,&scope));
    let coefficients=match result {Ok(value)=>value,Err(Stop::Py(error))=>return Err(error),Err(Stop::Decline)=>return Ok(None)};
    let Some(coefficients)=coefficients else{return Ok(None)};
    if scope.failed(){return Ok(None);}
    let zero=DomainCoeff::zero();
    let output=PyList::empty(py);
    for i in 0..=n {
        if i%64==0{py.check_signals()?;}
        output.append(PyList::new(py,[object(py,coefficients.get(&i).unwrap_or(&zero))?])?)?;
    }
    Ok(Some(output.into_any().unbind()))
}

/// `sdm_berk` over sparse rows already parsed by the caller's reader.
pub(crate) fn sparse_berk(rows:&symox::poly::sdm::SparseRows<DomainCoeff>,n:usize,scope:&DomainScope,py:Python<'_>)->PyResult<Option<PyObject>> {
    let result=symox::poly::sdm::berkowitz_checked(rows,n,&mut checker(py,scope));
    let coefficients=match result {Ok(value)=>value,Err(Stop::Py(error))=>return Err(error),Err(Stop::Decline)=>return Ok(None)};
    let Some(coefficients)=coefficients else{return Ok(None)};
    if scope.failed(){return Ok(None);}
    let out=PyDict::new(py);
    for (i,value) in coefficients{py.check_signals()?;out.set_item(i,object(py,&value)?)?;}
    Ok(Some(out.into_any().unbind()))
}

/// FFLU phase 1 (packed Bareiss elimination with row permutation).
pub(crate) fn fflu_eliminate(outer:&Bound<'_,PyList>,mut rows:Vec<Bound<'_,PyList>>,cells:&[Vec<Bound<'_,PyAny>>],
    perm:&Bound<'_,PyList>,n:usize,domain:ExactDomain,py:Python<'_>)->PyResult<bool> {
    let Some(scope)=scope(&domain) else{return Ok(false)};
    let Some(mut values)=values(cells,&domain,py)? else{return Ok(false)};
    let mut writes=Vec::new();
    let result=symox::poly::fflu::eliminate(&mut values,n,&mut checker(py,&scope),&mut |event|{
        use symox::poly::fflu::Eliminate;
        writes.push(match event {Eliminate::Swap(i,j)=>Write::Swap(i,j),Eliminate::Set(i,j,value)=>Write::Set(i,j,value.clone(),None),Eliminate::Retain(i,j)=>Write::Retain(i,j)});
        Ok(())
    });
    if settle(result,&scope)?.is_none(){return Ok(false);}
    replay(py,outer,&mut rows,writes,Some(perm))?;
    Ok(true)
}

/// FFLU phase 4 diagonal products D[i][i] = L[i-1][i-1]*L[i][i].
pub(crate) fn fflu_diagonal(d_rows:&[Bound<'_,PyList>],l_cells:&[Vec<Bound<'_,PyAny>>],d_cells:&[Vec<Bound<'_,PyAny>>],
    domain:ExactDomain,py:Python<'_>)->PyResult<bool> {
    let Some(scope)=scope(&domain) else{return Ok(false)};
    let (Some(lv),Some(mut dv))=(values(l_cells,&domain,py)?,values(d_cells,&domain,py)?) else{return Ok(false)};
    let mut writes=Vec::new();
    let result=symox::poly::fflu::diagonal(&lv,&mut dv,&mut checker(py,&scope),&mut |i,value|{writes.push((i,value.clone()));Ok(())});
    match result {Ok(true)=>{},Ok(false)|Err(Stop::Decline)=>return Ok(false),Err(Stop::Py(error))=>return Err(error)}
    if scope.failed(){return Ok(false);}
    for (i,value) in writes{py.check_signals()?;d_rows[i].set_item(i,object(py,&value)?)?;}
    Ok(true)
}
