//! Probability model metadata adapters; mathematical kernels live in Rust core.
use pyo3::prelude::*;
use symox::stats::{self,Law};
use crate::{to_expr,wrap,PyExpr};

fn order(e: &Bound<'_,PyAny>) -> PyResult<Option<usize>> {
    let e=to_expr(e)?;
    Ok(e.as_number().filter(|n|n.is_integer() && !n.is_negative()).and_then(|n|n.to_i64()).and_then(|n|usize::try_from(n).ok()))
}

#[pyfunction]
fn _stats_binomial_moment(n:&Bound<'_,PyAny>,p:&Bound<'_,PyAny>,success:&Bound<'_,PyAny>,failure:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>) -> PyResult<Option<PyExpr>> {
    let Some(k)=order(k)? else {return Ok(None)};
    Ok(Some(wrap(stats::binomial_moment(&to_expr(n)?,&to_expr(p)?,&to_expr(success)?,&to_expr(failure)?,k,&mut ||n.py().check_signals())?)))
}
#[pyfunction]
fn _stats_normal_moment(n:&Bound<'_,PyAny>,mean:&Bound<'_,PyAny>,sigma:&Bound<'_,PyAny>) -> PyResult<Option<PyExpr>> {
    let Some(k)=order(n)? else {return Ok(None)};
    Ok(Some(wrap(stats::normal_moment(&to_expr(mean)?,&to_expr(sigma)?,k,&mut ||n.py().check_signals())?)))
}

#[pyfunction]
fn _stats_poisson_moment(rate:&Bound<'_,PyAny>,n:&Bound<'_,PyAny>)->PyResult<Option<PyExpr>> {
    let Some(order)=order(n)? else{return Ok(None)};
    // An optimization budget only: the original distribution formula remains
    // available at every declined order, including symbolic/Float orders.
    if order>4096 {return Ok(None);}
    let value=to_expr(rate)?;
    if value.has_host()||value.has_matrix() {return Ok(None);}
    Ok(Some(wrap(stats::poisson_moment(&value,order,&mut ||rate.py().check_signals())?)))
}

#[pyfunction]
fn _stats_normal_raw(py:Python<'_>,powers:Vec<Bound<'_,PyAny>>,mean:&Bound<'_,PyAny>,covariance:&Bound<'_,PyAny>)->PyResult<Option<(PyExpr,bool)>> {
    let mut counts=Vec::with_capacity(powers.len());
    for power in powers {
        let Some(power)=order(&power)?.and_then(|power|u32::try_from(power).ok()) else{return Ok(None)};
        counts.push(power);
    }
    if counts.iter().all(|power|*power==0) {py.check_signals()?;return Ok(Some((wrap(symox::Expr::one()),true)));}
    if !stats::normal_raw_budget(&counts) {return Ok(None);}
    let mean=if let Ok(matrix)=mean.downcast::<crate::matrices::PyMatrix>() {
        matrix.borrow().m.data().to_vec()
    }else if mean.is_exact_instance_of::<pyo3::types::PyList>()||mean.is_exact_instance_of::<pyo3::types::PyTuple>() {
        let mut values=Vec::new();
        for value in mean.try_iter()? {
            let value=value?;
            if !value.is_instance_of::<crate::PyExprNode>()&&!value.is_exact_instance_of::<pyo3::types::PyInt>()&&!value.is_exact_instance_of::<pyo3::types::PyFloat>() {return Ok(None);}
            values.push(to_expr(&value)?);
        }
        values
    }else{return Ok(None)};
    let Ok(matrix)=covariance.downcast::<crate::matrices::PyMatrix>() else{return Ok(None)};
    let covariance={let matrix=matrix.borrow();if matrix.m.rows()!=counts.len()||matrix.m.cols()!=counts.len(){return Ok(None)};matrix.m.data().to_vec()};
    Ok(stats::normal_raw(&counts,&mean,&covariance,&mut ||py.check_signals())?.map(|(value,exact)|(wrap(value),exact)))
}

type Descriptor<'py>=(String,Vec<Bound<'py,PyAny>>,Vec<Bound<'py,PyAny>>,Vec<Bound<'py,PyAny>>);
#[pyfunction]
fn _stats_polynomial_query(e:&Bound<'_,PyAny>,atoms:Vec<Bound<'_,PyAny>>,models:Vec<Descriptor<'_>>,variance:bool) -> PyResult<Option<PyExpr>> {
    let atoms=atoms.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    let mut laws=Vec::new();
    for (name,params,values,weights) in models {
        let p=params.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
        laws.push(match (name.as_str(),p.as_slice()) {
            ("Binomial",[n,p,s,f])=>Law::Binomial{n:n.clone(),p:p.clone(),success:s.clone(),failure:f.clone()},
            ("Normal",[mean,sigma])=>Law::Normal{mean:mean.clone(),sigma:sigma.clone()},
            ("Poisson",[rate])=>Law::Poisson{rate:rate.clone()},
            ("Finite",_) if values.len()==weights.len()=>Law::Finite(values.iter().zip(&weights).map(|(v,w)|Ok((to_expr(v)?,to_expr(w)?))).collect::<PyResult<Vec<_>>>()?),
            _=>return Ok(None),
        });
    }
    Ok(stats::polynomial_query(&to_expr(e)?,&atoms,&laws,variance,&mut ||e.py().check_signals())?.map(wrap))
}

pub fn register(m:&Bound<'_,PyModule>)->PyResult<()> {
    m.add_function(wrap_pyfunction!(_stats_binomial_moment,m)?)?;
    m.add_function(wrap_pyfunction!(_stats_normal_moment,m)?)?;
    m.add_function(wrap_pyfunction!(_stats_poisson_moment,m)?)?;
    m.add_function(wrap_pyfunction!(_stats_normal_raw,m)?)?;
    m.add_function(wrap_pyfunction!(_stats_polynomial_query,m)?)?;
    Ok(())
}
