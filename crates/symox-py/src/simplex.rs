//! Protocol adapter for exact rational simplex inputs, without scalar hooks.
use pyo3::prelude::*;
use symox::{poly::Rat, simplex::{self, Bounds, Constraint, Failure, Relation}};
use crate::{PyExpr, wrap};

fn rational(value: &Bound<'_, PyAny>) -> Option<Rat> {
    let node = value.extract::<PyRef<'_, crate::PyExprNode>>().ok()?;
    let standard = crate::node_type(value.py(), &node.e)?;
    if !value.get_type().is(&standard) { return None; }
    Rat::from_number(node.e.as_number()?)
}

type ResultObject = (String, Option<(PyExpr, Vec<PyExpr>)>);

#[pyfunction]
fn _native_simplex(cost: Vec<Bound<'_, PyAny>>,
    constraints: Vec<(Vec<Bound<'_, PyAny>>, String, Bound<'_, PyAny>)>,
    py: Python<'_>) -> PyResult<Option<ResultObject>> {
    let Some(cost) = cost.iter().map(rational).collect::<Option<Vec<_>>>() else { return Ok(None) };
    let mut native = Vec::with_capacity(constraints.len());
    for (coefficients, relation, rhs) in constraints {
        let Some(coefficients) = coefficients.iter().map(rational).collect::<Option<Vec<_>>>() else { return Ok(None) };
        let Some(rhs) = rational(&rhs) else { return Ok(None) };
        let relation = match relation.as_str() { "<=" => Relation::Le, ">=" => Relation::Ge, "=" => Relation::Eq, _ => return Ok(None) };
        native.push(Constraint { coefficients, relation, rhs });
    }
    let result = simplex::solve(&cost, &native, &mut || py.check_signals())?;
    Ok(Some(match result {
        Ok((optimum, solution)) => ("ok".into(), Some((wrap(symox::Expr::num(optimum.to_number())), solution.into_iter().map(|v| wrap(symox::Expr::num(v.to_number()))).collect()))),
        Err(Failure::Infeasible) => ("infeasible".into(), None),
        Err(Failure::Unbounded) => ("unbounded".into(), None),
        Err(Failure::InvalidShape|Failure::InvalidBounds) => return Ok(None),
    }))
}

// Eligibility never calls an iterator, conversion, equality or indexing hook
// on a user-defined object. Unsupported inputs keep the complete general API.
fn sequence<'py>(value:&Bound<'py,PyAny>)->Option<Vec<Bound<'py,PyAny>>> {
    use pyo3::types::{PyList,PyTuple};
    if value.is_exact_instance_of::<PyList>() {Some(value.downcast::<PyList>().ok()?.iter().collect())}
    else if value.is_exact_instance_of::<PyTuple>() {Some(value.downcast::<PyTuple>().ok()?.iter().collect())}
    else {None}
}

fn scalar(value:&Bound<'_,PyAny>)->PyResult<Option<Rat>> {
    if value.is_exact_instance_of::<pyo3::types::PyInt>() {
        let value=crate::to_expr(value)?;
        return Ok(value.as_number().and_then(Rat::from_number));
    }
    Ok(rational(value))
}

fn vector(value:&Bound<'_,PyAny>)->PyResult<Option<Vec<Rat>>> {
    let Some(values)=sequence(value) else{return Ok(None)};
    let mut result=Vec::with_capacity(values.len());
    for (i,value) in values.into_iter().enumerate() {
        if i%64==0 {value.py().check_signals()?;}
        let Some(value)=scalar(&value)? else{return Ok(None)};
        result.push(value);
    }
    Ok(Some(result))
}

fn endpoint(value:&Bound<'_,PyAny>,lower:bool)->PyResult<Option<Option<Rat>>> {
    if value.is_none() {return Ok(Some(None));}
    if let Some(value)=scalar(value)? {return Ok(Some(Some(value)));}
    if let Ok(node)=value.extract::<PyRef<'_,crate::PyExprNode>>() {
        if crate::node_type(value.py(),&node.e).is_some_and(|standard|value.get_type().is(&standard))
            && matches!(node.e.kind(),symox::Kind::Const(constant) if *constant==if lower{symox::Constant::NegInfinity}else{symox::Constant::Infinity}) {
            return Ok(Some(None));
        }
    }
    Ok(None)
}

fn bounds(value:&Bound<'_,PyAny>,n:usize)->PyResult<Option<Vec<Bounds>>> {
    use pyo3::types::{PyDict,PyInt,PyList,PyTuple};
    if value.is_none() {return Ok(Some(vec![Bounds{lower:Some(Rat::zero()),upper:None};n]));}
    let pairs=if value.is_exact_instance_of::<PyDict>() {
        let dictionary=value.downcast::<PyDict>()?;
        // In particular, do not look up Python int keys in a dictionary that
        // could compare them to an equal custom key with an observable hook.
        if dictionary.keys().iter().any(|key|!key.is_exact_instance_of::<PyInt>()) {return Ok(None);}
        let mut pairs=vec![None;n];
        for (key,value) in dictionary.iter() {
            let key=crate::to_expr(&key)?;
            if let Some(i)=key.as_number().and_then(|key|key.to_i64()).and_then(|i|usize::try_from(i).ok()).filter(|i|*i<n) {pairs[i]=Some(value);}
        }
        pairs
    }else {
        let Some(values)=sequence(value) else{return Ok(None)};
        if values.len()==2&&values.iter().all(|v|v.is_none()||(!v.is_instance_of::<PyList>()&&!v.is_instance_of::<PyTuple>())) {
            vec![Some(value.clone());n]
        }else {
            if values.len()!=n {return Ok(None);}
            values.into_iter().map(Some).collect()
        }
    };
    let mut result=Vec::with_capacity(n);
    for (i,pair) in pairs.into_iter().enumerate() {
        if i%64==0 {value.py().check_signals()?;}
        let Some(pair)=pair else {result.push(Bounds{lower:Some(Rat::zero()),upper:None});continue};
        let Some(pair)=sequence(&pair).filter(|pair|pair.len()==2) else{return Ok(None)};
        let Some(lower)=endpoint(&pair[0],true)? else{return Ok(None)};
        let Some(upper)=endpoint(&pair[1],false)? else{return Ok(None)};
        result.push(Bounds{lower,upper});
    }
    Ok(Some(result))
}

#[pyfunction]
fn _native_linprog(c:&Bound<'_,PyAny>,a:&Bound<'_,PyAny>,b:&Bound<'_,PyAny>,
    a_eq:&Bound<'_,PyAny>,b_eq:&Bound<'_,PyAny>,bound:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<Option<ResultObject>> {
    py.check_signals()?;
    let Some(cost)=vector(c)? else{return Ok(None)};
    let n=cost.len();let mut constraints=Vec::new();
    for (matrix,rhs,relation) in [(a,b,Relation::Le),(a_eq,b_eq,Relation::Eq)] {
        let rows=if matrix.is_none(){Vec::new()}else {
            let Some(rows)=sequence(matrix) else{return Ok(None)};
            let mut result=Vec::with_capacity(rows.len());
            for row in rows {
                py.check_signals()?;
                let Some(row)=vector(&row)?.filter(|row|row.len()==n) else{return Ok(None)};
                result.push(row);
            }
            result
        };
        let rhs=if rhs.is_none(){Vec::new()}else {let Some(rhs)=vector(rhs)? else{return Ok(None)};rhs};
        if rows.len()!=rhs.len() {return Ok(None);}
        constraints.extend(rows.into_iter().zip(rhs).map(|(coefficients,rhs)|Constraint{coefficients,relation,rhs}));
    }
    let Some(bounds)=bounds(bound,n)? else{return Ok(None)};
    let result=simplex::linprog(&cost,&constraints,&bounds,&mut ||py.check_signals())?;
    py.check_signals()?;
    Ok(Some(match result {
        Ok((optimum,solution))=>("ok".into(),Some((wrap(symox::Expr::num(optimum.to_number())),solution.into_iter().map(|v|wrap(symox::Expr::num(v.to_number()))).collect()))),
        Err(Failure::Infeasible)=>("infeasible".into(),None),
        Err(Failure::Unbounded)=>("unbounded".into(),None),
        Err(Failure::InvalidBounds)=>("invalid_bounds".into(),None),
        Err(Failure::InvalidShape)=>return Ok(None),
    }))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_native_simplex, m)?)?;
    m.add_function(wrap_pyfunction!(_native_linprog, m)?)?;
    Ok(())
}
