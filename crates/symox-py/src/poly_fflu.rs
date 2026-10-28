//! Exact stock ZZ/QQ FFLU adapters; Python retains its real allocation barriers.
use pyo3::prelude::*;
use pyo3::exceptions::PyTypeError;
use pyo3::types::{PyDict,PyInt,PyList,PyModule,PyString,PyTuple};
use symox::poly::{Coeff,fflu,domains::ExactDomain};
use super::{Scalar,sequence,entries,integers,rationals,element_object};
use super::super::fraction_free::{domain,failure};

struct Matrix<'py>{outer:Bound<'py,PyList>,rows:Vec<Bound<'py,PyList>>,cells:Vec<Vec<Bound<'py,PyAny>>>}
fn matrix<'py>(value:&Bound<'py,PyAny>,stock:&Bound<'py,PyAny>,m:usize,n:usize)->Option<Matrix<'py>>{
    let rows=sequence(value,Some(stock),true)?;
    if rows.len()!=m{return None;}
    let mut seen=rustc_hash::FxHashSet::default();
    if rows.iter().any(|row|!row.is_exact_instance_of::<PyList>()||!seen.insert(row.as_ptr() as usize)){return None;}
    let cells=entries(&rows,true)?;
    if cells.iter().any(|row|row.len()!=n){return None;}
    Some(Matrix{outer:value.downcast::<PyList>().ok()?.clone(),rows:rows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect(),cells})
}
fn dimensions(value:&Bound<'_,PyAny>)->Option<(usize,usize)>{
    let value=value.downcast_exact::<PyTuple>().ok()?;
    if value.len()!=2{return None;}
    let m=value.get_item(0).ok()?;let n=value.get_item(1).ok()?;
    if !m.is_exact_instance_of::<PyInt>()||!n.is_exact_instance_of::<PyInt>(){return None;}
    Some((m.extract().ok()?,n.extract().ok()?))
}
fn disjoint(a:&Matrix<'_>,b:&Matrix<'_>)->bool{
    !a.outer.is(&b.outer)&&!a.rows.iter().any(|row|b.rows.iter().any(|other|row.is(other)))
}
fn indices<'py>(value:&Bound<'py,PyAny>,m:usize)->Option<Bound<'py,PyList>>{
    let value=value.downcast_exact::<PyList>().ok()?;
    if value.len()!=m{return None;}
    let mut seen=vec![false;m];
    for item in value.iter(){
        if !item.is_exact_instance_of::<PyInt>(){return None;}
        let j=item.extract::<usize>().ok()?;
        if j>=m||seen[j]{return None;}seen[j]=true;
    }
    Some(value.clone())
}

fn eliminate<C:Coeff>(mut a:Matrix<'_>,perm:&Bound<'_,PyList>,mut values:Vec<Vec<C>>,n:usize,k:ExactDomain)->PyResult<bool>{
    let py=perm.py();
    let result=fflu::eliminate(&mut values,n,&mut ||py.check_signals(),&mut |event|{
        match event{
            fflu::Eliminate::Swap(i,j)=>{
                a.outer.set_item(i,&a.rows[j])?;a.outer.set_item(j,&a.rows[i])?;a.rows.swap(i,j);
                let left=perm.get_item(i)?;let right=perm.get_item(j)?;perm.set_item(i,right)?;perm.set_item(j,left)?;Ok(())
            },
            fflu::Eliminate::Set(i,j,value)=>a.rows[i].set_item(j,element_object(py,value.to_rat(),k.clone())?),
            fflu::Eliminate::Retain(i,j)=>{let value=a.rows[i].get_item(j)?;a.rows[i].set_item(j,value)},
        }
    }).map_err(failure)?;
    Ok(result.is_some())
}
#[pyfunction]
fn _ddm_fflu(a:&Bound<'_,PyAny>,perm:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>,shape:&Bound<'_,PyAny>,stock:&Bound<'_,PyAny>)->PyResult<bool>{
    let py=a.py();py.check_signals()?;
    let Some(k)=domain(k) else{return Ok(false)};let Some((m,n))=dimensions(shape) else{return Ok(false)};
    let Some(a)=matrix(a,stock,m,n) else{return Ok(false)};let Some(perm)=indices(perm,m) else{return Ok(false)};
    match k{
        ExactDomain::ZZ=>match integers(&a.cells,Scalar::ZZ,py)?{Some(values)=>eliminate(a,&perm,values,n,k),None=>Ok(false)},
        ExactDomain::QQ=>match rationals(&a.cells,py)?{Some(values)=>eliminate(a,&perm,values,n,k),None=>Ok(false)},
        _=>Ok(false),
    }
}
#[pyfunction]
fn _ddm_fflu_p(p:&Bound<'_,PyAny>,perm:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>,shape:&Bound<'_,PyAny>,stock:&Bound<'_,PyAny>)->PyResult<bool>{
    let py=p.py();py.check_signals()?;
    let Some(k)=domain(k) else{return Ok(false)};let Some((m,_))=dimensions(shape) else{return Ok(false)};
    let Some(p)=matrix(p,stock,m,m) else{return Ok(false)};let Some(perm)=indices(perm,m) else{return Ok(false)};
    // Canonical zero matrices are admitted; values are parsed without calling
    // user scalar protocols even though only K.one is published by this phase.
    let scalar=if matches!(k,ExactDomain::ZZ){Scalar::ZZ}else{Scalar::QQ};
    if p.cells.iter().flatten().any(|v|super::scalar(v)!=Some(scalar)){return Ok(false);}
    let indices=perm.iter().map(|v|v.extract::<usize>()).collect::<PyResult<Vec<_>>>()?;
    fflu::permutation(&indices,m,&mut ||py.check_signals(),&mut |i,j|p.rows[i].set_item(j,element_object(py,symox::poly::Rat::one(),k.clone())?))
}
fn lower<C:Coeff>(u:Matrix<'_>,l:Matrix<'_>,mut uv:Vec<Vec<C>>,mut lv:Vec<Vec<C>>,n:usize,k:ExactDomain)->PyResult<bool>{
    let py=u.outer.py();
    fflu::lower(&mut uv,&mut lv,n,&mut ||py.check_signals(),&mut |event|{
        match event{
            fflu::Lower::Copy(i,j,r,c,_)=>{let value=u.rows[r].get_item(c)?;l.rows[i].set_item(j,value)},
            fflu::Lower::ZeroUpper(i,j)=>u.rows[i].set_item(j,element_object(py,symox::poly::Rat::zero(),k.clone())?),
            fflu::Lower::OneLower(i,j)=>l.rows[i].set_item(j,element_object(py,symox::poly::Rat::one(),k.clone())?),
        }
    })
}
#[pyfunction]
fn _ddm_fflu_l(u:&Bound<'_,PyAny>,l:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>,shape:&Bound<'_,PyAny>,stock:&Bound<'_,PyAny>)->PyResult<bool>{
    let py=u.py();py.check_signals()?;
    let Some(k)=domain(k) else{return Ok(false)};let Some((m,n))=dimensions(shape) else{return Ok(false)};
    let (Some(u),Some(l))=(matrix(u,stock,m,n),matrix(l,stock,m,m)) else{return Ok(false)};
    if !disjoint(&u,&l){return Ok(false);}
    match k{
        ExactDomain::ZZ=>match (integers(&u.cells,Scalar::ZZ,py)?,integers(&l.cells,Scalar::ZZ,py)?){(Some(uv),Some(lv))=>lower(u,l,uv,lv,n,k),_=>Ok(false)},
        ExactDomain::QQ=>match (rationals(&u.cells,py)?,rationals(&l.cells,py)?){(Some(uv),Some(lv))=>lower(u,l,uv,lv,n,k),_=>Ok(false)},
        _=>Ok(false),
    }
}
fn diagonal<C:Coeff>(d:Matrix<'_>,lv:Vec<Vec<C>>,mut dv:Vec<Vec<C>>,k:ExactDomain)->PyResult<bool>{
    let py=d.outer.py();fflu::diagonal(&lv,&mut dv,&mut ||py.check_signals(),&mut |i,v|d.rows[i].set_item(i,element_object(py,v.to_rat(),k.clone())?))
}
#[pyfunction]
fn _ddm_fflu_d(l:&Bound<'_,PyAny>,d:&Bound<'_,PyAny>,k:&Bound<'_,PyAny>,shape:&Bound<'_,PyAny>,stock:&Bound<'_,PyAny>)->PyResult<bool>{
    let py=l.py();py.check_signals()?;
    let Some(k)=domain(k) else{return Ok(false)};let Some((m,_))=dimensions(shape) else{return Ok(false)};
    let (Some(l),Some(d))=(matrix(l,stock,m,m),matrix(d,stock,m,m)) else{return Ok(false)};
    if !disjoint(&l,&d){return Ok(false);}
    match k{
        ExactDomain::ZZ=>match (integers(&l.cells,Scalar::ZZ,py)?,integers(&d.cells,Scalar::ZZ,py)?){(Some(lv),Some(dv))=>diagonal(d,lv,dv,k),_=>Ok(false)},
        ExactDomain::QQ=>match (rationals(&l.cells,py)?,rationals(&d.cells,py)?){(Some(lv),Some(dv))=>diagonal(d,lv,dv,k),_=>Ok(false)},
        _=>Ok(false),
    }
}
fn dictionary<'py>(value:&Bound<'py,PyAny>)->Option<Bound<'py,PyDict>>{
    let value=value.downcast_exact::<PyDict>().ok()?;
    if value.iter().any(|(key,_)|!key.is_exact_instance_of::<PyString>()){return None;}
    Some(value.clone())
}
fn lookup<'py>(g:&Bound<'py,PyDict>,b:&Bound<'py,PyDict>,name:&str)->PyResult<PyObject>{
    Ok(g.get_item(name)?.or(b.get_item(name)?).map_or_else(||g.py().None(),Bound::unbind))
}
fn dependencies(py:Python<'_>,function:&Bound<'_,PyAny>,phase:i32,metadata:Option<&PyObject>)->PyResult<Option<PyObject>>{
    let Some(metadata)=metadata else{return Ok(None)};
    let frame=unsafe{pyo3::ffi::PyEval_GetFrame()};if frame.is_null(){return Ok(None);}
    let frame=unsafe{Bound::from_borrowed_ptr(py,frame.cast())};
    let (g,b)=(frame.getattr(pyo3::intern!(py,"f_globals"))?,frame.getattr(pyo3::intern!(py,"f_builtins"))?);
    let (Some(g),Some(b))=(dictionary(&g),dictionary(&b)) else{return Ok(None)};
    let names:&[&str]=match phase{0=>&["range","all"],1=>&["enumerate"],2|3=>&["range"],_=>return Ok(None)};
    let mut out=Vec::new();for name in names{out.push(lookup(&g,&b,name)?);}
    let data=metadata.bind(py).call1((function,))?;
    let Ok(data)=data.downcast_exact::<PyTuple>() else{return Ok(None)};
    if data.len()!=5{return Ok(None);}
    let globals=data.get_item(1)?;let Some(globals)=dictionary(&globals) else{return Ok(None)};
    let Some(builtins)=globals.get_item("__builtins__")? else{return Ok(None)};
    let builtins=if let Some(builtins)=dictionary(&builtins){builtins}else if let Ok(module)=builtins.downcast_exact::<PyModule>(){let d=module.dict();let Some(d)=dictionary(d.as_any()) else{return Ok(None)};d}else{return Ok(None)};
    out.push(lookup(&globals,&builtins,"type")?);out.push(lookup(&globals,&builtins,"int")?);
    Ok(Some(PyTuple::new(py,out)?.into_any().unbind()))
}
pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()>{
    module.add_function(wrap_pyfunction!(_ddm_fflu,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_fflu_p,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_fflu_l,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_fflu_d,module)?)?;
    let metadata=if let Some(d)=dictionary(module.dict().as_any()){d.get_item("_ddm_function_metadata")?.map(Bound::unbind)}else{None};
    module.add("_ddm_fflu_dependencies",pyo3::types::PyCFunction::new_closure(module.py(),Some(c"_ddm_fflu_dependencies"),None,
        move |args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>|->PyResult<Option<PyObject>>{
            if kwargs.is_some()||args.len()!=2{return Err(PyTypeError::new_err("expected comparison helper and phase"));}
            let phase=args.get_item(1)?;if !phase.is_exact_instance_of::<PyInt>(){return Ok(None);}
            let Ok(phase)=phase.extract::<i32>() else{return Ok(None)};
            dependencies(args.py(),&args.get_item(0)?,phase,metadata.as_ref())
        })?)?;Ok(())
}
