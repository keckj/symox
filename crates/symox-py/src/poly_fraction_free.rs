//! Typed FF adapters; Python sets reproduce representation order, not algebra.
use pyo3::prelude::*;
use pyo3::types::{PyBool,PyDict,PyInt,PyList,PySet,PyTuple};
use symox::poly::{Coeff,Int,Rat,fraction_free::{self,OrderedMap,SetOrder},domains::ExactDomain};
use super::ddm::{self,Scalar};

pub(super) fn domain(value:&Bound<'_,PyAny>)->Option<ExactDomain> {
    if !value.is_exact_instance_of::<super::domains::PyDomain>(){return None;}
    let value=value.downcast::<super::domains::PyDomain>().ok()?.get();
    match value.inner {ExactDomain::ZZ=>Some(ExactDomain::ZZ),ExactDomain::QQ=>Some(ExactDomain::QQ),_=>None}
}
pub(super) fn failure<C:Coeff>(error:fraction_free::Failure<PyErr,C>)->PyErr {
    match error {
        fraction_free::Failure::Interrupted(error)=>error,
        fraction_free::Failure::NotExact(_,divisor) if divisor.is_zero()=>pyo3::exceptions::PyZeroDivisionError::new_err("division by zero"),
        fraction_free::Failure::NotExact(_,_)=>super::domains::ExactQuotientFailed::new_err("division is not exact in the selected ring"),
    }
}
fn integer_key(value:&Bound<'_,PyAny>)->PyResult<Option<usize>> {
    if !value.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let value=crate::to_expr(value)?;
    Ok(value.as_number().and_then(|v|v.to_i64()).and_then(|v|usize::try_from(v).ok()))
}
fn order(py:Python<'_>,operation:SetOrder,left:&[usize],right:&[usize])->PyResult<Vec<usize>> {
    fn dictionary<'py>(py:Python<'py>,keys:&[usize])->PyResult<Bound<'py,PyDict>> {
        let out=PyDict::new(py);for (i,&key) in keys.iter().enumerate(){if i%64==0{py.check_signals()?;}out.set_item(key,py.None())?;}Ok(out)
    }
    let a=dictionary(py,left)?;let b=dictionary(py,right)?;
    let (a,b,method)=match operation {
        SetOrder::ViewIntersection=>(a.call_method0("keys")?,b.call_method0("keys")?,"__and__"),
        SetOrder::Difference=>(py.get_type::<PySet>().call1((a,))?,py.get_type::<PySet>().call1((b,))?,"__sub__"),
        SetOrder::Intersection=>(py.get_type::<PySet>().call1((a,))?,py.get_type::<PySet>().call1((b,))?,"__and__"),
    };
    a.call_method1(method,(b,))?.try_iter()?.map(|key|key?.extract::<usize>()).collect()
}
fn build<C:Coeff>(keys:&[Vec<usize>],values:Vec<Vec<C>>)->Vec<OrderedMap<C>> {
    keys.iter().zip(values).map(|(keys,values)|{let mut row=OrderedMap::default();for (&j,value) in keys.iter().zip(values){row.insert(j,value);}row}).collect()
}
fn output<C:Coeff>(py:Python<'_>,rows:Vec<OrderedMap<C>>,denominator:C,pivots:Vec<usize>,domain:ExactDomain)->PyResult<(PyObject,PyObject,Vec<usize>)> {
    let den=ddm::element_object(py,denominator.to_rat(),domain.clone())?;
    let out=PyDict::new(py);
    for (i,(row,&pivot)) in rows.into_iter().zip(&pivots).enumerate() {
        py.check_signals()?;
        let values=PyDict::new(py);
        for (offset,(j,value)) in row.iter().enumerate() {
            if offset%64==0{py.check_signals()?;}
            if j==pivot{values.set_item(j,den.bind(py))?;}
            else{values.set_item(j,ddm::element_object(py,value.to_rat(),domain.clone())?)?;}
        }
        out.set_item(i,values)?;
    }
    py.check_signals()?;
    Ok((out.into_any().unbind(),den,pivots))
}

#[pyfunction]
#[pyo3(signature=(a,k,stock=None))]
fn _sdm_rref_den<'py>(a:&Bound<'py,PyAny>,k:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,py:Python<'py>)->PyResult<Option<(PyObject,PyObject,Vec<usize>)>> {
    py.check_signals()?;
    let Some(domain)=domain(k) else{return Ok(None)};
    if !a.is_exact_instance_of::<PyDict>()&&!stock.is_some_and(|stock|a.get_type().is(stock)){return Ok(None);}
    let Ok(a)=a.downcast::<PyDict>() else{return Ok(None)};
    // The original one-row branch returns a shallow copy and its existing
    // pivot object, so leave this inexpensive protocol-preserving branch intact.
    if a.len()==1{return Ok(None);}
    let mut rows=Vec::with_capacity(a.len());let mut aliases=rustc_hash::FxHashSet::default();
    for (i,row) in a.iter() {
        py.check_signals()?;
        let Some(i)=integer_key(&i)? else{return Ok(None)};
        if !row.is_exact_instance_of::<PyDict>()||!aliases.insert(row.as_ptr() as usize){return Ok(None);}
        let row=row.downcast::<PyDict>()?;
        if row.is_empty(){return Ok(None);}
        let mut keys=Vec::with_capacity(row.len());let mut values=Vec::with_capacity(row.len());
        for (offset,(key,value)) in row.iter().enumerate() {
            if offset%64==0{py.check_signals()?;}
            let Some(key)=integer_key(&key)? else{return Ok(None)};
            keys.push(key);values.push(value);
        }
        rows.push((i,keys,values));
    }
    rows.sort_by_key(|(i,_,_)|*i);
    let (keys,values):(Vec<_>,Vec<_>)=rows.into_iter().map(|(_,keys,values)|(keys,values)).unzip();
    let mut ordering=|operation,left:&[usize],right:&[usize]|order(py,operation,left,right);
    if domain==ExactDomain::ZZ {
        let Some(values)=ddm::integers(&values,Scalar::ZZ,py)? else{return Ok(None)};
        if values.iter().flatten().any(Coeff::is_zero){return Ok(None);}
        let Some((rows,den,pivots))=fraction_free::sparse(build(&keys,values),&mut ||py.check_signals(),&mut ordering).map_err(failure::<Int>)? else{return Ok(None)};
        Ok(Some(output(py,rows,den,pivots,domain)?))
    }else {
        let Some(values)=ddm::rationals(&values,py)? else{return Ok(None)};
        if values.iter().flatten().any(Coeff::is_zero){return Ok(None);}
        let Some((rows,den,pivots))=fraction_free::sparse(build(&keys,values),&mut ||py.check_signals(),&mut ordering).map_err(failure::<Rat>)? else{return Ok(None)};
        Ok(Some(output(py,rows,den,pivots,domain)?))
    }
}

// Parsing uses builtin containers and exact native QQ elements only. Any
// foreign protocol remains on the original Python path, without being probed.
#[pyfunction]
fn _dm_clear_denoms<'py>(a: &Bound<'py, PyAny>, shape: &Bound<'py, PyAny>,
    k: &Bound<'py, PyAny>, convert: &Bound<'py, PyAny>, rowwise: bool,
    stock: &Bound<'py, PyAny>, py: Python<'py>) -> PyResult<Option<(PyObject, PyObject, PyObject)>> {
    py.check_signals()?;
    if domain(k) != Some(ExactDomain::QQ) || !convert.is_exact_instance_of::<PyBool>()
        || !a.get_type().is(stock) || !shape.is_exact_instance_of::<PyTuple>() { return Ok(None); }
    let shape_tuple = shape.downcast::<PyTuple>()?;
    if shape_tuple.len() != 2 { return Ok(None); }
    let (Some(m), Some(n)) = (integer_key(&shape_tuple.get_item(0)?)?, integer_key(&shape_tuple.get_item(1)?)?) else { return Ok(None); };
    if rowwise && m > isize::MAX as usize { return Ok(None); }
    let convert = convert.extract::<bool>()?;
    let mut indices = Vec::new();
    let mut columns = Vec::new();
    let mut originals = Vec::new();
    let mut key_objects = Vec::new();
    if let Ok(a) = a.downcast::<PyList>() {
        if a.len() != m { return Ok(None); }
        for (i, row) in a.iter().enumerate() {
            py.check_signals()?;
            if !row.is_exact_instance_of::<PyList>() { return Ok(None); }
            let row = row.downcast::<PyList>()?;
            if row.len() != n { return Ok(None); }
            let mut keys = Vec::new(); let mut values = Vec::new();
            for (j, value) in row.iter().enumerate() {
                if j % 64 == 0 { py.check_signals()?; }
                if ddm::scalar(&value) != Some(Scalar::QQ) { return Ok(None); }
                let element = value.downcast::<super::domains::PyDomainElement>()?.get();
                let symox::poly::domains::ElementValue::Rational(number) = &element.value.value else { unreachable!() };
                if !number.is_zero() { keys.push(j); values.push(value); }
            }
            if !values.is_empty() { indices.push(i); columns.push(keys); originals.push(values); key_objects.push(None); }
        }
    } else if let Ok(a) = a.downcast::<PyDict>() {
        for (i, row) in a.iter() {
            py.check_signals()?;
            let row_key = i.clone();
            let Some(i) = integer_key(&i)? else { return Ok(None); };
            if i >= m || !row.is_exact_instance_of::<PyDict>() { return Ok(None); }
            let row = row.downcast::<PyDict>()?;
            if row.is_empty() && rowwise { return Ok(None); }
            let mut keys = Vec::new(); let mut values = Vec::new(); let mut objects = Vec::new();
            for (offset, (j, value)) in row.iter().enumerate() {
                if offset % 64 == 0 { py.check_signals()?; }
                objects.push(j.clone());
                let Some(j) = integer_key(&j)? else { return Ok(None); };
                if j >= n { return Ok(None); }
                keys.push(j); values.push(value);
            }
            indices.push(i); columns.push(keys); originals.push(values); key_objects.push(Some((row_key, objects)));
        }
    } else { return Ok(None); }
    let Some(values) = ddm::rationals(&originals, py)? else { return Ok(None); };
    // Explicit sparse zeros have historical dup_strip semantics (including
    // index shifts) and must retain that complete, observable fallback.
    if values.iter().flatten().any(Coeff::is_zero) { return Ok(None); }
    let (factors, numerators) = fraction_free::clear_denominators(&values, rowwise, &mut || py.check_signals())?;
    let target = if convert { ExactDomain::ZZ } else { ExactDomain::QQ };
    let dod = PyDict::new(py); let flat = PyList::empty(py); let flat_indices = PyList::empty(py);
    let diagonals = PyList::empty(py);
    if rowwise {
        let one = ddm::element_object(py, Rat::one(), ExactDomain::QQ)?;
        for i in 0..m { if i % 64 == 0 { py.check_signals()?; } diagonals.append(one.bind(py))?; }
    }
    for (offset, ((&i, keys), numbers)) in indices.iter().zip(&columns).zip(numerators).enumerate() {
        py.check_signals()?;
        let factor = &factors[if rowwise { offset } else { 0 }];
        if rowwise { diagonals.set_item(i, ddm::element_object(py, Rat::from(factor.clone()), ExactDomain::ZZ)?)?; }
        let row = PyDict::new(py);
        for (j, (&column, number)) in keys.iter().zip(numbers).enumerate() {
            if j % 64 == 0 { py.check_signals()?; }
            let value = if factor.is_one() && !convert { originals[offset][j].clone().unbind() }
                else { ddm::element_object(py, Rat::from(number), target.clone())? };
            if let Some((row_key, columns)) = &key_objects[offset] {
                if rowwise { row.set_item(&columns[j], value)?; }
                else { flat.append(value)?; flat_indices.append((row_key, &columns[j]))?; }
            } else if rowwise { row.set_item(column, value)?; }
            else { flat.append(value)?; flat_indices.append((i,column))?; }
        }
        if rowwise {
            if let Some((row_key, _)) = &key_objects[offset] { dod.set_item(row_key, row)?; }
            else { dod.set_item(i,row)?; }
        }
    }
    py.check_signals()?;
    if rowwise { Ok(Some((diagonals.into_any().unbind(), dod.into_any().unbind(), py.None()))) }
    else {
        let denominator = ddm::element_object(py, Rat::from(factors[0].clone()), ExactDomain::ZZ)?;
        let data = (PyTuple::new(py,flat_indices.iter())?,shape).into_pyobject(py)?.into_any().unbind();
        Ok(Some((denominator,flat.into_any().unbind(),data)))
    }
}


// AUTO preparation reads builtin containers and native exact scalars only.
// Sparse short-circuits intentionally do not inspect scalar values or keys:
// the old heuristic only measures dictionary lengths before those decisions.
#[pyfunction]
fn _dm_rref_auto<'py>(a:&Bound<'py,PyAny>,shape:&Bound<'py,PyAny>,k:&Bound<'py,PyAny>,
    stock:&Bound<'py,PyAny>,py:Python<'py>)->PyResult<Option<&'static str>> {
    use fraction_free::{AutoStats,AutoPrefix};
    py.check_signals()?;
    let Some(domain)=domain(k) else{return Ok(None)};
    if !a.get_type().is(stock)||!shape.is_exact_instance_of::<PyTuple>() {return Ok(None);}
    let shape=shape.downcast::<PyTuple>()?;
    if shape.len()!=2{return Ok(None);}
    let dimension=|value:&Bound<'_,PyAny>| {
        if value.is_exact_instance_of::<PyInt>() {value.extract::<usize>().ok()}else{None}
    };
    let (Some(m),Some(n))=(dimension(&shape.get_item(0)?),dimension(&shape.get_item(1)?)) else{return Ok(None)};
    let rational=domain==ExactDomain::QQ;
    let mut rows=0;let mut entries=0u128;let mut values=Vec::new();let mut sparse=None;
    if let Ok(a)=a.downcast::<PyDict>() {
        let mut refs=Vec::with_capacity(a.len());
        for (_,row) in a.iter() {
            py.check_signals()?;
            if !row.is_exact_instance_of::<PyDict>(){return Ok(None);}
            let row=row.downcast_into::<PyDict>()?;
            entries+=row.len() as u128;refs.push(row);
        }
        rows=a.len();sparse=Some((a,refs));
    }else if let Ok(a)=a.downcast::<PyList>() {
        if a.len()!=m{return Ok(None);}
        for row in a.iter() {
            py.check_signals()?;
            if !row.is_exact_instance_of::<PyList>(){return Ok(None);}
            let row=row.downcast::<PyList>()?;
            if row.len()!=n{return Ok(None);}
            let before=values.len();
            for (j,value) in row.iter().enumerate() {
                if j%64==0{py.check_signals()?;}
                let Some(kind)=ddm::scalar(&value) else{return Ok(None)};
                let zero=match kind {
                    Scalar::Integer=>return Ok(None),
                    Scalar::ZZ if !rational=>{
                        let element=value.downcast::<super::domains::PyDomainElement>()?.get();
                        let symox::poly::domains::ElementValue::Rational(number)=&element.value.value else{unreachable!()};number.is_zero()
                    },
                    Scalar::QQ if rational=>{
                        let element=value.downcast::<super::domains::PyDomainElement>()?.get();
                        let symox::poly::domains::ElementValue::Rational(number)=&element.value.value else{unreachable!()};number.is_zero()
                    },
                    _=>return Ok(None),
                };
                if !zero{values.push(value);}
            }
            let length=values.len()-before;
            if length!=0{rows+=1;entries+=length as u128;}
        }
    }else{return Ok(None);}
    let stats=AutoStats{rows,entries,columns:n};
    let step=fraction_free::auto_prefix(stats,rational);
    if let AutoPrefix::Done(method)=step{return Ok(Some(method.name()));}
    if let Some((a,refs))=sparse {
        for ((i,_),row) in a.iter().zip(refs) {
            py.check_signals()?;
            // to_dok hashes (row,column) pairs on the established path.
            // Foreign keys must retain those callbacks when values are read.
            if !i.is_exact_instance_of::<PyInt>(){return Ok(None);}
            for (offset,(j,value)) in row.iter().enumerate() {
                if offset%64==0{py.check_signals()?;}
                if !j.is_exact_instance_of::<PyInt>(){return Ok(None);}
                values.push(value);
            }
        }
    }
    if step==AutoPrefix::IntegerBits {
        let mut bits=None;
        for (j,value) in values.iter().enumerate() {
            if j%64==0{py.check_signals()?;}
            let number=match ddm::scalar(value) {
                Some(Scalar::Integer)=>Int::from(value.extract::<num_bigint::BigInt>()?),
                Some(Scalar::ZZ)=>{
                    let element=value.downcast::<super::domains::PyDomainElement>()?.get();
                    let symox::poly::domains::ElementValue::Rational(number)=&element.value.value else{unreachable!()};number.n.clone()
                },
                _=>return Ok(None),
            };
            bits=Some(bits.unwrap_or(0).max(number.bits()));
        }
        Ok(fraction_free::auto_integer_bits(stats,bits.unwrap_or(1)).map(|method|method.name()))
    }else {
        let mut fractions=Vec::with_capacity(values.len());
        for (j,value) in values.iter().enumerate() {
            if j%64==0{py.check_signals()?;}
            let value=match ddm::scalar(value) {
                Some(Scalar::Integer)=>Rat::from(Int::from(value.extract::<num_bigint::BigInt>()?)),
                Some(Scalar::QQ)=>{
                    let element=value.downcast::<super::domains::PyDomainElement>()?.get();
                    let symox::poly::domains::ElementValue::Rational(number)=&element.value.value else{unreachable!()};number.clone()
                },
                _=>return Ok(None),
            };
            fractions.push(value);
        }
        Ok(Some(fraction_free::auto_rational_denominators(&fractions,&mut ||py.check_signals())?.name()))
    }
}



// Persistent CPython sets are representation oracles only. The Rust kernel
// checks every membership against its independent mathematical set state.
// Retaining original key objects also preserves pivot/dictionary key aliases.
struct IrrefPythonSets<'py> {
    py:Python<'py>,sets:Vec<Option<Bound<'py,PySet>>>,free_sets:Vec<usize>,keys:Vec<Bound<'py,PyAny>>,
    origins:rustc_hash::FxHashMap<usize,usize>,
}
impl<'py> IrrefPythonSets<'py> {
    fn new(py:Python<'py>)->Self {Self{py,sets:vec![],free_sets:vec![],keys:vec![],origins:Default::default()}}
    fn key(&mut self,value:&Bound<'py,PyAny>)->PyResult<Option<symox::poly::sdm::IrrefKey>> {
        if !value.is_exact_instance_of::<PyInt>(){return Ok(None);}
        let Ok(index)=value.extract::<usize>() else{return Ok(None)};
        let pointer=value.as_ptr() as usize;
        let origin=if let Some(&origin)=self.origins.get(&pointer){origin}else{
            let origin=self.keys.len();self.keys.push(value.clone());self.origins.insert(pointer,origin);origin
        };
        Ok(Some(symox::poly::sdm::IrrefKey{index,origin}))
    }
    fn collect(&mut self,set:&Bound<'py,PySet>)->PyResult<Vec<symox::poly::sdm::IrrefKey>> {
        let references=set.iter().collect::<Vec<_>>();
        references.iter().map(|key|self.key(key)?.ok_or_else(||pyo3::exceptions::PyTypeError::new_err("invalid native sparse key"))).collect()
    }
}
impl symox::poly::sdm::IrrefSets<PyErr> for IrrefPythonSets<'_> {
    fn create(&mut self,keys:&[symox::poly::sdm::IrrefKey],from_dictionary:bool)->PyResult<usize> {
        self.py.check_signals()?;
        let set=if from_dictionary {
            // set(dict) uses CPython's dictionary merge/preallocation path;
            // set(list(dict)) can have a different table/iteration order.
            let dictionary=PyDict::new(self.py);
            for (i,key)in keys.iter().enumerate(){if i%64==0{self.py.check_signals()?;}dictionary.set_item(&self.keys[key.origin],self.py.None())?;}
            self.py.get_type::<PySet>().call1((dictionary,))?.downcast_into::<PySet>()?
        }else {
            let set=PySet::empty(self.py)?;
            for key in keys{set.add(&self.keys[key.origin])?;}set
        };
        let id=if let Some(id)=self.free_sets.pop(){self.sets[id]=Some(set);id}
            else{let id=self.sets.len();self.sets.push(Some(set));id};Ok(id)
    }
    fn release(&mut self,set:usize)->PyResult<()> {
        if self.sets[set].take().is_some(){self.free_sets.push(set);}Ok(())
    }
    fn insert(&mut self,set:usize,key:symox::poly::sdm::IrrefKey)->PyResult<()> {self.sets[set].as_ref().unwrap().add(&self.keys[key.origin])}
    fn remove(&mut self,set:usize,key:symox::poly::sdm::IrrefKey)->PyResult<()> {
        self.sets[set].as_ref().unwrap().call_method1("remove",(&self.keys[key.origin],))?;Ok(())
    }
    fn members(&mut self,set:usize)->PyResult<Vec<symox::poly::sdm::IrrefKey>> {
        let set=self.sets[set].as_ref().unwrap().clone();self.collect(&set)
    }
    fn binary(&mut self,op:symox::poly::sdm::IrrefSetOp,left:usize,right:usize)->PyResult<Vec<symox::poly::sdm::IrrefKey>> {
        let method=match op {symox::poly::sdm::IrrefSetOp::Difference=>"__sub__",symox::poly::sdm::IrrefSetOp::Intersection=>"__and__"};
        let result=self.sets[left].as_ref().unwrap().call_method1(method,(self.sets[right].as_ref().unwrap(),))?.downcast_into::<PySet>()?;
        self.collect(&result)
    }
}

#[pyfunction]
#[pyo3(signature=(a,stock=None))]
fn _sdm_irref<'py>(a:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,py:Python<'py>)
    ->PyResult<Option<(PyObject,PyObject,PyObject)>> {
    py.check_signals()?;
    if !a.is_exact_instance_of::<PyDict>()&&!stock.is_some_and(|s|a.get_type().is(s)){return Ok(None);}
    let Ok(a)=a.downcast::<PyDict>() else{return Ok(None)};
    // No arithmetic domain can be inferred from an empty free-function input.
    // Its cheap original branch also preserves customized container factories.
    if a.is_empty(){return Ok(None);}
    // The domain is read from the first stored value without any callback;
    // every other value must be a stock element of exactly that domain.
    let first=a.iter().next().and_then(|(_,row)|row.downcast_exact::<PyDict>().ok().and_then(|row|row.iter().next().map(|(_,v)|v)));
    let Some(first)=first else{return Ok(None)};
    if ddm::scalar(&first)==Some(Scalar::QQ) {
        return sdm_irref_typed::<Rat>(a,py,&|value|{
            if ddm::scalar(value)!=Some(Scalar::QQ){return None;}
            let element=value.downcast::<super::domains::PyDomainElement>().ok()?.get();
            let symox::poly::domains::ElementValue::Rational(value)=&element.value.value else{return None};
            Some(value.clone())
        },&|py,value|ddm::element_object(py,value,ExactDomain::QQ));
    }
    // GF(p), QQ<a>, QQ_I: the Python body uses only the native
    // pow/mul/sub/neg/bool operations authenticated by the caller.
    let Some((domain,kind))=ddm::field_domain(&first,true) else{return Ok(None)};
    let result=match kind {
        ddm::FieldKind::Finite=>ddm::with_field_scope(&domain,||sdm_irref_typed::<symox::poly::f5b::Fp>(a,py,
            &|value|ddm::field_value(value,&domain),&|py,value|ddm::field_object(py,value,&domain)))?,
        ddm::FieldKind::Algebraic=>ddm::with_field_scope(&domain,||sdm_irref_typed::<symox::poly::algebraic::AN>(a,py,
            &|value|ddm::field_value(value,&domain),&|py,value|ddm::field_object(py,value,&domain)))?,
    };
    Ok(result.flatten())
}

type IrrefOutput=Option<(PyObject,PyObject,PyObject)>;

fn sdm_irref_typed<'py,C:symox::poly::Field>(a:&Bound<'py,PyDict>,py:Python<'py>,
    extract:&dyn Fn(&Bound<'py,PyAny>)->Option<C>,wrap:&dyn Fn(Python<'py>,C)->PyResult<PyObject>)->PyResult<IrrefOutput> {
    let references=a.iter().map(|(_,row)|row).collect::<Vec<_>>();
    let mut representation=IrrefPythonSets::new(py);let mut rows=Vec::with_capacity(references.len());
    for row in references {
        py.check_signals()?;
        if !row.is_exact_instance_of::<PyDict>(){return Ok(None);}
        let row=row.downcast::<PyDict>()?;
        if row.is_empty(){return Ok(None);}
        let entries=row.iter().collect::<Vec<_>>();
        let mut values:symox::poly::sdm::IrrefRow<C>=Default::default();
        for (i,(key,value))in entries.iter().enumerate() {
            if i%64==0{py.check_signals()?;}
            let Some(key)=representation.key(key)? else{return Ok(None)};
            let Some(value)=extract(value) else{return Ok(None)};
            // Explicit stored zeros have historical errors/structural metadata.
            // Do not silently delete them or change which pivot is selected.
            if value.is_zero(){return Ok(None);}
            values.insert(key.index,(key.origin,value));
        }
        rows.push(values);
    }
    let Some(result)=symox::poly::sdm::irref(rows,&mut representation,&mut ||py.check_signals())? else{return Ok(None)};
    let output=PyDict::new(py);let pivots=PyList::empty(py);let columns=PyDict::new(py);
    // The old pivot2row numbering objects are shared across metadata sets;
    // rref's separate enumerate has distinct large integer key objects.
    let mut numbering=Vec::with_capacity(result.rows.len());
    for i in 0..result.rows.len(){numbering.push(i.into_pyobject(py)?.into_any());}
    for (i,row)in result.rows.into_iter().enumerate() {
        py.check_signals()?;let out=PyDict::new(py);
        for (step,(_,(origin,value)))in row.iter().enumerate(){
            if step%64==0{py.check_signals()?;}
            out.set_item(&representation.keys[*origin],wrap(py,value.clone())?)?;
        }
        output.set_item(i,out)?;
    }
    for key in result.pivots{pivots.append(&representation.keys[key.origin])?;}
    for (key,indices)in result.nonzero_columns {
        py.check_signals()?;let set=PySet::empty(py)?;
        for index in indices{set.add(&numbering[index])?;}
        columns.set_item(&representation.keys[key.origin],set)?;
    }
    py.check_signals()?;
    Ok(Some((output.into_any().unbind(),pivots.into_any().unbind(),columns.into_any().unbind())))
}

pub fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_sdm_rref_den,module)?)?;
    module.add_function(wrap_pyfunction!(_dm_clear_denoms,module)?)?;
    module.add_function(wrap_pyfunction!(_dm_rref_auto,module)?)?;
    module.add_function(wrap_pyfunction!(_sdm_irref,module)?)?;
    Ok(())
}
