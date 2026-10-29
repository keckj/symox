//! Live constructor/conversion proof for exact LLL admission.
use super::*;
use crate::poly::domains::{PyDomain,PyDomainElement};
use symox::poly::domains::ExactDomain;
use pyo3::types::PyCFunction;

const HELPER:&str="symox.polys.matrices._lll_native";
const METHODS:&str="symox.polys.domains._methods";
const DDM:&str="symox.polys.matrices.ddm";
const DENSE:&str="symox.polys.matrices._dense_poly";
const LLL:&str="symox.polys.matrices.lll";

fn class_method<'py>(c:&Context<'py>,value:&Bound<'py,PyAny>)->VResult<Bound<'py,PyAny>>{
    guard!(value.get_type().is(c.builtin("classmethod")?),"classmethod descriptor");
    Ok(value.getattr("__func__")?)
}

fn descriptor<'py>(c:&Context<'py>,value:&Bound<'py,PyAny>,owner:&Bound<'py,PyAny>,name:&Bound<'py,PyAny>)->VResult<()> {
    let native_type=c.py.get_type::<PyType>().into_any();let list=c.builtin("list")?;
    let getset=need(&c.raw(&native_type)?,"__dict__")?.get_type();
    let method=need(&c.raw(&list)?,"copy")?.get_type();
    let wrapper=need(&c.raw(&list)?,"__add__")?.get_type();
    guard!(value.get_type().is(&getset)||value.get_type().is(&method)||value.get_type().is(&wrapper),"native descriptor type");
    guard!(value.getattr("__objclass__")?.is(owner)&&value.getattr("__name__")?.eq(name)?,"native descriptor origin");
    Ok(())
}

fn named_descriptor<'py>(c:&Context<'py>,class:&Bound<'py,PyAny>,name:&str)->VResult<()> {
    descriptor(c,&need(&c.raw(class)?,name)?,class,PyString::new(c.py,name).as_any())
}

fn snapshot_function<'py>(c:&Context<'py>,function:&Bound<'py,PyAny>,code:&Bound<'py,PyAny>,defaults:&Bound<'py,PyAny>,kwdefaults:&Bound<'py,PyAny>,closure:&Bound<'py,PyAny>)->VResult<()> {
    guard!(c.api.is_function(function),"snapshot function type");
    guard!(c.api.code(function)?.is(code)&&c.api.defaults(function)?.is(defaults)&&c.api.kwdefaults(function)?.is(kwdefaults),"snapshot function metadata");
    let expected=tuple(closure.clone())?;let actual=c.api.closure(function)?;
    if expected.is_empty(){guard!(actual.is_none()||tuple(actual)?.is_empty(),"snapshot empty closure");}else{
        let actual=tuple(actual)?;guard!(actual.len()==expected.len(),"snapshot closure length");
        for (cell,value) in actual.iter().zip(expected.iter()){
            let current=match cell.getattr("cell_contents"){Ok(value)=>value,Err(error) if error.is_instance_of::<pyo3::exceptions::PyValueError>(c.py)=>return Err(Failure::Decline("snapshot empty cell")),Err(error)=>return Err(error.into())};
            guard!(current.is(&value),"snapshot closure contents");
        }
    }
    Ok(())
}

fn inherited_native<'py>(c:&Context<'py>,class:&Bound<'py,PyAny>,owner:&Bound<'py,PyAny>,names:&[&str])->VResult<()> {
    for name in names{guard!(c.effective(class,name)?.is(need(&c.raw(owner)?,name)?),"inherited native protocol");}Ok(())
}

fn slot_cache(namespace:&Bound<'_,PyAny>)->VResult<bool> {
    Ok(match item(namespace,"__slotnames__")?{Some(value)=>value.downcast_exact::<PyList>().is_ok_and(|list|list.is_empty()),None=>false})
}

fn rank_error<'py>(c:&Context<'py>)->VResult<()> {
    c.builtins(&["Exception","BaseException"])?;
    let module="symox.polys.matrices.exceptions";
    let rank=c.member(module,"DMRankError")?;let base=c.member(module,"DMError")?;
    guard!(c.member(LLL,"DMRankError")?.is(&rank)&&c.member(HELPER,"_RANK_ERROR")?.is(&rank),"rank exception aliases");
    let kind=c.py.get_type::<PyType>();
    guard!(c.state(&rank)?.meta.is(&kind)&&c.state(&base)?.meta.is(&kind)
        &&same_tuple(&c.state(&rank)?.bases,&[&base])&&same_tuple(&c.state(&base)?.bases,&[&c.builtin("Exception")?]),"rank exception hierarchy");
    for class in [&rank,&base]{
        let namespace=c.raw(class)?;
        for key in namespace.try_iter()?{
            let key=key?;let mut allowed=false;
            for name in ["__module__","__doc__","__weakref__"]{if key.eq(name)?{allowed=true;break;}}
            // copy/pickle store copyreg's __slotnames__ cache in the class dict.
            // No class of this verified hierarchy declares __slots__, so the
            // only interpreter-written value is an exact empty list.
            if !allowed&&key.eq("__slotnames__")?{allowed=slot_cache(&namespace)?;}
            guard!(allowed,"rank exception custom attribute");
        }
        let name=need(&namespace,"__module__")?;guard!(name.downcast_exact::<PyString>().is_ok()&&name.eq(module)?,"rank exception module");
        if let Some(value)=item(&namespace,"__weakref__")?{descriptor(c,&value,class,PyString::new(c.py,"__weakref__").as_any())?;}
    }
    let records=tuple(c.member(HELPER,"_RANK_CLASSES")?)?;
    guard!(records.len()==2,"rank exception captured classes");
    for record in records.iter(){
        let record=tuple(record)?;guard!(record.len()==2,"rank class record");let class=record.get_item(0)?;
        guard!(class.is(&rank)||class.is(&base),"rank class identity");let namespace=c.raw(&class)?;
        let attributes=tuple(record.get_item(1)?)?;
        let mut recorded=0usize;
        for attribute in attributes.iter(){let attribute=tuple(attribute)?;guard!(attribute.len()==2,"rank attribute record");let key=attribute.get_item(0)?;guard!(key.downcast_exact::<PyString>().is_ok(),"rank attribute name");if key.eq("__slotnames__")?{continue;}recorded+=1;guard!(namespace.get_item(key)?.is(attribute.get_item(1)?),"rank class attribute changed");}
        let cached=usize::from(item(&namespace,"__slotnames__")?.is_some());
        guard!(namespace.len()?==recorded+cached,"rank class namespace length");
    }
    let captured=c.member(HELPER,"_RANK_FUNCTIONS")?.downcast_into_exact::<PyList>().map_err(|_|Failure::Decline("rank function records"))?;
    guard!(captured.is_empty(),"rank custom constructor function");
    Ok(())
}

struct Protocol<'py>{ddm:Bound<'py,PyAny>,element:Bound<'py,PyAny>}
fn live_protocol<'py>(c:&Context<'py>)->VResult<Protocol<'py>> {
    for name in ["_CAPTURED","_BOOTSTRAP_STOCK"]{let flag=c.member(HELPER,name)?;guard!(flag.downcast_exact::<PyBool>().is_ok()&&flag.is_truthy()?,"helper bootstrap declined");}
    c.builtins(&["type","object","list","tuple","dict","str","bool","int","isinstance","issubclass","super","classmethod","staticmethod","sum","abs","len","all","any","zip","getattr","hasattr","range","min","max","enumerate","float","complex","bytes","bytearray","set","frozenset","ZeroDivisionError","AttributeError","TypeError"])?;
    let zero_error=c.ns(LLL)?.get_item("ZeroDivisionError")?.unwrap_or(c.builtin("ZeroDivisionError")?);
    guard!(zero_error.is(c.builtin("ZeroDivisionError")?),"rank division exception alias");
    let native_type=c.py.get_type::<PyType>().into_any();let object=c.builtin("object")?;let list=c.builtin("list")?;
    for (alias,name) in [("_TYPE","type"),("_ALL","all"),("_ANY","any"),("_LEN","len"),("_ZIP","zip"),("_DICT","dict"),("_STR","str"),("_OBJECT","object"),("_TUPLE","tuple"),("_LIST","list"),("_GETATTR","getattr")]{guard!(c.member(HELPER,alias)?.is(c.builtin(name)?),"helper builtin aliases");}
    for (alias,module) in [("_core","symox._core"),("_basic","symox._basic"),("_methods",METHODS),("_abc","abc"),("_abc_native","_abc"),("_builtins","builtins"),("_sys","sys"),("_LLL",LLL),("_DDM_MODULE",DDM),("_DENSE_MODULE",DENSE)]{guard!(c.member(HELPER,alias)?.is(c.module(module)?),"helper module namespace changed");}
    let ddm=c.member(DDM,"DDM")?;let domain=c.member("symox._core","Domain")?;let element=c.member("symox._core","DomainElement")?;
    guard!(domain.is(c.py.get_type::<PyDomain>())&&element.is(c.py.get_type::<PyDomainElement>()),"actual native domain classes");
    guard!(c.member(HELPER,"_DDM")?.is(&ddm)&&c.member(HELPER,"_DOMAIN")?.is(&domain)&&c.member(HELPER,"_ELEMENT")?.is(&element),"captured protocol classes");
    guard!(c.state(&ddm)?.meta.is(&native_type)&&same_tuple(&c.state(&ddm)?.bases,&[&list]),"DDM class hierarchy");
    inherited_native(c,&ddm,&list,&["__iter__","__getitem__","__setitem__","__len__","__new__","__getattribute__"])?;
    inherited_native(c,&ddm,&object,&["__setattr__","__delattr__"])?;
    for name in ["shape","rows","cols","domain","__getattr__","__del__","__class__"]{guard!(item(&c.raw(&ddm)?,name)?.is_none(),"DDM foreign attribute");}
    let dictionary_descriptor=need(&c.raw(&ddm)?,"__dict__")?;
    guard!(dictionary_descriptor.is(c.member(HELPER,"_DDM_STATE")?),"DDM dictionary descriptor changed");
    descriptor(c,&dictionary_descriptor,&ddm,PyString::new(c.py,"__dict__").as_any())?;
    for class in [&domain,&element]{
        guard!(same_tuple(&c.state(class)?.bases,&[&object]),"native scalar class hierarchy");
        for name in ["__getattribute__","__getattr__","__class__"]{guard!(item(&c.raw(class)?,name)?.is_none(),"native scalar lookup hook");}
    }
    guard!(item(&c.raw(&domain)?,"_symox_monogenic_extension")?.is_none(),"monogenic domain protocol");
    for name in ["__eq__","__ne__"]{named_descriptor(c,&domain,name)?;}
    // Captured containers are mutable Python bindings. Their shape and full
    // consumed name set must be anchored here, not chosen by those bindings.
    let domain_names=&["of_type","zero","one","is_EXRAW","is_EX","is_ZZ_I","is_QQ_I"][..];
    let element_names=&["__add__","__radd__","__sub__","__rsub__","__mul__","__rmul__",
        "__pow__","__abs__","__float__","__floor__","__lt__","__le__","__gt__","__ge__","__bool__","parent",
        "__iadd__","__isub__","__imul__","__itruediv__","__ipow__"][..];
    for (field,class,names) in [("_DOMAIN_DESCRIPTORS",&domain,domain_names),("_ELEMENT_DESCRIPTORS",&element,element_names)]{
        let records=tuple(c.member(HELPER,field)?)?;
        guard!(records.len()==names.len(),"descriptor record count");
        for (record,expected_name) in records.iter().zip(names){
            let record=tuple(record)?;guard!(record.len()==2,"descriptor record");let name=record.get_item(0)?;
            guard!(name.downcast_exact::<PyString>().is_ok()&&name.eq(*expected_name)?,"descriptor name");let value=record.get_item(1)?;
            let namespace=c.raw(class)?;let live=match namespace.get_item(&name){Ok(v)=>Some(v),Err(error) if error.is_instance_of::<pyo3::exceptions::PyKeyError>(c.py)=>None,Err(error)=>return Err(error.into())};
            if value.is_none(){
                guard!(expected_name.starts_with("__i")&&live.is_none(),"absent inplace descriptor");
            }else{
                guard!(live.is_some_and(|live|live.is(&value)),"native descriptor changed");descriptor(c,&value,class,&name)?;
            }
        }
    }
    rank_error(c)?;
    let zz=c.member(LLL,"ZZ")?;let qq=c.member(LLL,"QQ")?;
    guard!(zz.is_exact_instance_of::<PyDomain>()&&qq.is_exact_instance_of::<PyDomain>(),"native coefficient domains");
    guard!(zz.downcast::<PyDomain>().map_err(PyErr::from)?.get().inner==ExactDomain::ZZ&&qq.downcast::<PyDomain>().map_err(PyErr::from)?.get().inner==ExactDomain::QQ,"coefficient domain values");
    guard!(c.member(HELPER,"_ZZ")?.is(&zz)&&c.member(HELPER,"_QQ")?.is(&qq),"coefficient domain aliases");
    let floor=c.member(LLL,"mfloor")?;guard!(floor.is(c.member(HELPER,"_FLOOR")?)&&c.origin(&floor,"floor","math")?,"floor callback");
    for (field,name) in [("_DUMP","_get_dump"),("_TOKEN","get_cache_token")]{let value=c.member(HELPER,field)?;guard!(value.is(c.member("_abc",name)?)&&c.origin(&value,name,"_abc")?,"ABC getter alias");}
    let import=c.member(HELPER,"_IMPORT")?;guard!(import.is(c.builtin("__import__")?)&&c.origin(&import,"__import__","builtins")?,"import callback");
    // Validate every captured lookup, including the actually resolved module
    // shadow. No hashing before exact string-key validation.
    let lookup_names=["isinstance","getattr","hasattr","type","all","any","dict","set","frozenset","bool","float","complex",
        "str","bytes","bytearray","list","tuple","int","len","range","min","max","enumerate","sum","abs","super"];
    let lookup_modules=[LLL,METHODS,"symox._basic",DDM,DENSE];
    let lookups=tuple(c.member(HELPER,"_LOOKUPS")?)?;
    guard!(lookups.len()==lookup_modules.len()*lookup_names.len(),"lookup record count");
    for ((module,expected_name),record) in lookup_modules.iter().flat_map(|module|lookup_names.iter().map(move |name|(*module,*name))).zip(lookups.iter()){
        let record=tuple(record)?;guard!(record.len()==3,"lookup record");let namespace=dictionary(record.get_item(0)?)?;
        guard!(namespace.is(c.ns(module)?)&&str_keys(&namespace),"lookup namespace keys");let name=record.get_item(1)?;
        guard!(name.downcast_exact::<PyString>().is_ok()&&name.eq(expected_name)?,"lookup name");let value=record.get_item(2)?;
        guard!(value.is(c.builtin(expected_name)?)&&support::_lll_native_owner(&value,&name,c.module("builtins")?.as_any())?,"captured builtin origin");
        let live=namespace.get_item(&name)?.or(c.ns("builtins")?.get_item(&name)?).ok_or(Failure::Decline("lookup missing"))?;
        guard!(live.is(&value),"lookup value changed");
    }
    // All captured function metadata is reread through the public C getters.
    for record in tuple(c.member(HELPER,"_FUNCTIONS")?)?.iter(){
        let record=tuple(record)?;guard!(record.len()==7,"function record");let module=record.get_item(0)?.downcast_into_exact::<PyModule>().map_err(|_|Failure::Decline("function module"))?;
        let namespace=module.dict();guard!(str_keys(&namespace),"function module keys");let name=record.get_item(1)?;guard!(name.downcast_exact::<PyString>().is_ok(),"function name");
        let function=record.get_item(2)?;guard!(namespace.get_item(name)?.is_some_and(|live|live.is(&function)),"function binding changed");
        guard!(c.api.is_function(&function),"function type");
        guard!(c.api.code(&function)?.is(record.get_item(3)?)&&c.api.closure(&function)?.is(record.get_item(4)?)&&c.api.defaults(&function)?.is(record.get_item(5)?)&&c.api.kwdefaults(&function)?.is(record.get_item(6)?),"function metadata changed");
    }
    let constructors=c.member(HELPER,"_CLASS_FUNCTIONS")?.downcast_into_exact::<PyList>().map_err(|_|Failure::Decline("constructor records"))?;
    guard!(constructors.len()==4,"constructor record count");
    for record in constructors.iter(){
        let record=tuple(record)?;guard!(record.len()==7,"constructor record");let name=record.get_item(0)?;
        guard!(name.downcast_exact::<PyString>().is_ok(),"constructor name");let saved=record.get_item(1)?;
        guard!(c.raw(&ddm)?.get_item(&name)?.is(&saved),"constructor descriptor changed");
        let function=if saved.get_type().is(c.builtin("classmethod")?){class_method(c,&saved)?}else{saved};
        guard!(c.api.is_function(&function)&&c.api.globals(&function)?.is(c.ns(DDM)?),"constructor function namespace");
        guard!(c.api.closure(&function)?.is(record.get_item(5)?),"constructor closure changed");
        snapshot_function(c,&function,&record.get_item(2)?,&record.get_item(3)?,&record.get_item(4)?,&record.get_item(6)?)?;
    }
    // Immutable templates also reject callbacks captured before registration.
    let native_elements=c.member(DENSE,"_native_elements")?;
    guard!(c.member(DDM,"_native_elements")?.is(&native_elements),"constructor conversion binding");
    c.function_plain(&native_elements,"_native_elements",DENSE,&c.with_builtins(&["isinstance","list","dict","all","TypeError","AttributeError"],vec![])?)?;
    let ddm_namespace=c.raw(&ddm)?;
    c.function(&need(&ddm_namespace,"__init__")?,"DDM.__init__",DDM,&c.with_builtins(&["isinstance","list","all","type","len","any","super"],vec![("_native_elements",native_elements)])?,&[&ddm],Defaults::None,false)?;
    c.function_plain(&need(&ddm_namespace,"copy")?,"DDM.copy",DDM,&[("DDM",ddm.clone())])?;
    c.function_plain(&class_method(c,&need(&ddm_namespace,"zeros")?)?,"DDM.zeros",DDM,&c.with_builtins(&["range"],vec![("DDM",ddm.clone())])?)?;
    c.function_plain(&class_method(c,&need(&ddm_namespace,"eye")?)?,"DDM.eye",DDM,&c.with_builtins(&["isinstance","tuple","int","range","min"],vec![])?)?;
    let native_call=c.member(METHODS,"_native_call")?;
    let native_convert=c.member(METHODS,"_native_convert")?;
    descriptor(c,&native_call,&domain,PyString::new(c.py,"__call__").as_any())?;
    descriptor(c,&native_convert,&domain,PyString::new(c.py,"convert").as_any())?;
    guard!(native_call.is(c.member(HELPER,"_NATIVE_CALL")?)&&native_convert.is(c.member(HELPER,"_NATIVE_CONVERT")?),"native conversion delegates");
    let call=c.member(METHODS,"_call")?;let convert=c.member(METHODS,"convert")?;let convert_from=c.member(METHODS,"convert_from")?;
    guard!(need(&c.raw(&domain)?,"__call__")?.is(&call)&&need(&c.raw(&domain)?,"convert")?.is(&convert)&&need(&c.raw(&domain)?,"convert_from")?.is(&convert_from),"domain Python delegates");
    let inexact=c.member(METHODS,"_inexact_value")?;let not_coeff=c.member(METHODS,"_not_a_coeff")?;
    c.function_plain(&inexact,"_inexact_value",METHODS,&c.with_builtins(&["type","isinstance","float","complex"],vec![])?)?;
    c.function_plain(&call,"_call",METHODS,&c.with_builtins(&["len","isinstance","bool","list"],vec![("_inexact_value",inexact),("_native_call",native_call)])?)?;
    c.function_plain(&convert_from,"convert_from",METHODS,&c.with_builtins(&["getattr"],vec![])?)?;
    c.function(&convert,"convert",METHODS,&c.with_builtins(&["getattr"],vec![("_native_convert",native_convert),("_not_a_coeff",not_coeff.clone())])?,&[],Defaults::OneNone,false)?;
    c.function_plain(&not_coeff,"_not_a_coeff",METHODS,&c.with_builtins(&["isinstance","float","any"],vec![])?)?;
    c.function_plain(&c.member(METHODS,"_float_from_sympy")?,"_float_from_sympy",METHODS,&c.with_builtins(&["isinstance","getattr"],vec![])?)?;
    let division=need(&c.raw(&element)?,"__truediv__")?;
    guard!(division.is(c.member(HELPER,"_PY_DIV")?)&&c.api.is_function(&division)
        &&c.api.code(&division)?.is(c.member(HELPER,"_PY_DIV_CODE")?)&&c.api.closure(&division)?.is(c.member(HELPER,"_PY_DIV_CLOSURE")?),"division function changed");
    let cells=tuple(c.member(HELPER,"_PY_DIV_CELLS")?)?;guard!(cells.len()==3,"division closure shape");
    let native=cells.get_item(0)?;let operation=cells.get_item(1)?;let reflected=cells.get_item(2)?;
    descriptor(c,&native,&element,PyString::new(c.py,"__truediv__").as_any())?;
    guard!(operation.downcast_exact::<PyString>().is_ok()&&operation.eq("truediv")?&&reflected.downcast_exact::<PyBool>().is_ok()&&!reflected.is_truthy()?,"division closure values");
    c.function(&division,"_DomainElement.__truediv__","symox.polys.domains._algebraic_methods",&[],&[&native,&operation,&reflected],Defaults::None,false)?;
    for name in ["_protocol","_matrix","apply"]{
        let function=c.member(HELPER,name)?;let label=match name{"_protocol"=>"_lll_native._protocol","_matrix"=>"_lll_native._matrix","prefix"=>"_lll_native.prefix","ready"=>"_lll_native.ready",_=>"_lll_native.apply"};
        c.function_plain(&function,label,HELPER,&[])?;
    }
    // A diagnostic study can still use the immutable legacy Python guards.
    // The actual public route binds these two slots directly to native entries.
    for (name,native_name,label) in [("prefix","_lll_prefix_admit","_lll_native.prefix"),("ready","_lll_ready_admit","_lll_native.ready")]{
        let function=c.member(HELPER,name)?;
        if c.api.is_function(&function){c.function_plain(&function,label,HELPER,&[])?;}else{
            let native=c.member("symox._core",native_name)?;
            guard!(function.is(&native)&&c.origin(&native,native_name,"symox._core")?,"native admission entry point");
        }
    }
    Ok(Protocol{ddm,element})
}

fn matrix_inputs<'py>(c:&Context<'py>,p:&Protocol<'py>,x:&Bound<'py,PyAny>)->VResult<()> {
    use super::super::super::super::{scalar,Scalar};
    guard!(x.get_type().is(&p.ddm),"LLL matrix type");
    let rows=x.downcast::<PyList>().map_err(PyErr::from)?;
    let values=c.api.dictionary(x)?.ok_or(Failure::Decline("matrix dictionary"))?;
    let values=dictionary(values)?;guard!(str_keys(&values)&&values.len()==4,"matrix instance attributes");
    for name in ["shape","rows","cols","domain"]{guard!(values.contains(name)?,"matrix shape attributes");}
    let shape=tuple(values.get_item("shape")?.ok_or(Failure::Decline("matrix shape"))?)?;
    guard!(shape.len()==2,"matrix shape tuple");
    let m=shape.get_item(0)?;let n=shape.get_item(1)?;
    guard!(m.downcast_exact::<PyInt>().is_ok()&&n.downcast_exact::<PyInt>().is_ok(),"matrix dimensions");
    let (Ok(m),Ok(n))=(m.extract::<usize>(),n.extract::<usize>())else{return Err(Failure::Decline("matrix dimensions"));};
    guard!(m>0&&m<=n&&rows.len()==m,"matrix dimensions");
    for (name,dimension) in [("rows",m),("cols",n)]{let value=values.get_item(name)?.ok_or(Failure::Decline("matrix dimensions"))?;guard!(value.downcast_exact::<PyInt>().is_ok()&&value.eq(dimension)?,"matrix dimensions");}
    let domain=values.get_item("domain")?.ok_or(Failure::Decline("matrix coefficient domain"))?;
    guard!(domain.is_exact_instance_of::<PyDomain>()&&domain.downcast::<PyDomain>().map_err(PyErr::from)?.get().inner==ExactDomain::ZZ,"matrix coefficient domain");
    for row in rows.iter(){let row=row.downcast_into_exact::<PyList>().map_err(|_|Failure::Decline("matrix row type"))?;guard!(row.len()==n,"matrix row length");for value in row.iter(){guard!(scalar(&value)==Some(Scalar::ZZ),"matrix coefficient type");}}
    Ok(())
}

fn inputs<'py>(c:&Context<'py>,p:&Protocol<'py>,x:&Bound<'py,PyAny>,delta:&Bound<'py,PyAny>,transform:&Bound<'py,PyAny>)->VResult<()> {
    use super::super::super::super::{scalar,Scalar};
    guard!(transform.downcast_exact::<PyBool>().is_ok()&&scalar(delta)==Some(Scalar::QQ),"LLL argument types");
    matrix_inputs(c,p,x)?;
    let value=delta.downcast::<PyDomainElement>().map_err(PyErr::from)?.get();
    let symox::poly::domains::ElementValue::Rational(value)=&value.value.value else{return Err(Failure::Decline("delta scalar"));};
    let quarter=symox::poly::Rat::new(1.into(),4.into());
    guard!(value>&quarter&&value<&symox::poly::Rat::from(symox::poly::Int::one()),"delta range");
    Ok(())
}

fn python_observers<'py>(c:&Context<'py>)->VResult<()> {
    // This covers the Python sys.settrace/setprofile surface, not C tracing
    // installed with a NULL/None state object. It is explicitly a bound, not
    // a claim that gettrace/getprofile=None excludes every native observer.
    for name in ["gettrace","getprofile"]{
        let function=c.member("sys",name)?;guard!(function.downcast_exact::<PyCFunction>().is_ok()&&c.origin(&function,name,"sys")?,"trace getter origin");
        guard!(unsafe{pyo3::ffi::PyCFunction_GetFlags(function.as_ptr())}==pyo3::ffi::METH_NOARGS,"trace getter signature");
        let Some(pointer)=(unsafe{pyo3::ffi::PyCFunction_GetFunction(function.as_ptr())})else{return Err(Failure::Decline("trace getter symbol"));};
        // Direct documented METH_NOARGS signature, avoiding C_CALL profiling
        // events from invoking an otherwise stock getter through vectorcall.
        let value=unsafe{pointer(c.module("sys")?.as_ptr(),std::ptr::null_mut())};
        if value.is_null(){return Err(PyErr::fetch(c.py).into());}
        let value=unsafe{Bound::<PyAny>::from_owned_ptr(c.py,value)};
        guard!(value.is_none(),"active Python trace/profile");
    }
    Ok(())
}

fn checked_prefix<'py>(c:&Context<'py>,x:&Bound<'py,PyAny>,delta:&Bound<'py,PyAny>,transform:&Bound<'py,PyAny>)->VResult<(Protocol<'py>,usize)> {
    python_observers(c)?;
    let protocol=live_protocol(c)?;inputs(c,&protocol,x,delta,transform)?;
    let graph=super::graph(c)?;super::imports(c,&graph)?;
    Ok((protocol,graph.count))
}

fn negative_cache_ready<'py>(c:&Context<'py>,p:&Protocol<'py>)->VResult<()> {
    let basic=c.member("symox._basic","Basic")?;let dump=c.member("_abc","_get_dump")?;
    let state=c.abc(&basic,&dump)?;guard!(state.len()==4,"ABC ready state");
    let version=state.get_item(3)?;guard!(version.downcast_exact::<PyInt>().is_ok(),"ABC ready version");
    let token=c.member("_abc","get_cache_token")?.call0()?;
    guard!(token.downcast_exact::<PyInt>().is_ok()&&version.eq(token)?,"ABC ready token changed");
    let negative=state.get_item(2)?.downcast_into_exact::<PySet>().map_err(|_|Failure::Decline("ABC ready negative set"))?;
    for reference in negative.iter(){if support::_lll_weakref_target(&reference)?.bind(c.py).is(&p.element){return Ok(());}}
    Err(Failure::Decline("first real conversion has not classified the element"))
}

#[pyfunction]
fn _lll_prefix_contract<'py>(x:&Bound<'py,PyAny>,delta:&Bound<'py,PyAny>,transform:&Bound<'py,PyAny>,py:Python<'py>)->PyResult<(bool,Option<&'static str>,usize)>{
    let result=(||{
        let c=Context::full(py)?;let (_,count)=checked_prefix(&c,x,delta,transform)?;
        Ok::<_,Failure>(count)
    })();
    match result{Ok(count)=>Ok((true,None,count)),Err(Failure::Decline(reason))=>Ok((false,Some(reason),0)),Err(Failure::Python(error))=>Err(error)}
}

#[pyfunction]
fn _lll_ready_contract<'py>(x:&Bound<'py,PyAny>,y:&Bound<'py,PyAny>,t:&Bound<'py,PyAny>,delta:&Bound<'py,PyAny>,transform:&Bound<'py,PyAny>,py:Python<'py>)->PyResult<(bool,Option<&'static str>)>{
    let result=(||{
        // Deliberately reread the whole live proof in this initial private
        // implementation. No ABC state or verdict from entry is reused.
        let c=Context::full(py)?;let (protocol,_)=checked_prefix(&c,x,delta,transform)?;
        inputs(&c,&protocol,y,delta,transform)?;
        if transform.is_truthy()?{guard!(!t.is_none(),"missing transform");inputs(&c,&protocol,t,delta,transform)?;}else{guard!(t.is_none(),"unexpected transform");}
        negative_cache_ready(&c,&protocol)?;
        Ok::<_,Failure>(())
    })();
    match result{Ok(())=>Ok((true,None)),Err(Failure::Decline(reason))=>Ok((false,Some(reason))),Err(Failure::Python(error))=>Err(error)}
}
// The public helper signatures are unchanged. Delta/transform are certified
// before the original prefix; the unchanged kernel validates the actual delta
// and transform matrix again. No mutable admission verdict crosses this boundary.
/// This only reports immutable runtime/API availability. It contains no class,
/// function metadata, classifier state or cached admission verdict.
#[pyfunction]
fn _lll_first_supported()->bool { Getters::load().is_some() }

#[pyfunction]
fn _lll_prefix_admit<'py>(x:&Bound<'py,PyAny>,delta:&Bound<'py,PyAny>,transform:&Bound<'py,PyAny>,py:Python<'py>)->PyResult<bool>{
    let result=(||{let c=Context::full(py)?;checked_prefix(&c,x,delta,transform)?;Ok::<_,Failure>(())})();
    match result{Ok(())=>Ok(true),Err(Failure::Decline(_))=>Ok(false),Err(Failure::Python(error))=>Err(error)}
}

#[pyfunction]
fn _lll_ready_admit<'py>(x:&Bound<'py,PyAny>,y:&Bound<'py,PyAny>,t:&Bound<'py,PyAny>,py:Python<'py>)->PyResult<bool>{
    let result=(||{
        let c=Context::full(py)?;python_observers(&c)?;
        let protocol=live_protocol(&c)?;
        matrix_inputs(&c,&protocol,x)?;matrix_inputs(&c,&protocol,y)?;
        if !t.is_none(){matrix_inputs(&c,&protocol,t)?;}
        let graph=super::graph(&c)?;super::imports(&c,&graph)?;
        negative_cache_ready(&c,&protocol)?;
        Ok::<_,Failure>(())
    })();
    match result{Ok(())=>Ok(true),Err(Failure::Decline(_))=>Ok(false),Err(Failure::Python(error))=>Err(error)}
}
pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()>{
    module.add_function(wrap_pyfunction!(_lll_first_supported,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_prefix_admit,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_ready_admit,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_prefix_contract,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_ready_contract,module)?)
}
