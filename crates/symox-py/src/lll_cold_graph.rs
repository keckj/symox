//! Native class/import graph proof consumed by exact LLL admission.
//! No Python verifier, classification, serializer or persisted verdict is used.
use pyo3::prelude::*;
use pyo3::types::{PyBool,PyDict,PyInt,PyList,PyModule,PySet,PyString,PyTuple,PyType};
use std::collections::HashMap;
use super::{getters::Getters,reference,support};

enum Failure { Decline(&'static str), Python(PyErr) }
impl From<PyErr> for Failure {fn from(error:PyErr)->Self{Self::Python(error)}}
type VResult<T> = Result<T,Failure>;
macro_rules! guard {($condition:expr,$reason:literal)=>{if !$condition{return Err(Failure::Decline($reason));}};}
#[path="lll_cold_protocol.rs"]
mod protocol;
fn tuple(value:Bound<'_,PyAny>)->VResult<Bound<'_,PyTuple>>{value.downcast_into_exact::<PyTuple>().map_err(|_|Failure::Decline("tuple type"))}
fn dictionary(value:Bound<'_,PyAny>)->VResult<Bound<'_,PyDict>>{value.downcast_into_exact::<PyDict>().map_err(|_|Failure::Decline("dictionary type"))}
fn str_keys(value:&Bound<'_,PyDict>)->bool{value.iter().all(|(key,_)|key.downcast_exact::<PyString>().is_ok())}
fn item<'py>(mapping:&Bound<'py,PyAny>,name:&str)->VResult<Option<Bound<'py,PyAny>>>{
    match mapping.get_item(name){Ok(value)=>Ok(Some(value)),Err(error) if error.is_instance_of::<pyo3::exceptions::PyKeyError>(mapping.py())=>Ok(None),Err(error)=>Err(error.into())}
}
fn need<'py>(mapping:&Bound<'py,PyAny>,name:&str)->VResult<Bound<'py,PyAny>>{
    item(mapping,name)?.ok_or(Failure::Decline("missing binding"))
}
fn same_tuple(actual:&Bound<'_,PyTuple>,expected:&[&Bound<'_,PyAny>])->bool{
    actual.len()==expected.len() && actual.iter().zip(expected).all(|(a,b)|a.is(*b))
}

const MODULES:&[&str]=&["builtins","abc","_abc","weakref","symox","symox._core","symox._basic","symox._compat","symox._node_types","symox.logic.boolalg","symox.core.operations",
    "sys","_imp","_thread","_weakref","math","importlib","importlib._bootstrap","symox.tensor.tensor","symox.polys.domains._methods",
    "symox.polys.matrices._lll_native","symox.polys.matrices.lll","symox.polys.matrices.ddm","symox.polys.matrices._dense_poly","mpmath",
    "symox.polys.domains._algebraic_methods","symox.polys.matrices.exceptions"];

struct ClassState<'py>{namespace:Bound<'py,PyAny>,mro:Bound<'py,PyTuple>,bases:Bound<'py,PyTuple>,meta:Bound<'py,PyType>}
struct Context<'py>{py:Python<'py>,api:Getters,modules:HashMap<&'static str,(Bound<'py,PyModule>,Bound<'py,PyDict>)>}
#[derive(Clone,Copy)]
enum Defaults { None, OneNone, NoneZero }
impl<'py> Context<'py>{
    fn new(py:Python<'py>,imports:bool)->VResult<Self>{
        Self::selected(py,if imports{&MODULES[..20]}else{&MODULES[..11]})
    }
    fn full(py:Python<'py>)->VResult<Self>{Self::selected(py,MODULES)}
    fn selected(py:Python<'py>,names:&[&'static str])->VResult<Self>{
        let api=Getters::load().ok_or(Failure::Decline("unsupported optional public API"))?;
        let args=PyTuple::new(py,names)?;
        let values=support::_lll_loaded_modules(args.as_any(),py)?.ok_or(Failure::Decline("loaded module namespaces"))?;
        let mut modules=HashMap::new();
        for (name,value) in names.iter().zip(tuple(values.into_bound(py))?.iter()){
            let pair=tuple(value)?;
            let module=pair.get_item(0)?.downcast_into_exact::<PyModule>().map_err(|_|Failure::Decline("module type"))?;
            let namespace=dictionary(pair.get_item(1)?)?;
            modules.insert(*name,(module,namespace));
        }
        Ok(Self{py,api,modules})
    }
    fn module(&self,name:&str)->VResult<Bound<'py,PyModule>>{Ok(self.modules.get(name).ok_or(Failure::Decline("unloaded module"))?.0.clone())}
    fn ns(&self,name:&str)->VResult<Bound<'py,PyDict>>{Ok(self.modules.get(name).ok_or(Failure::Decline("unloaded module"))?.1.clone())}
    fn member(&self,module:&str,name:&str)->VResult<Bound<'py,PyAny>>{
        self.ns(module)?.get_item(name)?.ok_or(Failure::Decline("missing module binding"))
    }
    fn builtin(&self,name:&str)->VResult<Bound<'py,PyAny>>{self.member("builtins",name)}
    fn origin(&self,value:&Bound<'py,PyAny>,name:&str,module:&str)->VResult<bool>{
        Ok(support::_lll_native_owner(value,PyString::new(self.py,name).as_any(),self.module(module)?.as_any())?)
    }
    fn builtins(&self,names:&[&str])->VResult<()> {
        for name in names {guard!(self.origin(&self.builtin(name)?,name,"builtins")?,"builtin origin");}
        Ok(())
    }
    fn state(&self,class:&Bound<'py,PyAny>)->VResult<ClassState<'py>>{
        let (namespace,mro,bases,meta)=support::class_state(class)?.ok_or(Failure::Decline("class namespace"))?;
        Ok(ClassState{namespace,mro:tuple(mro)?,bases:tuple(bases)?,meta})
    }
    fn raw(&self,class:&Bound<'py,PyAny>)->VResult<Bound<'py,PyAny>>{Ok(self.state(class)?.namespace)}
    fn effective(&self,class:&Bound<'py,PyAny>,name:&str)->VResult<Bound<'py,PyAny>>{
        for parent in self.state(class)?.mro.iter(){if let Some(value)=item(&self.raw(&parent)?,name)?{return Ok(value);}}
        Err(Failure::Decline("missing effective attribute"))
    }
    fn function(&self,function:&Bound<'py,PyAny>,label:&str,module:&str,bindings:&[(&str,Bound<'py,PyAny>)],closure:&[&Bound<'py,PyAny>],defaults:Defaults,recursive:bool)->VResult<()> {
        guard!(reference::function_code_matches(function,label,&self.api)?,"canonical function code");
        let namespace=self.api.globals(function)?;
        guard!(namespace.is(self.ns(module)?.as_any()) && function.getattr("__builtins__")?.is(self.ns("builtins")?.as_any()),"function namespaces");
        let actual=self.api.defaults(function)?;
        match defaults {
            Defaults::None=>guard!(actual.is_none(),"function defaults"),
            Defaults::OneNone=>{let actual=tuple(actual)?;guard!(actual.len()==1&&actual.get_item(0)?.is_none(),"function defaults");},
            Defaults::NoneZero=>{let actual=tuple(actual)?;guard!(actual.len()==2&&actual.get_item(0)?.is_none(),"function defaults");let zero=actual.get_item(1)?;guard!(zero.downcast_exact::<PyInt>().is_ok()&&zero.eq(0i32)?,"function defaults");},
        }
        let kwargs=self.api.kwdefaults(function)?;
        if recursive {
            let kwargs=dictionary(kwargs)?;guard!(str_keys(&kwargs)&&kwargs.len()==1,"function keyword defaults");
            let value=kwargs.get_item("recursive")?.ok_or(Failure::Decline("function keyword defaults"))?;
            guard!(value.downcast_exact::<PyBool>().is_ok()&&!value.is_truthy()?,"function keyword defaults");
        } else {guard!(kwargs.is_none(),"function keyword defaults");}
        let cells=self.api.closure(function)?;
        if closure.is_empty(){guard!(cells.is_none(),"function closure");}else{
            let cells=tuple(cells)?;guard!(cells.len()==closure.len(),"function closure size");
            for (cell,expected) in cells.iter().zip(closure){
                let value=match cell.getattr("cell_contents"){Ok(value)=>value,Err(error) if error.is_instance_of::<pyo3::exceptions::PyValueError>(self.py)=>return Err(Failure::Decline("empty closure")),Err(error)=>return Err(error.into())};
                guard!(value.is(*expected),"function closure binding");
            }
        }
        let namespace=dictionary(namespace)?;
        for (name,expected) in bindings {
            let value=namespace.get_item(*name)?.or(self.ns("builtins")?.get_item(*name)?).ok_or(Failure::Decline("missing global"))?;
            guard!(value.is(expected),"function global binding");
        }
        Ok(())
    }
    fn function_plain(&self,function:&Bound<'py,PyAny>,label:&str,module:&str,bindings:&[(&str,Bound<'py,PyAny>)])->VResult<()> {
        self.function(function,label,module,bindings,&[],Defaults::None,false)
    }
    fn with_builtins(&self,names:&[&'static str],extra:Vec<(&'static str,Bound<'py,PyAny>)>)->VResult<Vec<(&'static str,Bound<'py,PyAny>)>>{
        let mut result=extra;for name in names{result.push((*name,self.builtin(name)?));}Ok(result)
    }
    fn abc(&self,class:&Bound<'py,PyAny>,dump:&Bound<'py,PyAny>)->VResult<Bound<'py,PyTuple>>{
        tuple(support::_lll_abc_state(class,dump,self.py)?.ok_or(Failure::Decline("opaque ABC state"))?.into_bound(self.py))
    }
}

struct Graph<'py>{count:usize,paths:Vec<Bound<'py,PyString>>,reference_type:Bound<'py,PyType>,basic:Bound<'py,PyAny>,element:Bound<'py,PyAny>,resolve:Bound<'py,PyAny>}
fn graph<'py>(c:&Context<'py>)->VResult<Graph<'py>>{
    c.builtins(&["type","object","dict","tuple","list","str","bool","int","isinstance","issubclass","getattr","any","super","hash","KeyError","TypeError"])?;
    let basic=c.member("symox._basic","Basic")?;let element=c.member("symox._core","DomainElement")?;
    let basic_meta=c.member("symox._basic","_BasicMeta")?;let abc_meta=c.member("abc","ABCMeta")?;
    let expr_meta=c.member("symox._compat","_ExprType")?;let bool_meta=c.member("symox.logic.boolalg","_BooleanType")?;
    let lattice_meta=c.member("symox.logic.boolalg","_BooleanLatticeType")?;let category_meta=c.member("symox.core.operations","_OperationCategory")?;
    let node_meta=c.member("symox._node_types","NodeType")?;let native_type=c.py.get_type::<PyType>().into_any();let object=c.builtin("object")?;
    guard!(c.state(&basic)?.meta.is(&basic_meta)&&c.state(&element)?.meta.is(&native_type),"classifier input metaclasses");
    guard!(same_tuple(&c.state(&element)?.mro,&[&element,&object]),"native element MRO");
    for name in ["__getattribute__","__getattr__","__class__"]{guard!(item(&c.raw(&element)?,name)?.is_none(),"native element lookup");}
    let module=need(&c.raw(&element)?,"__module__")?;guard!(module.downcast_exact::<PyString>().is_ok()&&module.eq("symox._core")?,"native element module");
    let matrix=c.member("symox._core","Matrix")?;guard!(c.state(&matrix)?.meta.is(&native_type)&&!matrix.is(&element)&&!matrix.is(&object),"matrix branch");
    let type_dict=c.raw(&native_type)?;
    let name_descriptor=need(&type_dict,"__name__")?;
    let element_name=name_descriptor.call_method1("__get__",(&element,))?;
    guard!(element_name.downcast_exact::<PyString>().is_ok()&&element_name.eq("DomainElement")?,"native element name");
    guard!(c.state(&node_meta)?.meta.is(&native_type)&&same_tuple(&c.state(&node_meta)?.bases,&[&native_type]),"nominal metaclass");
    let metas=[(&basic_meta,vec![&abc_meta]),(&abc_meta,vec![&native_type]),(&expr_meta,vec![&basic_meta]),(&category_meta,vec![&expr_meta]),(&bool_meta,vec![&basic_meta]),(&lattice_meta,vec![&bool_meta,&category_meta])];
    for (meta,bases) in &metas{
        let state=c.state(meta)?;guard!(state.meta.is(&native_type)&&same_tuple(&state.bases,bases),"metaclass bases");
        for name in ["__getattribute__","__dict__","__name__","__mro__","__subclasses__"]{guard!(c.effective(meta,name)?.is(need(&type_dict,name)?),"metaclass native lookup");}
        guard!(c.effective(meta,"__subclasshook__")?.is(need(&c.raw(&object)?,"__subclasshook__")?),"metaclass subclass hook");
        for parent in state.mro.iter(){guard!(item(&c.raw(&parent)?,"_abc_impl")?.is_none(),"metaclass ABC data descriptor");}
    }
    let mut delegates=Vec::new();for name in ["_abc_instancecheck","_abc_subclasscheck"]{
        let value=c.member("abc",name)?;guard!(value.is(c.member("_abc",name)?)&&c.origin(&value,name,"_abc")?,"ABC native delegate");delegates.push(value);
    }
    let instance=need(&c.raw(&basic_meta)?,"__instancecheck__")?;
    c.function(&instance,"_BasicMeta.__instancecheck__","symox._basic",&c.with_builtins(&["isinstance","super","type"],vec![("Basic",basic.clone()),("_core",c.module("symox._core")?.into_any())])?,&[&basic_meta],Defaults::None,false)?;
    c.function_plain(&need(&c.raw(&abc_meta)?,"__instancecheck__")?,"ABCMeta.__instancecheck__","abc",&[("_abc_instancecheck",delegates[0].clone())])?;
    let abc_check=need(&c.raw(&abc_meta)?,"__subclasscheck__")?;
    c.function_plain(&abc_check,"ABCMeta.__subclasscheck__","abc",&[("_abc_subclasscheck",delegates[1].clone())])?;
    let nominal_subclass=c.member("symox._compat","_nominal_subclass")?;
    let expr_check=need(&c.raw(&expr_meta)?,"__subclasscheck__")?;
    c.function_plain(&expr_check,"_ExprType.__subclasscheck__","symox._compat",&c.with_builtins(&["isinstance","type"],vec![("_nominal_subclass",nominal_subclass.clone())])?)?;
    let resolve=c.member("symox.core.operations","_resolve")?;let category_check=need(&c.raw(&category_meta)?,"__subclasscheck__")?;
    c.function_plain(&category_check,"_OperationCategory.__subclasscheck__","symox.core.operations",&c.with_builtins(&["any","isinstance","type"],vec![("_ExprType",expr_meta.clone()),("_resolve",resolve.clone())])?)?;
    let bool_check=need(&c.raw(&bool_meta)?,"__subclasscheck__")?;
    c.function(&bool_check,"_BooleanType.__subclasscheck__","symox.logic.boolalg",&c.with_builtins(&["isinstance","issubclass","super"],vec![])?,&[&bool_meta],Defaults::None,false)?;
    guard!(item(&c.raw(&basic_meta)?,"__subclasscheck__")?.is_none()&&item(&c.raw(&lattice_meta)?,"__subclasscheck__")?.is_none(),"hidden subclass checker");
    c.function_plain(&nominal_subclass,"_nominal_subclass","symox._compat",&c.with_builtins(&["isinstance","issubclass"],vec![("_node_types",c.module("symox._node_types")?.into_any())])?)?;
    let native_module=c.member("symox._node_types","native_module")?;let is_symox=c.member("symox._node_types","is_symox_class")?;
    let anchor=c.member("symox._node_types","nominal_anchor")?;let node_class=c.member("symox._node_types","node_class")?;
    let anchors=c.member("symox._node_types","_anchors")?;let classes=c.member("symox._node_types","_classes")?;let class_table=c.member("symox._node_types","CLASSES")?;
    c.function_plain(&native_module,"native_module","symox._node_types",&c.with_builtins(&["type"],vec![])?)?;
    c.function_plain(&is_symox,"is_symox_class","symox._node_types",&c.with_builtins(&["getattr","isinstance","str","type"],vec![("native_module",native_module)])?)?;
    c.function_plain(&anchor,"nominal_anchor","symox._node_types",&c.with_builtins(&["KeyError","TypeError"],vec![("_anchors",anchors.clone()),("is_symox_class",is_symox),("node_class",node_class.clone())])?)?;
    c.function_plain(&node_class,"node_class","symox._node_types",&[("_classes",classes.clone()),("CLASSES",class_table.clone())])?;
    let classes=dictionary(classes)?;let class_table=dictionary(class_table)?;guard!(str_keys(&classes)&&str_keys(&class_table),"nominal dictionary keys");
    guard!(classes.get_item(&element_name)?.is_none_or(|value|value.is_none())&&class_table.get_item(&element_name)?.is_none(),"element nominal class");
    let dump=c.member("_abc","_get_dump")?;guard!(c.origin(&dump,"_get_dump","_abc")?,"ABC dump origin");
    let root_dump=c.abc(&basic,&dump)?;let registry=root_dump.get_item(0)?.downcast_into_exact::<PySet>().map_err(|_|Failure::Decline("ABC registry set"))?;
    let witness=registry.iter().next().ok_or(Failure::Decline("native weakref witness"))?;let reference_type=witness.get_type();
    let weak_dictionary=c.member("weakref","WeakKeyDictionary")?;guard!(anchors.get_type().is(&weak_dictionary),"nominal weak dictionary class");
    let (anchor_state,generic)=support::_lll_instance_state(&anchors)?.ok_or(Failure::Decline("weak dictionary instance"))?;guard!(generic,"weak dictionary generic lookup");
    let anchor_state=dictionary(anchor_state.into_bound(c.py))?;guard!(str_keys(&anchor_state),"weak dictionary instance keys");
    for name in ["data","_remove"]{guard!(anchor_state.contains(name)?,"weak dictionary fields");for base in c.state(&weak_dictionary)?.mro.iter(){guard!(item(&c.raw(&base)?,name)?.is_none(),"weak dictionary field descriptor");}}
    for (name,label) in [("__getitem__","WeakKeyDictionary.__getitem__"),("__setitem__","WeakKeyDictionary.__setitem__")]{c.function_plain(&need(&c.raw(&weak_dictionary)?,name)?,label,"weakref",&[("ref",reference_type.clone().into_any())])?;}
    let anchor_data=dictionary(anchor_state.get_item("data")?.ok_or(Failure::Decline("weak dictionary data"))?)?;let query_hash=element.hash()?;
    for (reference,value) in anchor_data.iter(){
        guard!(reference.get_type().is(&reference_type),"nominal reference type");
        let target=support::_lll_weakref_target(&reference)?.into_bound(c.py);
        if target.is(&element){guard!(value.is_none(),"existing element nominal anchor");}else{guard!(reference.hash()?!=query_hash,"nominal reference collision");}
    }
    let native_check=need(&type_dict,"__subclasscheck__")?;let default_hook=need(&c.raw(&object)?,"__subclasshook__")?;
    let allowed_checks=[&native_check,&abc_check,&expr_check,&category_check,&bool_check];
    let mut todo=vec![(basic.clone(),false)];let mut seen:Vec<Bound<'py,PyAny>>=Vec::new();let mut active:Vec<Bound<'py,PyAny>>=Vec::new();let mut paths=Vec::new();
    while let Some((class,leaving))=todo.pop(){
        if leaving{guard!(active.pop().is_some_and(|value|value.is(&class)),"graph traversal");continue;}
        guard!(!active.iter().any(|ancestor|ancestor.is(&class)),"cyclic ABC classifier graph");
        if seen.iter().any(|old|old.is(&class)){continue;}guard!(seen.len()<128,"class graph budget");seen.push(class.clone());
        active.push(class.clone());todo.push((class.clone(),true));
        let state=c.state(&class)?;
        for parent in state.mro.iter(){let namespace=c.raw(&parent)?;for name in ["__subclasses__","__subclasscheck__"]{guard!(item(&namespace,name)?.is_none(),"class ABC dispatch override");}}
        guard!(state.meta.is(&native_type)||metas.iter().any(|(meta,_)|state.meta.is(*meta)),"custom graph metaclass");
        guard!(!class.is(&element)&&!class.is(&object),"positive native subtype");guard!(c.effective(&class,"__subclasshook__")?.is(&default_hook),"custom subclass hook");
        let checker=c.effective(state.meta.as_any(),"__subclasscheck__")?;guard!(allowed_checks.iter().any(|expected|checker.is(*expected)),"custom subclass checker");
        if checker.is(&native_check){continue;}
        if checker.is(&expr_check){guard!(item(&state.namespace,"_symox_expr_facade")?.is_none(),"expression facade");guard!(need(&state.namespace,"__module__")?.downcast_exact::<PyString>().is_ok(),"expression module");continue;}
        if checker.is(&category_check){
            if let Some(category)=item(&state.namespace,"_category")?{if !category.is_none(){let category=tuple(category)?;guard!(category.len()==2,"category shape");for path in tuple(category.get_item(1)?)?.iter(){paths.push(path.downcast_into_exact::<PyString>().map_err(|_|Failure::Decline("category path"))?);}}}
            continue;
        }
        if checker.is(&bool_check)&&state.meta.is(&lattice_meta){guard!(item(&state.namespace,"_category")?.is_none_or(|v|v.is_none()),"lattice category");}
        let abc=c.abc(&class,&dump)?;guard!(abc.len()==4&&abc.get_item(3)?.downcast_exact::<PyInt>().is_ok(),"ABC state shape");
        for index in 0..3 {
            let group=abc.get_item(index)?.downcast_into_exact::<PySet>().map_err(|_|Failure::Decline("ABC set type"))?;
            for reference in group.iter(){
                guard!(reference.get_type().is(&reference_type),"ABC reference type");let target=support::_lll_weakref_target(&reference)?.into_bound(c.py);
                if target.is(&element){guard!(index==2,"ABC positive element");}else{guard!(reference.hash()?!=query_hash,"ABC reference collision");}
                if index==0&&!target.is_none(){todo.push((target,false));}
            }
        }
        let children=support::_lll_type_subclasses(&class,c.py)?.ok_or(Failure::Decline("native subclasses"))?.into_bound(c.py).downcast_into_exact::<PyList>().map_err(|_|Failure::Decline("native subclasses list"))?;
        todo.extend(children.iter().map(|child|(child,false)));
    }
    Ok(Graph{count:seen.len(),paths,reference_type,basic,element,resolve})
}

fn plain_class<'py>(c:&Context<'py>,class:&Bound<'py,PyAny>,fields:&[&str])->VResult<()> {
    let object=c.builtin("object")?;let native_type=c.py.get_type::<PyType>();let state=c.state(class)?;
    guard!(state.meta.is(&native_type)&&same_tuple(&state.bases,&[&object]),"import helper class bases");
    let object_dict=c.raw(&object)?;
    for name in ["__getattribute__","__setattr__","__delattr__","__new__"]{guard!(c.effective(class,name)?.is(need(&object_dict,name)?),"import helper native object protocol");}
    for base in state.mro.iter(){
        let namespace=c.raw(&base)?;guard!(item(&namespace,"__getattr__")?.is_none()&&item(&namespace,"__del__")?.is_none(),"import helper custom hook");
        for name in fields{guard!(item(&namespace,name)?.is_none(),"import helper field descriptor");}
    }
    Ok(())
}

const CATEGORIES:&[(&str,&str,&str)]=&[("symox:Add","symox","Add"),("symox:Mul","symox","Mul"),("symox:MatAdd","symox","MatAdd"),("symox:MatMul","symox","MatMul"),("symox:Union","symox","Union"),("symox:Intersection","symox","Intersection"),("symox.tensor.tensor:TensAdd","symox.tensor.tensor","TensAdd"),("symox.tensor.tensor:TensMul","symox.tensor.tensor","TensMul")];
fn imports<'py>(c:&Context<'py>,g:&Graph<'py>)->VResult<()> {
    c.builtins(&["__import__","hasattr","set","RuntimeError","property"])?;
    // IMPORT_NAME/IMPORT_FROM and importlib use hasattr even for already
    // loaded modules (notably __path__). A module-level fallback can execute
    // arbitrary user code despite an exact ModuleType and a complete registry.
    for name in ["symox","symox._node_types","importlib","math","mpmath"]{
        if let Some((_,namespace))=c.modules.get(name){
            guard!(namespace.get_item("__getattr__")?.is_none(),"module attribute fallback");
        }
    }
    for name in ["symox","importlib"]{c.member(name,"__path__")?;}
    if c.modules.contains_key("mpmath"){
        c.member("mpmath","__path__")?;c.member("mpmath","libmp")?;
    }
    let registry=support::_lll_module_registry(c.py)?.ok_or(Failure::Decline("module registry"))?.into_bound(c.py);
    guard!(c.member("sys","modules")?.is(&registry),"authoritative module registry");
    let bootstrap="importlib._bootstrap";
    let locks=dictionary(c.member(bootstrap,"_module_locks")?)?;
    let blocking=dictionary(c.member(bootstrap,"_blocking_on")?)?;
    guard!(locks.is_empty()&&blocking.is_empty(),"import lock state is not empty");
    let needs_loading=c.member(bootstrap,"_NEEDS_LOADING")?;guard!(needs_loading.get_type().is(c.builtin("object")?),"module-missing sentinel");
    for name in ["acquire_lock","release_lock"]{guard!(c.origin(&c.member("_imp",name)?,name,"_imp")?,"native import lock function");}
    for name in ["allocate_lock","get_ident"]{guard!(c.origin(&c.member("_thread",name)?,name,"_thread")?,"native thread function");}
    guard!(c.member("_weakref","ref")?.is(&g.reference_type),"native import lock weakref factory");
    for module in ["_imp","_thread","_weakref"]{guard!(c.member(bootstrap,module)?.is(c.module(module)?),"import lock native modules");}
    let lock=c.member(bootstrap,"_ModuleLock")?;let manager=c.member(bootstrap,"_ModuleLockManager")?;let spec=c.member(bootstrap,"ModuleSpec")?;
    plain_class(c,&lock,&["lock","wakeup","name","owner","count","waiters"])?;
    plain_class(c,&manager,&["_name","_lock"])?;
    plain_class(c,&spec,&["name","submodule_search_locations","_initializing"])?;
    let get_lock=c.member(bootstrap,"_get_module_lock")?;let unlock=c.member(bootstrap,"_lock_unlock_module")?;let deadlock=c.member(bootstrap,"_DeadlockError")?;
    let state=c.state(&deadlock)?;guard!(state.meta.is(c.py.get_type::<PyType>())&&same_tuple(&state.bases,&[&c.builtin("RuntimeError")?]),"import deadlock exception bases");
    for key in state.namespace.try_iter()?{
        let key=key?;let mut allowed=false;for name in ["__module__","__weakref__","__doc__"]{if key.eq(name)?{allowed=true;break;}}
        guard!(allowed,"import deadlock exception protocol");
    }
    c.function_plain(&get_lock,"_frozen_importlib._get_module_lock",bootstrap,&c.with_builtins(&["KeyError"],vec![("_imp",c.module("_imp")?.into_any()),("_module_locks",locks.clone().into_any()),("_thread",c.module("_thread")?.into_any()),("_ModuleLock",lock.clone()),("_weakref",c.module("_weakref")?.into_any())])?)?;
    c.function_plain(&unlock,"_frozen_importlib._lock_unlock_module",bootstrap,&[("_get_module_lock",get_lock.clone()),("_DeadlockError",deadlock.clone())])?;
    let lock_ns=c.raw(&lock)?;let manager_ns=c.raw(&manager)?;
    c.function_plain(&need(&lock_ns,"__init__")?,"_frozen_importlib._ModuleLock.__init__",bootstrap,&[("_thread",c.module("_thread")?.into_any())])?;
    c.function_plain(&need(&lock_ns,"acquire")?,"_frozen_importlib._ModuleLock.acquire",bootstrap,&[("_thread",c.module("_thread")?.into_any()),("_blocking_on",blocking.clone().into_any()),("_DeadlockError",deadlock)])?;
    c.function_plain(&need(&lock_ns,"release")?,"_frozen_importlib._ModuleLock.release",bootstrap,&c.with_builtins(&["RuntimeError"],vec![("_thread",c.module("_thread")?.into_any())])?)?;
    c.function_plain(&need(&lock_ns,"has_deadlock")?,"_frozen_importlib._ModuleLock.has_deadlock",bootstrap,&c.with_builtins(&["set"],vec![("_thread",c.module("_thread")?.into_any()),("_blocking_on",blocking.into_any())])?)?;
    c.function_plain(&need(&manager_ns,"__init__")?,"_frozen_importlib._ModuleLockManager.__init__",bootstrap,&[])?;
    c.function_plain(&need(&manager_ns,"__enter__")?,"_frozen_importlib._ModuleLockManager.__enter__",bootstrap,&[("_get_module_lock",get_lock)])?;
    c.function_plain(&need(&manager_ns,"__exit__")?,"_frozen_importlib._ModuleLockManager.__exit__",bootstrap,&[])?;
    let find=c.member(bootstrap,"_find_and_load")?;let gcd=c.member(bootstrap,"_gcd_import")?;let sanity=c.member(bootstrap,"_sanity_check")?;
    c.function_plain(&sanity,"_frozen_importlib._sanity_check",bootstrap,&c.with_builtins(&["isinstance","str"],vec![])?)?;
    c.function_plain(&find,"_frozen_importlib._find_and_load",bootstrap,&[("_ModuleLockManager",manager),("sys",c.module("sys")?.into_any()),("_NEEDS_LOADING",needs_loading),("_lock_unlock_module",unlock)])?;
    c.function(&gcd,"_frozen_importlib._gcd_import",bootstrap,&[("_sanity_check",sanity),("_find_and_load",find),("_gcd_import",gcd.clone())],&[],Defaults::NoneZero,false)?;
    c.function(&c.member(bootstrap,"_handle_fromlist")?,"_frozen_importlib._handle_fromlist",bootstrap,&c.with_builtins(&["isinstance","str","hasattr"],vec![])?,&[],Defaults::None,true)?;
    c.function(&c.member("importlib","import_module")?,"import_module","importlib",&[("_bootstrap",c.module(bootstrap)?.into_any())],&[],Defaults::OneNone,false)?;
    c.function_plain(&g.resolve,"_resolve","symox.core.operations",&c.with_builtins(&["getattr","tuple"],vec![])?)?;
    let parent=need(&c.raw(&spec)?,"parent")?;guard!(parent.get_type().is(c.builtin("property")?),"spec parent descriptor");
    c.function_plain(&parent.getattr("fget")?,"_frozen_importlib.ModuleSpec.parent",bootstrap,&[])?;
    for module_name in ["symox","symox._node_types","symox.polys.domains._methods","symox.tensor.tensor","importlib","math","mpmath"]{
        // mpmath belongs to the complete constructor/rounding proof; the
        // independent classifier-only sub-proof does not consume it.
        if module_name=="mpmath"&&!c.modules.contains_key(module_name){continue;}
        let namespace=c.ns(module_name)?;
        let name=namespace.get_item("__name__")?.ok_or(Failure::Decline("module name"))?;
        let package=namespace.get_item("__package__")?.ok_or(Failure::Decline("module package"))?;
        guard!(name.downcast_exact::<PyString>().is_ok()&&name.eq(module_name)?&&package.downcast_exact::<PyString>().is_ok(),"import package strings");
        let value=namespace.get_item("__spec__")?.ok_or(Failure::Decline("module spec"))?;
        guard!(value.get_type().is(&spec),"module spec type");
        let (state,generic)=support::_lll_instance_state(&value)?.ok_or(Failure::Decline("module spec dictionary"))?;
        guard!(generic,"module spec generic lookup");let state=dictionary(state.into_bound(c.py))?;guard!(str_keys(&state),"module spec keys");
        if let Some(initializing)=state.get_item("_initializing")?{guard!(initializing.downcast_exact::<PyBool>().is_ok()&&!initializing.is_truthy()?,"module still initializing");}
        let spec_name=state.get_item("name")?.ok_or(Failure::Decline("module spec name"))?;
        guard!(spec_name.downcast_exact::<PyString>().is_ok()&&spec_name.eq(module_name)?,"module spec name");
        let search=state.get_item("submodule_search_locations")?.ok_or(Failure::Decline("module spec package field"))?;
        let expected_parent=if search.is_none(){module_name.rsplit_once('.').map_or("",|(parent,_)|parent)}else{module_name};
        guard!(package.eq(expected_parent)?,"relative import package mismatch could warn");
    }
    guard!(c.member("symox","Basic")?.is(&g.basic),"relative Basic binding");
    for name in ["S","CoercionFailed","Expr","Set"]{c.member("symox",name)?;}
    for path in &g.paths{
        let mut selected=None;
        for (expected,module,name) in CATEGORIES{if path.as_any().eq(*expected)?{selected=Some((*module,*name));break;}}
        let (module,name)=selected.ok_or(Failure::Decline("unsupported category path"))?;
        let target=c.member(module,name)?;c.state(&target)?;
        guard!(!target.is(&g.element)&&!target.is(c.builtin("object")?),"positive category subtype");
    }
    Ok(())
}

#[pyfunction]
fn _lll_classifier_contract(py:Python<'_>,imports:bool)->PyResult<(bool,Option<&'static str>,usize)> {
    let result=(||{let context=Context::new(py,imports)?;let value=graph(&context)?;if imports{self::imports(&context,&value)?;}Ok(value)})();
    match result{Ok(value)=>Ok((true,None,value.count)),Err(Failure::Decline(reason))=>Ok((false,Some(reason),0)),Err(Failure::Python(error))=>Err(error)}
}
pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()>{protocol::register(module)?;module.add_function(wrap_pyfunction!(_lll_classifier_contract,module)?)}
