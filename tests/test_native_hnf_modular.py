"""Modulo-D HNF: exact recurrence, whole lattices, aliases and full fallback."""
import ast
import copy
import importlib
import json
import types
from fractions import Fraction
from pathlib import Path
import pytest
import symox as sx
import sympy as sp

nf=importlib.import_module('symox.polys.matrices.normalforms')
public=importlib.import_module('symox.matrices.normalforms')
ref=importlib.import_module('sympy.polys.matrices.normalforms')
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices.ddm import DDM
from sympy.polys.matrices.domainmatrix import DomainMatrix as ReferenceMatrix
CASES=[json.loads(line) for line in (Path(__file__).parent/'data/native_hnf_modular.jsonl').read_text().splitlines()]


def matrix(rows,shape):return DomainMatrix.from_rep(DDM([[sx.ZZ(v) for v in row] for row in rows],shape,sx.ZZ))
def values(a):return [[int(v) for v in row] for row in a]


def original(namespace=None):
    tree=ast.parse(Path(nf.__file__).read_text())
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_hermite_normal_form_modulo_D')
    node.body=[n for n in node.body if not (isinstance(n,ast.If) and isinstance(n.test,ast.Call)
               and isinstance(n.test.func,ast.Name) and n.test.func.id=='_native_hnf_modular')]
    code=compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),nf.__file__,'exec')
    return types.FunctionType(next(c for c in code.co_consts if type(c)is types.CodeType),nf.__dict__ if namespace is None else namespace)


def multiply(a,b,n):return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(len(a))]
def determinant(a):
    a=[list(map(Fraction,row)) for row in a];d=Fraction(1)
    for j in range(len(a)):
        p=next((i for i in range(j,len(a)) if a[i][j]),None)
        if p is None:return 0
        if p!=j:a[p],a[j]=a[j],a[p];d=-d
        v=a[j][j];d*=v
        for i in range(j+1,len(a)):
            q=a[i][j]/v
            for k in range(j+1,len(a)):a[i][k]-=q*a[j][k]
    return d


def lattice_certificate(source,output,shape):
    m,n=shape
    work=DDM([[sx.ZZ(v) for v in row] for row in source],shape,sx.ZZ)
    k,operations=sx._core._ddm_hnf(work,shape,DDM,True)
    transform=[[int(i==j) for j in range(n)] for i in range(n)]
    for i,j,a,b,c,d in operations:
        a,b,c,d=map(int,(a,b,c,d));element=[[int(p==q) for q in range(n)] for p in range(n)]
        if i==j:
            assert (a,b,c,d)==(-1,0,-1,0);element[i][i]=-1
        else:
            assert a*d-b*c==1
            element[i][i],element[j][i],element[i][j],element[j][j]=a,b,c,d
        transform=multiply(transform,element,n)
    assert abs(determinant(transform))==1
    full=values(work)
    assert multiply(source,transform,n)==full
    assert all(v==0 for row in full for v in row[:k])
    assert n-k==m and [row[k:] for row in full]==output
    for i in range(m):
        assert output[i][i]>0
        assert all(output[j][i]==0 for j in range(i+1,m))
        assert all(0<=output[i][j]<output[i][i] for j in range(i+1,m))


@pytest.mark.parametrize('case',CASES,ids=lambda c:str(c['shape'])+'-D'+str(c['D']))
@pytest.mark.parametrize('fmt',['dense','sparse'])
def test_all_frozen_exact_cases_public_and_complete_lattice_certificates(case,fmt,monkeypatch):
    calls=[];native=sx._core._ddm_hnf_modular
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_hnf_modular',record)
    shape=tuple(case['shape']);a=matrix(case['input'],shape)
    if fmt=='sparse':a=a.to_sparse()
    before=a.to_list();out=nf.hermite_normal_form(a,D=sx.ZZ(case['D']))
    assert calls==[True] and out.shape==(shape[0],shape[0]) and out.domain is sx.ZZ
    assert out.rep.fmt=='dense' and values(out.to_list())==case['output'] and a.to_list()==before
    if fmt=='dense':
        r=ReferenceMatrix([[sp.ZZ(v) for v in row] for row in case['input']],shape,sp.ZZ)
        expected=ref.hermite_normal_form(r,D=sp.ZZ(case['D']))
        assert values(expected.to_list())==case['output']
        if case['valid_determinant_multiple']:lattice_certificate(case['input'],case['output'],shape)


@pytest.mark.parametrize('entry',['private','dense','sparse','matrix','immutable'])
@pytest.mark.parametrize('multiple',[1,2,7])
def test_first_public_routes_and_return_kinds(entry,multiple,monkeypatch):
    rows=[[12,6,4],[3,9,6],[2,16,14]];calls=[];native=sx._core._ddm_hnf_modular
    def record(*a):
        result=native(*a);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_hnf_modular',record)
    if entry in ('matrix','immutable'):
        a=(sx.Matrix if entry=='matrix' else sx.ImmutableMatrix)(rows)
        out=public.hermite_normal_form(a,D=300*multiple)
        assert out==sx.Matrix([[10,0,2],[0,15,3],[0,0,2]]) and a==sx.Matrix(rows)
    else:
        a=matrix(rows,(3,3))
        if entry=='sparse':a=a.to_sparse()
        out=(nf._hermite_normal_form_modulo_D(a,sx.ZZ(300*multiple)) if entry=='private'
             else nf.hermite_normal_form(a,D=sx.ZZ(300*multiple)))
        assert values(out.to_list())==[[10,0,2],[0,15,3],[0,0,2]] and values(a.to_list())==rows
    assert calls==[True]


def normalized_aliases(refs):
    seen={};result=[]
    for v in refs:
        key=id(v)
        if key not in seen:seen[key]=len(seen)
        result.append(seen[key])
    return result


def observe(function,rows,shape,modulus,constructor,error=False):
    d=sx.ZZ(modulus);work=[[sx.ZZ(v) for v in row] for row in rows]
    before=[v for row in work for v in row];initial_rows=list(work);captured=[]
    class Source:
        domain=sx.ZZ
        def to_list(self):return work
    a=Source();a.shape=shape
    def capture(w,result_shape,domain):
        refs=[d,work]+initial_rows+before+[w,w.default_factory]
        refs+=list(w)
        refs+=list(w.values())
        refs+=[key for row in w.values() for key in row]
        refs+=[value for row in w.values() for value in row.values()]
        refs+=[v for row in work for v in row]
        captured.append({'working':values(work),'rows':[list(row) for row in w.values()],
                         'outer_keys':list(w),'aliases':normalized_aliases(refs),
                         'output':[[int(v) for v in row.values()] for row in w.values()]})
        if error:raise RuntimeError('constructor after completed W')
        return constructor(w,result_shape,domain)
    nf.DomainMatrix=capture
    try:
        try:result=function(a,d)
        except RuntimeError as exc:
            assert error;out=(type(exc).__name__,str(exc))
        else:out=values(result.to_list());assert not error
    finally:nf.DomainMatrix=constructor
    return out,captured


@pytest.mark.parametrize('rows,d', [([[1,0]],1),([[0,1]],1),([[2,0]],2),([[2,0]],4),([[0,2]],4),
                                  ([[0,0]],3),([[0,0],[0,0]],6),([[2,0],[0,3]],6),
                                  ([[12,6,4],[3,9,6],[2,16,14]],600)])
@pytest.mark.parametrize('error',[False,True])
def test_complete_work_output_modulus_alias_graph_and_final_constructor_barrier(rows,d,error,monkeypatch):
    native=sx._core._ddm_hnf_modular;calls=[]
    def record(*a):
        result=native(*a);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_hnf_modular',record)
    shape=(len(rows),len(rows[0]));constructor=nf.DomainMatrix
    expected=observe(original(),rows,shape,d,constructor,error)
    actual=observe(nf._hermite_normal_form_modulo_D,rows,shape,d,constructor,error)
    assert actual==expected and calls==[True]


def test_keys_above_small_int_cache_share_columns_and_keep_outer_key_provenance(monkeypatch):
    # Independently reconstruct only Python's range/dictionary key provenance.
    # Numeric reduction never deletes a key, so its later assignments cannot
    # replace the original equal key object stored by these exact dictionaries.
    from collections import defaultdict
    m=258;expected=defaultdict(dict)
    for i in range(m-1,-1,-1):
        for ii in range(m):expected[ii][i]=None
    def keys(w):return normalized_aliases(list(w)+[key for row in w.values() for key in row])
    expected_keys=keys(expected);observed=[];constructor=nf.DomainMatrix
    def capture(w,shape,domain):
        observed.append((list(w),[list(row) for row in w.values()],keys(w)))
        return constructor(w,shape,domain)
    monkeypatch.setattr(nf,'DomainMatrix',capture)
    calls=[];native=sx._core._ddm_hnf_modular
    def record(*a):
        r=native(*a);calls.append(r);return r
    monkeypatch.setattr(sx._core,'_ddm_hnf_modular',record)
    # The unchecked zero matrix with D=1 follows the original recurrence's
    # identity output. It is an explicit outside-hypothesis control, not HNF.
    a=matrix([[0]*m for _ in range(m)],(m,m))
    out=nf._hermite_normal_form_modulo_D(a,sx.ZZ(1))
    assert values(out.to_list())==[[int(i==j) for j in range(m)] for i in range(m)]
    assert calls==[True] and observed==[(list(expected),[list(row) for row in expected.values()],expected_keys)]


@pytest.mark.parametrize('name',['_gcdex','add_columns','symmetric_residue','range','len'])
def test_actual_helper_callbacks_are_preserved(name,monkeypatch):
    events=[];original_value=nf.__dict__.get(name,getattr(__import__('builtins'),name,None))
    def hook(*a):events.append(name);return original_value(*a)
    monkeypatch.setattr(nf,name,hook,raising=False)
    a=matrix([[12,6,4],[3,9,6],[2,16,14]],(3,3))
    expected=original()(a,sx.ZZ(300));before=events[:];events.clear()
    actual=nf._hermite_normal_form_modulo_D(a,sx.ZZ(300))
    assert actual==expected and events==before and events


@pytest.mark.parametrize('name',['__lt__','__le__','__neg__','__mul__','__rmul__','__add__','__sub__','__floordiv__','__mod__','__eq__','__ne__','parent'])
def test_actual_coefficient_callbacks_keep_full_original_body(name,monkeypatch):
    events=[];method=getattr(sx._core.DomainElement,name)
    def hook(*a):events.append(name);return method(*a)
    monkeypatch.setattr(sx._core.DomainElement,name,hook)
    a=matrix([[12,6,4],[3,9,6],[2,16,14]],(3,3));d=sx.ZZ(300)
    expected=original()(a,d);before=events[:];events.clear()
    actual=nf._hermite_normal_form_modulo_D(a,d);after=events[:]
    assert values(actual.to_list())==values(expected.to_list()) and after==before


def test_inplace_modulus_callback_is_not_bypassed(monkeypatch):
    events=[];divide=sx._core.DomainElement.__floordiv__
    def hook(self,other):events.append('ifloordiv');return divide(self,other)
    monkeypatch.setattr(sx._core.DomainElement,'__ifloordiv__',hook,raising=False)
    a=matrix([[2,0],[0,3]],(2,2));d=sx.ZZ(6)
    expected=original()(a,d);before=events[:];events.clear()
    actual=nf._hermite_normal_form_modulo_D(a,d)
    assert values(actual.to_list())==values(expected.to_list()) and events==before==['ifloordiv','ifloordiv']


@pytest.mark.parametrize('kind',['subclass','foreign-factory','prefix-wrapper'])
def test_real_defaultdict_prefix_and_factory_callbacks(kind,monkeypatch):
    from collections import defaultdict
    events=[]
    class Subclass(defaultdict):
        def __getitem__(self,key):events.append(('getitem',key));return super().__getitem__(key)
    def create(factory):
        events.append(('create',factory is dict))
        if kind=='subclass':return Subclass(factory)
        if kind=='foreign-factory':
            def custom():events.append(('factory',));return {}
            return defaultdict(custom)
        return defaultdict(factory)
    monkeypatch.setattr(nf,'defaultdict',create)
    a=matrix([[2,0],[0,3]],(2,2));d=sx.ZZ(6)
    expected=original()(a,d);before=events[:];events.clear();actual=nf._hermite_normal_form_modulo_D(a,d)
    assert values(actual.to_list())==values(expected.to_list()) and events==before


@pytest.mark.parametrize('case',[('domain',None),('D-int',3),('D-zero',0),('D-negative',-3),('shape',3)])
def test_original_error_order_and_actual_prefixes(case,monkeypatch):
    kind,d=case;events=[]
    class Source:
        domain=sx.QQ if kind=='domain' else sx.ZZ
        shape=(2,1) if kind=='shape' else (1,2)
        def to_list(self):events.append('to_list');return [[sx.ZZ(1),sx.ZZ(0)]]
    constructor=nf.defaultdict
    def make(*a):events.append('defaultdict');return constructor(*a)
    monkeypatch.setattr(nf,'defaultdict',make)
    d=sx.ZZ(d) if kind in ('D-zero','D-negative','shape') else d
    def execute(f):
        events.clear()
        try:f(Source(),d)
        except Exception as e:return type(e).__name__,str(e),events[:]
        raise AssertionError('expected original error')
    assert execute(nf._hermite_normal_form_modulo_D)==execute(original())


def test_rank_check_declines_modular_for_deficient_rank_and_keeps_ordinary(monkeypatch):
    calls=[];native=sx._core._ddm_hnf_modular
    def record(*a):calls.append(1);return native(*a)
    monkeypatch.setattr(sx._core,'_ddm_hnf_modular',record)
    a=matrix([[1,2],[2,4]],(2,2));out=nf.hermite_normal_form(a,D=sx.ZZ(2),check_rank=True)
    assert out==nf.hermite_normal_form(a) and calls==[]
    b=matrix([[1,0],[0,2]],(2,2));assert nf.hermite_normal_form(b,D=sx.ZZ(2),check_rank=True)==b and calls==[1]


@pytest.mark.parametrize('digits',[200,768])
def test_large_modulus_and_coefficients_are_exact(digits,monkeypatch):
    big=2**digits;rows=[[big,-3],[7,-big]];a=matrix(rows,(2,2));d=sx.ZZ(big*big-21)
    calls=[];native=sx._core._ddm_hnf_modular
    def record(*args):
        r=native(*args);calls.append(r);return r
    monkeypatch.setattr(sx._core,'_ddm_hnf_modular',record)
    out=nf.hermite_normal_form(a,D=d)
    assert out==original()(a,d) and calls==[True]
    lattice_certificate(rows,values(out.to_list()),(2,2))


PROTOCOL_SCRIPT="import ast,builtins,copy,importlib,importlib.abc,json,sys,types\nfrom pathlib import Path\nclass NoRef(importlib.abc.MetaPathFinder):\n def find_spec(self,fullname,path=None,target=None):\n  if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('reference import')\nsys.meta_path.insert(0,NoRef())\nimport symox as sx\nnf=importlib.import_module('symox.polys.matrices.normalforms')\nprinting=importlib.import_module('symox.polys.domains._element_print')\nalgebraic=importlib.import_module('symox.polys.domains._algebraic_methods')\nfrom symox.polys.matrices.domainmatrix import DomainMatrix\nfrom symox.polys.matrices.ddm import DDM\nmode,target=sys.argv[1:3]\ntree=ast.parse(Path(nf.__file__).read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_hermite_normal_form_modulo_D')\nnode.body=[n for n in node.body if not (isinstance(n,ast.If) and isinstance(n.test,ast.Call) and isinstance(n.test.func,ast.Name) and n.test.func.id=='_native_hnf_modular')]\ncode=compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),nf.__file__,'exec');pure=types.FunctionType(next(c for c in code.co_consts if type(c)is types.CodeType),nf.__dict__)\na=DomainMatrix.from_rep(DDM([[sx.ZZ(v) for v in row] for row in [[12,6,4],[3,9,6],[2,16,14]]],(3,3),sx.ZZ));d=sx.ZZ(300)\ncore=sx._core;events=[];calls=[];audit=[];native=core._ddm_hnf_modular\n\ndef record(*a):\n r=native(*a);calls.append(r);return r\ncore._ddm_hnf_modular=record\nclass TupleSubclass(tuple):\n def __len__(self):events.append('tuple-len');raise AssertionError('tuple-len')\n def __iter__(self):events.append('tuple-iter');raise AssertionError('tuple-iter')\n def __getitem__(self,k):events.append('tuple-getitem');raise AssertionError('tuple-getitem')\n def __eq__(self,other):events.append('tuple-eq');raise AssertionError('tuple-eq')\nclass DictSubclass(dict):\n def __getitem__(self,k):events.append('dict-getitem');raise AssertionError('dict-getitem')\n def __iter__(self):events.append('dict-iter');raise AssertionError('dict-iter')\n def __len__(self):events.append('dict-len');raise AssertionError('dict-len')\nowner=printing if target in ('_int_congruent','_eq','_ne') else core.DomainElement if target.startswith('__') else nf\noriginal=owner.__dict__.get(target)\ndef replace(function):\n if isinstance(owner,types.ModuleType):owner.__dict__[target]=function\n else:setattr(owner,target,function)\n if target in ('_eq','_ne'):setattr(core.DomainElement,'__'+target.strip('_')+'__',function)\nif mode in ('closure','defaults','kwdefaults','consts','cached-builtins','changed-before-import'):\n code,namespace,name,defaults,closure=original.__code__,original.__globals__,original.__name__,original.__defaults__,original.__closure__\n if mode=='closure':closure=TupleSubclass(closure)\n if mode=='defaults':defaults=TupleSubclass((None,))\n if mode=='consts':code=code.replace(co_consts=TupleSubclass(code.co_consts))\n if mode=='changed-before-import':code=code.replace(co_consts=code.co_consts+('not the original semantics',))\n if mode=='cached-builtins':\n  old=namespace.get('__builtins__');foreign=dict(builtins.__dict__)\n  for key in ('range','type','int','len'):\n   def hook(*args,_name=key,_native=getattr(builtins,key)):\n    events.append('cached-'+_name);return _native(*args)\n   foreign[key]=hook\n  namespace['__builtins__']=foreign\n function=types.FunctionType(code,namespace,name,defaults,closure)\n if mode=='cached-builtins':namespace['__builtins__']=old\n if mode=='kwdefaults':function.__kwdefaults__=DictSubclass()\n replace(function)\nif mode=='outer-builtins':\n old=nf.__dict__.get('__builtins__');foreign=dict(builtins.__dict__)\n def foreign_range(*args):events.append('outer-range');return range(*args)\n foreign['range']=foreign_range;nf.__dict__['__builtins__']=foreign\n actual=types.FunctionType(nf._hermite_normal_form_modulo_D.__code__,nf.__dict__)\n pure=types.FunctionType(pure.__code__,nf.__dict__);nf.__dict__['__builtins__']=old\nelse:actual=nf._hermite_normal_form_modulo_D\nif mode=='metadata-none':core._ddm_function_metadata=None\nif mode=='metadata-object':\n class Foreign:\n  def __getattribute__(self,name):events.append('foreign-'+name);raise RuntimeError(name)\n core._ddm_function_metadata=Foreign()\nrestore_core=False\nif mode=='core-module':\n class ForeignModule(types.ModuleType):\n  def __getattribute__(self,name):events.append('foreign-'+name);raise RuntimeError(name)\n core.__class__=ForeignModule;restore_core=True\nif mode in ('core-key','namespace-key'):\n class ForeignKey:\n  def __hash__(self):return hash('_ddm_function_metadata')\n  def __eq__(self,other):events.append('equality');raise RuntimeError('equality')\n namespace=core.__dict__ if mode=='core-key' else nf.__dict__\n namespace.pop('_ddm_function_metadata',None);namespace[ForeignKey()]=1\nif mode=='changed-before-import' or target=='preimport':\n module_name='symox.polys.matrices._hnf_modular_native';sys.modules.pop(module_name,None)\n guard=importlib.import_module(module_name);guard.register(nf.__dict__);nf._native_hnf_modular=guard.try_hnf_modular\nif mode=='audit':\n def audit_hook(event,args):\n  if event=='object.__getattr__':audit.append(args[1]);raise AssertionError('new audit event '+args[1])\n sys.addaudithook(audit_hook)\ndef execute(f):\n events.clear();calls.clear()\n try:out=f(a,d)\n except Exception as e:return {'error':type(e).__name__,'message':str(e)},events[:],calls[:]\n observed=events[:];selected=calls[:]\n return {'values':[[str(v) for v in row] for row in out.to_list()],'shape':out.shape},observed,selected\nactual_result,actual_events,actual_calls=execute(actual)\nexpected_result,expected_events,expected_calls=execute(pure)\nif restore_core:core.__class__=types.ModuleType\nassert actual_result==expected_result,(mode,target,actual_result,expected_result)\nassert actual_events==expected_events,(mode,target,actual_events,expected_events)\nassert expected_calls==[] and audit==[]\nassert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)\nif mode=='audit':assert actual_calls==[True]\nelse:assert actual_calls==[],(mode,target,actual_calls)\nprint(json.dumps({'mode':mode,'target':target,'actual':actual_result,'actual_events':actual_events,'expected_events':expected_events,'native_calls':actual_calls,'audit_events':audit,'NoRef':True}))\n"


@pytest.mark.parametrize('mode,target',[('audit','call'),('outer-builtins','call'),('metadata-none','call')]
    + [('closure',name) for name in ('__floordiv__','__mod__')]
    + [(mode,name) for mode in ('defaults','kwdefaults','consts','cached-builtins','changed-before-import')
                  for name in ('_int_congruent','_eq','_ne','_gcdex','add_columns','symmetric_residue')]
    + [(mode,phase) for mode in ('metadata-object','core-key','namespace-key','core-module') for phase in ('call','preimport')])
def test_fresh_protocol_fallback_and_audit_without_reference_import(mode,target):
    import os,subprocess,sys
    r=subprocess.run([sys.executable,'-c',PROTOCOL_SCRIPT,mode,target],capture_output=True,text=True,
                     env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).parent.parent)))
    assert r.returncode==0,r.stdout+r.stderr
    result=json.loads(r.stdout)
    assert result['NoRef'] and result['audit_events']==[] and result['actual_events']==result['expected_events']
    assert result['native_calls']==([True] if mode=='audit' else [])
