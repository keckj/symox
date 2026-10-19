"""Check scoped source equivalence and closed standard-library dependencies.

These certificates cover signatures, control flow and standard-library calls.
They do not certify the protocols of every symbolic object supplied as input.
Reference code is inspected only by this development audit, never at runtime.
"""
import ast
import builtins
import dis
import hashlib
import inspect
import itertools
import importlib
import json
from pathlib import Path
import textwrap
import types

FUNCTIONS = ('dict_merge', 'group', 'take', 'prefixes', 'postfixes')
RENAMES = {'dict_merge': {'mapping': 'dict'}, 'group': {'values': 'g', 'key': 'k'}}
EXTENDED=(('sympy','prod'),('sympy.ntheory.modular','symmetric_residue'))+tuple(
    ('sympy.utilities.iterables',name) for name in (
        'rotate_left','rotate_right','common_prefix','common_suffix',
        'is_palindromic','roundrobin','filter_symbols','necklaces','bracelets',
        'multiset','ordered_partitions','runs','uniq','binary_partitions'))
STDLIB={'math','itertools','collections','functools','_functools','operator','_operator','builtins'}
AMBIENT={'eval','exec','compile','globals','locals','vars','dir','__import__',
         'open','input','print'}


def normalized(function, renames):
    tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
    node = tree.body[0]
    if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
        node.body.pop(0)
    class Names(ast.NodeTransformer):
        def visit_Name(self, value):
            value.id = renames.get(value.id, value.id)
            return value
        def visit_Attribute(self, value):
            if isinstance(value.value, ast.Name) and value.value.id == 'itertools' and value.attr == 'groupby':
                return ast.Name(id='groupby', ctx=ast.Load())
            return self.generic_visit(value)
    return ast.dump(Names().visit(tree), include_attributes=False)


def _global_names(code):
    # co_names also contains attribute spellings; LOAD_GLOBAL distinguishes
    # operator.mul from an unrelated global named mul.
    result={instruction.argval for instruction in dis.get_instructions(code)
            if instruction.opname=='LOAD_GLOBAL'}
    for constant in code.co_consts:
        if isinstance(constant,types.CodeType):result.update(_global_names(constant))
    return result


def verify_dependencies(reference,native,seen=None):
    """Fail closed on non-equivalent helpers, ambient state and opaque globals."""
    seen=set() if seen is None else seen
    pair=id(reference),id(native)
    if pair in seen:return []
    if not inspect.isfunction(reference) or not inspect.isfunction(native):
        raise ValueError('dependency is not a plain Python function')
    if inspect.signature(reference)!=inspect.signature(native):
        raise ValueError('dependency signatures differ')
    if reference.__code__.co_freevars or native.__code__.co_freevars:
        raise ValueError('closure dependencies require a separate review')
    a,b=normalized(reference,{}),normalized(native,{})
    if a!=b:raise ValueError('dependency source differs')
    for function in [reference,native]:
        for parameter in inspect.signature(function).parameters.values():
            value=parameter.default
            if value is inspect.Parameter.empty:continue
            if type(value) in (str,int,float,bool,type(None)):continue
            if getattr(value,'__module__',None) in STDLIB:continue
            raise ValueError('mutable or opaque default requires a separate review')
        tree=ast.parse(textwrap.dedent(inspect.getsource(function)))
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom) and (node.level or node.module not in STDLIB):
                raise ValueError('non-standard-library local dependency')
            if isinstance(node,ast.Import) and any(alias.name not in STDLIB for alias in node.names):
                raise ValueError('non-standard-library local dependency')
    seen.add(pair)
    names=_global_names(reference.__code__)
    if names!=_global_names(native.__code__):raise ValueError('global dependency names differ')
    def resolve(function,name):
        if name in function.__globals__:return function.__globals__[name]
        namespace=function.__globals__.get('__builtins__',builtins)
        return namespace[name] if isinstance(namespace,dict) else getattr(namespace,name)
    evidence=[]
    for name in sorted(names):
        if name in AMBIENT:raise ValueError('ambient-state dependency '+name)
        x,y=resolve(reference,name),resolve(native,name)
        if inspect.isfunction(x) and inspect.isfunction(y):
            evidence.extend(verify_dependencies(x,y,seen))
        elif x is y and (getattr(x,'__module__',None) in STDLIB or
                         isinstance(x,types.ModuleType) and x.__name__ in STDLIB):
            evidence.append({'name':name,'kind':'identical_standard_library_object',
                             'origin':getattr(x,'__module__',None) or x.__name__})
        elif type(x) in (str,int,float,bool,type(None)) and type(x) is type(y) and x==y:
            evidence.append({'name':name,'kind':'equal_immutable_constant','value':repr(x)})
        else:raise ValueError('unverified dependency '+name)
    evidence.append({'name':reference.__qualname__,'kind':'equivalent_python_helper',
        'normalized_ast_sha256':hashlib.sha256(a.encode()).hexdigest(),
        'reference_source_sha256':hashlib.sha256(inspect.getsource(reference).encode()).hexdigest(),
        'native_source_sha256':hashlib.sha256(inspect.getsource(native).encode()).hexdigest()})
    return evidence


def certificates(extended=False):
    import sympy
    from sympy.utilities import iterables as reference
    from symox.utilities import iterables as native
    if sympy.__version__ != '1.14.0':
        raise ValueError('These certificates target SymPy 1.14.0')
    records = []
    for name in FUNCTIONS:
        a, b = getattr(reference, name), getattr(native, name)
        assert str(inspect.signature(a)) == str(inspect.signature(b)), name
        source_a, source_b = normalized(a, {}), normalized(b, RENAMES.get(name, {}))
        assert source_a == source_b, name
        # Every external operation in these small function bodies is either a
        # verified builtin or the same itertools.groupby implementation.
        for global_name in ('range', 'zip', 'len', 'list'):
            target = getattr(builtins, global_name)
            assert a.__globals__.get(global_name, target) is target
            assert b.__globals__.get(global_name, target) is target
        assert reference.groupby is native.itertools.groupby is itertools.groupby
        records.append({
            'id': 'sympy.utilities.iterables.' + name,
            'reference_version': sympy.__version__,
            'status': 'source_equivalent',
            'signature': str(inspect.signature(a)),
            'normalized_ast_sha256': hashlib.sha256(source_a.encode()).hexdigest(),
            'reference_source_sha256': hashlib.sha256(inspect.getsource(a).encode()).hexdigest(),
            'native_source_sha256': hashlib.sha256(inspect.getsource(b).encode()).hexdigest(),
            'native_file': str(Path('python/symox')/Path(inspect.getsourcefile(b)).relative_to(Path(native.__file__).parents[1])),
            'local_alpha_renames': RENAMES.get(name, {}),
            'scope': 'Same Python runtime and input protocols: matching signatures, control flow, exact object preservation, error propagation and asymptotic time/space. Symbolic input-type protocols need their own certificates.',
            'global_sympy_parity_certified': False,
        })
    if extended:
        for module,name in EXTENDED:
            a=getattr(importlib.import_module(module),name)
            b=getattr(importlib.import_module(module.replace('sympy','symox',1)),name)
            dependencies=verify_dependencies(a,b)
            records.append({'id':module+'.'+name,'reference_version':sympy.__version__,
                'status':'source_equivalent','signature':str(inspect.signature(a)),
                'normalized_ast_sha256':hashlib.sha256(normalized(a,{}).encode()).hexdigest(),
                'reference_source_sha256':hashlib.sha256(inspect.getsource(a).encode()).hexdigest(),
                'native_source_sha256':hashlib.sha256(inspect.getsource(b).encode()).hexdigest(),
                'native_file':str(Path(inspect.getsourcefile(b)).resolve().relative_to(Path(__file__).resolve().parents[1])),
                'dependency_evidence':dependencies,
                'scope':'Same Python runtime and compatible input protocols; identical control flow, signatures, recursively verified helper bodies and shared standard-library operations. Preserves exact objects, exception operations, laziness, state transitions and asymptotic time/space. Does not certify symbolic operand protocols or inherited class APIs.',
                'global_sympy_parity_certified':False})
    from sympy_parity import native_fingerprint
    return {'schema': 1,'native_source_sha256':native_fingerprint(),'records': records}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    report=certificates(extended=True)
    Path(args.out).write_text(json.dumps(report, indent=2) + '\n')
    print(f"{len(report['records'])} scoped source-equivalence certificates checked")
