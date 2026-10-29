"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'lambdify')

import symox as _S
_implemented={}

def implemented_function(symfunc,implementation):
    if not callable(implementation):raise TypeError('implementation must be callable')
    if isinstance(symfunc,str):symfunc=_S.Function(symfunc)
    if not isinstance(symfunc,_S.Function):raise ValueError('implementation requires an undefined function')
    symfunc._imp_=implementation;_implemented[symfunc.name]=implementation
    return symfunc

def _recursive_to_string(doprint,arg):
    # SymPy's lambdify._recursive_to_string.
    from ..utilities.iterables import iterable
    if isinstance(arg,_S.Basic) or getattr(arg,'is_Matrix',False):return doprint(arg)
    if iterable(arg):
        if isinstance(arg,list):left,right='[',']'
        elif isinstance(arg,tuple):
            if not arg:return '()'
            left,right='(',',)'
        else:raise NotImplementedError('unhandled type: %s, %s'%(type(arg),arg))
        return left+', '.join(_recursive_to_string(doprint,e) for e in arg)+right
    if isinstance(arg,str):return arg
    return doprint(arg)

def lambdastr(args,expr,printer=None,dummify=None):
    # Port of SymPy 1.14 lambdify.lambdastr.
    import inspect
    from ..utilities.iterables import flatten,iterable,NotIterable
    DeferredVector=_S.DeferredVector
    if printer is not None:
        if inspect.isfunction(printer):lambdarepr=printer
        elif inspect.isclass(printer):lambdarepr=lambda expr:printer().doprint(expr)
        else:lambdarepr=lambda expr:printer.doprint(expr)
    else:
        from ..printing.lambdarepr import lambdarepr
    def sub_args(args,dummies_dict):
        if isinstance(args,str):return args
        elif isinstance(args,DeferredVector):return str(args)
        elif iterable(args):
            dummies=flatten([sub_args(a,dummies_dict) for a in args])
            return ','.join(str(a) for a in dummies)
        elif isinstance(args,(_S.Function,_S.Symbol,_S.Derivative)) or getattr(args,'is_Function',False):
            dummies=_S.Dummy()
            dummies_dict.update({args:dummies})
            return str(dummies)
        else:return str(args)
    def sub_expr(expr,dummies_dict):
        expr=_S.sympify(expr)
        if isinstance(expr,_S.Basic):expr=expr.xreplace(dummies_dict)
        elif isinstance(expr,list):expr=[sub_expr(a,dummies_dict) for a in expr]
        return expr
    def isiter(l):return iterable(l,exclude=(str,DeferredVector,NotIterable))
    def flat_indexes(iterable):
        n=0
        for el in iterable:
            if isiter(el):
                for ndeep in flat_indexes(el):yield (n,)+ndeep
            else:yield (n,)
            n+=1
    if dummify is None:
        dummify=any(isinstance(a,_S.Basic) and a.atoms(_S.Function,_S.Derivative)
                    for a in (args if isiter(args) else [args]))
    if isiter(args) and any(isiter(i) for i in args):
        dum_args=[str(_S.Dummy(str(i))) for i in range(len(args))]
        indexed_args=','.join([dum_args[ind[0]]+''.join(['[%s]'%k for k in ind[1:]])
                               for ind in flat_indexes(args)])
        lstr=lambdastr(flatten(args),expr,printer=printer,dummify=dummify)
        return 'lambda %s: (%s)(%s)'%(','.join(dum_args),lstr,indexed_args)
    dummies_dict={}
    if dummify:args=sub_args(args,dummies_dict)
    elif isinstance(args,str):pass
    elif iterable(args,exclude=DeferredVector):args=','.join(str(a) for a in args)
    if dummify and not isinstance(expr,str):expr=sub_expr(expr,dummies_dict)
    expr=_recursive_to_string(lambdarepr,expr)
    return 'lambda %s: (%s)'%(args,expr)

def _user_namespace(modules):
    # SymPy's lambdify namespace: name mappings of the given modules, the
    # first module taking precedence.
    table={}
    candidates=modules if isinstance(modules,(list,tuple)) else [modules]
    for module in reversed(candidates):
        if isinstance(module,dict):
            table.update((name,value) for name,value in module.items() if isinstance(name,str) and callable(value))
    return table

def _module_present(name,modules):
    # SymPy's lambdify._module_present over the module list.
    candidates=modules if isinstance(modules,(list,tuple)) else [modules]
    for module in candidates:
        if module==name if isinstance(module,str) else getattr(module,'__name__',None)==name:return True
    return False


def _normalized_sinc(modules):
    # SymPy's NumPy-family printers (selected unless mpmath is present) print
    # sinc(x) as a call of the normalized sinc on x/pi; a namespace entry
    # named ``sinc`` receives that argument.
    if _module_present('mpmath',modules):return False
    return any(_module_present(name,modules) for name in ('scipy','numpy','cupy','jax'))


def _integer_kernel(parameters,expressions,modules,native,vectors={}):
    """Wrap ``native`` so ``int`` arguments follow the Python integer
    arithmetic of SymPy's printed source, as ``_ExactIntegerLambda`` does.

    ``expressions`` is one expression or a list; ``vectors`` maps the
    parameters standing for a DeferredVector to its name, its entries being
    the symbols ``name[i]``. Returns ``native`` unless every expression is a polynomial
    with Integer coefficients in the parameters."""
    if modules is not None:
        names=modules if isinstance(modules,(list,tuple)) else [modules]
        if not all(isinstance(m,str) and m in ('math','mpmath','sympy','numpy','scipy') for m in names):return native
    import re
    names={p:'_a%d'%i for i,p in enumerate(parameters)}
    by_label={label:names[p] for p,label in dict(vectors).items()}
    def source(e):
        if e in names:return names[e]
        if getattr(e,'is_Symbol',False):
            match=re.fullmatch(r'(.+)\[(\d+)\]',e.name)
            if match and match.group(1) in by_label:return '%s[%s]'%(by_label[match.group(1)],match.group(2))
            return None
        if getattr(e,'is_Integer',False):return '(%d)'%int(e)
        if e.is_Add or e.is_Mul:
            parts=[source(a) for a in e.args]
            if None in parts:return None
            return '(%s)'%(' + ' if e.is_Add else '*').join(parts)
        if e.is_Pow and e.exp.is_Integer and e.exp>=0:
            base=source(e.base)
            return None if base is None else '(%s**%d)'%(base,int(e.exp))
        return None
    single=not isinstance(expressions,(list,tuple))
    try:bodies=[source(_S.sympify(e)) for e in ([expressions] if single else expressions)]
    except (AttributeError,TypeError):return native
    if not bodies or None in bodies:return native
    body=bodies[0] if single else '[%s]'%', '.join(bodies)
    exact=eval('lambda %s: %s'%(', '.join(names[p] for p in parameters),body),{})
    def integral(value):
        if type(value) is int:return True
        return isinstance(value,(list,tuple)) and any(type(v) is int for v in value)
    def call(*values,**kw):
        if not kw and any(integral(v) for v in values):return exact(*values)
        return native(*values,**kw)
    return call


def _compile_implemented(args,expr,modules,options):
    from .._host import is_host_function_class
    namespace=_user_namespace(modules)
    normalized_sinc=_normalized_sinc(modules)
    root=namespace.get('sqrt')
    def square_root(node):
        # SymPy's _hprint_Pow prints b**(1/2) as sqrt(b) and, for a
        # commutative power, b**(-1/2) as 1/sqrt(b); the unqualified name
        # resolves to a namespace entry called ``sqrt``.
        if root is None or not getattr(node,'is_Pow',False):return 0
        exponent=node.args[1]
        if exponent==_S.Rational(1,2):return 1
        if exponent==_S.Rational(-1,2) and node.is_commutative:return -1
        return 0
    def reciprocal_root(value):return 1/root(value)
    def implementation(node):
        kind=square_root(node)
        if kind:return root if kind>0 else reciprocal_root
        if is_host_function_class(type(node)):
            # A Function subclass: a namespace entry, else its _imp_.
            found=namespace.get(type(node).__name__)
            return found if found is not None else getattr(type(node),'_imp_',None)
        if not getattr(node,'is_Function',False):
            return None
        name=getattr(node,'class_name','')
        found=namespace.get(name)
        if found is None and getattr(node,'is_UndefinedFunction',False):
            found=_implemented.get(name)
        return found
    if not namespace and not _implemented and not _S._host.host_classes_defined():return None
    arguments=[args] if isinstance(args,_S.Expr) else list(args)
    def matches(node):return implementation(node) is not None
    def calls(value):
        # Visit each distinct subexpression once: shared subgraphs would
        # otherwise make the walk exponential in the expression depth.
        stack=[value];result=[];seen=set()
        while stack:
            node=stack.pop()
            try:
                if node in seen:continue
                seen.add(node)
            except TypeError:pass
            if matches(node):result.append(node)
            else:stack.extend(getattr(node,'args',()))
        return list(dict.fromkeys(result))
    if isinstance(expr,(list,tuple)):
        plans=[_compile_implemented(args,e,modules,options) or _integer_kernel(arguments,e,modules,_S._native_lambdify(args,e,modules=modules,**options)) for e in expr]
        return lambda *values:type(expr)(plan(*values) for plan in plans)
    found=calls(expr)
    if not found:return None
    def plan(value):
        outer=calls(value)
        if not outer:return _integer_kernel(arguments,value,modules,_S._native_lambdify(arguments,value,modules=modules,**options))
        temporary=[_S.Dummy('implemented') for _ in outer]
        body=_S.sympify(value).xreplace(dict(zip(outer,temporary)))
        evaluator=_integer_kernel(arguments+temporary,body,modules,_S._native_lambdify(arguments+temporary,body,modules=modules,**options))
        def operands(node):
            if square_root(node):return [node.args[0]]
            if normalized_sinc and not is_host_function_class(type(node)) and getattr(node,'class_name','')=='sinc':
                return [node.args[0]/_S.pi]
            return node.args
        subplans=[(implementation(node),[plan(a) for a in operands(node)]) for node in outer]
        def evaluate(*values):return evaluator(*values,*(function(*(p(*values) for p in children)) for function,children in subplans))
        return evaluate
    return plan(expr)
__all__ += ['implemented_function','lambdastr']


def _compile_object(args, expr):
    """Preserve Python/NumPy operand protocols outside the float64 fast path."""
    import numpy
    from ..printing.numpy import NumPyPrinter
    source = lambdastr(args, expr, printer=NumPyPrinter, dummify=True)
    return eval(source, {'numpy': numpy})


def _numpy_modules(modules):
    """Whether ``modules`` selects SymPy's NumPy-based printers only."""
    candidates=modules if isinstance(modules,(list,tuple)) else [modules]
    if not candidates:return False
    for module in candidates:
        if module is None or isinstance(module,dict):continue
        name=module if isinstance(module,str) else getattr(module,'__name__',None)
        if name not in ('numpy','scipy'):return False
    return True


def _compile_structured(args,expr,modules,options):
    """Flatten native kernels once, then restore argument/output structure."""
    from .._cse_matrices import is_concrete_matrix
    from ..tensor.array import NDimArray
    def sequence(value):return isinstance(value,(list,tuple)) or getattr(value,'is_Tuple',False)
    def contains_array(value):
        if isinstance(value,NDimArray):return True
        return sequence(value) and any(contains_array(child) for child in value)
    arrays=contains_array(expr) or contains_array(args)
    if not _numpy_modules(modules):
        if arrays:
            # SymPy prints arrays only with its NumPy-based printers.
            raise NotImplementedError('Unsupported by the %r lambdify modules: N-dim arrays require NumPy'%(modules,))
        return None
    roots=list(args) if sequence(args) else [args]
    if not (arrays or any(sequence(arg) for arg in roots) or sequence(expr) or is_concrete_matrix(expr)):
        return None
    arguments=[];vectors=[]
    def argument_plan(value):
        # SymPy unpacks an array argument like the nested list it iterates as.
        if sequence(value) or isinstance(value,NDimArray):return tuple(argument_plan(child) for child in value)
        arguments.append(_S.sympify(value))
        if isinstance(value,_S.DeferredVector):vectors.append(arguments[-1])
        return None
    input_plans=[argument_plan(arg) for arg in roots]
    expressions=[]
    def output_plan(value):
        if sequence(value):
            plans=[output_plan(child) for child in value]
            constructor=list if isinstance(value,list) else tuple
            return lambda values:constructor(plan(values) for plan in plans)
        if isinstance(value,NDimArray):
            if value._mutable:
                raise NotImplementedError('unhandled type: %s'%type(value))
            # SymPy's NumPy printer emits numpy.array of the nested entries,
            # or numpy.zeros(shape) for an array without entries.
            import numpy
            value._materialize_limit()
            shape=tuple(value.shape)
            indices=[len(expressions)+i for i in range(len(value))]
            expressions.extend(_S.sympify(entry) for entry in value._flat_values())
            def array(values):
                if not indices:return numpy.zeros(shape)
                entries=numpy.asarray([values[index] for index in indices])
                return entries.reshape(shape+entries.shape[1:])
            return array
        if is_concrete_matrix(value):
            import numpy
            indices=[len(expressions)+i for i in range(len(value))]
            expressions.extend(value)
            def matrix(values):
                if not indices:return numpy.empty(value.shape)
                entries=numpy.asarray([values[index] for index in indices])
                return entries.reshape(value.shape+entries.shape[1:])
            return matrix
        index=len(expressions);expressions.append(_S.sympify(value))
        return lambda values:values[index]
    restore=output_plan(expr)
    parameters=[_S.Dummy('argument') for _ in arguments]
    replacements=dict(zip(arguments,parameters))
    expressions=[expression.xreplace(replacements) for expression in expressions]
    evaluator=(_compile_implemented(parameters,expressions,modules,options)
               or _integer_kernel(parameters,expressions,modules,
                                  _S._native_lambdify(parameters,expressions,modules=modules,**options),
                                  {replacements[v]:v.name for v in vectors}))
    def call(*values):
        if len(values)!=len(roots):raise TypeError('expected %d argument(s), got %d'%(len(roots),len(values)))
        flattened=[]
        def unpack(plan,value):
            if plan is None:flattened.append(value);return
            children=tuple(value)
            if len(children)!=len(plan):raise ValueError('nested numerical argument has incorrect length')
            for child,entry in zip(plan,children):unpack(child,entry)
        for plan,value in zip(input_plans,values):unpack(plan,value)
        return restore(evaluator(*flattened))
    return call


# --- mpmath -----------------------------------------------------------------
# SymPy's lambdify with the 'mpmath' module prints the expression with the
# MpmathPrinter and evaluates it in a namespace built by
# ``from mpmath import *`` and MPMATH_TRANSLATIONS, so the callable computes
# with mpf numbers at the current mpmath precision.

_MPMATH_TRANSLATIONS = {
    'Abs': 'fabs', 'elliptic_k': 'ellipk', 'elliptic_f': 'ellipf', 'elliptic_e': 'ellipe',
    'elliptic_pi': 'ellippi', 'ceiling': 'ceil', 'chebyshevt': 'chebyt', 'chebyshevu': 'chebyu',
    'assoc_legendre': 'legenp', 'E': 'e', 'I': 'j', 'ln': 'log', 'oo': 'inf',
    'LambertW': 'lambertw', 'MutableDenseMatrix': 'matrix', 'ImmutableDenseMatrix': 'matrix',
    'conjugate': 'conj', 'dirichlet_eta': 'altzeta', 'Ei': 'ei', 'Shi': 'shi', 'Chi': 'chi',
    'Si': 'si', 'Ci': 'ci', 'RisingFactorial': 'rf', 'FallingFactorial': 'ff',
    'betainc_regularized': 'betainc',
}
_MPMATH_CONSTANTS = {'Pi': 'pi', 'Exp1': 'e', 'ImaginaryUnit': '1j', 'Infinity': 'inf',
                     'NegativeInfinity': '(-inf)', 'NaN': 'nan', 'ComplexInfinity': 'zoo',
                     'EulerGamma': 'euler', 'Catalan': 'catalan', 'GoldenRatio': 'phi'}
# Printed names of SymPy's PythonCodePrinter (builtins or Piecewise forms).
_MPMATH_PRINTED = {'Abs': 'abs', 'Min': 'min', 'Max': 'max'}
_MPMATH_EXCLUDED = ('numpy', 'scipy', 'cupy', 'jax', 'tensorflow', 'numexpr', 'math', 'cmath', 'sympy', 'symox')


class _Unprintable(Exception):
    pass


def _mpmath_only(modules):
    """Whether ``modules`` selects the mpmath namespace and printer only."""
    candidates=modules if isinstance(modules,(list,tuple)) else [modules]
    if not _module_present('mpmath',candidates):return False
    for module in candidates:
        if isinstance(module,dict):continue
        name=module if isinstance(module,str) else getattr(module,'__name__',None)
        if name!='mpmath':return False
    return True


class _MpmathSource:
    def __init__(self,symbol_names):
        self.names=symbol_names;self.functions=set()
    def number(self,value):
        if value.is_Integer:return str(int(value)) if value>=0 else '(%d)'%int(value)
        if value.is_Rational:return '(mpf(%d)/mpf(%d))'%(int(value.p),int(value.q))
        if value.is_Float:
            from .._host import _float_mpf
            return 'mpf(%s)'%(tuple(map(int,_float_mpf(value))),)
        raise _Unprintable(value)
    def call(self,name,args):
        self.functions.add(name)
        return name+'('+', '.join(args)+')'
    def __call__(self,e):
        S=_S
        if isinstance(e,(list,tuple)):
            body=', '.join(self(v) for v in e)
            return ('['+body+']') if isinstance(e,list) else '('+body+(',' if len(e)==1 else '')+')'
        if isinstance(e,bool):return str(e)
        e=S.sympify(e)
        if e in self.names:return self.names[e]
        name=getattr(e,'class_name',type(e).__name__)
        if getattr(e,'is_Matrix',False) and hasattr(e,'tolist'):
            return self.call('ImmutableDenseMatrix',['['+', '.join('['+', '.join(self(v) for v in row)+']' for row in e.tolist())+']'])
        if name in _MPMATH_CONSTANTS:
            value=_MPMATH_CONSTANTS[name]
            if value[0].isalpha():self.functions.add(value)
            elif 'inf' in value:self.functions.add('inf')
            return value
        if e is S.true or name=='BooleanTrue':return 'True'
        if e is S.false or name=='BooleanFalse':return 'False'
        if getattr(e,'is_Symbol',False):
            text=str(e)
            if not text.isidentifier():raise _Unprintable(e)
            return text
        if getattr(e,'is_Number',False):return self.number(e)
        args=e.args
        if getattr(e,'is_Add',False):return '('+' + '.join(self(a) for a in args)+')'
        if getattr(e,'is_Mul',False):return '('+'*'.join(self(a) for a in args)+')'
        if getattr(e,'is_Pow',False):
            base,exp=args
            if exp==S.Rational(1,2):return self.call('sqrt',[self(base)])
            if exp==S.Rational(-1,2) and e.is_commutative:return '(1/'+self.call('sqrt',[self(base)])+')'
            return '('+self(base)+'**'+self(exp)+')'
        relations={'Equality':'==','Unequality':'!=','LessThan':'<=','StrictLessThan':'<',
                   'GreaterThan':'>=','StrictGreaterThan':'>'}
        if name in relations:return '('+self(args[0])+' '+relations[name]+' '+self(args[1])+')'
        if name in ('And','Or'):return '('+(' and ' if name=='And' else ' or ').join(self(a) for a in args)+')'
        if name=='Not':return '(not '+self(args[0])+')'
        if name=='Piecewise':
            from ..printing.mathml import _pairs
            tail=None
            for value,condition in reversed(_pairs(e)):
                if str(condition) in ('True','true'):tail=self(value)
                else:tail='(('+self(value)+') if ('+self(condition)+') else ('+(tail if tail is not None else 'None')+'))'
            return '('+tail+')'
        if name=='sinc':
            a=self(args[0]);self.functions.add('sin')
            return '(((sin('+a+')/'+a+') if ('+a+' != 0) else (1)))'
        if name=='Heaviside':
            a=self(args[0]);h0=self(args[1]) if len(args)>1 else 'mpf(1)/mpf(2)'
            return '(((0) if ('+a+' < 0) else ('+h0+') if ('+a+' == 0) else (1)))'
        if name=='log1p':return self.call('log1p',[self(args[0])])
        if name=='expm1':return self.call('expm1',[self(args[0])])
        if name=='log2':return '('+self.call('log',[self(args[0])])+'/log(2))'
        if name=='uppergamma':return self.call('gammainc',[self(args[0]),self(args[1]),'inf'])
        if name=='lowergamma':return self.call('gammainc',[self(args[0]),'0',self(args[1])])
        if getattr(e,'is_Function',False) and name.isidentifier():
            return self.call(_MPMATH_PRINTED.get(name,name),[self(a) for a in args])
        raise _Unprintable(e)


def _mpmath_namespace(modules,expr):
    import builtins
    import mpmath
    namespace={}
    candidates=modules if isinstance(modules,(list,tuple)) else [modules]
    for module in reversed(candidates):
        if isinstance(module,dict):namespace.update(module);continue
        exec('from mpmath import *',namespace)
        for name,target in _MPMATH_TRANSLATIONS.items():namespace[name]=namespace[target]
    namespace.setdefault('Abs',abs)
    namespace.update({'builtins':builtins,'range':range,'mpmath':mpmath})
    return namespace


def _compile_mpmath(args,expr,modules,options):
    """SymPy's mpmath lambdify, or None for forms the printer does not cover."""
    if not _mpmath_only(modules) or options.get('printer') is not None or options.get('cse'):return None
    from .._host import is_host_function_class
    if is_host_function_class(type(expr)) or (isinstance(expr,(list,tuple)) and any(is_host_function_class(type(v)) for v in expr)):
        return None
    roots=[args] if isinstance(args,(_S.Expr,str)) else list(args)
    mapping={};parameters=[]
    def descend(value,path):
        if isinstance(value,(tuple,list)) or getattr(value,'is_Tuple',False):
            for i,child in enumerate(value):descend(child,path+'['+str(i)+']')
        else:mapping[_S.sympify(value)]=path
    try:
        for i,arg in enumerate(roots):
            name=getattr(arg,'name',str(arg))
            if options.get('dummify') or not name.isidentifier() or isinstance(arg,(tuple,list)) or getattr(arg,'is_Tuple',False):
                name='_arg'+str(i)
            parameters.append(name);descend(arg,name)
        printer=_MpmathSource(mapping)
        body=printer(expr)
    except _Unprintable:
        return None
    except Exception:
        return None
    namespace=_mpmath_namespace(modules,expr)
    import builtins
    for name in printer.functions:
        if name not in namespace and not hasattr(builtins,name):
            if name in _implemented:namespace[name]=_implemented[name]
            else:return None
    source='def _lambdifygenerated('+', '.join(parameters)+'):\n    return '+body+'\n'
    import linecache,textwrap,weakref
    global _lambdify_generated_counter
    filename='<lambdifygenerated-%s>'%_lambdify_generated_counter
    _lambdify_generated_counter+=1
    exec(compile(source,filename,'exec'),namespace)
    linecache.cache[filename]=(len(source),None,source.splitlines(True),filename)
    function=namespace['_lambdifygenerated']
    weakref.finalize(function,linecache.cache.pop,filename,None)
    # SymPy's docstring (sig, str(expr) cut at 78 characters, source).
    sig=textwrap.fill('func('+', '.join(parameters)+')',subsequent_indent=' '*8)
    limit=options.get('docstring_limit',1000)
    if isinstance(expr,_S.Basic) and _too_large_for_docstring(expr,limit):
        expr_str="EXPRESSION REDACTED DUE TO LENGTH, (see lambdify's `docstring_limit`)"
        src_str="SOURCE CODE REDACTED DUE TO LENGTH, (see lambdify's `docstring_limit`)"
    else:
        expr_str=str(expr)
        if len(expr_str)>78:expr_str=textwrap.wrap(expr_str,75)[0]+'...'
        src_str=source
    function.__doc__=('Created with lambdify. Signature:\n\n'+sig+'\n\nExpression:\n\n'+expr_str
                      +'\n\nSource code:\n\n'+src_str+'\n\nImported modules:\n\n')
    return function


# SymPy's lambdify module namespaces (sympy.utilities.lambdify.MODULES):
# the default entries, the import commands and the name translations.
_NAMESPACE_MODULES = {
    'math': ({}, ('from math import *',), {'ceiling': 'ceil', 'E': 'e', 'ln': 'log'}),
    'cmath': ({}, ('import cmath; from cmath import *',), {}),
    'mpmath': ({}, ('from mpmath import *',), {
        'Abs': 'fabs', 'elliptic_k': 'ellipk', 'elliptic_f': 'ellipf', 'elliptic_e': 'ellipe',
        'elliptic_pi': 'ellippi', 'ceiling': 'ceil', 'chebyshevt': 'chebyt', 'chebyshevu': 'chebyu',
        'assoc_legendre': 'legenp', 'E': 'e', 'I': 'j', 'ln': 'log', 'oo': 'inf',
        'LambertW': 'lambertw', 'MutableDenseMatrix': 'matrix', 'ImmutableDenseMatrix': 'matrix',
        'conjugate': 'conj', 'dirichlet_eta': 'altzeta', 'Ei': 'ei', 'Shi': 'shi', 'Chi': 'chi',
        'Si': 'si', 'Ci': 'ci', 'RisingFactorial': 'rf', 'FallingFactorial': 'ff',
        'betainc_regularized': 'betainc'}),
    'numpy': ({'I': 1j}, ('import numpy; from numpy import *; from numpy.linalg import *',),
              {'Heaviside': 'heaviside'}),
    'scipy': ({'I': 1j}, ('import scipy; import numpy; from scipy.special import *',),
              {'jn': 'spherical_jn', 'yn': 'spherical_yn'}),
    'cupy': ({'I': 1j}, ('import cupy',), {}),
    'jax': ({'I': 1j}, ('import jax',), {}),
    'tensorflow': ({}, ('import tensorflow',), {}),
    'torch': ({'I': 1j}, ('import torch',), {}),
    'sympy': ({}, ('from symox.functions import *', 'from symox.matrices import *',
                   'from symox import Integral, pi, oo, nan, zoo, E, I'), {}),
    'numexpr': ({}, ('import numexpr',), {}),
}
_NAMESPACE_MODULES['symox'] = _NAMESPACE_MODULES['sympy']
_NAMESPACE_CACHE = {}


def _module_namespace(name):
    """SymPy's ``lambdify._import``: the translation dictionary of a module
    name, built once."""
    if name in _NAMESPACE_CACHE:
        return _NAMESPACE_CACHE[name]
    try:
        default, commands, translations = _NAMESPACE_MODULES[name]
    except KeyError:
        raise NameError("'%s' module cannot be used for lambdification" % name) from None
    namespace = dict(default)
    for command in commands:
        try:
            exec(command, {}, namespace)
        except ImportError:
            raise ImportError("Cannot import '%s' with '%s' command" % (name, command)) from None
    for symbol_name, target in translations.items():
        namespace[symbol_name] = namespace[target]
    namespace.setdefault('Abs', abs)
    _NAMESPACE_CACHE[name] = namespace
    return namespace


def _lambdify_namespace(modules, exprs):
    """The ``__globals__`` of SymPy's generated function: the module
    namespaces (the first module taking precedence), the symbols of the
    expression, ``builtins`` and ``range``."""
    import builtins
    from symox import Symbol
    if modules is None:
        try:
            _module_namespace('scipy')
        except ImportError:
            try:
                _module_namespace('numpy')
            except ImportError:
                modules = ['math', 'mpmath', 'sympy']
            else:
                modules = ['numpy']
        else:
            modules = ['numpy', 'scipy']
    if isinstance(modules, (dict, str)) or not hasattr(modules, '__iter__'):
        candidates = [modules]
    else:
        candidates = list(modules)
    namespace = {}
    for module in reversed(candidates):
        if isinstance(module, str):
            namespace.update(_module_namespace(module))
        elif isinstance(module, dict):
            namespace.update(module)
        elif hasattr(module, '__dict__'):
            namespace.update(module.__dict__)
        else:
            raise TypeError('Argument must be either a string, dict or module but it is: %s' % module)
    for expr in exprs:
        if hasattr(expr, 'atoms'):
            namespace.update({str(term): term for term in expr.atoms(Symbol)})
    namespace.update({'builtins': builtins, 'range': range, '__builtins__': builtins.__dict__})
    return namespace


# ----------------------------------------------------------------------
# SymPy's generated function (source, docstring, signature, result types)
# ----------------------------------------------------------------------

_lambdify_generated_counter = 1
_KEYWORDS = None


def _module_named(name, modlist):
    # SymPy's lambdify._module_present.
    if name in modlist:
        return True
    return any(getattr(m, '__name__', None) == name for m in modlist)


def _reference_flavor(modules, options):
    """The printer SymPy's lambdify selects ('math', 'numpy', 'scipy'), with
    the user functions of dict modules, or None for printers whose source is
    not reproduced natively (mpmath, sympy, cmath, cupy, jax, ...)."""
    if options.get('printer') is not None or options.get('cse') or options.get('dummify'):
        return None
    if set(options) - {'docstring_limit', 'use_imps'}:
        return None
    if modules is None:
        try:
            _module_namespace('scipy')
        except ImportError:
            try:
                _module_namespace('numpy')
            except ImportError:
                return None
            namespaces = ['numpy']
        else:
            namespaces = ['numpy', 'scipy']
    elif isinstance(modules, (dict, str)) or not hasattr(modules, '__iter__'):
        namespaces = [modules]
    else:
        namespaces = list(modules)
    for name in ('mpmath', 'cupy', 'jax', 'numexpr', 'tensorflow', 'torch', 'sympy', 'symox', 'cmath'):
        if _module_named(name, namespaces):
            return None
    if _module_named('scipy', namespaces):
        flavor = 'scipy'
    elif _module_named('numpy', namespaces):
        flavor = 'numpy'
    else:
        flavor = 'math'
    user_functions = {}
    for m in namespaces[::-1]:
        if isinstance(m, dict):
            for k in m:
                if not isinstance(k, str):
                    return None
                user_functions[k] = k
    return flavor, namespaces, user_functions


_REFERENCE_NODE_LIMIT = 20000


def _reference_plan(args, expr, modules, options, shown=None):
    """``(build, kind, math)``: a builder of SymPy's generated
    ``_lambdifygenerated`` function from the natively printed source, its
    result type code for float arguments and the math flavor; None when the
    generated source is not reproduced."""
    import keyword
    from .. import _core
    if not isinstance(expr, _S.Basic) or getattr(expr, 'is_Matrix', False):
        return None
    plan = _reference_flavor(modules, options)
    if plan is None:
        return None
    flavor, namespaces, user_functions = plan
    iterable_args = (args,) if isinstance(args, _S.Expr) else args
    try:
        iterable_args = list(iterable_args)
    except TypeError:
        return None
    names = []
    for arg in iterable_args:
        if not isinstance(arg, _S.Symbol) or type(arg).__name__ != 'Symbol':
            return None
        name = arg.name
        if not name or not name.isidentifier() or keyword.iskeyword(name):
            return None
        names.append(name)
    if len(set(names)) != len(names):
        return None
    # The printed source is a tree: shared subexpressions (a DAG such as 40
    # nested sin(e) + cos(e)) would be exponential, keep the native result.
    from .. import preorder_traversal
    for count, _ in enumerate(preorder_traversal(expr)):
        if count > _REFERENCE_NODE_LIMIT:
            return None
    # The native path raises for free symbols that are not arguments and for
    # undefined functions without an implementation: keep that behavior.
    if not expr.free_symbols <= set(iterable_args):
        return None
    from ..core.function import AppliedUndef
    if any(str(f.func) not in (user_functions or {}) for f in expr.atoms(AppliedUndef)):
        return None
    printed = _core.lambdify_source(expr, flavor, user_functions)
    if printed is None:
        return None
    body, imports = printed
    if not _names_resolved(body, names, imports, namespaces):
        return None
    kind = _core.lambdify_result_kind(expr, flavor)
    limit = options.get('docstring_limit', 1000)
    shown = expr if shown is None else shown

    def build():
        return _generated_function(names, body, imports, namespaces, expr, shown, limit)
    return build, kind, flavor == 'math'


def _names_resolved(body, names, imports, namespaces):
    """Every global name of the printed return expression is an argument,
    an import, a builtin or a name of the module namespaces (otherwise the
    printed form is not SymPy's, e.g. a held ``UnevaluatedExpr`` or a
    ``Catalan`` constant that SymPy prints as a Float)."""
    import builtins, keyword, re
    known = set(names) | {k for _, k in imports} | set(keyword.kwlist)
    for name in re.findall(r'(?<![\w.])[A-Za-z_]\w*(?!\w|\s*=(?!=))', re.sub(r"'[^']*'|\"[^\"]*\"", '', body)):
        if name in known or hasattr(builtins, name):
            continue
        for m in namespaces:
            space = _module_namespace(m) if isinstance(m, str) else m if isinstance(m, dict) else getattr(m, '__dict__', {})
            if name in space:
                break
        else:
            return False
    return True


def _attach_reference(native, args, expr, modules, options, shown=None):
    """Attach SymPy's generated function to a native ``CompiledFunction``
    (built on first need), or leave it unchanged."""
    if not hasattr(native, '_attach_reference'):
        return native
    plan = _reference_plan(args, expr, modules, options, shown)
    if plan is not None:
        native._attach_reference(*plan)
    return native


def _too_large_for_docstring(expr, limit):
    if limit is None:
        return False
    from .. import postorder_traversal
    i = 0
    for _ in postorder_traversal(expr):
        i += 1
        if i > limit:
            return True
    return False


def _generated_function(names, body, imports, namespaces, expr, shown, limit):
    """SymPy's lambdify assembly: the namespace, the module imports, the
    source registered in ``linecache`` and the docstring."""
    import builtins, linecache, textwrap, weakref
    global _lambdify_generated_counter
    namespace = {}
    for m in namespaces[::-1]:
        if isinstance(m, str):
            namespace.update(_module_namespace(m))
        elif isinstance(m, dict):
            namespace.update(m)
        elif hasattr(m, '__dict__'):
            namespace.update(m.__dict__)
        else:
            return None
    for term in expr.atoms(_S.Symbol):
        namespace[str(term)] = term
    funcstr = 'def _lambdifygenerated(%s):\n    return %s\n' % (', '.join(names), body)
    imp_mod_lines = []
    for mod, k in imports:
        if k not in namespace:
            ln = 'from %s import %s' % (mod, k)
            try:
                exec(ln, {}, namespace)
            except ImportError:
                ln = '%s = %s.%s' % (k, mod, k)
                try:
                    exec(ln, {}, namespace)
                except Exception:
                    return None
            imp_mod_lines.append(ln)
    namespace.update({'builtins': builtins, 'range': range})
    funclocals = {}
    filename = '<lambdifygenerated-%s>' % _lambdify_generated_counter
    _lambdify_generated_counter += 1
    exec(compile(funcstr, filename, 'exec'), namespace, funclocals)
    linecache.cache[filename] = (len(funcstr), None, funcstr.splitlines(True), filename)

    def cleanup_linecache(filename):
        def _cleanup():
            if filename in linecache.cache:
                del linecache.cache[filename]
        return _cleanup
    func = funclocals['_lambdifygenerated']
    weakref.finalize(func, cleanup_linecache(filename))
    sig = textwrap.fill('func({})'.format(', '.join(names)), subsequent_indent=' ' * 8)
    if _too_large_for_docstring(shown, limit):
        expr_str = "EXPRESSION REDACTED DUE TO LENGTH, (see lambdify's `docstring_limit`)"
        src_str = "SOURCE CODE REDACTED DUE TO LENGTH, (see lambdify's `docstring_limit`)"
    else:
        expr_str = str(shown)
        if len(expr_str) > 78:
            expr_str = textwrap.wrap(expr_str, 75)[0] + '...'
        src_str = funcstr
    func.__doc__ = (
        "Created with lambdify. Signature:\n\n"
        "{sig}\n\n"
        "Expression:\n\n"
        "{expr}\n\n"
        "Source code:\n\n"
        "{src}\n\n"
        "Imported modules:\n\n"
        "{imp_mods}"
        ).format(sig=sig, expr=expr_str, src=src_str, imp_mods='\n'.join(imp_mod_lines))
    return func
