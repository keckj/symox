"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'lambdify')

import symox as _S
_implemented={}

def implemented_function(symfunc,implementation):
    # SymPy creates a new undefined function holding the implementation,
    # with the assumptions of a given undefined function.
    if isinstance(symfunc,_S.Function):
        symfunc=_S.Function(symfunc.name,nargs=symfunc._nargs,**symfunc._assumptions)
    elif isinstance(symfunc,str):symfunc=_S.Function(symfunc)
    else:
        raise ValueError('\nsymfunc should be either a string or an UndefinedFunction instance.')
    symfunc._imp_=implementation;_implemented[symfunc.name]=implementation
    # evalf consults the implementation (SymPy's Function._eval_evalf).
    register=getattr(_S._core,'_register_implemented_evalf',None)
    if register is not None:register()
    return symfunc

def _imp_evalf(expr,prec):
    """SymPy's ``Function._eval_evalf`` fallback of an applied undefined
    function: ``Float(imp(*[a.evalf(prec) for a in args]), prec)``, or None
    when there is no implementation or its value is not a real number."""
    imp=_implemented.get(getattr(expr.func,'name',None))
    if imp is None:return None
    try:
        return _S.Float(imp(*[a.evalf(prec) for a in expr.args]),prec)
    except (TypeError,ValueError):
        return None

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
        if _is_codegen_node(node):
            return _codegen_implementation(node,modules)
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
    import functools
    return eval(source, {'numpy': numpy, 'functools': functools})


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
        # sympy.codegen nodes MpmathPrinter prints under their class name
        # (minimum, fma, Cbrt, ...): without a namespace entry the generated
        # function raises NameError when called, as SymPy's does.
        self.codegen=set()
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
        if _is_codegen_node(e):
            if name=='Sqrt':return self.call('sqrt',[self(args[0])])
            self.codegen.add(name)
            return self.call(name,[self(a) for a in args])
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
    def host(value):return is_host_function_class(type(value)) and not _is_codegen_node(value)
    if host(expr) or (isinstance(expr,(list,tuple)) and any(host(v) for v in expr)):
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
            elif name not in printer.codegen:return None
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
    # SymPy's SciPyPrinter adds ``from numpy import logaddexp`` (and
    # logaddexp2) through its module_imports.
    'scipy': ({'I': 1j}, ('import scipy; import numpy; from scipy.special import *',
                          'from numpy import logaddexp, logaddexp2'),
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


_REFERENCE_FLAVORS = {}


def _reference_flavor(modules, options):
    """The printer SymPy's lambdify selects ('math', 'numpy', 'scipy'), with
    the user functions of dict modules, or None for printers whose source is
    not reproduced natively (mpmath, sympy, cmath, cupy, jax, ...)."""
    if options and (options.get('printer') is not None or options.get('cse') or options.get('dummify')
                    or set(options) - {'docstring_limit', 'use_imps'}):
        return None
    # A module name (or the default) always selects the same printer.
    cacheable = modules is None or type(modules) is str
    if cacheable and modules in _REFERENCE_FLAVORS:
        return _REFERENCE_FLAVORS[modules]
    plan = _select_reference_flavor(modules)
    if cacheable:
        _REFERENCE_FLAVORS[modules] = plan
    return plan


def _select_reference_flavor(modules):
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


def _reference_plan(args, expr, modules, options, shown=None, introspection=False, codegen=False):
    """``(build, kind, flavor[, body, names, resolve[, evaluate]])``: a builder of
    SymPy's generated ``_lambdifygenerated`` function from the natively
    printed source, its result type code for float arguments, the printer
    flavor ('math', 'numpy' or 'scipy') and, without user functions, the
    printed return expression, the argument names and the resolver of its
    global names (for its native evaluation); None when the generated source
    is not reproduced."""
    from .. import _core
    if not isinstance(expr, _core.Expr) and not isinstance(expr, _S.Basic) or getattr(expr, 'is_Matrix', False):
        return None
    plan = _reference_flavor(modules, options)
    if plan is None:
        return None
    flavor, namespaces, user_functions = plan
    iterable_args = (args,) if isinstance(args, _core.Expr) or isinstance(args, _S.Expr) else args
    try:
        iterable_args = list(iterable_args)
    except TypeError:
        return None
    names = []
    for arg in iterable_args:
        if type(arg).__name__ != 'Symbol' or not isinstance(arg, _S.Symbol):
            return None
        name = arg.name
        if not name or not name.isidentifier() or _iskeyword(name):
            return None
        names.append(name)
    if len(set(names)) != len(names):
        return None
    # The printed source is a tree: shared subexpressions (a DAG such as 40
    # nested sin(e) + cos(e)) would be exponential, keep the native result
    # (the native printer counts the tree up to the limit).
    printed = _core.lambdify_reference(expr, flavor, user_functions, names, _REFERENCE_NODE_LIMIT)
    if printed is None:
        return None
    body, imports, kind, undefined, unknown = printed
    # The native path raises for free symbols that are not arguments and for
    # undefined functions without an implementation: keep that behavior.
    if not expr.free_symbols <= set(iterable_args):
        return None
    if undefined:
        from ..core.function import AppliedUndef
        if any(str(f.func) not in (user_functions or {}) for f in expr.atoms(AppliedUndef)):
            return None
    undefined = bool(unknown) and not _names_resolved(unknown, namespaces)
    if undefined and not ((introspection or codegen) and _undefined_functions(unknown, expr)):
        return None
    limit = options.get('docstring_limit', 1000)
    shown = expr if shown is None else shown
    if undefined:
        # SymPy prints a function it does not support under its own name
        # and its generated function raises NameError; the native function
        # computes the values and keeps the generated one for introspection
        # and for the arguments it does not evaluate (complex, objects).
        def build_undefined():
            return _generated_function(names, body, imports, namespaces, expr, shown, limit)
        return build_undefined, kind, flavor, None, None, None, False

    def build():
        return _generated_function(names, body, imports, namespaces, expr, shown, limit)
    if user_functions:
        return build, kind, flavor

    def resolve(used):
        return _resolve_globals(used, imports, namespaces)
    return build, kind, flavor, body, names, resolve


# Functions that SymPy 1.14's NumPy and math lambdify printers print under
# their own name without a definition (their generated function raises
# NameError), checked against SymPy with and without SciPy. Other functions
# are rewritten (beta, binomial, catalan, erf2, RisingFactorial, ...) or
# refused (fresnels, lowergamma, ...) by SymPy and keep the previous path.
_UNDEFINED_PRINTED = frozenset((
    'Chi', 'Ci', 'Ei', 'LambertW', 'Li', 'Shi', 'Si', 'airyai', 'airyaiprime', 'airybi',
    'arg', 'besseli', 'besselj', 'besselk', 'bessely', 'chebyshevt', 'conjugate',
    'dirichlet_eta', 'elliptic_e', 'elliptic_k', 'erfcinv', 'erfi', 'erfinv', 'expint',
    'factorial', 'fibonacci', 'hermite', 'im', 'legendre', 'li', 'lucas', 'polygamma',
    're', 'zeta'))


def _undefined_functions(unknown, expr):
    """Whether every unresolved name of the printed source is a function
    applied in ``expr`` that SymPy prints under its own name without a
    definition (``_UNDEFINED_PRINTED``)."""
    from symox import Function
    functions = expr.atoms(Function)
    applied = {type(f).__name__ for f in functions}
    # The ``sympy.codegen`` classes the native printer prints under their
    # class name are the ones SymPy's printers do not know (fma, Cbrt,
    # cosm1 without SciPy, minimum with the math printer, ...).
    codegen = {type(f).__name__ for f in functions if _is_codegen_node(f)}
    return all(name in codegen or (name in applied and name in _UNDEFINED_PRINTED) for name in unknown)


# Modules of the ``sympy.codegen`` function classes whose lambdify printing
# follows SymPy's printers (class-name dispatch).
_CODEGEN_NODE_MODULES = ('symox.codegen.numpy_nodes', 'symox.codegen.scipy_nodes',
                         'symox.codegen.cfunctions')


def _is_codegen_node(node):
    return type(node).__module__ in _CODEGEN_NODE_MODULES


def _codegen_flavor(modules):
    """The SymPy printer lambdify selects for ``modules`` ('mpmath', 'scipy',
    'numpy' or 'math'), or None for the other printers."""
    if modules is None:
        try:
            _module_namespace('scipy')
            return 'scipy'
        except ImportError:
            try:
                _module_namespace('numpy')
                return 'numpy'
            except ImportError:
                return 'mpmath'
    candidates = [modules] if isinstance(modules, (dict, str)) or not hasattr(modules, '__iter__') else list(modules)
    for name in ('mpmath', 'scipy', 'numpy'):
        if _module_named(name, candidates):
            return name
    for name in ('cupy', 'jax', 'numexpr', 'tensorflow', 'torch', 'sympy', 'symox', 'cmath'):
        if _module_named(name, candidates):
            return None
    return 'math'


# The ``sympy.codegen`` functions each SymPy printer prints as a module
# function (``numpy.expm1``, ``math.hypot``, ``mpmath.log10``); the others it
# prints under their class name (``fma(x, y, z)``).
_CODEGEN_KNOWN = {
    'math': frozenset(('expm1', 'log1p', 'log2', 'log10', 'hypot', 'isnan', 'isinf', 'Sqrt')),
    'mpmath': frozenset(('expm1', 'log1p', 'log2', 'log10', 'hypot', 'isnan', 'isinf', 'Sqrt')),
    'numpy': frozenset(('expm1', 'log1p', 'log2', 'log10', 'hypot', 'isnan', 'isinf', 'Sqrt',
                        'exp2', 'logaddexp', 'logaddexp2', 'minimum', 'maximum')),
}
_CODEGEN_KNOWN['scipy'] = _CODEGEN_KNOWN['numpy'] | {'cosm1', 'powm1'}


def _codegen_implementation(node, modules):
    """The callable SymPy's generated function applies for a
    ``sympy.codegen`` node: the printed name looked up in the lambdify
    namespace (``reduce`` over the binary ufunc for numpy ``minimum`` and
    ``maximum``); a name the namespace lacks raises NameError when called,
    like SymPy's generated function. None for printers not reproduced."""
    import functools, importlib
    flavor = _codegen_flavor(modules)
    if flavor is None:
        return None
    name = type(node).__name__
    printed = 'sqrt' if name == 'Sqrt' else name
    known = name in _CODEGEN_KNOWN[flavor]
    namespace = _lambdify_namespace(modules, ())
    module = {'math': 'math', 'mpmath': 'mpmath', 'numpy': 'numpy', 'scipy': 'numpy'}[flavor]
    if flavor == 'scipy' and name in ('cosm1', 'powm1'):
        module = 'scipy.special'

    def resolve(short):
        found = namespace.get(short)
        if found is not None:
            return found
        if known:
            # The printer's ``from module import name`` line.
            try:
                return getattr(importlib.import_module(module), short)
            except (ImportError, AttributeError):
                pass

        def undefined(*values):
            raise NameError("name '%s' is not defined" % short)
        return undefined
    if flavor in ('numpy', 'scipy') and name in ('minimum', 'maximum'):
        if len(node.args) == 1:
            return lambda value: value
        ufunc = resolve(name)
        return lambda *values: functools.reduce(ufunc, list(values))
    if flavor == 'mpmath' and name == 'log2':
        log = resolve('log')
        return lambda value: log(value)/log(2)
    return resolve(printed)


def _codegen_nodes(expr):
    """The ``sympy.codegen`` function applications of ``expr``."""
    import sys
    if not any(name in sys.modules for name in _CODEGEN_NODE_MODULES):
        return []
    from symox import Function
    try:
        return [f for f in expr.atoms(Function) if _is_codegen_node(f)]
    except (AttributeError, TypeError):
        return []


def _compile_codegen(args, expr, modules, options, shown):
    """SymPy's generated function for an expression with ``sympy.codegen``
    nodes (numpy_nodes, scipy_nodes, cfunctions) or an ``amin``/``amax``
    reduction, or None to keep the general path. SymPy prints these through
    its printers' class-name dispatch: ``minimum(x, y, z)`` as
    ``reduce(minimum, [x, y, z])``, ``expm1`` as ``numpy.expm1``, and the
    nodes a printer does not know under their own name (the generated
    function raises NameError when the namespace has no such name)."""
    import sys
    nodes = sys.modules.get('symox.codegen.numpy_nodes')
    if nodes is not None and isinstance(expr, nodes._AxisReduction):
        return _compile_axis_reduction(args, expr, modules, options, shown)
    from .. import _core
    if not (isinstance(expr, _core.Expr) or isinstance(expr, _S.Basic)) or _mpmath_only(modules):
        return None
    if not _codegen_nodes(expr):
        return None
    plan = _reference_plan(args, expr, modules, options, shown, codegen=True)
    return plan[0]() if plan is not None else None


def _argument_names(args):
    """The parameter names of SymPy's generated function for plain symbol
    arguments, or None."""
    from .. import _core
    iterable_args = (args,) if isinstance(args, _core.Expr) or isinstance(args, _S.Expr) else args
    try:
        iterable_args = list(iterable_args)
    except TypeError:
        return None
    names = []
    for arg in iterable_args:
        if type(arg).__name__ != 'Symbol' or not isinstance(arg, _S.Symbol):
            return None
        name = arg.name
        if not name or not name.isidentifier() or _iskeyword(name):
            return None
        names.append(name)
    return names if len(set(names)) == len(names) else None


def _axis_reduction_call(args, expr, modules, options):
    """``amin``/``amax`` when the generated source is not reproduced (cse,
    dummify, nested arguments, ...): the array lambdified on its own, then
    the namespace's ``amin``/``amax`` (numpy's by default) along ``axis``."""
    from ..codegen.ast import none
    head = type(expr).__name__
    reduction = _lambdify_namespace(modules, ()).get(head)
    if reduction is None:
        import numpy
        reduction = getattr(numpy, head)
    array = expr.array
    if isinstance(array, tuple):
        array = list(array)
    inner = _S.lambdify(args, array, modules=modules, **options)
    axis = expr.axis
    if axis is None or axis is none:
        axis = None
    elif getattr(axis, 'is_Integer', False):
        axis = int(axis)
    elif getattr(axis, 'is_Tuple', False) and all(getattr(a, 'is_Integer', False) for a in axis):
        axis = tuple(int(a) for a in axis)
    else:
        return None

    def call(*values, **kwargs):
        return reduction(inner(*values, **kwargs), axis=axis)
    return call


def _compile_axis_reduction(args, expr, modules, options, shown):
    """``amin``/``amax``: NumPyPrinter's ``numpy.amin(array, axis=...)``;
    the math and mpmath printers refuse the node."""
    flavor = _codegen_flavor(modules)
    if flavor == 'mpmath':
        from ..printing.pycode import MpmathPrinter
        MpmathPrinter({'allow_unknown_functions': True}).doprint(expr)
        return None
    if flavor == 'math':
        from ..printing.pycode import PythonCodePrinter
        PythonCodePrinter({'allow_unknown_functions': True}).doprint(expr)
        return None
    if flavor is None:
        return None
    plan = _reference_flavor(modules, options)
    names = _argument_names(args)
    if plan is None or names is None:
        return _axis_reduction_call(args, expr, modules, options)
    flavor, namespaces, user_functions = plan
    if names is None:
        return None
    from .. import _core
    from ..codegen.ast import none

    def printed(value):
        printed = _core.lambdify_reference(value, flavor, user_functions, names, _REFERENCE_NODE_LIMIT)
        if printed is None:
            raise _Unprintable(value)
        body, imports, _, _, unknown = printed
        if unknown and not _names_resolved(unknown, namespaces):
            raise _Unprintable(value)
        return body, imports
    imports = []
    try:
        array = expr.array
        if isinstance(array, (list, tuple)):
            parts = []
            for value in array:
                body, more = printed(_S.sympify(value))
                parts.append(body)
                imports += more
            body = '[%s]' % ', '.join(parts) if isinstance(array, list) else '(%s)' % (
                parts[0] + ',' if len(parts) == 1 else ', '.join(parts))
        else:
            body, imports = printed(_S.sympify(array))
        axis = expr.axis
        if axis is None or axis is none:
            axis_text = 'None'
        else:
            axis_text, more = printed(_S.sympify(axis))
            imports += more
    except (_Unprintable, TypeError, ValueError, AttributeError):
        return _axis_reduction_call(args, expr, modules, options)
    head = type(expr).__name__
    imports = list(dict.fromkeys(imports + [('numpy', head)]))
    source = '%s(%s, axis=%s)' % (head, body, axis_text)
    atoms = _S.Tuple(*(array if isinstance(array, (list, tuple)) else [array]))
    return _generated_function(names, source, imports, namespaces, atoms, expr, options.get('docstring_limit', 1000))


def _resolve_globals(used, imports, namespaces):
    """The objects the global names ``used`` by the generated function are
    bound to (the namespace of ``_generated_function``), or None."""
    import builtins, importlib
    out = []
    for name in used:
        if name in ('builtins', 'range'):
            return None
        for m in namespaces:
            space = _module_namespace(m) if isinstance(m, str) else m if isinstance(m, dict) else getattr(m, '__dict__', None)
            if space is None:
                return None
            if name in space:
                out.append(space[name])
                break
        else:
            mod = next((mod for mod, k in imports if k == name), None)
            if mod is not None:
                try:
                    out.append(getattr(importlib.import_module(mod), name))
                except (ImportError, AttributeError):
                    return None
            elif hasattr(builtins, name):
                out.append(getattr(builtins, name))
            else:
                return None
    return out


def _iskeyword(name):
    import keyword
    return keyword.iskeyword(name)


def _names_resolved(unknown, namespaces):
    """Every global name of the printed return expression that is not an
    argument, an import or a keyword is a builtin or a name of the module
    namespaces (otherwise the printed form is not SymPy's, e.g. a held
    ``UnevaluatedExpr`` or a ``Catalan`` constant that SymPy prints as a
    Float)."""
    import builtins
    for name in unknown:
        if hasattr(builtins, name):
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
    plan = _reference_plan(args, expr, modules, options, shown, introspection=True)
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
