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

def lambdastr(args,expr,printer=None,dummify=None):
    from ..printing.lambdarepr import LambdaPrinter
    roots=[args] if isinstance(args,(_S.Expr,str)) else list(args)
    mapping={};parameters=[]
    def descend(value,path):
        if isinstance(value,(tuple,list)) or getattr(value,'is_Tuple',False):
            for i,child in enumerate(value):descend(child,path+'['+str(i)+']')
        else:mapping[_S.sympify(value)]=path
    for i,arg in enumerate(roots):
        name=getattr(arg,'name',str(arg))
        if dummify or not name.isidentifier() or isinstance(arg,(tuple,list)):name='_arg'+str(i)
        parameters.append(name);descend(arg,name)
    settings={'symbol_names':mapping}
    emitter=LambdaPrinter(settings) if printer is None else printer(settings) if isinstance(printer,type) else printer
    if callable(emitter) and not hasattr(emitter,'doprint'):body=emitter(expr)
    else:
        if hasattr(emitter,'_settings'):emitter._settings.setdefault('symbol_names',{}).update(mapping)
        body=emitter.doprint(expr)
    return 'lambda '+', '.join(parameters)+': ('+body+')'

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
        stack=[value];result=[]
        while stack:
            node=stack.pop()
            if matches(node):result.append(node)
            else:stack.extend(getattr(node,'args',()))
        return list(dict.fromkeys(result))
    if isinstance(expr,(list,tuple)):
        plans=[_compile_implemented(args,e,modules,options) or _S._native_lambdify(args,e,modules=modules,**options) for e in expr]
        return lambda *values:type(expr)(plan(*values) for plan in plans)
    found=calls(expr)
    if not found:return None
    def plan(value):
        outer=calls(value)
        if not outer:return _S._native_lambdify(arguments,value,modules=modules,**options)
        temporary=[_S.Dummy('implemented') for _ in outer]
        body=_S.sympify(value).xreplace(dict(zip(outer,temporary)))
        evaluator=_S._native_lambdify(arguments+temporary,body,modules=modules,**options)
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
    arguments=[]
    def argument_plan(value):
        # SymPy unpacks an array argument like the nested list it iterates as.
        if sequence(value) or isinstance(value,NDimArray):return tuple(argument_plan(child) for child in value)
        arguments.append(_S.sympify(value));return None
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
               or _S._native_lambdify(parameters,expressions,modules=modules,**options))
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
