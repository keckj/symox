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

def _compile_implemented(args,expr,modules,options):
    if not _implemented:return None
    arguments=[args] if isinstance(args,_S.Expr) else list(args)
    def matches(node):return (getattr(node,'is_UndefinedFunction',False) and getattr(node,'class_name','') in _implemented)
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
        subplans=[( _implemented[node.class_name],[plan(a) for a in node.args]) for node in outer]
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
