"""Newton iterations as executable structured program nodes."""
import symox as S
from .ast import (Scope,CodeBlock,Variable,Declaration,Assignment,While,Return,
 FunctionDefinition,real,integer,Print)


def newtons_method(expr,wrt,atol=1e-12,delta=None,*,rtol=4e-16,debug=False,itermax=None,counter=None,
                   delta_fn=lambda e,x:-e/S.diff(e,x),cse=False,handle_nan=None,bounds=None):
    expr=S.sympify(expr);wrt=S.sympify(wrt);delta=S.Symbol('delta') if delta is None else S.sympify(delta)
    declarations=[Declaration(Variable(delta,real,S.oo))];condition=S.Abs(delta)>S.sympify(atol)+S.sympify(rtol)*S.Abs(wrt)
    update=delta_fn(expr,wrt);body=[]
    if cse:
        replacements,reduced=S.cse(update)
        for symbol,value in replacements:body.append(Declaration(Variable(symbol,real,value)))
        update=reduced[0]
    body.append(Assignment(delta,update));candidate=wrt+delta
    if bounds is not None:
        low,high=map(S.sympify,bounds);candidate=S.Min(S.Max(candidate,low),high)
    body.append(Assignment(wrt,candidate))
    if debug:body.append(Print((wrt,delta)))
    if itermax is not None:
        counter=S.Symbol('counter',integer=True) if counter is None else S.sympify(counter)
        declarations.append(Declaration(Variable(counter,integer,0)))
        condition=S.And(condition,counter<S.sympify(itermax));body.append(Assignment(counter,counter+1))
    if handle_nan is not None:
        # A NaN delta fails the ordered self-comparison; run the supplied
        # recovery block after the loop through a one-iteration conditional.
        recovery=CodeBlock(*handle_nan) if isinstance(handle_nan,(tuple,list)) else handle_nan
        declarations.append(While(S.Ne(delta,delta,evaluate=False),CodeBlock(recovery,Assignment(delta,0))))
    loop=While(condition,CodeBlock(*body))
    if handle_nan is not None:return Scope(CodeBlock(*declarations[:-1],loop,declarations[-1]))
    return Scope(CodeBlock(*declarations,loop))


def newtons_method_function(expr,wrt,params=None,func_name='newton',attrs=(),*,delta=None,**kwargs):
    expr=S.sympify(expr);params=list(params) if params is not None else [wrt]+sorted(expr.free_symbols-{wrt},key=str)
    missing=expr.free_symbols-set(params)
    if missing:raise ValueError('Newton function parameters omit '+', '.join(map(str,missing)))
    scope=newtons_method(expr,wrt,delta=delta,**kwargs)
    return FunctionDefinition(real,func_name,tuple(Variable(p,real) for p in params),CodeBlock(scope,Return(wrt)),attrs)
__all__=['newtons_method','newtons_method_function']
