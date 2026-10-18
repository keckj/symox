"""Public utility and assumption interfaces over the native expression engine."""
from contextlib import redirect_stdout
from functools import wraps
from io import StringIO
from itertools import permutations
import importlib
import inspect
import math
import sys
import textwrap
import timeit


def capture(func):
    output = StringIO()
    with redirect_stdout(output):
        func()
    return output.getvalue()


def filldedent(s,w=70,**kwargs):
    return '\n'+textwrap.fill(textwrap.dedent(str(s)).strip(),width=w,**kwargs)


def public(obj):
    namespace = vars(sys.modules[obj.__module__])
    names = namespace.setdefault('__all__',[])
    if obj.__name__ not in names:
        if not hasattr(names,'append'):
            names = namespace['__all__'] = list(names)
        names.append(obj.__name__)
    return obj


def _threaded(func,use_add):
    @wraps(func)
    def operation(expr,*args,**kwargs):
        from . import Expr,Add,sympify
        if hasattr(expr,'applyfunc') and not isinstance(expr,Expr):
            return expr.applyfunc(lambda e:func(e,*args,**kwargs))
        if isinstance(expr,(list,tuple,set,frozenset)):
            return type(expr)(func(e,*args,**kwargs) for e in expr)
        if isinstance(expr,dict):
            return {key:func(value,*args,**kwargs) for key,value in expr.items()}
        expr = sympify(expr)
        if expr.is_Tuple:
            from . import Tuple
            return Tuple(*(func(e,*args,**kwargs) for e in expr.args))
        if use_add and expr.is_Add:
            return Add(*(func(e,*args,**kwargs) for e in expr.args))
        if expr.is_Relational:
            return expr.func(func(expr.lhs,*args,**kwargs),func(expr.rhs,*args,**kwargs))
        return func(expr,*args,**kwargs)
    return operation


def threaded(func):
    return _threaded(func,True)


def xthreaded(func):
    return _threaded(func,False)


class vectorize:
    def __init__(self,*mdargs):
        if any(not isinstance(a,(int,str)) for a in mdargs):
            raise TypeError('vectorized argument selectors must be integers or names')
        self.mdargs = mdargs

    def __call__(self,func):
        @wraps(func)
        def operation(*args,**kwargs):
            selectors = self.mdargs or tuple(range(len(args)))+tuple(kwargs)
            def apply(position,arguments,options):
                if position == len(selectors):
                    return func(*arguments,**options)
                key = selectors[position]
                if isinstance(key,int):
                    if not -len(arguments) <= key < len(arguments):
                        return apply(position+1,arguments,options)
                    value = arguments[key]
                else:
                    if key not in options:
                        return apply(position+1,arguments,options)
                    value = options[key]
                def descend(item):
                    if getattr(item,'is_Tuple',False):
                        item = item.args
                    if isinstance(item,(list,tuple)) or type(item).__module__.startswith('numpy') and hasattr(item,'__iter__'):
                        return [descend(v) for v in item]
                    changed,updated = list(arguments),dict(options)
                    if isinstance(key,int):
                        changed[key] = item
                    else:
                        updated[key] = item
                    return apply(position+1,changed,updated)
                return descend(value)
            return apply(0,list(args),dict(kwargs))
        return operation


def timed(func,setup='pass',limit=None):
    if limit is not None and (not isinstance(limit,int) or limit < 1):
        raise ValueError('timing limit must be a positive integer')
    timer = timeit.Timer(func,setup=setup)
    number = 1
    while number < (limit or 10**9) and timer.timeit(number) < 0.2:
        number = min(number*10,limit or 10**9)
    elapsed = min(timer.repeat(3,number))/number
    order = min(max(-int(math.floor(math.log10(elapsed))//3),0),3) if elapsed > 0 else 3
    return number,elapsed,elapsed*10**(3*order),('s','ms','μs','ns')[order]


def arity(cls):
    counts = getattr(cls,'_nargs',None)
    if counts is not None:
        return counts[0] if len(counts) == 1 else tuple(counts)
    try:
        parameters = list(inspect.signature(cls,follow_wrapped=False).parameters.values())
    except (ValueError,TypeError):
        return None
    if not any(p.kind == p.VAR_POSITIONAL for p in parameters):
        positional = [p for p in parameters if p.kind in (p.POSITIONAL_ONLY,p.POSITIONAL_OR_KEYWORD)]
        required = sum(p.default is p.empty for p in positional)
        return required if required == len(positional) else tuple(range(required,len(positional)+1))
    head = getattr(cls,'_symox_head_name',None)
    if head:
        from . import _core,Dummy
        argument,accepted = Dummy('arity'),[]
        for count in range(9):
            try:_core._apply(head,*([argument]*count),evaluate=False)
            except (TypeError,ValueError):continue
            accepted.append(count)
        if accepted and accepted != list(range(9)):
            return accepted[0] if len(accepted) == 1 else tuple(accepted)
    return None


def comp(z1,z2,tol=None):
    from . import sympify,Float,Abs
    if isinstance(z2,str) or tol == '':
        return str(z1) == str(z2)
    a,b = sympify(z1),sympify(z2)
    if not a.is_number or not b.is_number:
        raise ValueError('comparison requires numeric expressions')
    if tol is None:
        floats = a.atoms(Float) | b.atoms(Float)
        if not floats:
            return a == b
        bits = min(f._prec for f in floats)
        tol = Float(5)*2**(-bits+1)
    tol = sympify(tol)
    if tol < 0:
        raise ValueError('tolerance must be nonnegative')
    error = Abs(a-b)
    if b != 0 and Abs(a) > 1:
        error /= Abs(a)
    return bool(error.evalf(30) <= tol.evalf(30))


def failing_assumptions(expr,**assumptions):
    from . import sympify,ask,Q
    expr = sympify(expr)
    failures = {}
    for name,wanted in assumptions.items():
        if wanted not in (True,False,None):
            raise ValueError('assumption values must be True, False or None')
        predicate = getattr(Q,name,None)
        if predicate is None:
            raise ValueError('unknown assumption: '+name)
        actual = ask(predicate(expr))
        if actual != wanted:
            failures[name] = actual
    return failures


def check_assumptions(expr,against=None,**assume):
    if against is not None:
        if assume:
            raise ValueError('use either against or explicit assumptions')
        from . import sympify
        assume = dict(sympify(against).assumptions0)
    failures = failing_assumptions(expr,**assume)
    if any(value is not None for value in failures.values()):
        return False
    return None if failures else True


_handlers = {}


def register_handler(key,handler):
    from . import Q,Predicate
    name = key if isinstance(key,str) else key.name
    if isinstance(handler,str):
        module,attribute = handler.rsplit('.',1)
        handler = getattr(importlib.import_module(module),attribute)
    if not callable(handler):
        raise TypeError('assumption handler must be a callable or handler class')
    _handlers.setdefault(name,[]).append(handler)
    if not hasattr(Q,name):
        setattr(Q,name,Predicate(name))
    return getattr(Q,name)


def remove_handler(key,handler):
    name = key if isinstance(key,str) else key.name
    if isinstance(handler,str):
        module,attribute = handler.rsplit('.',1)
        handler = getattr(importlib.import_module(module),attribute)
    values = _handlers.get(name,[])
    if handler in values:
        values.remove(handler)


def ask(proposition,assumptions=True,context=None):
    from .assumptions.assume import global_assumptions
    from . import _core
    result = _core._ask(proposition,assumptions,global_assumptions if context is None else context)
    if result is not None:
        return result
    if getattr(proposition,'class_name','') not in ('AppliedPredicate','AppliedBinaryRelation'):
        return result
    predicate,*arguments = proposition.args
    values = []
    for handler in _handlers.get(predicate.name,()):
        if inspect.isclass(handler):
            if len(arguments) != 1:
                continue
            arg = arguments[0]
            method = getattr(handler,arg.class_name,None) or getattr(handler,'Expr',None) or getattr(handler,'Basic',None)
            value = None if method is None else method(arg,assumptions)
        else:
            value = handler(*arguments,assumptions=assumptions)
        if value is not None:
            values.append(bool(value))
    if len(set(values)) > 1:
        raise ValueError('inconsistent custom assumption handlers')
    return values[0] if values else None


def bool_map(bool1,bool2):
    from . import simplify_logic,satisfiable,Xor
    a,b = simplify_logic(bool1),simplify_logic(bool2)
    left,right = sorted(a.free_symbols,key=str),sorted(b.free_symbols,key=str)
    if len(left) != len(right):
        return False
    for ordering in permutations(right):
        mapping = dict(zip(left,ordering))
        if satisfiable(Xor(a.xreplace(mapping),b)) is False:
            return a,mapping
    return False


def expand_multinomial(expr,deep=True):
    from . import sympify
    return sympify(expr).expand(deep=deep,mul=False,multinomial=True,power_base=False,
                                power_exp=False,log=False,basic=False)
