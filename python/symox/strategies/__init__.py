"""Composable deterministic expression rewrite strategies."""
from functools import wraps


def identity(expr):return expr


def chain(*rules):
    def apply(expr):
        for rule in rules:expr = rule(expr)
        return expr
    return apply


def do_one(*rules):
    def apply(expr):
        for rule in rules:
            result = rule(expr)
            if result != expr:return result
        return expr
    return apply


def exhaust(rule):
    def apply(expr):
        seen = set()
        while True:
            try:
                if expr in seen:raise ValueError('rewrite strategy entered a cycle')
                seen.add(expr)
            except TypeError:pass
            result = rule(expr)
            if result == expr:return result
            expr = result
    return apply


def condition(predicate,rule):return lambda expr:rule(expr) if predicate(expr) else expr


def null_safe(rule):
    def apply(expr):
        result = rule(expr)
        return expr if result is None else result
    return apply


def tryit(rule,exception):
    def apply(expr):
        try:return rule(expr)
        except exception:return expr
    return apply


def minimize(*rules,**kwargs):
    objective = kwargs.pop('objective',lambda expr:expr)
    if kwargs:raise TypeError('unknown strategy options')
    return lambda expr:min((rule(expr) for rule in rules),key=objective)


def debug(rule,file=None):
    def apply(expr):
        import sys
        result = rule(expr)
        if result != expr:print('%s: %s -> %s' % (getattr(rule,'__name__','rule'),expr,result),file=file or sys.stdout)
        return result
    return apply


def _rebuild(expr,args):
    if isinstance(expr,tuple):return tuple(args)
    if isinstance(expr,list):return list(args)
    return expr.func(*args)


def rm_id(isid,new=None):
    def apply(expr):
        args = getattr(expr,'args',())
        remaining = [a for a in args if not isid(a)]
        if len(remaining) == len(args):return expr
        if not remaining:remaining = list(args[:1])
        return (new or _rebuild)(expr,remaining)
    return apply


def glom(key,count,combine):
    def apply(expr):
        groups = {}
        for arg in expr.args:groups.setdefault(key(arg),[]).append(arg)
        return _rebuild(expr,[combine(sum(count(a) for a in group),k) for k,group in groups.items()])
    return apply


def flatten(expr):
    return _rebuild(expr,[child for arg in expr.args for child in (arg.args if type(arg) is type(expr) and getattr(arg,'func',None)==expr.func else (arg,))])


def unpack(expr):return expr.args[0] if len(expr.args) == 1 else expr


def sort(key):return lambda expr:_rebuild(expr,sorted(expr.args,key=key))


def canon(*rules):return exhaust(do_one(*rules))


__all__ = ['identity','chain','do_one','exhaust','condition','null_safe','tryit','minimize','debug',
           'rm_id','glom','flatten','unpack','sort','canon']
