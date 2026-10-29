"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Dummy Symbol Wild symbols var')
from symox import Atom as _Atom, _host
_host_defined=_host._defined[0]


class Str(_Atom):
    """An immutable label which does not contribute free mathematical symbols.

    As in SymPy it is an ``Atom``: it can be stored in native expressions
    (``Eq(Y[1], Str('a'))``, ``FiniteSet(Str('a'), Str('b'))``) as a host atom."""
    # Created as a host atom class like a user Atom subclass (the metaclass
    # only does so outside the symox package), then given its public module.
    __module__='symox_host'
    __slots__=('name',)
    is_Atom=True
    is_Symbol=False
    def __new__(cls,name,**kwargs):
        if not isinstance(name,str):raise TypeError('name should be a string, not %s' % repr(type(name)))
        # Host paths are consulted once a host instance exists (see below).
        _host._defined[0]=True
        obj=_Atom.__new__(cls,**kwargs);object.__setattr__(obj,'name',name);return obj
    def __getnewargs__(self):return (self.name,)
    def _hashable_content(self):return (self.name,)
    def _sympystr(self,printer):return self.name
    def _sympyrepr(self,printer):return "%s(%s)" % (type(self).__name__,printer._print(self.name))
    def __str__(self):return self.name
    def __repr__(self):return self.name
    def subs(self,*args,**kwargs):return self
    def xreplace(self,rule):return rule.get(self,self)


Str.__module__=__name__
# Defining the class at import is not a user host class: the host paths
# stay disabled until an instance (or a user host class) exists.
_host._defined[0]=_host_defined


def uniquely_named_symbol(xname,exprs=(),compare=str,modify=None,**assumptions):
    import re
    from symox import Symbol
    if isinstance(xname,(list,tuple)):xname=xname[0]
    expressions=[exprs] if hasattr(exprs,'free_symbols') else exprs
    names={compare(s.name) for expr in expressions for s in getattr(expr,'free_symbols',())}
    name=str(getattr(xname,'name',xname))
    def increment(text):
        match=re.search(r'(\d+)$',text)
        return text[:match.start()]+str(int(match.group())+1) if match else text+'0'
    modifier=modify or increment
    while compare(name) in names:name=modifier(name)
    return Symbol(name,**assumptions)


def disambiguate(*iter):
    # Verbatim port of the reference disambiguate.
    from symox import Tuple,Symbol,ordered
    new_iter = Tuple(*iter)
    key = lambda x:tuple(sorted(x.assumptions0.items()))
    syms = ordered(new_iter.free_symbols, keys=key)
    mapping = {}
    for s in syms:
        mapping.setdefault(str(s).lstrip('_'), []).append(s)
    reps = {}
    for k in mapping:
        mapk0 = Symbol("%s" % (k), **mapping[k][0].assumptions0)
        if mapping[k][0] != mapk0:
            reps[mapping[k][0]] = mapk0
        skip = 0
        for i in range(1, len(mapping[k])):
            while True:
                name = "%s_%i" % (k, i + skip)
                if name not in mapping:
                    break
                skip += 1
            ki = mapping[k][i]
            reps[ki] = Symbol(name, **ki.assumptions0)
    return new_iter.xreplace(reps)
__all__ += ['Str','uniquely_named_symbol','disambiguate']
