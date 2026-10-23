"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Dummy Symbol Wild symbols var')


class Str:
    """An immutable label which does not contribute free mathematical symbols."""
    __slots__=('name',)
    def __init__(self,name):object.__setattr__(self,'name',str(name))
    def __setattr__(self,name,value):raise AttributeError('string atoms are immutable')
    args=()
    free_symbols=frozenset()
    is_Atom=True
    def __str__(self):return self.name
    __repr__=__str__
    def __hash__(self):return hash((type(self),self.name))
    def __eq__(self,other):return isinstance(other,Str) and self.name==other.name
    def subs(self,*args,**kwargs):return self
    def xreplace(self,rule):return rule.get(self,self)


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
