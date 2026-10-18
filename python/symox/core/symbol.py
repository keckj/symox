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


def disambiguate(*expressions):
    from symox import sympify,Tuple,Symbol
    expressions=[sympify(expr) for expr in expressions]
    symbols=set().union(*(expr.free_symbols for expr in expressions));groups={}
    for symbol in sorted(symbols,key=lambda s:(s.name,repr(sorted(s.assumptions0.items())),getattr(s,'dummy_index',-1))):groups.setdefault(symbol.name.lstrip('_'),[]).append(symbol)
    occupied=set(groups);mapping={}
    for name,group in groups.items():
        for i,symbol in enumerate(group):
            candidate=name if i==0 else name+'_'+str(i)
            while i and candidate in occupied:i+=1;candidate=name+'_'+str(i)
            occupied.add(candidate);mapping[symbol]=Symbol(candidate,**symbol.assumptions0)
    return Tuple(*(expr.xreplace(mapping) for expr in expressions))
__all__ += ['Str','uniquely_named_symbol','disambiguate']
