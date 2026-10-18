"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox._public_tools', 'check_assumptions failing_assumptions')

# Expose native assumption properties through the standard introspection API.
_assume_defined = frozenset(('algebraic antihermitian commutative complex composite '
    'even extended_negative extended_nonnegative extended_nonpositive '
    'extended_nonzero extended_positive extended_real finite hermitian imaginary '
    'infinite integer irrational negative noninteger nonnegative nonpositive '
    'nonzero odd polar positive prime rational real transcendental zero').split())


def assumptions(expr, _check=None):
    from symox import sympify
    expr = sympify(expr)
    if expr.is_Symbol:
        facts = expr.assumptions0
        return facts if _check is None else {key:facts[key] for key in _check if key in facts}
    names = _assume_defined if _check is None else _check
    result = {}
    for name in names:
        value = getattr(expr, 'is_'+name, None)
        if value is not None:result[name] = value
    return result


def common_assumptions(exprs, check=None):
    facts = [assumptions(expr, _check=check) for expr in exprs]
    if not facts:return {}
    return {name:value for name,value in facts[0].items()
            if all(name in other and other[name] == value for other in facts[1:])}

__all__ += ['assumptions', 'common_assumptions']


class StdFactKB(dict):
    def __init__(self,facts=None):
        self._generator=dict(getattr(facts,'generator',facts) or {})
        super().__init__();self.deduce_all_facts(self._generator)
    @property
    def generator(self):return dict(self._generator)
    def deduce_all_facts(self,facts):
        from symox import Symbol
        supplied=dict(facts);known={key:value for key,value in self.items() if value is not None}
        for key,value in supplied.items():
            if key in known and value is not None and known[key]!=bool(value):raise ValueError('inconsistent assumption '+key)
            if value is not None:known[key]=bool(value)
        if known:
            inferred=Symbol('_assumption_probe',**known).assumptions0
            dict.update(self,inferred)
        for key,value in supplied.items():
            if key not in self:dict.__setitem__(self,key,value)
    def copy(self):
        result=type(self)(self);result._generator=dict(self._generator);return result


def as_property(fact):return 'is_'+fact

def make_property(fact):
    def getter(self):
        knowledge=getattr(self,'_assumptions',None)
        if knowledge is None or knowledge is getattr(type(self),'default_assumptions',None):
            knowledge=getattr(type(self),'default_assumptions',StdFactKB()).copy();self._assumptions=knowledge
        if fact not in knowledge:
            handler=getattr(self,'_eval_is_'+fact,None);value=handler() if handler else None
            knowledge.deduce_all_facts({fact:value})
            if value is None and not getattr(self,'_assumptions_query_active',False):
                self._assumptions_query_active=True
                try:
                    for other in _assume_defined:
                        if knowledge.get(other) is not None:continue
                        handler=getattr(self,'_eval_is_'+other,None)
                        if handler:knowledge.deduce_all_facts({other:handler()})
                        if knowledge.get(fact) is not None:break
                finally:self._assumptions_query_active=False
        return knowledge.get(fact)
    return property(getter)

class ManagedProperties(type):
    def __new__(mcls,name,bases,namespace,**kwargs):
        cls=super().__new__(mcls,name,bases,namespace,**kwargs)
        facts={}
        for base in reversed(bases):facts.update(getattr(base,'_explicit_class_assumptions',{}))
        facts.update({key[3:]:value for key,value in namespace.items() if key.startswith('is_') and key[3:] in _assume_defined and value in (True,False,None)})
        cls._explicit_class_assumptions=facts;cls.default_assumptions=StdFactKB(facts)
        for fact in _assume_defined:
            if fact in cls.default_assumptions:setattr(cls,'is_'+fact,cls.default_assumptions[fact])
            elif 'is_'+fact not in namespace:setattr(cls,'is_'+fact,make_property(fact))
        return cls

__all__ += ['StdFactKB','ManagedProperties','as_property','make_property']
