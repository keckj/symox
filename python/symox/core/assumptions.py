"""Compatible imports backed by existing native implementations.

The rule engine objects (``_assume_rules``, ``StdFactKB``, ``_ask`` and the
class-assumption preparation) follow sympy/core/assumptions.py (SymPy 1.14).
Native expressions answer ``is_<fact>`` queries in the Rust core; these
objects serve Python classes and downstream code that inspect the rules.
"""
from symox._api_exports import export as _export
_export(globals(), 'symox._public_tools', 'check_assumptions failing_assumptions')

from symox.utilities.exceptions import sympy_deprecation_warning
from .facts import FactRules, FactKB
from .random import _assumptions_shuffle as shuffle
from .assumptions_generated import generated_assumptions as _assumptions


def _load_pre_generated_assumption_rules() -> FactRules:
    """Load the assumption rules from pre-generated data."""
    return FactRules._from_python(_assumptions)


def _generate_assumption_rules():
    """Generate the default assumption rules (source of the generated tables)."""
    return FactRules([
    'integer        ->  rational',
    'rational       ->  real',
    'rational       ->  algebraic',
    'algebraic      ->  complex',
    'transcendental ==  complex & !algebraic',
    'real           ->  hermitian',
    'imaginary      ->  complex',
    'imaginary      ->  antihermitian',
    'extended_real  ->  commutative',
    'complex        ->  commutative',
    'complex        ->  finite',

    'odd            ==  integer & !even',
    'even           ==  integer & !odd',

    'real           ->  complex',
    'extended_real  ->  real | infinite',
    'real           ==  extended_real & finite',

    'extended_real        ==  extended_negative | zero | extended_positive',
    'extended_negative    ==  extended_nonpositive & extended_nonzero',
    'extended_positive    ==  extended_nonnegative & extended_nonzero',

    'extended_nonpositive ==  extended_real & !extended_positive',
    'extended_nonnegative ==  extended_real & !extended_negative',

    'real           ==  negative | zero | positive',
    'negative       ==  nonpositive & nonzero',
    'positive       ==  nonnegative & nonzero',

    'nonpositive    ==  real & !positive',
    'nonnegative    ==  real & !negative',

    'positive       ==  extended_positive & finite',
    'negative       ==  extended_negative & finite',
    'nonpositive    ==  extended_nonpositive & finite',
    'nonnegative    ==  extended_nonnegative & finite',
    'nonzero        ==  extended_nonzero & finite',

    'zero           ->  even & finite',
    'zero           ==  extended_nonnegative & extended_nonpositive',
    'zero           ==  nonnegative & nonpositive',
    'nonzero        ->  real',

    'prime          ->  integer & positive',
    'composite      ->  integer & positive & !prime',
    '!composite     ->  !positive | !even | prime',

    'irrational     ==  real & !rational',

    'imaginary      ->  !extended_real',

    'infinite       ==  !finite',
    'noninteger     ==  extended_real & !integer',
    'extended_nonzero == extended_real & !zero',
    ])


_assume_rules = _load_pre_generated_assumption_rules()
_assume_defined = _assume_rules.defined_facts.copy()
_assume_defined.add('polar')
_assume_defined = frozenset(_assume_defined)


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


class StdFactKB(FactKB):
    """A FactKB specialized for the built-in rules."""
    def __init__(self, facts=None):
        super().__init__(_assume_rules)
        if not facts:
            self._generator = {}
        elif not isinstance(facts, FactKB):
            self._generator = facts.copy()
        else:
            self._generator = facts.generator
        if facts:
            self.deduce_all_facts(facts)

    def copy(self):
        return self.__class__(self)

    @property
    def generator(self):
        return self._generator.copy()


def as_property(fact):
    """Convert a fact name to the name of the corresponding property"""
    return 'is_%s' % fact


def make_property(fact):
    """Create the automagic property corresponding to a fact."""
    def getit(self):
        try:
            return self._assumptions[fact]
        except KeyError:
            if self._assumptions is self.default_assumptions:
                self._assumptions = self.default_assumptions.copy()
            return _ask(fact, self)

    getit.func_name = as_property(fact)
    return property(getit)


def _ask(fact, obj):
    """Find the truth value of ``fact`` for ``obj`` from its handlers and rules.

    Handlers ``_eval_is_<fact>`` are tried first; prerequisite facts are then
    queued breadth first (in shuffled order) and every value found is deduced
    into ``obj._assumptions``. An unanswerable query caches None.
    """
    assumptions = obj._assumptions
    handler_map = obj._prop_handler
    facts_to_check = [fact]
    facts_queued = {fact}
    for fact_i in facts_to_check:
        if fact_i in assumptions:
            continue
        fact_i_value = None
        handler_i = handler_map.get(fact_i)
        if handler_i is not None:
            fact_i_value = handler_i(obj)
        if fact_i_value is not None:
            assumptions.deduce_all_facts(((fact_i, fact_i_value),))
        fact_value = assumptions.get(fact)
        if fact_value is not None:
            return fact_value
        new_facts_to_check = list(_assume_rules.prereq[fact_i] - facts_queued)
        shuffle(new_facts_to_check)
        facts_to_check.extend(new_facts_to_check)
        facts_queued.update(new_facts_to_check)
    if fact in assumptions:
        return assumptions[fact]
    assumptions._tell(fact, None)
    return None


def _prepare_class_assumptions(cls):
    """Precompute class level assumptions and generate handlers."""
    local_defs = {}
    for k in _assume_defined:
        attrname = as_property(k)
        v = cls.__dict__.get(attrname, '')
        if isinstance(v, (bool, int, type(None))):
            if v is not None:
                v = bool(v)
            local_defs[k] = v

    defs = {}
    for base in reversed(cls.__bases__):
        assumptions = getattr(base, '_explicit_class_assumptions', None)
        if assumptions is not None:
            defs.update(assumptions)
    defs.update(local_defs)

    cls._explicit_class_assumptions = defs
    cls.default_assumptions = StdFactKB(defs)

    cls._prop_handler = {}
    for k in _assume_defined:
        eval_is_meth = getattr(cls, '_eval_is_%s' % k, None)
        if eval_is_meth is not None:
            cls._prop_handler[k] = eval_is_meth

    for k, v in cls.default_assumptions.items():
        setattr(cls, as_property(k), v)

    derived_from_bases = set()
    for base in cls.__bases__:
        default_assumptions = getattr(base, 'default_assumptions', None)
        if default_assumptions is not None:
            derived_from_bases.update(default_assumptions)

    for fact in derived_from_bases - set(cls.default_assumptions):
        pname = as_property(fact)
        if pname not in cls.__dict__:
            setattr(cls, pname, make_property(fact))

    for fact in _assume_defined:
        pname = as_property(fact)
        if not hasattr(cls, pname):
            setattr(cls, pname, make_property(fact))


class ManagedProperties(type):
    """Deprecated metaclass that prepares class assumptions (SymPy 1.12+)."""
    def __init__(cls, *args, **kwargs):
        msg = ("The ManagedProperties metaclass. "
               "Basic does not use metaclasses any more")
        sympy_deprecation_warning(msg,
            deprecated_since_version="1.12",
            active_deprecations_target='managedproperties')
        _prepare_class_assumptions(cls)

__all__ += ['StdFactKB','ManagedProperties','as_property','make_property']
