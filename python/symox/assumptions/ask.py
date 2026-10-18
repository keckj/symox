"""Three-valued queries and the built-in predicate namespace."""
from .. import _core
from .assume import global_assumptions


class AssumptionKeys:
    """Namespace of built-in predicates, such as ``Q.positive``."""
    def __init__(self):
        self.__dict__.update((name, _core._predicate(name))
                             for name in _core._predicate_names())


Q = AssumptionKeys()


def ask(proposition, assumptions=True, context=global_assumptions):
    """Return True, False or None when the proposition is undecidable."""
    from .._public_tools import ask as public_ask
    return public_ask(proposition, assumptions, context)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox._public_tools', 'register_handler remove_handler')


def _fact_predicate(name):
    if hasattr(Q,name):return getattr(Q,name)
    return _core._predicate(name,undefined=True)


def get_all_known_facts():
    from ._facts import CLAUSES
    from .cnf import Literal
    return {frozenset(Literal(_fact_predicate(name),negated) for name,negated in clause) for clause in CLAUSES}


def get_known_facts_dict():
    from ._facts import KNOWN
    return {_fact_predicate(name):({*map(_fact_predicate,positive)},{*map(_fact_predicate,negative)}) for name,(positive,negative) in KNOWN.items()}

__all__ += ['get_all_known_facts','get_known_facts_dict']
