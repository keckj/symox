"""Bounded response reuse for exact integer prime-count/partition queries."""
from functools import lru_cache, wraps


def _cached(factory, evaluation):
    @lru_cache(maxsize=32)
    def integer_result(value):
        return factory(value, evaluate=True)

    @wraps(factory)
    def call(*args, **kwargs):
        # Restrict only cache eligibility. Every other domain, call form,
        # option and evaluation context retains its original factory path.
        if (not kwargs and len(args) == 1 and type(args[0]) is int
                and -1_000_000 <= args[0] <= 1_000_000 and evaluation.get()):
            return integer_result(args[0])
        return factory(*args, **kwargs)

    call._symox_integer_results = integer_result
    return call


def install(namespace, evaluation):
    from ._function_heads import BuiltinHeadType, FunctionHead, _own
    for name in ('primepi', 'partition'):
        head = namespace[name]
        # Keep the existing public object and aliases, signatures and pickle
        # restoration: only the head's own factory is wrapped, so user
        # subclasses and other mathematical heads keep their paths.
        if isinstance(head, BuiltinHeadType):
            factory = _own(head).get('_symox_factory')
            if factory is not None and not hasattr(factory, '_symox_integer_results'):
                type.__setattr__(head, '_symox_factory', _cached(factory, evaluation))
        elif isinstance(head, FunctionHead):
            factory = head.__dict__['_factory']
            if not hasattr(factory, '_symox_integer_results'):
                object.__setattr__(head, '_factory', _cached(factory, evaluation))
