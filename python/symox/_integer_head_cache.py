"""Bounded response reuse for exact integer prime-count/partition queries."""
from functools import lru_cache


def install(namespace, evaluation):
    base = namespace['_BuiltinFunctionHead']

    class CachedIntegerHead(base):
        def __call__(self,*args,**kwargs):
            # Restrict only cache eligibility. Every other domain, call form,
            # option and evaluation context retains its original factory path.
            if (not kwargs and len(args)==1 and type(args[0]) is int
                    and -1_000_000 <= args[0] <= 1_000_000 and evaluation.get()):
                return self._integer_results(args[0])
            return self._factory(*args,**kwargs)

    for name in ('primepi','partition'):
        head = namespace[name]
        if type(head) is not base:
            continue
        factory = head._factory
        def integer_result(value,_factory=factory):
            return _factory(value,evaluate=True)
        head._integer_results = lru_cache(maxsize=32)(integer_result)
        # Keep the existing public object and aliases, signatures and pickle
        # restoration. Other mathematical heads do not pay a cache-dispatch tax.
        head.__class__ = CachedIntegerHead
