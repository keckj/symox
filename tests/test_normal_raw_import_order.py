"""Stock optional dispatchers do not disable native Gaussian moments."""
import subprocess
import sys

import pytest


@pytest.mark.parametrize('order',['stats-first','diffgeom-first'])
def test_stock_diffgeom_import_order_and_real_overrides(order):
    source=r'''
import functools,math
import symox as S
if ORDER=='stats-first':
    from symox import stats
    import symox.diffgeom
else:
    import symox.diffgeom
    from symox import stats
import symox.stats._joint as joint
cache=joint.lru_cache
native=S._core._stats_normal_raw

def forbidden(*args,**kwargs):
    raise AssertionError('native Gaussian route fell back to Python recurrence')
joint.lru_cache=forbidden
assert joint._normal_raw((2,2),[S.Integer(0)]*2,S.ImmutableMatrix([[2,1],[1,3]]))==8
value=joint._normal_raw((1200,),[S.Integer(0)],S.ImmutableMatrix([[1]]))
assert int(value)==math.prod(range(1,1200,2))
joint.lru_cache=cache

# Copying stock wrapper metadata does not turn a real callback into a trusted
# dispatcher. Both simple replacements and functools.wraps stay observable.
original=S.simplify
for decorated in (False,True):
    calls=[]
    def override(value,*args,**kwargs):
        calls.append(value)
        raise KeyboardInterrupt
    if decorated:override=functools.wraps(original)(override)
    S.simplify=override
    S._core._stats_normal_raw=forbidden
    try:
        joint._normal_raw((2,),[S.Integer(0)],S.ImmutableMatrix([[1]]))
    except KeyboardInterrupt:
        assert len(calls)==1
    else:
        raise AssertionError('a real final simplifier override was bypassed')
    S._core._stats_normal_raw=native
    S.simplify=original
'''.replace('ORDER',repr(order))
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
