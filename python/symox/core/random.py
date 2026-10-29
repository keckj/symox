"""Reproducible numeric probes with an isolated random generator."""
import random as _random
import symox as S
rng=_random.Random()
choice=rng.choice;random=rng.random;randint=rng.randint;randrange=rng.randrange
sample=rng.sample;shuffle=rng.shuffle;uniform=rng.uniform
# Assumption queries shuffle prerequisite facts with their own seeded generator.
_assumptions_rng=_random.Random()
_assumptions_shuffle=_assumptions_rng.shuffle

def seed(a=None,version=2):
    rng.seed(a,version=version);_assumptions_rng.seed(a,version=version)
def random_complex_number(a=2,b=-1,c=3,d=1,rational=False,tolerance=None):
    real,imag=uniform(float(a),float(c)),uniform(float(b),float(d))
    if rational:
        return S.nsimplify(real,rational=True,tolerance=tolerance)+S.I*S.nsimplify(imag,rational=True,tolerance=tolerance)
    return S.Float(real)+S.I*S.Float(imag)
def verify_numerically(f,g,z=None,tol=1e-6,a=2,b=-1,c=3,d=1):
    # SymPy 1.14: both sides are evaluated separately at the default
    # precision and compared with comp; evaluating the (often identically
    # zero) difference instead made evalf raise the precision to its limit.
    from .numbers import comp
    f,g=(S.sympify(i) for i in (f,g))
    if z is None:
        z=f.free_symbols|g.free_symbols
    elif isinstance(z,S.Symbol):
        z=[z]
    reps=list(zip(z,[random_complex_number(a,b,c,d) for _ in z]))
    z1=f.subs(reps).n()
    z2=g.subs(reps).n()
    return comp(z1,z2,tol)
def test_derivative_numerically(f,z,tol=1e-6,a=2,b=-1,c=3,d=1):
    point=random_complex_number(a,b,c,d)
    point=S.Float(str(S.re(point).evalf(50)),60)+S.I*S.Float(str(S.im(point).evalf(50)),60)
    step=S.Rational(1,10**8)
    finite=((f.subs(z,point+step)-f.subs(z,point-step))/(2*step)).evalf(50)
    exact=S.diff(f,z).subs(z,point).evalf(50)
    return abs(complex(finite-exact))<tol

test_derivative_numerically.__test__=False
__all__=['seed','random_complex_number','verify_numerically','test_derivative_numerically']


def _randrange(seed=None):
    """Return a randrange generator seeded by an int or by a value sequence."""
    if seed is None:
        return randrange
    if isinstance(seed,int):
        rng.seed(seed);return randrange
    from ..utilities.iterables import is_sequence
    if is_sequence(seed):
        seed=list(seed);seed.reverse()
        def give(a,b=None,seq=seed):
            from ..utilities.misc import as_int
            if b is None:a,b=0,a
            a,b=as_int(a),as_int(b)
            if b-a<1:raise ValueError('_randrange got empty range')
            try:x=seq.pop()
            except IndexError:raise ValueError('_randrange sequence was too short')
            if a<=x<b:return x
            return give(a,b,seq)
        return give
    raise ValueError('_randrange got an unexpected seed')
