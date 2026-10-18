"""Reproducible numeric probes with an isolated random generator."""
import random as _random
import symox as S
rng=_random.Random()
choice=rng.choice;random=rng.random;randint=rng.randint;randrange=rng.randrange
sample=rng.sample;shuffle=rng.shuffle;uniform=rng.uniform

def seed(a=None,version=2):rng.seed(a,version=version)
def random_complex_number(a=2,b=-1,c=3,d=1,rational=False,tolerance=None):
    real,imag=uniform(float(a),float(c)),uniform(float(b),float(d))
    if rational:
        return S.nsimplify(real,rational=True,tolerance=tolerance)+S.I*S.nsimplify(imag,rational=True,tolerance=tolerance)
    return S.Float(real)+S.I*S.Float(imag)
def verify_numerically(f,g,z=None,tol=1e-6,a=2,b=-1,c=3,d=1):
    f,g=S.sympify(f),S.sympify(g)
    variables=f.free_symbols|g.free_symbols if z is None else [z] if getattr(z,'is_Symbol',False) else z
    mapping={v:random_complex_number(a,b,c,d) for v in variables}
    try:return abs(complex((f-g).subs(mapping,simultaneous=True).evalf(30)))<tol
    except (TypeError,ValueError):return False
def test_derivative_numerically(f,z,tol=1e-6,a=2,b=-1,c=3,d=1):
    point=random_complex_number(a,b,c,d)
    point=S.Float(str(S.re(point).evalf(50)),60)+S.I*S.Float(str(S.im(point).evalf(50)),60)
    step=S.Rational(1,10**8)
    finite=((f.subs(z,point+step)-f.subs(z,point-step))/(2*step)).evalf(50)
    exact=S.diff(f,z).subs(z,point).evalf(50)
    return abs(complex(finite-exact))<tol

test_derivative_numerically.__test__=False
__all__=['seed','random_complex_number','verify_numerically','test_derivative_numerically']
