import sys, inspect, json
lib = sys.argv[1]
if lib == 'sympy':
    import sympy as S
else:
    import symox as S
from_ = S
x, y, z = S.symbols('x y z')
f = S.Function('f')
exprs = [
 ('sin', S.sin(x)), ('poly', x**2 + 2*x + 1), ('rat', S.Rational(1,3)*x), ('sqrt', S.sqrt(x)), ('invsqrt', 1/S.sqrt(x)),
 ('pi', S.pi*x), ('E', S.E**x), ('exp', S.exp(x) + S.cos(y)), ('I', S.I*x), ('abs', S.Abs(x)), ('float', 1.5*x + 0.1),
 ('pw', S.Piecewise((x, x > 0), (-x, True))), ('rel', x > 0), ('and', S.And(x > 0, y < 1)), ('max', S.Max(x, y)),
 ('log', S.log(x) + S.atan2(y, x)), ('gamma', S.gamma(x)), ('besselj', S.besselj(1, x)), ('erf', S.erf(x)),
 ('pow', x**y), ('inv', 1/x), ('ceil', S.ceiling(x) + S.floor(y)), ('mat', S.Matrix([[x, 1], [y, 2]])),
 ('list', [x, y]), ('tuple', (x, x*y)), ('acos', S.acos(x)), ('sinh', S.sinh(x)), ('sign', S.sign(x)),
 ('oo', x + S.oo), ('negpow', x**-2), ('re', S.re(x)), ('conj', S.conjugate(x)), ('Heaviside', S.Heaviside(x)),
 ('sinc', S.sinc(x)), ('tan', S.tan(x)**2), ('half', x/2), ('log10', S.log(x, 10)), ('Mod', S.Mod(x, 3)),
 ('factorial', S.factorial(x)), ('zeta', S.zeta(x)), ('loggamma', S.loggamma(x)), ('LambertW', S.LambertW(x)),
]
out = {}
for mod in [None, 'numpy', 'math', 'mpmath', 'sympy', ['math', 'mpmath'], 'scipy']:
    for name, e in exprs:
        key = '%s|%s' % (mod, name)
        args = [x, y] if not isinstance(e, (list, tuple)) and len(getattr(e, 'free_symbols', [x, y])) > 1 or isinstance(e, (list, tuple)) else x
        try:
            fn = S.lambdify(args, e, mod)
        except Exception as err:
            out[key] = 'LAMBDIFY-ERR %s: %s' % (type(err).__name__, str(err)[:80]); continue
        try:
            src = inspect.getsource(fn)
        except Exception as err:
            src = 'SRC-ERR %s' % type(err).__name__
        try:
            doc = fn.__doc__
        except Exception as err:
            doc = 'DOC-ERR'
        try:
            sig = str(inspect.signature(fn))
        except Exception as err:
            sig = 'SIG-ERR'
        out[key] = {'src': src, 'doc': doc, 'sig': sig, 'name': getattr(fn, '__name__', None)}
json.dump(out, open(sys.argv[2], 'w'), indent=1, sort_keys=True)
