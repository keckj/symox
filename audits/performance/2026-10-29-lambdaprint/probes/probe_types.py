import sys, json, warnings
warnings.simplefilter('ignore')
lib = sys.argv[1]
if lib == 'sympy':
    import sympy as S
else:
    import symox as S
import numpy as np
x, y = S.symbols('x y')
exprs = [
 ('sin', S.sin(x)), ('poly', x**2 + 2*x + 1), ('rat', S.Rational(1,3)*x), ('sqrt', S.sqrt(x)), ('invsqrt', 1/S.sqrt(x)),
 ('pi', S.pi*x), ('E', S.E**x), ('exp', S.exp(x) + S.cos(x)), ('I', S.I*x), ('abs', S.Abs(x)), ('float', 1.5*x + 0.1),
 ('pw', S.Piecewise((x, x > 0), (-x, True))), ('rel', x > 0), ('and', S.And(x > 0, x < 1)), ('max', S.Max(x, 1)),
 ('log', S.log(x) + S.atan2(1, x)), ('gamma', S.gamma(x)), ('erf', S.erf(x)),
 ('pow', x**x), ('inv', 1/x), ('ceil', S.ceiling(x) + S.floor(x)), ('acos', S.acos(x/4)), ('sinh', S.sinh(x)), ('sign', S.sign(x)),
 ('oo', x + S.oo), ('negpow', x**-2), ('re', S.re(x)), ('im', S.im(x)), ('conj', S.conjugate(x)), ('Heaviside', S.Heaviside(x)),
 ('sinc', S.sinc(x)), ('tan', S.tan(x)**2), ('half', x/2), ('log10', S.log(x, 10)), ('Mod', S.Mod(x, 3)), ('const', S.Integer(3)),
 ('fconst', S.Float(2.5)), ('piconst', S.pi), ('sinpi', S.sin(x) + S.pi), ('expI', S.exp(S.I*x)), ('sqrtneg', S.sqrt(x - 5)),
 ('list', [S.sin(x), x]), ('arg', S.arg(x)), ('cosh2', S.cosh(x)*x + 1), ('min', S.Min(x, 2)), ('exp2', 2**x), ('atan', S.atan(x)),
]
out = {}
for mod in [None, 'numpy', 'math', 'mpmath', 'scipy']:
    for name, e in exprs:
        try:
            fn = S.lambdify(x, e, mod)
        except Exception as err:
            out['%s|%s' % (mod, name)] = 'LAMBDIFY-ERR %s' % type(err).__name__; continue
        for vname, v in [('f', 0.5), ('i', 2), ('npf', np.float64(0.5)), ('npi', np.int64(2)), ('arr', np.array([0.5, 2.0]))]:
            key = '%s|%s|%s' % (mod, name, vname)
            try:
                r = fn(v)
                t = type(r).__module__ + '.' + type(r).__name__
                if isinstance(r, list): t += '[' + ','.join(type(i).__name__ for i in r) + ']'
                if isinstance(r, np.ndarray): t += '/' + str(r.dtype) + '/' + str(r.shape)
                out[key] = t + ' ' + repr(r)[:60]
            except Exception as err:
                out[key] = 'ERR %s: %s' % (type(err).__name__, str(err)[:60])
json.dump(out, open(sys.argv[2], 'w'), indent=0, sort_keys=True)
