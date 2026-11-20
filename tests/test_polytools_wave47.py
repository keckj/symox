"""Structural cancellation, polynomial options and nonfinite generators."""
import pytest
import sympy
import symox


class FailingTruth:
    def __bool__(self):
        raise RuntimeError('truth-value failure')


def cases(M):
    x, y = M.symbols('x y')
    A = M.Symbol('A', commutative=False)
    r = (x**2 - 1)/(x - 1)
    p = M.Piecewise((r, x < 0), ((x**2 - 4)/(x - 2), True), evaluate=False)
    q = M.Piecewise((1, y < 0), (2, True), evaluate=False)
    z = M.zoo
    tau = (-2*x*y + z + y/x**2 + (x**2 + 1/x)*(z*x**3 + z)/x**2
           + 2*(z*x**3 + z)/x**3)/(x**2*y + y/x + (z*x**3 + z)/x**2)
    inputs = {
        'rational': r,
        'piece': p,
        'piece sum': p + r,
        'piece product': p*r,
        'piece two factors': p*q*r,
        'piece function': M.sin(p),
        'piece power': p**2,
        'piece reciprocal': 1/p,
        'piece condition': M.Piecewise((r, r < 0), (0, True), evaluate=False),
        'piece nested': M.Piecewise((p*r, y < 0), (r, True), evaluate=False),
        'function ordinary': M.sin(r),
        'derivative': M.Derivative(M.Function('f')(r), x, evaluate=False),
        'noncommutative': (x**2*A - A)/(x - 1),
        'noncommutative sum': r + A,
        'noncommutative product': r*A,
        'noncommutative function': M.Function('f')(r*A),
        'noncommutative radical': (M.sqrt(3) + M.sqrt(15)*A)/(M.sqrt(2) + M.sqrt(10)*A),
        'zoo sum': z + x,
        'zoo product': z*x,
        'zoo expanded': z*(x**2 + 1)/x,
        'zoo ratio': (z*x + 1)/(z*x - 1),
        'zoo unrelated ratio': (z*x + y)/(z*x + 1),
        'zoo mixed ratio': (x + z*y)/(z*x + y),
        'zoo liegroup tau': tau,
        'infinity sum': M.oo + x,
        'infinity ratio': (M.oo*x + 1)/(M.oo*x - 1),
        'negative infinity': (-M.oo*x + 1)/(-M.oo*x - 1),
        'held zero zoo': M.Mul(0, z, evaluate=False),
        'held zoo cancellation': M.Mul(y, M.Pow(z, -1, evaluate=False), z, evaluate=False),
    }
    result = {
        f'{name} {mode}': (lambda e=e, kw=kw: M.cancel(e, **kw))
        for name, e in inputs.items()
        for mode, kw in [('default', {}), ('without signsimp', {'_signsimp': False})]
    }
    result.update({
        'signsimp None': lambda: M.cancel(M.Mul(0, z, evaluate=False), _signsimp=None),
        'signsimp zero': lambda: M.cancel(M.Mul(0, z, evaluate=False), _signsimp=0),
        'signsimp truth error': lambda: M.cancel(r, _signsimp=FailingTruth()),
        'piece explicit generator': lambda: M.cancel(p*r, x),
        'piece explicit QQ': lambda: M.cancel(p*r, x, domain='QQ'),
        'piece expressions tuple': lambda: M.cancel((p*r, x)),
        'piece expressions tuple polys': lambda: M.cancel((p*r, x), polys=True),
        'piece algebraic': lambda: M.cancel(M.Piecewise(((x**2 - 2)/(x + M.sqrt(2)), x < 0), (r, True))*r),
        'rational modulus control': lambda: M.cancel((x**2 + 1)/(x + 1), x, modulus=2),
        'algebraic extension control': lambda: M.cancel((x**2 - 2)/(x + M.sqrt(2)), extension=M.sqrt(2)),
        'gaussian control': lambda: M.cancel((x**2 + 1)/(x - M.I)),
        'high precision control': lambda: M.cancel((M.Float('0.125', 80)*x**2 - M.Float('0.125', 80))/(x - 1)),
    })
    for name, e in {'integer': M.S(2), 'relation': M.Eq(x, 1), 'boolean': M.true,
                    'infinity': M.oo, 'piece': p, 'expression': r,
                    'tuple': (x**2 - 1, x - 1)}.items():
        for mode, kw in [('deep', {'deep': True}), ('flag', {'auto': False}),
                         ('unknown', {'unknown': True})]:
            result[f'option {name} {mode}'] = lambda e=e, kw=kw: M.cancel(e, **kw)
    return result


def snapshot(value):
    if isinstance(value, (tuple, sympy.Tuple, symox.Tuple)):
        return tuple(snapshot(v) for v in value)
    if isinstance(value, (sympy.Poly, symox.Poly)):
        return ('Poly', str(value.as_expr()), tuple(map(str, value.gens)), str(value.domain))
    return (str(value), getattr(value, '_prec', None) if getattr(value, 'is_Float', False) else None)


def run(M, name):
    try:
        return ('value', snapshot(cases(M)[name]()))
    except Exception as exc:
        return ('error', type(exc).__name__)


@pytest.mark.parametrize('name', sorted(cases(sympy)))
def test_reference_contract(name):
    assert run(symox, name) == run(sympy, name)


def test_piecewise_branch_certificates():
    x = symox.Symbol('x')
    r = (x**2 - 1)/(x - 1)
    p = symox.Piecewise((r, x < 0), ((x**2 - 4)/(x - 2), True), evaluate=False)
    actual = symox.cancel(p*r)
    for value in (-3, -2, 0, 3):
        assert symox.cancel(actual.subs(x, value) - (p*r).subs(x, value)) == 0
