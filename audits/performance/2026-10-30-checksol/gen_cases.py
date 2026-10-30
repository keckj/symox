"""Generate the differential checksol cases: equations x candidate
solutions x every combination of the five flags (plus `rational`, which
SymPy's checksol ignores), in the pair, dict and list forms."""
import itertools
import json
import sys

# (equation(s), solution) as Python source; `S:` prefix = dict form,
# a list equation = several equations.
BASE = [
    ("x**2 - 1", "x", "1"), ("x**2 - 1", "x", "-1"), ("x**2 - 1", "x", "0"), ("x**2 - 2", "x", "sqrt(2)"),
    ("x**2 - 2", "x", "-sqrt(2)"), ("x**2 - 2", "x", "1.41421356237310"), ("x**2 - 2", "x", "1.4142"),
    ("x**3 - 2", "x", "2**(S(1)/3)"), ("x**3 - 2", "x", "-2**(S(1)/3)/2 + 2**(S(1)/3)*sqrt(3)*I/2"),
    ("x - pi", "x", "pi"), ("x - pi", "x", "3.14159265358979"), ("x - pi", "x", "S(22)/7"),
    ("sin(x)", "x", "pi"), ("sin(x)", "x", "2*pi"), ("sin(x)", "x", "pi/2"), ("cos(x) - S(1)/2", "x", "pi/3"),
    ("exp(x) - 2", "x", "log(2)"), ("exp(x) + 1", "x", "I*pi"), ("exp(x) - y", "x", "log(y)"),
    ("x*exp(x) - 1", "x", "LambertW(1)"), ("x*exp(x) - y", "x", "LambertW(y)"), ("x - y", "x", "y"),
    ("x - y", "x", "y + 1"), ("x + 1", "x", "x"), ("x - x**2", "x", "x"), ("x - y", "y", "x + 1"),
    ("sqrt(x) + 1", "x", "1"), ("sqrt(x) - 1", "x", "1"), ("sqrt(x) + x - 2", "x", "1"), ("sqrt(x) + x - 2", "x", "4"),
    ("sqrt(x - y) - 1", "x", "y + 1"), ("sqrt(x**2) - x", "x", "y"), ("Abs(x) - x", "x", "y"), ("Abs(x) - x", "x", "-1"),
    ("Abs(x) - 1", "x", "-1"), ("Abs(x) - 1", "x", "I"), ("x**y - 1", "x", "1"), ("z**x - y", "x", "log(y)/log(z)"),
    ("y**(1/x) - z", "x", "log(y)/log(z)"), ("log(x) - 1", "x", "E"), ("log(x) - 1", "x", "exp(1 + 2*pi*I)"),
    ("1/x", "x", "0"), ("1/x - 1", "x", "1"), ("(x**2 - 1)/(x - 1)", "x", "1"), ("(x**2 - 1)/(x - 1) - 2", "x", "1"),
    ("x - 1", "x", "oo"), ("x - 1", "x", "zoo"), ("x - 1", "x", "nan"), ("x + oo", "x", "1"), ("nan", "x", "1"),
    ("sin(x)**2 + cos(x)**2 - 1", "x", "y"), ("sin(x)**2 + cos(x)**2 - 1", "x", "x"), ("tan(x) - 1", "x", "pi/4"),
    ("x**4 - 1", "x", "I"), ("x**4 - 1", "x", "-I"), ("x**4 - 1", "x", "1 + I"), ("x**2 + 1", "x", "sqrt(-1)"),
    ("x**2 + y**2 - 25", "S:", "{x: 3, y: 4}"), ("x**2 + y**2 - 25", "S:", "{x: 3, y: 5}"), ("x**2 + y**2 - 25", "S:", "{x: 3}"),
    ("x + y - 1", "S:", "{x: 1 - y}"), ("x*y - 1", "S:", "{x: 1/y}"), ("x*y - 1", "S:", "{x: 1/y, y: 2}"),
    ("x - y", "S:", "{x: 1}"), ("x - y", "S:", "{x: y}"), ("x**2 - y", "S:", "{x: sqrt(y)}"), ("x**2 - y", "S:", "{x: -sqrt(y)}"),
    ("x - 1", "S:", "{}"), ("x**2 + x - x*(x + 1)", "S:", "{}"), ("sin(x)**2 + cos(x)**2 - 1", "S:", "{}"), ("S(0)", "S:", "{}"),
    ("f(x) - 1", "f(x)", "1"), ("f(x) - x", "S:", "{f(x): x}"), ("f(x)**2 - x", "S:", "{f(x): sqrt(x)}"),
    ("x + y - 2", "S:", "{x: 1, y: 1}"), ("x + y - 2", "S:", "{y: 1, x: 1}"), ("x*y + z - 3", "S:", "{x: 1, y: 2, z: 1}"),
    ("[x - 1, y - 2]", "S:", "{x: 1, y: 2}"), ("[x - 1, y - 2]", "S:", "{x: 1, y: 3}"), ("[x - 1, y - 2]", "S:", "{x: 1}"),
    ("[x + y - 3, x - y - 1]", "S:", "{x: 2, y: 1}"), ("[x**2 - 2, y - x]", "S:", "{x: sqrt(2), y: sqrt(2)}"),
    ("[exp(x) - y, x - 1]", "S:", "{x: 1, y: E}"), ("[x - 1, sqrt(x) + 1]", "x", "1"), ("(x - 1, x**2 - 1)", "x", "1"),
    ("[x*exp(x) - 1, x - LambertW(1)]", "x", "LambertW(1)"), ("[sin(x) - y*cos(x), y - 1]", "S:", "{x: atan(y), y: 1}"),
    ("[]", "x", "1"), ("Matrix([x - 1, y - 2])", "S:", "{x: 1, y: 2}"),
    ("Eq(x, 1)", "x", "1"), ("Eq(x, 1)", "x", "2"), ("Eq(x**2, 4)", "x", "-2"), ("Eq(sin(x), 0)", "x", "pi"),
    ("Eq(x, y)", "S:", "{x: y}"), ("Eq(x, y)", "x", "1"), ("Eq(x, True)", "x", "True"), ("Eq(True, x)", "x", "True"),
    ("Eq(x, True)", "x", "1"), ("Eq(x, False)", "x", "False"), ("Eq(x, True)", "x", "y"),
    ("Ne(x, 1)", "x", "2"), ("Ne(x, 1)", "x", "1"), ("Ne(x, y)", "S:", "{x: 1}"), ("x < 1", "x", "0"), ("x < y", "x", "2"),
    ("Ge(x, 0)", "x", "y"), ("And(x > 0, x < 2)", "x", "1"), ("Or(Eq(x, 1), Eq(x, 2))", "x", "2"),
    ("True", "x", "1"), ("False", "x", "1"), ("S.true", "S:", "{}"), ("0", "S:", "{}"), ("1", "x", "1"),
    ("Poly(x**2 - 1, x)", "x", "1"), ("Poly(x**2 - 1, x)", "x", "2"), ("x - 1", "'x'", "1"), ("x - 1", "S:", "{'x': 1}"),
    ("Piecewise((x, x > 0), (0, True)) - 1", "x", "1"), ("Piecewise((x, x > 0), (0, True))", "x", "-1"),
    ("Piecewise((1, x > 0), (2, True)) - 1", "x", "y"), ("x - Piecewise((1, y > 0), (2, True))", "x", "1"),
    ("erf(x) - S(1)/2", "x", "erfinv(S(1)/2)"), ("x**(S(1)/3) - 2", "x", "8"), ("x**(S(1)/3) + 2", "x", "-8"),
    ("x**(S(1)/3) - y", "x", "y**3"), ("(x + 1)**2 - x**2 - 2*x - 1", "x", "y"), ("x - sqrt(2)*sqrt(3)", "x", "sqrt(6)"),
    ("x - sqrt(5 + 2*sqrt(6))", "x", "sqrt(2) + sqrt(3)"), ("x - cos(pi/5)", "x", "(1 + sqrt(5))/4"),
    ("x - 2*cos(2*pi/7)", "x", "CRootOf(x**3 + x**2 - 2*x - 1, 2)"), ("x**3 + x**2 - 2*x - 1", "x", "CRootOf(x**3 + x**2 - 2*x - 1, 0)"),
    ("x**5 - x - 1", "x", "CRootOf(x**5 - x - 1, 0)"), ("x**5 - x - 1", "x", "1.16730397826142"),
    ("x - y*z", "x", "y*z*(1 + 1e-12)"), ("x - 1", "x", "1 + 1e-10"), ("x - 1", "x", "1 + 1e-8"),
    ("a*x + b", "x", "-b/a"), ("a*x**2 + b*x + c", "x", "(-b + sqrt(b**2 - 4*a*c))/(2*a)"),
    ("a*x**2 + b*x + c", "x", "(-b - sqrt(b**2 - 4*a*c))/(2*a)"), ("x**2 - a", "x", "sqrt(a)"), ("x**2 - a**2", "x", "-a"),
    ("p*x - 1", "x", "1/p"), ("x**2 - p", "x", "-sqrt(p)"), ("sqrt(x) - p", "x", "p**2"), ("sqrt(x) + p", "x", "p**2"),
    ("x - log(exp(y))", "x", "y"), ("exp(x) - exp(y)", "x", "y + 2*pi*I"), ("x - asin(sin(y))", "x", "y"),
    ("a**(3*x + 5) + 3*log(a**(3*x + 5))", "x", "log(3*LambertW(S(1)/3)/a**5)/(3*log(a))"),
    ("x**2 - 4*x + 4", "x", "2"), ("(x - 1)**3", "x", "1"), ("x*(x - 1) - x**2 + x", "x", "5"),
    ("x - 1", "x", "Integral(1, (y, 0, 1))"), ("x - Sum(1/2**n, (n, 1, oo))", "x", "1"),
    ("x - 1", "x", "1*m"), ("Derivative(f(x), x) - 1", "f(x)", "x"),
]

FLAGNAMES = ('numerical', 'minimal', 'warn', 'simplify', 'force')


def main(out):
    cases = []
    for fe, sym, sol in BASE:
        for bits in itertools.product((None, True, False), repeat=5) if False else itertools.product((True, False), repeat=5):
            flags = dict(zip(FLAGNAMES, bits))
            cases.append({'f': fe, 'sym': sym, 'sol': sol, 'flags': flags})
        # defaults (no flags) and the ignored `rational` flag
        cases.append({'f': fe, 'sym': sym, 'sol': sol, 'flags': {}})
        cases.append({'f': fe, 'sym': sym, 'sol': sol, 'flags': {'rational': True}})
    for i, c in enumerate(cases):
        c['id'] = i
    json.dump(cases, open(out, 'w'))
    print(len(cases), 'cases')


main(sys.argv[1])
