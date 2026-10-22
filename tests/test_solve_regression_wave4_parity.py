"""Polynomial systems whose elimination polynomial decomposes into quartic
radicals.  `roots_binomial` expands `sqrt` of each quartic root with
`expand(complex=True)`, and the system check evaluates the exact residuals.
Before the fix the solve did not finish within 200 s (SymPy 1.14: about 12 s).

Each system runs in a subprocess with a time limit; the solution set is
compared at 45 digits with the values of SymPy 1.14.
"""
import json
import os
import subprocess
import sys

import pytest

import symox as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SCRIPT = r"""
import json, sys, time
import symox as sp
x, y = sp.symbols('x y')
eqs = list(sp.sympify('[' + sys.argv[1] + ']'))
start = time.perf_counter()
result = sp.solve(eqs, [x, y])
elapsed = time.perf_counter() - start
rows = []
for row in result:
    values = []
    for value in row:
        number = sp.N(value, 60)
        values.append([str(sp.N(sp.re(number), 45)), str(sp.N(sp.im(number), 45))])
    rows.append(values)
exact = all(not value.atoms(sp.Float) for row in result for value in row)
print(json.dumps(dict(time=elapsed, rows=rows, exact=exact)))
"""

# SymPy 1.14 solutions, evaluated at 45 digits: ((re x, im x), (re y, im y)).
CASES = {
    'x**2 + y**2 - 4, y - x**4 + 1': [
        (('-1.26369678126881728554213547428086705508960980', '0'), ('1.55018400359790545804540128930874953234181707', '0')),
        (('1.26369678126881728554213547428086705508960980', '0'), ('1.55018400359790545804540128930874953234181707', '0')),
        (('-0.752989911121086840448506105991089475102791621', '-0.664082754259080106435166274140051897041210935'), ('-1.98431751410846420932585884411247677028105846', '0.251999798697164713872548391590551987107633865')),
        (('0.752989911121086840448506105991089475102791621', '0.664082754259080106435166274140051897041210935'), ('-1.98431751410846420932585884411247677028105846', '0.251999798697164713872548391590551987107633865')),
        (('-0.752989911121086840448506105991089475102791621', '0.664082754259080106435166274140051897041210935'), ('-1.98431751410846420932585884411247677028105846', '-0.251999798697164713872548391590551987107633865')),
        (('0.752989911121086840448506105991089475102791621', '-0.664082754259080106435166274140051897041210935'), ('-1.98431751410846420932585884411247677028105846', '-0.251999798697164713872548391590551987107633865')),
        (('0', '-1.35974459310592664047743122161793598009894360'), ('2.41845102461902296060631639891620400822029985', '0')),
        (('0', '1.35974459310592664047743122161793598009894360'), ('2.41845102461902296060631639891620400822029985', '0')),
    ],
    'x**2 - y**2 - 1, y - x**4 + 2': [
        (('-0.220593985066903876318669477755510181673164495', '1.22684065668885156231324774650610025452185849'), ('-0.171647147024426874961790802874892066111974602', '1.57668609232740438575833275724943518263931118')),
        (('0.220593985066903876318669477755510181673164495', '-1.22684065668885156231324774650610025452185849'), ('-0.171647147024426874961790802874892066111974602', '1.57668609232740438575833275724943518263931118')),
        (('-0.220593985066903876318669477755510181673164495', '-1.22684065668885156231324774650610025452185849'), ('-0.171647147024426874961790802874892066111974602', '-1.57668609232740438575833275724943518263931118')),
        (('0.220593985066903876318669477755510181673164495', '1.22684065668885156231324774650610025452185849'), ('-0.171647147024426874961790802874892066111974602', '-1.57668609232740438575833275724943518263931118')),
        (('-1.11002083852196480592287589772291987523074404', '0'), ('-0.481815589155234647049148920022207174536217316', '0')),
        (('1.11002083852196480592287589772291987523074404', '0'), ('-0.481815589155234647049148920022207174536217316', '0')),
        (('-1.29645914681530339216434781861021287208523455', '0'), ('0.825109883204088396972730525771991306760166520', '0')),
        (('1.29645914681530339216434781861021287208523455', '0'), ('0.825109883204088396972730525771991306760166520', '0')),
    ],
}


def _close(a, b):
    return all(abs(float(sp.Rational(u) - sp.Rational(v))) < 1e-40 for u, v in zip(a, b))


def _same_set(got, expected):
    got = [tuple(tuple(part) for part in row) for row in got]
    unused = list(range(len(got)))
    for row in expected:
        match = next((k for k in unused
                      if all(_close(gv, ev) for gv, ev in zip(got[k], row))), None)
        if match is None:
            return False
        unused.remove(match)
    return not unused


@pytest.mark.parametrize('system', sorted(CASES))
def test_quartic_radical_systems_finish_with_the_reference_solutions(system):
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'python'))
    process = subprocess.run([sys.executable, '-c', SCRIPT, system], env=env, cwd=ROOT,
                             capture_output=True, text=True, timeout=60)
    assert process.returncode == 0, process.stderr[-2000:]
    result = json.loads(process.stdout.strip().splitlines()[-1])
    assert result['exact']
    assert len(result['rows']) == 8
    assert _same_set(result['rows'], CASES[system])
    assert result['time'] < 40


def test_sqrt_of_a_positive_radical_stays_whole_under_complex_expansion():
    # SymPy's Pow.as_real_imag keeps sqrt(a) for a nonnegative number a,
    # decided numerically; it was expanded through atan2 and Piecewise.
    c = sp.sympify('-421/432 + sqrt(34257)/144')**sp.Rational(1, 3)
    s = -16/(9*c) + sp.Rational(4, 3) + 2*c
    assert sp.sqrt(s).expand(complex=True) == sp.sqrt(s)
    assert (1/sp.sqrt(s)).expand(complex=True) == 1/sp.sqrt(s)
    assert sp.re(sp.sqrt(s)) == sp.sqrt(s)
    assert sp.im(sp.sqrt(s)) == 0
    assert sp.sqrt(-s).expand(complex=True) == sp.I*sp.sqrt(s)


def test_quartic_decomposition_roots_keep_their_values():
    # The elimination polynomial of the first system: t = x**2 gives a quartic,
    # then roots_binomial takes the square roots of its exact roots.
    x = sp.Symbol('x')
    roots = sp.solve(x**8 - 2*x**4 + x**2 - 3, x)
    assert all(not r.atoms(sp.Float) for r in roots)
    got = [[(str(sp.re(sp.N(r, 50))), str(sp.im(sp.N(r, 50))))] for r in roots]
    expected = [row[:1] for row in CASES['x**2 + y**2 - 4, y - x**4 + 1']]
    assert len(got) == 8
    assert _same_set(got, expected)


def test_arg_of_radical_complex_numbers_follows_the_numeric_split():
    # SymPy 1.14 outputs: the split keeps sqrt(c) whole (c > 0 numerically),
    # so arg evaluates through atan2 of the reference real and imaginary parts.
    c = -1 + 10*sp.Integer(10)**sp.Rational(1, 3)
    b = 1 + sp.I*sp.sqrt(c)
    expected = [
        (b, 'atan(sqrt(-1 + 10*10**(1/3)))'),
        (b**2, 'atan(2*sqrt(-1 + 10*10**(1/3))/(2 - 10*10**(1/3))) + pi'),
        (b**3, '-pi + atan((-(-1 + 10*10**(1/3))**(3/2) + 3*sqrt(-1 + 10*10**(1/3)))/(4 - 30*10**(1/3)))'),
        (1/b, '-atan(sqrt(-1 + 10*10**(1/3)))'),
        (b**-3, 'atan((-3*sqrt(-1 + 10*10**(1/3))/10000 + (-1 + 10*10**(1/3))**(3/2)/10000)'
                '/(1/2500 - 3*10**(1/3)/1000)) + pi'),
        (sp.sqrt(c), '0'),
    ]
    for value, text in expected:
        assert str(sp.arg(value)) == text


def test_phase_margin_of_a_third_order_lag():
    from symox.physics.control import TransferFunction, phase_margin
    s = sp.Symbol('s')
    pm = phase_margin(TransferFunction(100, (s + 1)**3, s))
    assert str(pm) == ('-180 + 180*(atan((-3*sqrt(-1 + 10*10**(1/3))/10000 + (-1 + 10*10**(1/3))**(3/2)/10000)'
                       '/(1/2500 - 3*10**(1/3)/1000)) + pi)/pi')
    assert abs(float(pm.n(15)) + 52.6753621532920) < 1e-9
