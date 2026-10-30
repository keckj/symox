"""Native constant_system over rational function fields (prdeconst).

The generated cases and the expected outputs (srepr of the entries of
to_Matrix(), domain and generators of both PolyMatrix results) are
hard-coded from SymPy 1.14 (sympy.integrals.prde.constant_system); the
fraction-field kernel must answer each case (_prde_native.STATS).
"""
from symox import (Symbol, Integer, Rational, Add, Mul, Pow, Poly, QQ, srepr,
                   symbols)
from symox.polys.polymatrix import PolyMatrix as Matrix
from symox.integrals.risch import DifferentialExtension
from symox.integrals.prde import (constant_system, is_deriv_k,
                                  is_log_deriv_k_t_radical_in_field)
from symox.integrals import _prde_native

x, t, t0, t1 = symbols('x t t0 t1')
_NS = dict(Symbol=Symbol, Integer=Integer, Rational=Rational, Add=Add, Mul=Mul, Pow=Pow)

# SymPy 1.14 reference: kind, D, shapes, items (srepr), input domains, output
CASES = [
    {"kind":"base","D":["1"],"rows":3,"cols":1,"A":["Mul(Rational(1, 3), Symbol('x'))","Integer(0)","Integer(-1)"],"u":["Integer(0)","Integer(0)","Mul(Integer(-1), Rational(1, 3), Pow(Symbol('x'), Integer(-1)), Add(Mul(Rational(1, 2), Symbol('x')), Integer(-1)))"],"inA":"ZZ(x)","inu":"ZZ(x)","out":[{"shape":[3,1],"domain":"ZZ(x)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)"]},{"shape":[3,1],"domain":"ZZ(x)","gens":["t1"],"items":["Integer(0)","Integer(1)","Integer(0)"]}]},
    {"kind":"base","D":["1"],"rows":2,"cols":2,"A":["Integer(-1)","Rational(-1, 6)","Integer(0)","Mul(Integer(4), Symbol('x'))"],"u":["Integer(0)","Mul(Integer(-1), Symbol('x'))"],"inA":"ZZ(x)","inu":"ZZ(x)","out":[{"shape":[2,2],"domain":"ZZ(x)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[2,1],"domain":"ZZ(x)","gens":["t1"],"items":["Rational(1, 24)","Rational(-1, 4)"]}]},
    {"kind":"base","D":["1"],"rows":1,"cols":1,"A":["Add(Mul(Integer(2), Symbol('x')), Integer(1))"],"u":["Integer(0)"],"inA":"ZZ(x)","inu":"QQ","out":[{"shape":[1,1],"domain":"QQ(x)","gens":["t1"],"items":["Integer(1)"]},{"shape":[1,1],"domain":"QQ(x)","gens":["t1"],"items":["Integer(0)"]}]},
    {"kind":"base","D":["1"],"rows":1,"cols":2,"A":["Integer(0)","Integer(2)"],"u":["Integer(0)"],"inA":"QQ","inu":"QQ","out":[{"shape":[1,2],"domain":"QQ","gens":["t1"],"items":["Integer(0)","Integer(1)"]},{"shape":[1,1],"domain":"QQ","gens":["t1"],"items":["Integer(0)"]}]},
    {"kind":"base","D":["1"],"rows":3,"cols":2,"A":["Integer(1)","Integer(-1)","Mul(Integer(-1), Rational(1, 4), Symbol('x'))","Integer(0)","Mul(Integer(-1), Rational(1, 3), Symbol('x'))","Mul(Symbol('x'), Pow(Add(Mul(Integer(-1), Symbol('x')), Integer(-3)), Integer(-1)))"],"u":["Integer(0)","Integer(0)","Integer(0)"],"inA":"ZZ(x)","inu":"QQ","out":[{"shape":[3,2],"domain":"QQ(x)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)"]},{"shape":[3,1],"domain":"QQ(x)","gens":["t1"],"items":["Integer(0)","Integer(0)","Integer(0)"]}]},
    {"kind":"base","D":["1"],"rows":2,"cols":3,"A":["Mul(Integer(2), Pow(Symbol('x'), Integer(-1)), Add(Mul(Integer(2), Symbol('x')), Integer(1)))","Mul(Integer(-1), Integer(4), Symbol('x'))","Integer(0)","Mul(Rational(1, 2), Symbol('x'))","Integer(1)","Integer(-1)"],"u":["Integer(0)","Mul(Integer(-1), Integer(3), Pow(Add(Integer(2), Mul(Integer(-1), Integer(3), Symbol('x'))), Integer(-1)))"],"inA":"ZZ(x)","inu":"ZZ(x)","out":[{"shape":[3,3],"domain":"QQ(x)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[3,1],"domain":"ZZ(x)","gens":["t1"],"items":["Mul(Integer(-1), Integer(18), Pow(Symbol('x'), Integer(3)), Pow(Add(Mul(Integer(9), Pow(Symbol('x'), Integer(5))), Mul(Integer(-1), Integer(12), Pow(Symbol('x'), Integer(4))), Mul(Integer(-1), Integer(14), Pow(Symbol('x'), Integer(3))), Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(16), Symbol('x')), Integer(-8)), Integer(-1)))","Mul(Add(Mul(Integer(-1), Integer(18), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(9), Symbol('x'))), Pow(Add(Mul(Integer(9), Pow(Symbol('x'), Integer(5))), Mul(Integer(-1), Integer(12), Pow(Symbol('x'), Integer(4))), Mul(Integer(-1), Integer(14), Pow(Symbol('x'), Integer(3))), Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(16), Symbol('x')), Integer(-8)), Integer(-1)))","Mul(Add(Mul(Integer(-1), Integer(18), Pow(Symbol('x'), Integer(4))), Mul(Integer(6), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(3), Symbol('x')), Integer(-12)), Pow(Add(Mul(Integer(9), Pow(Symbol('x'), Integer(5))), Mul(Integer(-1), Integer(12), Pow(Symbol('x'), Integer(4))), Mul(Integer(-1), Integer(14), Pow(Symbol('x'), Integer(3))), Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(16), Symbol('x')), Integer(-8)), Integer(-1)))"]}]},
    {"kind":"base","D":["1"],"rows":3,"cols":3,"A":["Rational(-1, 2)","Integer(0)","Symbol('x')","Integer(0)","Mul(Integer(-1), Rational(3, 2), Pow(Symbol('x'), Integer(-1)), Add(Mul(Rational(1, 2), Symbol('x')), Integer(-1)))","Mul(Rational(4, 9), Symbol('x'))","Mul(Integer(-1), Pow(Symbol('x'), Integer(-1)))","Integer(0)","Mul(Integer(-1), Pow(Symbol('x'), Integer(-1)), Add(Mul(Integer(-1), Integer(3), Symbol('x')), Rational(-2, 3)))"],"u":["Integer(0)","Integer(0)","Mul(Integer(-1), Integer(3), Symbol('x'), Pow(Add(Symbol('x'), Rational(1, 2)), Integer(-1)))"],"inA":"ZZ(x)","inu":"ZZ(x)","out":[{"shape":[3,3],"domain":"ZZ(x)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[3,1],"domain":"ZZ(x)","gens":["t1"],"items":["Mul(Integer(-1), Integer(36), Pow(Symbol('x'), Integer(3)), Pow(Add(Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(7), Symbol('x')), Integer(2)), Integer(-1)))","Mul(Integer(-1), Integer(32), Pow(Symbol('x'), Integer(4)), Pow(Add(Mul(Integer(18), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(15), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(36), Symbol('x')), Integer(-12)), Integer(-1)))","Mul(Integer(-1), Integer(18), Pow(Symbol('x'), Integer(2)), Pow(Add(Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(7), Symbol('x')), Integer(2)), Integer(-1)))"]}]},
    {"kind":"base","D":["1"],"rows":1,"cols":3,"A":["Integer(2)","Rational(1, 3)","Add(Mul(Integer(-1), Rational(2, 3), Symbol('x')), Integer(-1))"],"u":["Integer(0)"],"inA":"ZZ(x)","inu":"QQ","out":[{"shape":[2,3],"domain":"QQ(x)","gens":["t1"],"items":["Integer(1)","Rational(1, 6)","Integer(0)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[2,1],"domain":"QQ(x)","gens":["t1"],"items":["Integer(0)","Integer(0)"]}]},
    {"kind":"exp","D":["1","t0"],"rows":3,"cols":2,"A":["Mul(Integer(-1), Rational(2, 3), Symbol('t0'), Pow(Add(Symbol('t0'), Mul(Integer(-1), Rational(2, 3), Symbol('x'))), Integer(-1)))","Mul(Rational(1, 2), Pow(Symbol('t0'), Integer(-1)), Pow(Symbol('x'), Integer(-1)), Add(Mul(Symbol('t0'), Symbol('x')), Integer(2)))","Mul(Integer(-1), Rational(2, 3), Symbol('t0'), Pow(Add(Mul(Symbol('t0'), Symbol('x')), Mul(Integer(-1), Integer(3), Symbol('x'))), Integer(-1)))","Pow(Symbol('x'), Integer(-1))","Integer(0)","Rational(1, 3)"],"u":["Integer(0)","Integer(0)","Mul(Integer(-1), Integer(3), Pow(Add(Symbol('t0'), Integer(2)), Integer(-1)))"],"inA":"ZZ(x,t0)","inu":"ZZ(t0)","out":[{"shape":[3,2],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)"]},{"shape":[3,1],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(0)","Integer(0)","Integer(1)"]}]},
    {"kind":"exp","D":["1","t0"],"rows":2,"cols":3,"A":["Rational(1, 3)","Mul(Integer(-1), Integer(6), Symbol('t0'))","Integer(0)","Mul(Rational(1, 2), Pow(Symbol('t0'), Integer(-1)), Add(Mul(Integer(2), Symbol('x')), Rational(1, 2)))","Mul(Rational(4, 3), Symbol('x'))","Mul(Rational(3, 2), Symbol('x'), Pow(Add(Mul(Integer(-1), Integer(3), Symbol('x')), Rational(-2, 3)), Integer(-1)))"],"u":["Integer(0)","Integer(0)"],"inA":"ZZ(x,t0)","inu":"QQ","out":[{"shape":[3,3],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[3,1],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(0)","Integer(0)","Integer(0)"]}]},
    {"kind":"exp","D":["1","t0"],"rows":3,"cols":3,"A":["Mul(Integer(-1), Symbol('x'))","Add(Mul(Integer(-1), Integer(3), Symbol('t0')), Mul(Integer(-1), Rational(2, 3), Symbol('x')))","Integer(1)","Mul(Integer(-1), Symbol('t0'), Pow(Add(Mul(Integer(-1), Rational(2, 3), Symbol('t0'), Symbol('x')), Rational(-2, 3)), Integer(-1)))","Integer(0)","Integer(0)","Integer(0)","Mul(Rational(1, 18), Symbol('x'))","Mul(Integer(-1), Rational(2, 3), Symbol('x'), Pow(Add(Mul(Integer(2), Symbol('t0'), Symbol('x')), Mul(Integer(-1), Integer(3), Symbol('x'))), Integer(-1)))"],"u":["Integer(0)","Integer(0)","Add(Integer(1), Mul(Integer(-1), Symbol('t0')))"],"inA":"ZZ(x,t0)","inu":"ZZ(t0)","out":[{"shape":[3,3],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[3,1],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(0)","Mul(Add(Mul(Integer(-1), Integer(36), Pow(Symbol('t0'), Integer(2))), Mul(Integer(90), Symbol('t0')), Integer(-54)), Pow(Add(Mul(Integer(2), Symbol('t0'), Symbol('x')), Mul(Integer(-1), Integer(36), Symbol('t0')), Mul(Integer(-1), Integer(11), Symbol('x'))), Integer(-1)))","Mul(Pow(Add(Mul(Integer(2), Symbol('t0'), Symbol('x')), Mul(Integer(-1), Integer(36), Symbol('t0')), Mul(Integer(-1), Integer(11), Symbol('x'))), Integer(-1)), Add(Mul(Integer(-1), Integer(108), Pow(Symbol('t0'), Integer(3))), Mul(Integer(-1), Integer(24), Pow(Symbol('t0'), Integer(2)), Symbol('x')), Mul(Integer(270), Pow(Symbol('t0'), Integer(2))), Mul(Integer(60), Symbol('t0'), Symbol('x')), Mul(Integer(-1), Integer(162), Symbol('t0')), Mul(Integer(-1), Integer(36), Symbol('x'))))"]}]},
    {"kind":"exp","D":["1","t0"],"rows":2,"cols":3,"A":["Integer(0)","Mul(Integer(-1), Rational(1, 6), Symbol('t0'), Symbol('x'))","Integer(0)","Integer(2)","Integer(0)","Symbol('x')"],"u":["Mul(Add(Mul(Rational(1, 2), Symbol('t0')), Integer(-1)), Pow(Add(Mul(Integer(2), Symbol('x')), Integer(2)), Integer(-1)))","Integer(0)"],"inA":"ZZ(x,t0)","inu":"ZZ(x,t0)","out":[{"shape":[3,3],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)","Integer(0)","Integer(1)"]},{"shape":[3,1],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(0)","Mul(Add(Integer(6), Mul(Integer(-1), Integer(3), Symbol('t0'))), Pow(Add(Mul(Integer(2), Symbol('t0'), Pow(Symbol('x'), Integer(2))), Mul(Integer(2), Symbol('t0'), Symbol('x'))), Integer(-1)))","Integer(0)"]}]},
    {"kind":"exp","D":["1","t0"],"rows":1,"cols":1,"A":["Mul(Rational(1, 2), Pow(Symbol('t0'), Integer(-1)), Pow(Symbol('x'), Integer(-1)), Add(Mul(Rational(1, 2), Symbol('x')), Integer(-1)))"],"u":["Integer(0)"],"inA":"ZZ(x,t0)","inu":"QQ","out":[{"shape":[1,1],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(1)"]},{"shape":[1,1],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(0)"]}]},
    {"kind":"exp","D":["1","t0"],"rows":1,"cols":1,"A":["Integer(-2)"],"u":["Integer(0)"],"inA":"QQ","inu":"QQ","out":[{"shape":[1,1],"domain":"QQ","gens":["t1"],"items":["Integer(1)"]},{"shape":[1,1],"domain":"QQ","gens":["t1"],"items":["Integer(0)"]}]},
    {"kind":"primitive","D":["1","1/x"],"rows":3,"cols":2,"A":["Mul(Symbol('t0'), Symbol('x'))","Mul(Rational(1, 2), Symbol('t0'), Symbol('x'))","Rational(-6, 5)","Rational(1, 3)","Rational(-2, 3)","Mul(Integer(2), Symbol('t0'), Symbol('x'))"],"u":["Integer(0)","Integer(0)","Integer(0)"],"inA":"ZZ(x,t0)","inu":"QQ","out":[{"shape":[3,2],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)"]},{"shape":[3,1],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(0)","Integer(0)","Integer(0)"]}]},
    {"kind":"primitive","D":["1","1/x"],"rows":1,"cols":1,"A":["Rational(1, 3)"],"u":["Integer(0)"],"inA":"QQ","inu":"QQ","out":[{"shape":[1,1],"domain":"QQ","gens":["t1"],"items":["Integer(1)"]},{"shape":[1,1],"domain":"QQ","gens":["t1"],"items":["Integer(0)"]}]},
    {"kind":"primitive","D":["1","1/x"],"rows":3,"cols":1,"A":["Mul(Integer(-1), Rational(5, 2), Symbol('x'))","Mul(Integer(-1), Integer(3), Symbol('t0'), Symbol('x'))","Mul(Integer(2), Pow(Symbol('t0'), Integer(-1)), Add(Mul(Integer(2), Symbol('t0'), Symbol('x')), Symbol('x')))"],"u":["Integer(0)","Add(Mul(Integer(-1), Integer(3), Symbol('t0'), Symbol('x')), Mul(Integer(2), Symbol('x')))","Integer(1)"],"inA":"ZZ(x,t0)","inu":"ZZ(x,t0)","out":[{"shape":[3,1],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)"]},{"shape":[3,1],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(0)","Integer(1)","Integer(0)"]}]},
    {"kind":"primitive","D":["1","1/x"],"rows":2,"cols":1,"A":["Mul(Integer(-1), Rational(1, 3), Pow(Symbol('t0'), Integer(-1)), Add(Mul(Integer(2), Symbol('t0'), Symbol('x')), Integer(2)))","Mul(Integer(-1), Rational(2, 3), Symbol('x'))"],"u":["Mul(Integer(-1), Rational(2, 3), Symbol('x'), Pow(Add(Rational(1, 2), Mul(Integer(-1), Integer(3), Symbol('x'))), Integer(-1)))","Integer(0)"],"inA":"ZZ(x,t0)","inu":"ZZ(x)","out":[{"shape":[2,1],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)"]},{"shape":[2,1],"domain":"ZZ(x,t0)","gens":["t1"],"items":["Integer(0)","Integer(1)"]}]},
    {"kind":"primitive","D":["1","1/x"],"rows":1,"cols":3,"A":["Integer(-1)","Integer(2)","Integer(3)"],"u":["Integer(0)"],"inA":"QQ","inu":"QQ","out":[{"shape":[1,3],"domain":"QQ","gens":["t1"],"items":["Integer(1)","Integer(-2)","Integer(-3)"]},{"shape":[1,1],"domain":"QQ","gens":["t1"],"items":["Integer(0)"]}]},
    {"kind":"primitive","D":["1","1/x"],"rows":3,"cols":2,"A":["Symbol('t0')","Mul(Integer(-1), Rational(1, 2), Symbol('t0'))","Add(Mul(Rational(1, 3), Symbol('t0'), Symbol('x')), Rational(2, 9))","Mul(Rational(1, 2), Symbol('x'), Pow(Add(Mul(Integer(-1), Integer(3), Symbol('t0'), Symbol('x')), Integer(2)), Integer(-1)))","Mul(Integer(-1), Integer(3), Symbol('t0'))","Symbol('t0')"],"u":["Integer(0)","Integer(0)","Integer(0)"],"inA":"ZZ(x,t0)","inu":"QQ","out":[{"shape":[3,2],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(1)","Integer(0)","Integer(0)","Integer(1)","Integer(0)","Integer(0)"]},{"shape":[3,1],"domain":"QQ(x,t0)","gens":["t1"],"items":["Integer(0)","Integer(0)","Integer(0)"]}]},
]


def _expr(s):
    return eval(s, dict(_NS))


def _de(case):
    if case['kind'] == 'base':
        return DifferentialExtension(extension={'D': [Poly(1, x)]})
    d = {'t0': t0, '1/x': 1/x}[case['D'][1]]
    return DifferentialExtension(extension={'D': [Poly(1, x), Poly(d, t0)]})


def _dom(s):
    # Symox orders the generators of an inferred fraction field by
    # default_sort_key (ZZ(t0,x)) where SymPy uses _sort_gens (ZZ(x,t0)): a
    # pre-existing difference of the Poly constructor, shared by the stock
    # body, recorded in docs/algorithmic-parity.md (prdeconst).
    if '(' not in s:
        return s
    ground, rest = s[:-1].split('(')
    return ground, sorted(rest.split(','))


def _check(M, ref):
    assert list(M.shape) == ref['shape']
    assert _dom(str(M.domain)) == _dom(ref['domain'])
    assert [str(g) for g in M.gens] == ref['gens']
    assert [srepr(e) for e in M.to_Matrix()] == ref['items']


def _key(M):
    return (M.shape, M.gens, str(M.domain), type(M.domain), [srepr(e) for e in M._items])


def _inputs(case):
    A = Matrix(case['rows'], case['cols'], [_expr(a) for a in case['A']], t1)
    u = Matrix(case['rows'], 1, [_expr(b) for b in case['u']], t1)
    return A, u


def test_generated_cases_match_sympy_and_use_native():
    for case in CASES:
        A, u = _inputs(case)
        assert (_dom(str(A.domain)), _dom(str(u.domain))) == (_dom(case['inA']), _dom(case['inu']))
        DE = _de(case)
        before = _prde_native.STATS['constant_system_field']
        a, v = constant_system(A, u, DE)
        _check(a, case['out'][0])
        _check(v, case['out'][1])
        if case['inA'] != 'QQ' or case['inu'] != 'QQ':
            assert _prde_native.STATS['constant_system_field'] == before + 1


def test_generated_cases_match_stock_body(monkeypatch):
    native = [constant_system(*_inputs(case), _de(case)) for case in CASES]
    monkeypatch.setattr(_prde_native, 'constant_system', lambda A, u, DE: None)
    for case, (a, v) in zip(CASES, native):
        b, w = constant_system(*_inputs(case), _de(case))
        assert (_key(a), _key(v)) == (_key(b), _key(w))


def test_sympy_constant_system():
    # sympy/integrals/tests/test_prde.py::test_constant_system
    A = Matrix([[-(x + 3)/(x - 1), (x + 1)/(x - 1), 1],
                [-x - 3, x + 1, x - 1],
                [2*(x + 3)/(x - 1), 0, 0]], t)
    u = Matrix([[(x + 1)/(x - 1)], [x + 1], [0]], t)
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    R = QQ.frac_field(x)[t]
    before = _prde_native.STATS['constant_system_field']
    assert constant_system(A, u, DE) ==         (Matrix([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 0],
                 [0, 0, 1]], ring=R), Matrix([0, 1, 0, 0], ring=R))
    assert _prde_native.STATS['constant_system_field'] == before + 1
    a, v = constant_system(A, u, DE)
    assert (str(a.domain), str(v.domain), a.gens, v.gens) == ('QQ(x)', 'QQ(x)', (t,), (t,))


def test_sympy_tower_callers():
    # sympy/integrals/tests/test_prde.py::test_is_deriv_k and
    # ::test_is_log_deriv_k_t_radical_in_field: their constant systems over
    # ZZ(x), QQ(x, t1), ... are answered by the fraction-field kernel.
    from symox import S
    t2 = Symbol('t2')
    before = _prde_native.STATS['constant_system_field']
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)],
        'exts': [None, 'log', 'log'], 'extargs': [None, x, x + 1]})
    assert is_deriv_k(Poly(2*x**2 + 2*x, t2), Poly(1, t2), DE) == \
        ([(t1, 1), (t2, 1)], t1 + t2, 2)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(t2, t2)],
        'exts': [None, 'log', 'exp'], 'extargs': [None, x, x]})
    assert is_deriv_k(Poly(x**2*t2**3, t2), Poly(1, t2), DE) == \
        ([(x, 3), (t1, 2)], 2*t1 + 3*x, 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/x, t1)],
        'exts': [None, 'log'], 'extargs': [None, x**2]})
    assert is_deriv_k(Poly(x, t1), Poly(1, t1), DE) == \
        ([(t1, S.Half)], t1/2, 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-1/x, t)],
        'exts': [None, 'log'], 'extargs': [None, 1/x]})
    assert is_deriv_k(Poly(1, t), Poly(x, t), DE) == ([(t, 1)], t, 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(5*t + 1, t), Poly(2*t*x, t), DE) == \
        (2, t*x**5)
    assert is_log_deriv_k_t_radical_in_field(Poly(2 + 3*t, t), Poly(5*x*t, t), DE) == \
        (5, x**3*t**2)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-t/x**2, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(-(1 + 2*t), t),
        Poly(2*x**2 + 2*x**2*t, t), DE) == (2, t + t**2)
    assert is_log_deriv_k_t_radical_in_field(Poly(-1, t), Poly(x**2, t), DE) == (1, t)
    assert is_log_deriv_k_t_radical_in_field(Poly(1, t), Poly(2*x**2, t), DE) == (2, 1/t)
    assert _prde_native.STATS['constant_system_field'] > before
