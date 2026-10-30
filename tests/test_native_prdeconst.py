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
    # ZZ(x), QQ(x, t1), ... are answered by the fraction-field kernel
    # (inside the native is_deriv_k).
    from symox import S
    t2 = Symbol('t2')
    before = _prde_native.STATS['is_deriv_k']
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
    assert _prde_native.STATS['is_deriv_k'] > before


# derivation(p, DE, coefficientD) with basic=False: SymPy 1.14 reference
# (case, p, input domain or None for the inferred one, coefficientD, srepr of
# the result expression, result domain); the generators are always (t,).
x_, a_ = x, Symbol('a')
_DERIV_DE = {
    'exp': lambda: DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]}),
    'prim': lambda: DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]}),
    'tower': lambda: DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t0), Poly(t, t)]}),
    'tower2': lambda: DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0, t0), Poly(1/(x + 1), t)]}),
}
DERIV_CASES = [
    ('exp', 'x*t**2 + t + 1', None, False, "Add(Mul(Pow(Symbol('t'), Integer(2)), Add(Mul(Integer(2), Symbol('x')), Integer(1))), Symbol('t'))", 'ZZ[x]'),
    ('exp', 'x*t**2 + t + 1', None, True, "Pow(Symbol('t'), Integer(2))", 'ZZ'),
    ('exp', 't**3/x + x**2*t', 'ZZ(x)', False, "Add(Mul(Pow(Symbol('t'), Integer(3)), Pow(Symbol('x'), Integer(-2)), Add(Mul(Integer(3), Symbol('x')), Integer(-1))), Mul(Symbol('t'), Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')))))", 'ZZ(x)'),
    ('exp', 't**3/x + x**2*t', 'ZZ(x)', True, "Add(Mul(Integer(-1), Pow(Symbol('t'), Integer(3)), Pow(Symbol('x'), Integer(-2))), Mul(Integer(2), Symbol('t'), Symbol('x')))", 'ZZ(x)'),
    ('exp', '(x**2 + x)/(x + 1)*t + 1/x', 'ZZ(x)', False, "Add(Mul(Symbol('t'), Add(Symbol('x'), Integer(1))), Mul(Integer(-1), Pow(Symbol('x'), Integer(-2))))", 'ZZ(x)'),
    ('exp', '(x**2 + x)/(x + 1)*t + 1/x', 'ZZ(x)', True, "Add(Symbol('t'), Mul(Integer(-1), Pow(Symbol('x'), Integer(-2))))", 'ZZ(x)'),
    ('exp', 'x**2*t/2 + 3', 'QQ[x]', False, "Mul(Symbol('t'), Add(Mul(Rational(1, 2), Pow(Symbol('x'), Integer(2))), Symbol('x')))", 'QQ[x]'),
    ('exp', 'x**2*t/2 + 3', 'QQ[x]', True, "Mul(Symbol('t'), Symbol('x'))", 'ZZ[x]'),
    ('exp', 'a*x*t + a', None, False, "Mul(Symbol('t'), Add(Mul(Symbol('a'), Symbol('x')), Symbol('a')))", 'ZZ[x,a]'),
    ('exp', 'a*x*t + a', None, True, "Mul(Symbol('a'), Symbol('t'))", 'ZZ[a]'),
    ('prim', 'x*t**2 + t/x', 'ZZ(x)', False, "Add(Pow(Symbol('t'), Integer(2)), Mul(Symbol('t'), Pow(Symbol('x'), Integer(-2)), Add(Mul(Integer(2), Pow(Symbol('x'), Integer(2))), Integer(-1))), Pow(Symbol('x'), Integer(-2)))", 'ZZ(x)'),
    ('prim', 'x*t**2 + t/x', 'ZZ(x)', True, "Add(Pow(Symbol('t'), Integer(2)), Mul(Integer(-1), Symbol('t'), Pow(Symbol('x'), Integer(-2))))", 'ZZ(x)'),
    ('prim', 't**2 + 1', None, False, "Mul(Integer(2), Symbol('t'), Pow(Symbol('x'), Integer(-1)))", 'ZZ(x)'),
    ('prim', 't**2 + 1', None, True, 'Integer(0)', 'ZZ'),
    ('prim', '(x + 1)*t', 'ZZ[x]', False, "Add(Symbol('t'), Mul(Pow(Symbol('x'), Integer(-1)), Add(Symbol('x'), Integer(1))))", 'ZZ(x)'),
    ('prim', '(x + 1)*t', 'ZZ[x]', True, "Symbol('t')", 'ZZ'),
    ('prim', 'x**2*t/3', 'QQ(x)', False, "Add(Mul(Rational(2, 3), Symbol('t'), Symbol('x')), Mul(Rational(1, 3), Symbol('x')))", 'QQ[x]'),
    ('prim', 'x**2*t/3', 'QQ(x)', True, "Mul(Rational(2, 3), Symbol('t'), Symbol('x'))", 'QQ[x]'),
    ('tower', 'x*t**2 + t0*t + 1/x', 'ZZ(x,t0)', False, "Add(Mul(Pow(Symbol('t'), Integer(2)), Add(Mul(Integer(2), Symbol('x')), Integer(1))), Mul(Symbol('t'), Pow(Symbol('x'), Integer(-1)), Add(Mul(Symbol('t0'), Symbol('x')), Integer(1))), Mul(Integer(-1), Pow(Symbol('x'), Integer(-2))))", 'ZZ(x,t0)'),
    ('tower', 'x*t**2 + t0*t + 1/x', 'ZZ(x,t0)', True, "Add(Pow(Symbol('t'), Integer(2)), Mul(Symbol('t'), Pow(Symbol('x'), Integer(-1))), Mul(Integer(-1), Pow(Symbol('x'), Integer(-2))))", 'ZZ(x)'),
    ('tower', 't0**2*t + x', 'ZZ[x,t0]', False, "Add(Mul(Symbol('t'), Pow(Symbol('x'), Integer(-1)), Add(Mul(Pow(Symbol('t0'), Integer(2)), Symbol('x')), Mul(Integer(2), Symbol('t0')))), Integer(1))", 'ZZ(x,t0)'),
    ('tower', 't0**2*t + x', 'ZZ[x,t0]', True, "Add(Mul(Integer(2), Symbol('t'), Symbol('t0'), Pow(Symbol('x'), Integer(-1))), Integer(1))", 'ZZ(x,t0)'),
    ('tower', 't0*t/x', 'ZZ(x,t0)', False, "Mul(Symbol('t'), Pow(Symbol('x'), Integer(-2)), Add(Mul(Symbol('t0'), Symbol('x')), Mul(Integer(-1), Symbol('t0')), Integer(1)))", 'ZZ(x,t0)'),
    ('tower', 't0*t/x', 'ZZ(x,t0)', True, "Mul(Symbol('t'), Pow(Symbol('x'), Integer(-2)), Add(Integer(1), Mul(Integer(-1), Symbol('t0'))))", 'ZZ(x,t0)'),
    ('tower2', 't0*t**2 + x*t', 'ZZ[x,t0]', False, "Add(Mul(Pow(Symbol('t'), Integer(2)), Symbol('t0')), Mul(Symbol('t'), Pow(Add(Symbol('x'), Integer(1)), Integer(-1)), Add(Mul(Integer(2), Symbol('t0')), Symbol('x'), Integer(1))), Mul(Symbol('x'), Pow(Add(Symbol('x'), Integer(1)), Integer(-1))))", 'ZZ(x,t0)'),
    ('tower2', 't0*t**2 + x*t', 'ZZ[x,t0]', True, "Add(Mul(Pow(Symbol('t'), Integer(2)), Symbol('t0')), Symbol('t'))", 'ZZ[t0]'),
    ('tower2', 't/(x + 1) + t0', 'ZZ(x,t0)', False, "Add(Mul(Integer(-1), Symbol('t'), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1)), Integer(-1))), Mul(Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1)), Integer(-1)), Add(Mul(Symbol('t0'), Pow(Symbol('x'), Integer(2))), Mul(Integer(2), Symbol('t0'), Symbol('x')), Symbol('t0'), Integer(1))))", 'ZZ(x,t0)'),
    ('tower2', 't/(x + 1) + t0', 'ZZ(x,t0)', True, "Add(Mul(Integer(-1), Symbol('t'), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1)), Integer(-1))), Symbol('t0'))", 'ZZ(x,t0)'),
]


def _pdom(s):
    # Same generator-order normalisation as _dom, for rings and fields.
    for o, c in (('(', ')'), ('[', ']')):
        if o in s:
            ground, rest = s[:-1].split(o)
            return ground, o, sorted(rest.split(','))
    return s


def _deriv_input(case, e, dom):
    from symox import sympify
    DE = _DERIV_DE[case]()
    p = sympify(e, locals={'x': x, 't': t, 't0': t0, 'a': a_})
    return (Poly(p, DE.t, domain=dom) if dom else Poly(p, DE.t)), DE


def test_derivation_matches_sympy_and_uses_native():
    from symox.integrals.risch import derivation
    before = _prde_native.STATS['derivation']
    for case, e, dom, cd, out, odom in DERIV_CASES:
        p, DE = _deriv_input(case, e, dom)
        level = DE.level
        r = derivation(p, DE, coefficientD=cd)
        assert DE.level == level
        assert srepr(r.as_expr()) == out, (case, e, cd)
        assert _pdom(str(r.domain)) == _pdom(odom), (case, e, cd, r.domain)
        assert r.gens == (t,)
    assert _prde_native.STATS['derivation'] == before + len(DERIV_CASES)


def test_derivation_matches_stock_body(monkeypatch):
    from symox.integrals.risch import derivation
    native = []
    for case, e, dom, cd, out, odom in DERIV_CASES:
        p, DE = _deriv_input(case, e, dom)
        native.append(derivation(p, DE, coefficientD=cd))
    monkeypatch.setattr(_prde_native, 'derivation', lambda p, DE, coefficientD: None)
    for (case, e, dom, cd, out, odom), r in zip(DERIV_CASES, native):
        p, DE = _deriv_input(case, e, dom)
        s = derivation(p, DE, coefficientD=cd)
        assert (srepr(r.as_expr()), r.gens, r.domain, str(r.domain)) == \
            (srepr(s.as_expr()), s.gens, s.domain, str(s.domain))


def test_derivation_basic_and_declined_inputs():
    from symox import cancel, exp, log
    from symox.integrals.risch import derivation
    DE = _DERIV_DE['prim']()
    before = _prde_native.STATS['derivation']
    # basic=True stays on the stock body (an expression, canceled).
    assert derivation(x*t**2 + 1/x, DE, basic=True) == cancel(t**2 + 2*t - 1/x**2)
    # EX coefficients are declined: stock body.
    p = Poly(exp(x)*t, t, domain='EX')
    assert derivation(p, DE) == Poly(exp(x)*t + exp(x)/x, t, domain='EX')
    assert _prde_native.STATS['derivation'] == before


# is_deriv_k(fa, fd, DE): SymPy 1.14 reference (srepr of ans, u and const,
# None, or the raised NotImplementedError), SymPy test_is_deriv_k inputs plus
# generated logarithmic, exponential, mixed and non-rational cases.
_a2, _t2 = symbols('a t2')


def _dk_de(D, exts=None, extargs=None):
    ext = {'D': D}
    if exts is not None:
        ext['exts'] = exts
        ext['extargs'] = extargs
    return DifferentialExtension(extension=ext)


def _dk_show(r):
    if r is None:
        return 'None'
    ans, u, c = r
    return repr(([(srepr(i), srepr(j)) for i, j in ans], srepr(u), srepr(c)))


DK_CASES = []


def _dk_case(name, f):
    DK_CASES.append((name, f))


def _dk_cases():
    x, a, t, t0, t1, t2 = symbols('x a t t0 t1 t2')
    DE = _dk_de

    L2 = lambda: DE([Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)], [None, 'log', 'log'], [None, x, x + 1])
    _dk_case('sym1', lambda: is_deriv_k(Poly(2*x**2 + 2*x, t2), Poly(1, t2), L2()))
    LE = lambda: DE([Poly(1, x), Poly(1/x, t1), Poly(t2, t2)], [None, 'log', 'exp'], [None, x, x])
    _dk_case('sym2', lambda: is_deriv_k(Poly(x**2*t2**3, t2), Poly(1, t2), LE()))
    _dk_case('sym3', lambda: is_deriv_k(Poly(x, t1), Poly(1, t1), DE([Poly(1, x), Poly(2/x, t1)], [None, 'log'], [None, x**2])))
    _dk_case('sym4', lambda: is_deriv_k(Poly(1 + x, t0), Poly(1, t0), DE([Poly(1, x), Poly(2/(1 + x), t0)], [None, 'log'], [None, x**2 + 2*x + 1])))
    _dk_case('sym5', lambda: is_deriv_k(Poly(1, t), Poly(x, t), DE([Poly(1, x), Poly(-1/x, t)], [None, 'log'], [None, 1/x])))
    LX = lambda: DE([Poly(1, x), Poly(1/x, t)], [None, 'log'], [None, x])
    _dk_case('log_x3', lambda: is_deriv_k(Poly(x**3, t), Poly(1, t), LX()))
    _dk_case('log_2x', lambda: is_deriv_k(Poly(2*x, t), Poly(1, t), LX()))
    _dk_case('log_t', lambda: is_deriv_k(Poly(t, t), Poly(1, t), LX()))
    _dk_case('log_xp1', lambda: is_deriv_k(Poly(x + 1, t), Poly(1, t), LX()))
    _dk_case('log_const', lambda: is_deriv_k(Poly(3, t), Poly(1, t), LX()))
    _dk_case('log_ax', lambda: is_deriv_k(Poly(a*x, t), Poly(1, t), LX()))
    _dk_case('log_inv2', lambda: is_deriv_k(Poly(1, t), Poly(x**2, t), LX()))
    _dk_case('log_frac', lambda: is_deriv_k(Poly(5*x**3, t), Poly(7*x, t), LX()))
    EX = lambda: DE([Poly(1, x), Poly(t, t)], [None, 'exp'], [None, x])
    _dk_case('exp_t', lambda: is_deriv_k(Poly(t, t), Poly(1, t), EX()))
    _dk_case('exp_t2x', lambda: is_deriv_k(Poly(x*t**2, t), Poly(1, t), EX()))
    _dk_case('exp_inv', lambda: is_deriv_k(Poly(3, t), Poly(t**2, t), EX()))
    _dk_case('exp_2x', lambda: is_deriv_k(Poly(t, t), Poly(1, t), DE([Poly(1, x), Poly(2*t, t)], [None, 'exp'], [None, 2*x])))
    _dk_case('exp_inv_x', lambda: is_deriv_k(Poly(t, t), Poly(1, t), DE([Poly(1, x), Poly(t/x**2, t)], [None, 'exp'], [None, -1/x])))
    _dk_case('mix1', lambda: is_deriv_k(Poly(x**2*t2, t2), Poly(1, t2), LE()))
    _dk_case('mix2', lambda: is_deriv_k(Poly(t2, t2), Poly(x, t2), LE()))
    _dk_case('mix3', lambda: is_deriv_k(Poly(t2**2, t2), Poly(x**3, t2), LE()))
    _dk_case('two_logs', lambda: is_deriv_k(Poly(x**2*(x + 1)**3, t2), Poly(1, t2), L2()))
    _dk_case('two_logs_q', lambda: is_deriv_k(Poly(x**2, t2), Poly((x + 1)**5, t2), L2()))
    _dk_case('two_logs_none', lambda: is_deriv_k(Poly(x + 2, t2), Poly(1, t2), L2()))
    _dk_case('log_x2_q', lambda: is_deriv_k(Poly(x**3, t1), Poly(1, t1), DE([Poly(1, x), Poly(2/x, t1)], [None, 'log'], [None, x**2])))
    _dk_case('nonrational', lambda: is_deriv_k(Poly(x, t), Poly(1, t), DE([Poly(1, x), Poly(a/x, t)], [None, 'log'], [None, x**a])))
    _dk_case('rational_x2', lambda: is_deriv_k(Poly(x**2 + 1, t), Poly(1, t), DE([Poly(1, x), Poly(2*x/(x**2 + 1), t)], [None, 'log'], [None, x**2 + 1])))
    _dk_case('rational_sq', lambda: is_deriv_k(Poly(x**2 + 1, t), Poly(1, t), DE([Poly(1, x), Poly(4*x/(x**2 + 1), t)], [None, 'log'], [None, (x**2 + 1)**2])))


_dk_cases()
DK_REF = {
    'sym1': '([("Symbol(\'t1\')", \'Integer(1)\'), ("Symbol(\'t2\')", \'Integer(1)\')], "Add(Symbol(\'t1\'), Symbol(\'t2\'))", \'Integer(2)\')',
    'sym2': '([("Symbol(\'x\')", \'Integer(3)\'), ("Symbol(\'t1\')", \'Integer(2)\')], "Add(Mul(Integer(2), Symbol(\'t1\')), Mul(Integer(3), Symbol(\'x\')))", \'Integer(1)\')',
    'sym3': '([("Symbol(\'t1\')", \'Rational(1, 2)\')], "Mul(Rational(1, 2), Symbol(\'t1\'))", \'Integer(1)\')',
    'sym4': '([("Symbol(\'t0\')", \'Rational(1, 2)\')], "Mul(Rational(1, 2), Symbol(\'t0\'))", \'Integer(1)\')',
    'sym5': '([("Symbol(\'t\')", \'Integer(1)\')], "Symbol(\'t\')", \'Integer(1)\')',
    'log_x3': '([("Symbol(\'t\')", \'Integer(3)\')], "Mul(Integer(3), Symbol(\'t\'))", \'Integer(1)\')',
    'log_2x': '([("Symbol(\'t\')", \'Integer(1)\')], "Symbol(\'t\')", \'Integer(2)\')',
    'log_t': 'None',
    'log_xp1': 'None',
    'log_const': '([("Symbol(\'t\')", \'Integer(0)\')], \'Integer(0)\', \'Integer(3)\')',
    'log_ax': '([("Symbol(\'t\')", \'Integer(1)\')], "Symbol(\'t\')", "Symbol(\'a\')")',
    'log_inv2': '([("Symbol(\'t\')", \'Integer(-2)\')], "Mul(Integer(-1), Integer(2), Symbol(\'t\'))", \'Integer(1)\')',
    'log_frac': '([("Symbol(\'t\')", \'Integer(2)\')], "Mul(Integer(2), Symbol(\'t\'))", \'Rational(5, 7)\')',
    'exp_t': '([("Symbol(\'x\')", \'Integer(1)\')], "Symbol(\'x\')", \'Integer(1)\')',
    'exp_t2x': 'None',
    'exp_inv': '([("Symbol(\'x\')", \'Integer(-2)\')], "Mul(Integer(-1), Integer(2), Symbol(\'x\'))", \'Integer(3)\')',
    'exp_2x': '([("Mul(Integer(2), Symbol(\'x\'))", \'Integer(1)\')], "Mul(Integer(2), Symbol(\'x\'))", \'Integer(1)\')',
    'exp_inv_x': '([("Mul(Integer(-1), Pow(Symbol(\'x\'), Integer(-1)))", \'Integer(1)\')], "Mul(Integer(-1), Pow(Symbol(\'x\'), Integer(-1)))", \'Integer(1)\')',
    'mix1': '([("Symbol(\'x\')", \'Integer(1)\'), ("Symbol(\'t1\')", \'Integer(2)\')], "Add(Mul(Integer(2), Symbol(\'t1\')), Symbol(\'x\'))", \'Integer(1)\')',
    'mix2': '([("Symbol(\'x\')", \'Integer(1)\'), ("Symbol(\'t1\')", \'Integer(-1)\')], "Add(Mul(Integer(-1), Symbol(\'t1\')), Symbol(\'x\'))", \'Integer(1)\')',
    'mix3': '([("Symbol(\'x\')", \'Integer(2)\'), ("Symbol(\'t1\')", \'Integer(-3)\')], "Add(Mul(Integer(-1), Integer(3), Symbol(\'t1\')), Mul(Integer(2), Symbol(\'x\')))", \'Integer(1)\')',
    'two_logs': '([("Symbol(\'t1\')", \'Integer(2)\'), ("Symbol(\'t2\')", \'Integer(3)\')], "Add(Mul(Integer(2), Symbol(\'t1\')), Mul(Integer(3), Symbol(\'t2\')))", \'Integer(1)\')',
    'two_logs_q': '([("Symbol(\'t1\')", \'Integer(2)\'), ("Symbol(\'t2\')", \'Integer(-5)\')], "Add(Mul(Integer(2), Symbol(\'t1\')), Mul(Integer(-1), Integer(5), Symbol(\'t2\')))", \'Integer(1)\')',
    'two_logs_none': 'None',
    'log_x2_q': '([("Symbol(\'t1\')", \'Rational(3, 2)\')], "Mul(Rational(3, 2), Symbol(\'t1\'))", \'Integer(1)\')',
    'nonrational': 'RAISE NotImplementedError Cannot work with non-rational coefficients in this case.',
    'rational_x2': '([("Symbol(\'t\')", \'Integer(1)\')], "Symbol(\'t\')", \'Integer(1)\')',
    'rational_sq': '([("Symbol(\'t\')", \'Rational(1, 2)\')], "Mul(Rational(1, 2), Symbol(\'t\'))", \'Integer(1)\')',
}
# DifferentialExtension(expr, x): (newf, T, D, exts, extargs) in SymPy 1.14;
# building these towers calls is_deriv_k from _log_part and _exp_part.
DK_DE_REF = {
    'log(x) + log(x**2)': '("Mul(Integer(3), Symbol(\'t0\'))", ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')"], [None, \'log\'], [\'None\', "Symbol(\'x\')"])',
    'log(x)*log(2*x)': '("Mul(Symbol(\'t0\'), Add(Symbol(\'t0\'), log(Integer(2))))", ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')"], [None, \'log\'], [\'None\', "Symbol(\'x\')"])',
    'log(x**2 + 2*x + 1)/log(x + 1)': '(\'Integer(2)\', ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/(x + 1), t0, domain=\'ZZ(x)\')"], [None, \'log\'], [\'None\', "Add(Symbol(\'x\'), Integer(1))"])',
    'exp(x)*exp(2*x) + log(x)': '("Add(Symbol(\'t0\'), Symbol(\'t1\'))", ["Symbol(\'x\')", "Symbol(\'t0\')", "Symbol(\'t1\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')", "Poly(3*t1, t1, domain=\'ZZ\')"], [None, \'log\', \'exp\'], [\'None\', "Symbol(\'x\')", "Mul(Integer(3), Symbol(\'x\'))"])',
    'x**x': '("Symbol(\'t1\')", ["Symbol(\'x\')", "Symbol(\'t0\')", "Symbol(\'t1\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')", "Poly((t0 + 1)*t1, t1, domain=\'ZZ[t0]\')"], [None, \'log\', \'exp\'], [\'None\', "Symbol(\'x\')", "Mul(Symbol(\'t0\'), Symbol(\'x\'))"])',
    'log(exp(x))': '("Symbol(\'x\')", ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(t0, t0, domain=\'ZZ\')"], [None, \'exp\'], [\'None\', "Symbol(\'x\')"])',
    'log(x*exp(x))': '("Symbol(\'t1\')", ["Symbol(\'x\')", "Symbol(\'t0\')", "Symbol(\'t1\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(t0, t0, domain=\'ZZ\')", "Poly((x + 1)/x, t1, domain=\'ZZ(x)\')"], [None, \'exp\', \'log\'], [\'None\', "Symbol(\'x\')", "Mul(Symbol(\'t0\'), Symbol(\'x\'))"])',
    'log(x)*log(x**3)/x': '("Mul(Integer(3), Pow(Symbol(\'t0\'), Integer(2)), Pow(Symbol(\'x\'), Integer(-1)))", ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')"], [None, \'log\'], [\'None\', "Symbol(\'x\')"])',
    'exp(log(x)/2)': "RAISE NotImplementedError Couldn't find an elementary transcendental extension for sqrt(x).  Try using a manual extension with the extension flag.",
    'log(log(x))*log(log(x**2))': '("Mul(Symbol(\'t1\'), Add(Symbol(\'t1\'), log(Integer(2))))", ["Symbol(\'x\')", "Symbol(\'t0\')", "Symbol(\'t1\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')", "Poly(1/(t0*x), t1, domain=\'ZZ(x,t0)\')"], [None, \'log\', \'log\'], [\'None\', "Symbol(\'x\')", "Symbol(\'t0\')"])',
    '2**x*log(2)': '("Mul(Symbol(\'t0\'), log(Integer(2)))", ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(log(2)*t0, t0, domain=\'ZZ[log(2)]\')"], [None, \'exp\'], [\'None\', "Mul(Symbol(\'x\'), log(Integer(2)))"])',
    'log(3*x)*exp(x)': '("Mul(Symbol(\'t0\'), Symbol(\'t1\'))", ["Symbol(\'x\')", "Symbol(\'t0\')", "Symbol(\'t1\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')", "Poly(t1, t1, domain=\'ZZ\')"], [None, \'log\', \'exp\'], [\'None\', "Mul(Integer(3), Symbol(\'x\'))", "Symbol(\'x\')"])',
    'log(x)/(log(x/2) + 1)': '("Mul(Symbol(\'t0\'), Pow(Add(Symbol(\'t0\'), Mul(Integer(-1), log(Integer(2))), Integer(1)), Integer(-1)))", ["Symbol(\'x\')", "Symbol(\'t0\')"], ["Poly(1, x, domain=\'ZZ\')", "Poly(1/x, t0, domain=\'ZZ(x)\')"], [None, \'log\'], [\'None\', "Symbol(\'x\')"])',
}


def _dk_run(f):
    try:
        return _dk_show(f())
    except NotImplementedError as e:
        return 'RAISE NotImplementedError %s' % e


def test_is_deriv_k_matches_sympy_and_uses_native():
    before = _prde_native.STATS['is_deriv_k']
    for name, f in DK_CASES:
        assert _dk_run(f) == DK_REF[name], name
    # Every case but the NotImplementedError one returns through the kernel.
    assert _prde_native.STATS['is_deriv_k'] == before + len(DK_CASES) - 1


def test_is_deriv_k_matches_stock_body(monkeypatch):
    native = [_dk_run(f) for _, f in DK_CASES]
    monkeypatch.setattr(_prde_native, 'is_deriv_k', lambda fa, fd, DE: None)
    assert [_dk_run(f) for _, f in DK_CASES] == native


def test_is_deriv_k_tower_construction_matches_sympy():
    from symox import sympify
    x = symbols('x')
    before = _prde_native.STATS['is_deriv_k']
    for e, ref in DK_DE_REF.items():
        try:
            d = DifferentialExtension(sympify(e), x)
            out = repr((srepr(d.newf), [srepr(i) for i in d.T], [str(i) for i in d.D], list(d.exts),
                        [srepr(i) for i in d.extargs]))
        except NotImplementedError as exc:
            out = 'RAISE NotImplementedError %s' % exc
        assert out == ref, e
    assert _prde_native.STATS['is_deriv_k'] > before


def test_is_deriv_k_declined_and_unsupported_inputs():
    x, t = symbols('x t')
    before = _prde_native.STATS['is_deriv_k']
    # len(DE.exts) != len(DE.D): the stock body raises.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)], 'exts': [None], 'extargs': [None]})
    try:
        is_deriv_k(Poly(x, t), Poly(1, t), DE)
    except NotImplementedError:
        pass
    else:
        raise AssertionError
    # EX coefficients are declined: stock body.
    from symox import exp
    DE = _dk_de([Poly(1, x), Poly(1/x, t)], [None, 'log'], [None, x])
    assert is_deriv_k(Poly(exp(x)*x, t, domain='EX'), Poly(1, t, domain='EX'), DE) is None
    assert _prde_native.STATS['is_deriv_k'] == before



def test_is_deriv_k_with_the_stock_bounds_fold(monkeypatch):
    # Once an AccumBounds has been created, Add and Mul consult the stock
    # coefficient fold; the kernel stays in use with it, not with another
    # fold nor with a rebound fold list.
    from symox import _compat
    from symox.calculus.accumulationbounds import _fold_coefficient
    name, f = DK_CASES[0]
    saved = _compat._accum_fold[0]
    try:
        _compat._accum_fold[0] = _fold_coefficient
        before = _prde_native.STATS['is_deriv_k']
        assert _dk_run(f) == DK_REF[name]
        assert _prde_native.STATS['is_deriv_k'] == before + 1
        _compat._accum_fold[0] = lambda op, args: _fold_coefficient(op, args)
        assert _dk_run(f) == DK_REF[name]
        assert _prde_native.STATS['is_deriv_k'] == before + 1
    finally:
        _compat._accum_fold[0] = saved
    monkeypatch.setattr(_compat, '_accum_fold', [None])
    assert _dk_run(f) == DK_REF[name]
    assert _prde_native.STATS['is_deriv_k'] == before + 1
