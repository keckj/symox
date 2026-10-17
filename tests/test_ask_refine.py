"""Native assumption queries, contexts, reconstruction and scalar refinement."""
import os
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


@pytest.mark.parametrize("predicate,expression,assumption", [
    ("real", "x", "Q.positive(x)"),
    ("positive", "x", "Q.real(x)"),
    ("negative", "x", "Q.positive(x)"),
    ("integer", "x", "Q.even(x) | Q.odd(x)"),
    ("even", "x + y", "Q.odd(x) & Q.odd(y)"),
    ("odd", "x*y", "Q.odd(x) & Q.odd(y)"),
    ("positive", "x*y", "Q.negative(x) & Q.negative(y)"),
    ("positive", "x + y", "Q.positive(x) & Q.positive(y)"),
    ("real", "sin(x)", "Q.real(x)"),
    ("real", "exp(x)", "Q.real(x)"),
    ("positive", "exp(x)", "Q.real(x)"),
    ("nonnegative", "x**2", "Q.real(x)"),
    ("rational", "x + y", "Q.rational(x) & Q.rational(y)"),
    ("irrational", "x + y", "Q.rational(x) & Q.irrational(y)"),
    ("zero", "x*y", "Q.zero(x)"),
    ("imaginary", "I*x", "Q.real(x) & Q.nonzero(x)"),
])
def test_queries_match_sympy(predicate, expression, assumption):
    results = []
    for backend in (sp, sm):
        x, y = backend.symbols("x y")
        ns = {"x": x, "y": y, "Q": backend.Q, "I": backend.I,
              "sin": backend.sin, "exp": backend.exp}
        results.append(backend.ask(getattr(backend.Q, predicate)(eval(expression, ns)),
                                   eval(assumption, ns)))
    assert results[1] is results[0]


@pytest.mark.parametrize("predicate", ["integer", "rational", "real", "complex",
    "positive", "negative", "zero", "nonzero", "finite", "infinite", "even",
    "odd", "prime", "composite", "algebraic", "imaginary", "extended_real"])
@pytest.mark.parametrize("value", ["0", "1", "-3", "Rational(1, 2)", "pi", "I", "oo", "nan"])
def test_numeric_queries(predicate, value):
    reference = sp.ask(getattr(sp.Q, predicate)(sp.sympify(value)))
    actual = sm.ask(getattr(sm.Q, predicate)(sm.sympify(value)))
    assert actual is reference


@pytest.mark.parametrize("expr,assumptions", [
    ("Abs(x)", "Q.positive(x)"),
    ("Abs(x)", "Q.negative(x)"),
    ("Abs(x)", "Q.real(x)"),
    ("sqrt(x**2)", "Q.real(x)"),
    ("sqrt(x**2)", "Q.positive(x)"),
    ("sqrt(x**2)", "Q.negative(x)"),
    ("sqrt(x**4)", "Q.real(x)"),
    ("Abs(x*y)", "Q.positive(x)"),
    ("(-1)**x", "Q.even(x)"),
    ("(-1)**x", "Q.odd(x)"),
    ("(-1)**(x+y)", "Q.even(x)"),
    ("re(x)", "Q.real(x)"),
    ("im(x)", "Q.imaginary(x)"),
    ("arg(x)", "Q.negative(x)"),
    ("atan2(y,x)", "Q.real(y) & Q.positive(x)"),
    ("Q.real(x)", "Q.positive(x)"),
    ("Q.positive(x)", "Q.real(x)"),
])
def test_refine_matches_sympy(expr, assumptions):
    values = []
    for b in (sp, sm):
        x,y=b.symbols("x y")
        ns={name:getattr(b,name) for name in ("Q","Abs","sqrt","re","im","arg","atan2")}
        ns.update(x=x,y=y)
        value = b.refine(eval(expr,ns),eval(assumptions,ns))
        values.append(value if b is sp else value.to_sympy())
    assert values[0] == values[1]


def test_context_restore_and_explicit_context():
    x=sm.Symbol("x")
    sm.global_assumptions.clear()
    assert sm.ask(sm.Q.real(x)) is None
    with sm.assuming(sm.Q.positive(x)):
        assert sm.ask(sm.Q.real(x)) is True
        assert sm.refine(sm.Abs(x)) == x
        assert sm.Abs(x).refine() == x
        assert sm.ask(sm.Q.real(x),context=sm.AssumptionsContext()) is None
        with pytest.raises(RuntimeError):
            with sm.assuming(sm.Q.integer(x)):
                assert sm.ask(sm.Q.integer(x)) is True
                raise RuntimeError("unwind")
        assert sm.ask(sm.Q.integer(x)) is None
    assert not sm.global_assumptions
    assert sm.ask(sm.Q.real(x)) is None


def test_contradictions_and_invalid_inputs():
    x=sm.Symbol("x")
    with pytest.raises(ValueError,match="[Ii]nconsistent"):
        sm.ask(sm.Q.real(x), sm.Q.positive(x) & sm.Q.negative(x))
    with pytest.raises(TypeError,match="proposition"):
        sm.ask(x)
    with pytest.raises(TypeError,match="assumptions"):
        sm.ask(sm.Q.real(x),x)
    for args in [(),(x,), (x,x,x)]:
        with pytest.raises(ValueError,match="two arguments"):
            sm.Q.eq(*args)
    with pytest.raises(TypeError,match="Predicate"):
        sm.AppliedPredicate(x,x)


def test_predicate_structure_roundtrips():
    x=sm.Dummy("x")
    for p in (sm.Q.positive,sm.Q.eq,sm.Predicate("custom")):
        assert isinstance(p,sm.Predicate)
        assert not p.is_number
        applied=p(x,0) if p==sm.Q.eq else p(x)
        assert isinstance(applied,sm.AppliedPredicate)
        assert applied.args[0] == p
        assert applied.function == p
        assert applied.func(*applied.args) == applied
        assert sm.sympify(applied.to_sympy()) == applied
        assert sm.sympify(sm.srepr(applied)) == applied
        assert pickle.loads(pickle.dumps(applied)) == applied
        assert pickle.loads(pickle.dumps(p)) == p
        assert applied.free_symbols == {x}
        assert applied.subs(x,2) == (p(2,0) if p==sm.Q.eq else p(2))
    assert sm.Predicate("positive") != sm.Q.positive
    assert sm.sympify("Q.positive(x)") == sm.Q.positive(sm.Symbol("x"))


def test_matrix_queries():
    for b in (sm,sp):
        a=b.MatrixSymbol("A",3,3)
        assert b.ask(b.Q.square(a)) is True
        assert b.ask(b.Q.invertible(a),b.Q.positive_definite(a)) is True
        assert b.ask(b.Q.symmetric(a.T),b.Q.symmetric(a)) is True
        assert b.ask(b.Q.square(b.MatrixSymbol("B",2,3))) is False


def test_no_sympy_at_runtime():
    script = '''
import sys
class BlockSympy:
    def find_spec(self, fullname, *args):
        if fullname == "sympy" or fullname.startswith("sympy."):
            raise AssertionError("SymPy runtime dependency: " + fullname)
sys.meta_path.insert(0, BlockSympy())
from symox import Q, ask, refine, Symbol, Abs, assuming
x=Symbol("x")
assert ask(Q.real(x), Q.positive(x)) is True
assert refine(Abs(x), Q.negative(x)) == -x
with assuming(Q.positive(x)):
    assert ask(Q.real(x)) is True
assert "sympy" not in sys.modules
'''
    subprocess.run([sys.executable,"-c",script],check=True,env=os.environ.copy(),timeout=30)


def test_predicate_arguments_and_printing():
    x=sm.Symbol("x")
    assert sm.Q.positive(x).arg == x
    assert sm.Q.positive(x).arguments == (x,)
    assert sm.Q.eq(x,0).arguments == (x,0)
    with pytest.raises(TypeError,match="unary"):
        _=sm.Q.eq(x,0).arg
    for e in [sm.Q.positive,sm.Q.positive(x),sm.Q.eq(x,0),sm.Predicate("test_name")]:
        assert sm.latex(e) == sp.latex(e.to_sympy())
        assert sm.pretty(e,use_unicode=False) == sp.pretty(e.to_sympy(),use_unicode=False)


@pytest.mark.parametrize("n",[4,17,18])
def test_xor_assumptions_retain_parity_for_many_operands(n):
    predicates=[sm.Predicate("p")(sm.Symbol(f"x{i}")) for i in range(n)]
    # Fix all but the final operand. XOR then determines the last one.
    assumptions=sm.And(sm.Xor(*predicates),*predicates[:-1])
    assert sm.ask(predicates[-1],assumptions) is (n%2==1)


def test_matrix_predicates_preserve_typed_parents():
    a=sm.MatrixSymbol("A",3,3)
    for applied in (sm.Q.invertible(a),sm.Q.symmetric(a.T),
                    sm.AppliedPredicate(sm.Q.invertible,sm.Matrix([[1,0],[0,1]]))):
        assert applied.func(*applied.args) == applied
        assert pickle.loads(pickle.dumps(applied)) == applied
        assert sm.sympify(applied.to_sympy()) == applied
    assert sm.ask(sm.Q.invertible(sm.eye(2))) is True


@pytest.mark.parametrize("name,operands",[
    ("eq",("x","x")),("ne",("x","x")),("gt",("2","1")),
    ("lt",("x","y")),("ge",("0","0")),("le",("oo","0")),
])
def test_binary_predicates(name,operands):
    reference=sp.ask(getattr(sp.Q,name)(*(sp.sympify(v) for v in operands)))
    assert sm.ask(getattr(sm.Q,name)(*(sm.sympify(v) for v in operands))) is reference


def test_assumptions_submodule_imports_share_public_api():
    from symox.assumptions.ask import Q, ask
    from symox.assumptions.refine import refine
    from symox.assumptions.assume import global_assumptions
    from symox.assumptions import ask as package_ask, refine as package_refine
    assert ask is package_ask is sm.ask
    assert refine is package_refine is sm.refine
    assert Q is sm.Q
    assert global_assumptions is sm.global_assumptions
