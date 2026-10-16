"""Run the consumer's real symbolic modules with either symbolic backend.

Set TRISPLYNE_PATH to a checkout; a sibling checkout is detected automatically.
The source is loaded into isolated namespaces, changing only import names.
No files in the consumer checkout are changed and no CUDA runtime is needed.
"""

import logging
import os
from pathlib import Path
import pickle
import subprocess
import sys
import types

import numpy as np
import pytest

import symox


@pytest.fixture(scope="module")
def consumers():
    root = Path(os.environ.get("TRISPLYNE_PATH", Path(__file__).resolve().parents[2] / "trisplyne"))
    source_root = root / "trisplyne"
    if not (source_root / "symbolic" / "triangle.py").exists():
        pytest.skip("set TRISPLYNE_PATH to run the real consumer integration tests")
    sympy = pytest.importorskip("sympy")
    registered = []
    result = []
    try:
        for backend in (symox, sympy):
            prefix = "_trisplyne_" + backend.__name__
            for name in (prefix, prefix + ".symbolic", prefix + ".utils", prefix + ".codegen"):
                mod = types.ModuleType(name)
                mod.__path__ = []
                sys.modules[name] = mod
                registered.append(name)
            logger = types.ModuleType(prefix + ".utils.logger")
            logger.logger = logging.getLogger(prefix)
            sys.modules[logger.__name__] = logger
            registered.append(logger.__name__)
            loaded = {}
            for relative in (
                "utils.types", "utils.math", "utils.bidirectional_map",
                "symbolic.indexing", "symbolic.utils", "symbolic.horner",
                "symbolic.triangle", "symbolic.bform", "symbolic.variables", "codegen.printers",
            ):
                path = source_root.joinpath(*relative.split(".")).with_suffix(".py")
                source = path.read_text().replace("trisplyne.", prefix + ".")
                for original in ("sympy", "symox"):
                    source = source.replace("import " + original + " as sm", "import " + backend.__name__ + " as sm")
                    source = source.replace("from " + original + ".", "from " + backend.__name__ + ".")
                name = prefix + "." + relative
                module = types.ModuleType(name)
                module.__file__ = str(path)
                sys.modules[name] = module
                registered.append(name)
                exec(compile(source, str(path), "exec"), module.__dict__)
                loaded[relative] = module
            result.append((backend, loaded))
        yield result
    finally:
        for name in reversed(registered):
            sys.modules.pop(name, None)


def test_real_horner_and_printer(consumers):
    outputs = []
    for sm, modules in consumers:
        a, b, c, d, e = sm.symbols("a b c d e")
        cases = [
            a, a + b, a*b + c*d, (a+b)**2,
            a**3*b**2 + c**2*a*b + d*e,
            a**2*b**2*c**2 + a*b*c + d*e,
            a**3 + b**3 + c**3 - 3*a*b*c,
            (2*c + 3*d)*a**2 + (c-d)*a*b - (c+d)*b,
        ]
        printer = modules["codegen.printers"].CXXCodePrinter()
        strings = []
        for expr in cases:
            horner = modules["symbolic.horner"].greedy_multivar_horner(expr, [a, b])
            evaluated = modules["symbolic.utils"].remove_unevaluated_expr(horner)
            assert sm.expand(evaluated - expr) == 0
            restored = pickle.loads(pickle.dumps(horner))
            assert restored == horner
            strings.append(printer.doprint(horner))
        outputs.append(strings)
    assert outputs[0] == outputs[1]


def test_real_psubs(consumers):
    outputs = []
    for sm, modules in consumers:
        a, b, c, d = sm.symbols("a b c d")
        psubs = modules["symbolic.utils"].psubs
        cases = [
            (a**3*b**2 + a*b*c + d, {a*b: c}),
            ((a+b)**3 + d, {a+b: -c}),
            (a**2 + a*b, {a+b: -c}),
        ]
        outputs.append([str(psubs(expr, mapping).expand()) for expr, mapping in cases])
    assert outputs[0] == outputs[1]


def test_real_triangle(consumers):
    outputs = []
    for sm, modules in consumers:
        triangle = modules["symbolic.triangle"].SymbolicTriangle()
        printer = modules["codegen.printers"].CXXCodePrinter()
        strings = []
        for vertex in range(3):
            triangle.configure_substitutions(vertex)
            for vector in (triangle.Bx, triangle.By):
                assert all(isinstance(x, sm.Expr) for x in vector)
                for expr in vector:
                    assert pickle.loads(pickle.dumps(expr)) == expr
                    strings.append(printer.doprint(expr))
        for expr in triangle.edge_normal_vcoords.ravel():
            strings.append(printer.doprint(triangle.variables[expr]))
        utils = modules["symbolic.utils"]
        a = triangle.edge_normal_vcoords
        definitions = {var: triangle.variables[var] for var in a.ravel()}
        edges = triangle.edges.ravel()
        edge_pairs = {
            edges[i] * edges[j]: sm.Symbol(f"ee_{i}{j}")
            for i in range(len(edges)) for j in range(i, len(edges))
        }
        # The consumer relies on the simplified nine ratios being binary sums
        # or products before substituting edge products in each part.
        for i in range(3):
            lhs, rhs = utils.split_two_terms(definitions[a[i, i]])
            if utils.expr_has_negative_sign(rhs):
                definitions[a[i, i]] = sm.Mul(1 / (-(1 / lhs)), -rhs)
        for i in range(3):
            for j in range(3):
                expr = (1 if i == j else definitions[a[i, j]]) / definitions[a[i, i]]
                lhs, rhs = utils.split_two_terms(expr.simplify())
                strings.append(printer.doprint(lhs.subs(edge_pairs) * rhs.subs(edge_pairs)))
        outputs.append(strings)
    assert outputs[0] == outputs[1]


def test_numpy_object_arrays_and_symbol_cache():
    x, y, z = symox.symbols("x y z", real=True)
    a = np.array([x, y, z], dtype=object)
    b = np.array([y, z, x], dtype=object)
    assert np.dot(a, b) == x*y + y*z + z*x
    assert (a * b).tolist() == [x*y, y*z, z*x]
    expr = ((x + y + z)**12).expand()
    symbols = expr.free_symbols
    symbols.clear()
    # Each Python set is independent even though the Rust discovery is cached.
    assert expr.free_symbols == {x, y, z}
    replaced = expr.xreplace({x: 1, y: 2, z: 3})
    assert replaced == 6**12


def test_real_casteljau_and_variable_declarations(consumers):
    outputs = []
    for sm, modules in consumers:
        b = sm.symbols("b:3", real=True)
        bform = modules["symbolic.bform"]
        coeffs = bform.SymbolicBFormCoefficients(4)
        # eval_casteljau only needs degree, coefficients and index tables.
        form = object.__new__(bform.SymbolicBForm)
        form.coeffs = coeffs
        form.casteljau_indices = modules["symbolic.indexing"].generate_casteljau_indices(4)
        expr = form.eval_casteljau([b] * 4).expand()
        assert len(expr.args) == 15
        derivatives = [expr.diff(x).expand() for x in b]
        defs = modules["symbolic.variables"].VariableDefinitions()
        group = defs.new_variable_group("value and gradient")
        result = sm.Symbol("result")
        group[result] = expr
        for i, derivative in enumerate(derivatives):
            group[sm.Symbol("grad" + str(i))] = derivative
        group.horner(list(b))
        printer = modules["codegen.printers"].CXXCodePrinter()
        outputs.append("\n".join(printer.doprint(block) for block in defs.as_declarations()))
    assert outputs[0] == outputs[1]


def test_thousands_of_polynomial_terms():
    x = symox.symbols("x:6")
    expr = (sum(x)**10).expand()
    assert len(expr.args) == 3003
    assert expr.free_symbols == set(x)
    assert pickle.loads(pickle.dumps(expr)) == expr
    assert expr.xreplace(dict.fromkeys(x, 1)) == 6**10
    assert expr.diff(x[0]).xreplace(dict.fromkeys(x, 1)) == 10 * 6**9


def test_compatibility_without_sympy_at_runtime():
    script = """
import importlib.abc
import pickle
import sys
class BlockSympy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('runtime imported SymPy: ' + fullname)
sys.meta_path.insert(0, BlockSympy())
import symox as sm
from symox.printing.cxx import CXX17CodePrinter
from symox.codegen.ast import Assignment
x, y = sm.symbols('x y', real=True)
held = x * sm.UnevaluatedExpr(y + 1)
assert pickle.loads(pickle.dumps(held)) == held
assert sm.simplify(held - x*(y+1)) == 0
assert sm.reduced(x*x + x*y, [x+y]) == ([x], sm.Integer(0))
assert sm.lambdify((x, y), held, 'numpy')(2, 3) == 8
assert CXX17CodePrinter().doprint(Assignment(x, y+1)) == 'x = y + 1;'
"""
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
