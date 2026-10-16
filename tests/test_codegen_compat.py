"""Exercise the subclassing interface used by trisplyne's actual printer."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess

import pytest

import symox as sm
from symox.codegen.ast import Assignment, CodeBlock, Comment, Indexed, IndexedBase
from symox.printing.cxx import CXX17CodePrinter
from symox.printing.precedence import PRECEDENCE, precedence
from symox.printing.pretty import PrettyPrinter
from symox.printing.str import StrPrinter


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def trisplyne_printers(tmp_path_factory):
    source = Path(os.environ.get(
        "TRISPLYNE_PRINTERS", Path(__file__).resolve().parents[2] /
        "trisplyne/trisplyne/codegen/printers.py"))
    if not source.exists():
        pytest.skip("Set TRISPLYNE_PRINTERS to run the trisplyne integration checks")
    pytest.importorskip("sympy")
    reference = _load_module(source, "_trisplyne_sympy_printers")
    target = tmp_path_factory.mktemp("trisplyne_printer") / "printers.py"
    target.write_text(source.read_text().replace("import sympy", "import symox").replace(
        "from sympy.", "from symox."))
    compatible = _load_module(target, "_trisplyne_symox_printers")
    return reference.CXXCodePrinter, compatible.CXXCodePrinter


CASES = [
    ("x**2", "x*x"),
    ("x**3", "x*x*x"),
    ("2*x**2*y", "2*x*x*y"),
    ("(x+y)**2", "(x+y)*(x+y)"),
    ("x**4", "pow(x, T(4))"),
    ("1/x", "T(1)/x"),
    ("x**-2", "T(1)/(x*x)"),
    ("y/x", "y/x"),
    ("y/x**2", "y/(x*x)"),
    ("y/(x**2*z)", "y/(x*x*z)"),
    ("(x+y)/x**2", "(x+y)/(x*x)"),
    ("1/(x+y)", "T(1)/(x+y)"),
    ("sqrt(x)", "sqrt(x)"),
    ("1/sqrt(x)", "T(1)/sqrt(x)"),
    ("x**Rational(1,3)", f"pow(x, T({(1/3).hex()}))"),
    ("x**a", "pow(x, a)"),
    ("-UnevaluatedExpr(x+y)", "-(x+y)"),
    ("-UnevaluatedExpr(-x)", "-(-x)"),
    ("Mul(-1,x+y,evaluate=False)", "-(x+y)"),
    ("-2*x*y", "-2*x*y"),
    ("-x/y", "-x/y"),
    ("Rational(3,4)*x", f"({(3/4).hex()})*x"),
    ("Float(0.5)*x", f"{(0.5).hex()}*x"),
    ("Rational(1,2)", (0.5).hex()),
    ("Rational(-1,2)", (-0.5).hex()),
]


@pytest.mark.parametrize("source,expected", CASES)
def test_trisplyne_exact_code(trisplyne_printers, source, expected):
    import sympy
    reference, compatible = trisplyne_printers
    outputs = []
    for backend, printer in [(sympy, reference), (sm, compatible)]:
        namespace = vars(backend).copy()
        namespace.update(zip("xyza", backend.symbols("x y z a")))
        outputs.append(printer().doprint(eval(source, namespace)))
    assert outputs == [expected, expected]


def test_trisplyne_assignments_and_symbols(trisplyne_printers):
    _, compatible = trisplyne_printers
    x, y = sm.symbols("x y")
    printer = compatible(symbol_table={x: "values[0]"}, ftype="double")
    assert printer.doprint(Assignment(y, x**2)) == "const double y = values[0]*values[0];"
    assert printer.doprint(Assignment(sm.Symbol("out[0]"), x)) == "out[0] = values[0];"


def test_cpp_statements_and_dispatch():
    x, y = sm.symbols("x y")

    class Custom(CXX17CodePrinter):
        def _print_Symbol(self, symbol):
            return "v_" + super()._print_Symbol(symbol)

    printer = Custom()
    assert printer.doprint(CodeBlock(Comment("compute"), Assignment(y, 2*x))) == (
        "/* compute */;\nv_y = 2*v_x;")
    assert printer.doprint(x/y) == "v_x/v_y"
    assert printer.doprint(x+1, assign_to=y) == "v_y = v_x + 1;"
    assert printer.parenthesize(x+y, PRECEDENCE["Mul"]) == "(v_x + v_y)"


def test_indexed_references():
    x = sm.Symbol("x")
    a = IndexedBase("a", shape=(2, 3))
    assert isinstance(a[1, 2], Indexed)
    assert a[1, 2] == Indexed(a, 1, 2)
    assert CXX17CodePrinter().doprint(a[1, 2]) == "a[5]"
    assert CXX17CodePrinter().doprint(Assignment(IndexedBase("out")[0], x)) == "out[0] = x;"
    with pytest.raises(ValueError, match="Rank"):
        a[0]


def test_precedence_and_debug_printers():
    x, y = sm.symbols("x y")
    assert precedence(x+y) == PRECEDENCE["Add"]
    assert precedence(x*y) == PRECEDENCE["Mul"]
    assert precedence(x**2) == PRECEDENCE["Pow"]
    assert precedence(sm.Rational(-1, 2)) == PRECEDENCE["Add"]
    assert StrPrinter({"order": "none"}).doprint(x+y) == str(x+y)
    assert PrettyPrinter({"order": "none"}).doprint(x) == "x"


def test_trisplyne_generated_cpp_values(trisplyne_printers, tmp_path):
    """Compile real generated C++ and compare its values to SymPy's oracle."""
    compiler = shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        pytest.skip("C++ compiler not installed")
    import sympy
    _, compatible = trisplyne_printers
    sources = [
        "sqrt(x) + x**Rational(1,3) - y/x**2",
        "(x+y)**2/(x+2)**3 - UnevaluatedExpr(-y)*x**4",
        "(x+y)/x**2 + Rational(3,4)*y",
        "-UnevaluatedExpr(x+y) + UnevaluatedExpr(2*(x+y))/(x+1)",
    ]
    expressions = []
    oracles = []
    for source in sources:
        for backend in (sm, sympy):
            namespace = vars(backend).copy()
            x, y = backend.symbols("x y")
            namespace.update(x=x, y=y)
            expr = eval(source, namespace)
            if backend is sm:
                expressions.append(compatible().doprint(expr))
            else:
                oracles.append(sympy.lambdify((x, y), expr, "math"))
    cpp = tmp_path / "printer.cpp"
    executable = tmp_path / "printer"
    cpp.write_text("""
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
using std::pow;
using std::sqrt;
template <typename T>
std::array<T, 4> evaluate(T x, T y) {
    return {""" + ",\n".join(expressions) + """};
}
int main() {
    double x, y;
    std::cout << std::setprecision(17);
    while (std::cin >> x >> y) {
        for (double value : evaluate(x, y)) std::cout << value << ' ';
        std::cout << '\\n';
    }
}
""")
    subprocess.run([compiler, "-std=c++17", "-O2", str(cpp), "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    points = [(0.5, -1.0), (0.75, -0.25), (1.0, 0.0), (1.5, 0.5), (2.0, 1.5)]
    output = subprocess.run([str(executable)],
                            input="".join(f"{x} {y}\n" for x, y in points),
                            check=True, capture_output=True, text=True)
    rows = output.stdout.splitlines()
    assert len(rows) == len(points)
    for (x, y), row in zip(points, rows):
        assert list(map(float, row.split())) == pytest.approx(
            [oracle(x, y) for oracle in oracles], rel=1e-13, abs=1e-13)
