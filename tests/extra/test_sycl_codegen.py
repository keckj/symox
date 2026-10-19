"""SYCL source export and execution through a real optional AdaptiveCpp SDK."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import symox as sm


def test_source_export_needs_no_optional_runtime():
    subprocess.run([sys.executable, "-c", "import sys; import symox as s; x=s.Symbol('x'); assert 'sycl::sin' in s.extra.syclcode(s.sin(x)); assert s.extra.adaptivecppcode(x)==s.extra.syclcode(x); assert not {'slangpy', 'torch', 'jax', 'sympy'}.intersection(sys.modules)"],
                   check=True, capture_output=True, text=True)


@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_native_printer_math_and_literals(dtype):
    x = sm.Symbol("x")
    source = sm.extra.syclcode(sm.sin(x)+sm.sqrt(x)+sm.Rational(1, 3), dtype=dtype)
    assert "sycl::sin" in source and "sycl::sqrt" in source
    assert "std::" not in source and "1/3" not in source
    assert "float('" not in source and "torch" not in source
    assert "numeric_limits" in sm.extra.syclcode(sm.oo, dtype=dtype)
    assert "replacement" in sm.extra.syclcode(x, symbol_names={x: "replacement"})
    assert "sycl::cos" in sm.extra.syclcode(sm.Function("custom")(x), user_functions={"custom": "sycl::cos"})
    with pytest.raises((ValueError, NotImplementedError)):
        sm.extra.syclcode(sm.I*x)
    with pytest.raises((ValueError, NotImplementedError)):
        sm.extra.syclcode(sm.Function("unknown")(x))


@pytest.mark.parametrize("name", ["class", "bad-name", "_sx_hidden", "__reserved", "_Upper", "a::b"])
def test_invalid_export_identifiers(name):
    x = sm.Symbol("x")
    with pytest.raises((ValueError, NotImplementedError)):
        sm.extra.syclcode_function(name, x, x)


def test_export_cse_types_identifiers_and_asynchronous_contract():
    x, y = sm.Symbol("lambda"), sm.Symbol("α")
    shared = sm.sin(x+y) + sm.exp(x*y)
    source = sm.extra.syclcode_function("exported", (x, y), (shared**2, shared**3, x > y))
    assert source == sm.extra.adaptivecppcode_function("exported", (x, y), (shared**2, shared**3, x > y))
    assert source.count("sycl::sin") == 1
    assert "bool _sx_o_2" in source
    assert "lambda" not in source and "α" not in source
    kernel = sm.extra.sycl_kernel("exported", (x, y), shared)
    for expected in ["sycl::event", "sycl::queue&", "depends_on", "parallel_for", "aspect::fp64", "single_task"]:
        assert expected in kernel
    assert ".wait(" not in kernel and ".wait_and_throw(" not in kernel
    assert "malloc" not in kernel
    assert kernel == sm.extra.adaptivecpp_kernel("exported", (x, y), shared)
    assert "aspect::fp64" not in sm.extra.sycl_kernel("float_kernel", x, x, dtype="float32")
    for args, expr in [((x, x), x), ((x,), x+y)]:
        with pytest.raises(ValueError):
            sm.extra.syclcode_function("f", args, expr)


def test_printer_nested_subclass_hooks():
    class Custom(sm.extra.SYCLCodePrinter):
        def _print_Symbol(self, expr):
            return "v_" + expr.name

        def _print_sin(self, expr):
            return "sycl::cos(" + self._print(expr.args[0]) + ")"

    source = Custom().doprint(sm.sin(sm.Symbol("x")) + sm.Symbol("x")**2)
    assert "v_x" in source and "sycl::cos" in source


def test_cse_boolean_and_piecewise_domains():
    x = sm.Symbol("x")
    condition = x > 0
    source = sm.extra.syclcode_function("condition", x, (condition, sm.And(condition, x < 1), condition))
    assert "bool _sx_t_" in source
    expr = sm.Piecewise((sm.log(x)+1/x, x > 0), (x*x, True))
    assert sm.extra.syclcode_function("guarded", x, expr, cse=True) == sm.extra.syclcode_function("guarded", x, expr, cse=False)


def _compiler():
    compiler = os.environ.get("SYMOX_ACPP") or shutil.which("acpp")
    if not compiler:
        candidate = Path(sys.executable).parent / "acpp"
        if candidate.is_file():
            compiler = str(candidate)
    if not compiler:
        if os.environ.get("SYMOX_ACPP_REQUIRE_RUNTIME") == "1":
            pytest.fail("Required AdaptiveCpp compiler is unavailable")
        pytest.skip("AdaptiveCpp SDK not installed; set SYMOX_ACPP to test compilation")
    return compiler


def test_compile_and_execute_scalar_and_asynchronous_kernels(tmp_path):
    compiler = _compiler()
    x, y = sm.symbols("x y")
    expressions = [sm.Piecewise((sm.log(x)+1/x, x > 0), (x*x, True)),
                   sm.sin(x)+sm.sqrt(1+y*y), x > 0, x**3, sm.Rational(1, 3)]
    source = sm.extra.sycl_kernel("evaluate", (x, y), expressions)
    source += sm.extra.syclcode_function("single_float", (x, y), sm.sin(x)+sm.sqrt(1+y*y), dtype="float32")
    a, b = sm.Symbol("same", positive=True), sm.Symbol("same", negative=True)
    source += sm.extra.syclcode_function("assumptions", (a, b), a+2*b)
    source += sm.extra.syclcode_function("empty", (), sm.Matrix(0, 0, []))
    source += r'''
#include <iostream>
#include <iomanip>
#include <cmath>
int main() {
    sycl::queue q{sycl::cpu_selector_v};
    constexpr std::size_t n=7;
    auto x=sycl::malloc_shared<double>(n,q);
    auto y=sycl::malloc_shared<double>(n,q);
    auto a=sycl::malloc_shared<double>(n,q);
    auto b=sycl::malloc_shared<double>(n,q);
    auto c=sycl::malloc_shared<bool>(n,q);
    auto d=sycl::malloc_shared<double>(n,q);
    auto e=sycl::malloc_shared<double>(n,q);
    auto flag=sycl::malloc_shared<int>(1,q);
    *flag=0;
    for(std::size_t i=0;i<n;++i){x[i]=double(i)-3.;y[i]=.25*double(i);}
    auto prepare=q.submit([&](sycl::handler& h){h.single_task([=](){*flag=41;});});
    auto event=evaluate(q,n,x,y,a,b,c,d,e,{prepare});
    event.wait_and_throw();
    if(*flag!=41)return 3;
    auto pending=q.submit([&](sycl::handler& h){h.single_task([=](){*flag=42;});});
    auto zero=evaluate(q,0,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,{pending});
    zero.wait_and_throw();
    if(*flag!=42)return 6;
    std::cout<<std::setprecision(17);
    for(std::size_t i=0;i<n;++i)
        std::cout<<a[i]<<' '<<b[i]<<' '<<c[i]<<' '<<d[i]<<' '<<e[i]<<'\n';
    if(std::abs(single_float(.2f,.3f)-(std::sin(.2f)+std::sqrt(1.09f)))>1.e-6f)return 4;
    if(assumptions(2.,-3.)!=-4.)return 5;
    (void)empty();
    sycl::free(x,q);sycl::free(y,q);sycl::free(a,q);sycl::free(b,q);
    sycl::free(c,q);sycl::free(d,q);sycl::free(e,q);sycl::free(flag,q);
}
'''
    source_file, program = tmp_path / "generated.cpp", tmp_path / "generated"
    source_file.write_text(source)
    target = os.environ.get("SYMOX_ACPP_TARGETS", "omp")
    subprocess.run([compiler, "--acpp-targets="+target, "-std=c++17", "-O2", str(source_file), "-o", str(program)],
                   check=True, capture_output=True, text=True, timeout=180)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    result = subprocess.run([str(program)], check=True, capture_output=True, text=True, env=env, timeout=60)
    actual = np.asarray([[float(value) for value in line.split()] for line in result.stdout.splitlines()])
    tx, ty = np.arange(7.)-3, .25*np.arange(7.)
    safe = np.where(tx > 0, tx, 1.)
    expected = np.stack([np.where(tx > 0, np.log(safe)+1/safe, tx*tx),
                         np.sin(tx)+np.sqrt(1+ty*ty), tx > 0, tx**3, np.full(7, 1/3)], axis=-1)
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)
