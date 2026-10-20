"""A fast wrong, incomplete or unevaluated answer must never win a benchmark."""
from pathlib import Path
import sys

import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benches"))
from performance_cases import Case,catalog
from performance_validation import encode,verify,numeric_equal,equivalent
from performance_report import summaries


def test_catalog_is_deterministic_and_has_unique_names():
    assert catalog() == catalog()
    assert len({case.id for case in catalog()}) == len(catalog())


def test_machine_precision_roots_do_not_pass_thirty_digit_validation():
    actual=[sp.Float(2.0**0.5)]
    reference=[sp.sqrt(2).evalf(50)]
    case=Case("test","polynomials","nroots",1,"numeric_roots",30)
    assert verify(case,encode(actual),encode(reference))["status"] == "numerical_difference"
    assert not numeric_equal(sp.S.Zero,sp.Rational(1,10**1000),50)


def test_missing_root_and_held_work_are_excluded():
    case=Case("test","polynomials","nroots",2,"numeric_roots",30)
    assert verify(case,encode([1]),encode([1,2]))["status"] == "root_count_difference"
    x=sp.Symbol("x");held=sp.Integral(sp.sin(x),x)
    case=Case("test","integration","test")
    assert verify(case,encode(held),encode(held))["status"] == "unevaluated"


def test_exactness_and_factorization_quality_are_separate_requirements():
    case=Case("test","test","test")
    assert verify(case,encode(sp.Float(1)),encode(sp.Integer(1)))["status"] == "exactness_difference"
    x=sp.Symbol("x")
    case=Case("test","polynomials","factor_cyclotomic",2)
    assert verify(case,encode(x*x-1),encode((x-1)*(x+1)))["status"] == "factorization_quality_difference"


def test_bound_dummy_names_do_not_mask_lost_solution_branches():
    a,b=sp.Dummy("n",integer=True),sp.Dummy("n",integer=True)
    left=sp.ImageSet(sp.Lambda(a,2*sp.pi*a),sp.S.Integers)
    right=sp.ImageSet(sp.Lambda(b,2*sp.pi*b),sp.S.Integers)
    assert equivalent(left,right)
    assert not equivalent(sp.FiniteSet(0),left)
    z,w=sp.Dummy("free"),sp.Dummy("free")
    assert not equivalent(sp.FiniteSet(z),sp.FiniteSet(w))


def test_numeric_transport_preserves_tiny_values_and_python_integer_results():
    import mpmath
    mp=mpmath.mp.clone();mp.dps=100
    value=mp.mpf("1e-1000")
    assert numeric_equal(sp.sympify(encode(value)["value"]),sp.Rational(1,10**1000),90)
    assert encode(7)==encode(sp.Integer(7))


def test_unverified_or_unstable_results_have_no_speedup():
    trial=dict(backends={name:dict(first_seconds=1,warm_seconds=1) for name in ("symox","sympy")},validation_key="v")
    report=dict(environment_unchanged=True,results={"case":dict(trials=[trial],validations={"v":dict(status="numerical_difference")})})
    assert summaries(report)[0]["comparable"] is False
    report["results"]["case"]["validations"]["v"]["status"]="verified"
    report["environment_unchanged"]=False
    assert summaries(report)[0]["comparable"] is False
