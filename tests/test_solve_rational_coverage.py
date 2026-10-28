"""A complete QQ-root certificate avoids solving the factors a second time.

The certificate never supplies or removes roots. Unsupported or incomplete
coverage retains the existing postprocessing, domain and option behavior.
"""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


x = sx.Symbol('x')


def foreign(value):
    if isinstance(value, dict):
        return {foreign(key): foreign(item) for key, item in value.items()}
    if isinstance(value, list):
        return [foreign(item) for item in value]
    if isinstance(value, tuple):
        return tuple(foreign(item) for item in value)
    if isinstance(value, set):
        return {foreign(item) for item in value}
    return sx.to_sympy(value)


def certificate(expression, values, target=x):
    candidates = [sx.sympify(value) for value in values]
    return sx._core._solve_complete_rational_roots(expression, target, candidates)


FLAGS = [{}, {'simplify': True}, {'simplify': False}, {'check': False},
         {'dict': True}, {'set': True}, {'dict': True, 'set': True},
         {'cubics': False}, {'quartics': False}, {'incomplete': False}]


@pytest.mark.parametrize('degree', [4, 7])
@pytest.mark.parametrize('expanded', [False, True])
@pytest.mark.parametrize('options', FLAGS)
def test_complete_public_solution_lists_and_options_match_reference(degree, expanded, options):
    X = sp.Symbol('x')
    expression = sx.prod(x-index for index in range(1, degree+1))
    reference = sp.prod(X-index for index in range(1, degree+1))
    if expanded:
        expression, reference = sx.expand(expression), sp.expand(reference)
    actual = foreign(sx.solve(expression, x, **options))
    expected = sp.solve(reference, X, **options)
    if options.get('dict') and options.get('set') and expanded:
        # Pre-port 345253 already orders these mappings increasingly, whereas
        # SymPy's combined flags leave expanded-polynomial roots descending.
        # Preserve the established native order and compare every assignment.
        assert actual == [{X: sp.Integer(index)} for index in range(1, degree+1)]
        assert sorted(actual, key=lambda row: row[X]) == sorted(expected, key=lambda row: row[X])
    else:
        assert actual == expected


@pytest.mark.parametrize('source,values', [
    ('(2*x-1)**3*(3*x+2)**2*(x-4)**4', ['-2/3', '1/2', '4']),
    ('-7*(x-1/2)*(x+2/3)*(x-4)/3', ['-2/3', '1/2', '4']),
    ('x**5*(x-1)**3', ['0', '1']),
    ('x**3-6*x**2+11*x-6', ['1', '2', '3']),
    ('x**3-3*x**2+3*x-1', ['1']),
    ('(x-1)**1000000', ['1']),
    ('x**1000000', ['0']),
    ('x**100000001-x**100000000', ['0', '1']),
])
def test_native_certificate_removes_every_multiplicity_without_factoring_again(source, values):
    expression = sx.sympify(source, locals={'x': x})
    assert certificate(expression, values) is True


def test_exact_coefficient_and_root_sizes_do_not_become_machine_precision():
    small = sx.Rational(1, 10**80)
    large = sx.Integer(10**50)
    expression = (x-small)**3*(x-large)**2
    assert certificate(expression, [small, large]) is True
    actual = sx.solve(expression, x)
    assert actual == [small, large]
    assert not any(root.has(sx.Float) for root in actual)


@pytest.mark.skipif(sys.platform != 'linux', reason='the allocation regression uses Linux address-space accounting')
def test_huge_sparse_nonzero_root_refuses_before_any_dense_allocation():
    # A regressed dense conversion needs several GB. Bound the child process
    # after imports so NumPy's existing mappings do not affect the allowance.
    # An allocation error or a timeout fails this test without exhausting CI.
    script = '''
import os, resource
import symox as s
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
pages = int(open('/proc/self/statm').read().split()[0])
ceiling = pages*os.sysconf('SC_PAGE_SIZE') + 128*1024*1024
resource.setrlimit(resource.RLIMIT_AS, (ceiling, ceiling))
x = s.Symbol('budget_x')
assert s._core._solve_complete_rational_roots(x**100000000-1, x, [s.Integer(1)]) is False
assert s._core._solve_complete_rational_roots((x+1)**100000000-1, x, [s.Integer(0)]) is False
assert s._core._solve_complete_rational_roots(x**4294967297-1, x, [s.Integer(1)]) is False
assert s._core._solve_complete_rational_roots(x**100000001-x**100000000, x, [s.Integer(0), s.Integer(1)]) is True
assert s._core._solve_complete_rational_roots((x-1)**100000000, x, [s.Integer(1)]) is True
'''
    result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('source,values', [
    ('(x-1)*(x-2)', ['1']),
    ('(x-1)*(x-2)', ['1', '2', '3']),
    ('(x-1)*(x**2-2)', ['1']),
    ('x**2-2', ['-sqrt(2)', 'sqrt(2)']),
    ('a*(x-1)', ['1']),
    ('f(1)*(x-1)', ['1']),
    ('sin(1)*(x-1)', ['1']),
    ('1.0*x-1', ['1']),
    ('x-1', ['1.0']),
    ('(x**2-1)/(x-1)', ['-1']),
    ('sin(x)', ['0']),
    ('exp(x)-1', ['0']),
    ('sqrt(x)-1', ['1']),
    ('x**3', ['1']),
    ('0', ['1']),
    ('1', ['1']),
    ('x-1', []),
])
def test_native_refuses_partial_extraneous_inexact_or_nonqq_coverage(source, values):
    expression = sx.sympify(source, locals={'x': x})
    assert certificate(expression, values) is False


def test_duplicates_and_dictionary_containers_preserve_the_old_normalization():
    expression = (x-1)*(x-2)
    assert certificate(expression, [1, 1, 2]) is True
    assert sx._core._solve_complete_rational_roots(expression, x, [{x: sx.Integer(1)}, {x: sx.Integer(2)}]) is True
    other = sx.Symbol('other')
    assert sx._core._solve_complete_rational_roots(expression, x, [{other: sx.Integer(1)}]) is False
    assert sx._core._solve_complete_rational_roots(expression, x, [{x: sx.Integer(1), other: sx.Integer(2)}]) is False


@pytest.mark.parametrize('assumptions', [{'positive': True}, {'negative': True}, {'nonzero': True}])
def test_assumption_filtered_root_lists_decline_completeness_and_keep_checks(assumptions):
    target = sx.Symbol('target', **assumptions)
    reference_target = sp.Symbol('target', **assumptions)
    expression = target*(target-1)*(target+2)
    reference = reference_target*(reference_target-1)*(reference_target+2)
    filtered = sx.solve(expression, target)
    assert sx._core._solve_complete_rational_roots(expression, target, filtered) is False
    assert foreign(filtered) == sp.solve(reference, reference_target)
    assert foreign(sx.solve(expression, target, check=False)) == sp.solve(reference, reference_target, check=False)


def test_imaginary_assumptions_with_no_admissible_rational_roots_do_not_certify():
    target = sx.Symbol('target', imaginary=True)
    reference_target = sp.Symbol('target', imaginary=True)
    expression = target*(target-1)*(target+2)
    actual = sx.solve(expression, target)
    assert actual == sp.solve(reference_target*(reference_target-1)*(reference_target+2), reference_target) == []
    assert sx._core._solve_complete_rational_roots(expression, target, actual) is False


@pytest.mark.parametrize('options', [{}, {'check': False}, {'dict': True}, {'set': True}])
def test_variable_denominators_keep_existing_pole_checks_and_result_shapes(options):
    X = sp.Symbol('x')
    expression = (x**2-1)/(x-1)
    assert foreign(sx.solve(expression, x, **options)) == sp.solve((X**2-1)/(X-1), X, **options)


@pytest.mark.parametrize('source', ['x**2-2', '(x-1)*(x**2-2)', 'a*(x-1)*(x-2)', 'sin(1)*(x-1)*(x-2)'])
def test_other_domains_retain_the_complete_existing_solution_lists(source):
    actual = sx.solve(sx.sympify(source), x)
    assert foreign(actual) == sp.solve(sp.sympify(source), sp.Symbol('x'))


def test_float_rationalization_options_and_declared_precision_are_preserved():
    coefficient = sx.Float('1.25', 80)
    equation = coefficient*x-2*coefficient
    assert certificate(equation, [2]) is False
    X = sp.Symbol('x')
    reference_coefficient = sp.Float('1.25', 80)
    for rational in (False, True, None):
        result = sx.solve(equation, x, rational=rational)
        expected = sp.solve(reference_coefficient*X-2*reference_coefficient, X, rational=rational)
        assert foreign(result) == expected
        if rational is False:
            assert all(root._prec >= coefficient._prec for root in result)


def test_equations_systems_and_expression_targets_keep_their_existing_routes():
    assert sx.solve(sx.Eq((x-1)*(x-2), 0), x) == [1, 2]
    assert sx.solve([(x-1)*(x-2)], x) == [(1,), (2,)]
    assert sx.solve(x**4-3*x**2+2, x**2) == [1, 2]


@pytest.mark.parametrize('domain', [sx.Reals, sx.Integers, sx.Complexes, sx.FiniteSet(1)])
def test_the_ignored_solve_domain_flag_does_not_restrict_the_qq_root_list(domain):
    assert sx.solve((x-1)*(x-2), x, domain=domain) == [1, 2]


def test_native_certificate_declines_expression_targets():
    assert sx._core._solve_complete_rational_roots(x**4-3*x**2+2, x**2, [sx.Integer(1), sx.Integer(2)]) is False


@pytest.mark.parametrize('target', [sx.Symbol('nc', commutative=False), sx.Dummy('dummy'), sx.Wild('wild')])
def test_native_certificate_declines_nonordinary_symbols(target):
    assert sx._core._solve_complete_rational_roots(target, target, [sx.Integer(0)]) is False


def test_adapted_indexed_scalar_identity_is_outside_the_ordinary_symbol_proof():
    indexed = sx.IndexedBase('coverage_array')[0]._symox_()
    assert indexed._scalar_identity is not None
    assert sx._core._solve_complete_rational_roots(indexed, indexed, [sx.Integer(0)]) is False


def test_host_symbol_subclasses_keep_the_existing_solve_linear_protocol(monkeypatch):
    class HookedSymbol(sx.Symbol):
        def _eval_subs(self, old, new):
            return None

    target = HookedSymbol('coverage_hosted')
    equation = (target-1)*(target-2)*(target-3)
    assert sx._core._solve_complete_rational_roots(equation, target, [sx.Integer(i) for i in (1, 2, 3)]) is False
    original = sx.solve_linear
    calls = []

    def traced(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)

    monkeypatch.setattr(sx, 'solve_linear', traced)
    assert sx.solve(equation, target) == [1, 2, 3]
    assert len(calls) >= 3


def test_root_objects_and_symbolic_incomplete_results_keep_the_previous_fallback():
    polynomial = x**5-x-1
    root = sx.CRootOf(polynomial, 0)
    assert certificate(polynomial, [root]) is False
    parameter = sx.Symbol('parameter')
    incomplete = x**5+parameter*x+1
    assert sx.solve(incomplete, x) == []
    with pytest.raises(NotImplementedError):
        sx.solve(incomplete, x, incomplete=False)


def test_complete_native_proof_avoids_the_repeated_high_order_analysis(monkeypatch):
    import symox._solve_relations as adapter

    def repeated_analysis(*args):
        raise AssertionError('a certified complete QQ list must not solve the factors again')

    monkeypatch.setattr(adapter, '_high_order_roots', repeated_analysis)
    assert sx.solve(sx.prod(x-index for index in range(1, 8)), x) == list(range(1, 8))


def test_failed_native_proof_keeps_the_previous_high_order_fallback(monkeypatch):
    import symox._solve_relations as adapter

    original = adapter._high_order_roots
    calls = []

    def traced(*args):
        calls.append(args)
        return original(*args)

    monkeypatch.setattr(adapter, '_high_order_roots', traced)
    assert sx.solve(x**2-2, x) == [-sx.sqrt(2), sx.sqrt(2)]
    assert calls


@pytest.mark.parametrize('options', [{'simplify': True}, {'simplify': False}, {'set': True}])
def test_explicit_simplify_and_set_options_keep_their_previous_early_routes(monkeypatch, options):
    def unnecessary_proof(*args):
        raise AssertionError('this option does not need the default high-order decision')

    monkeypatch.setattr(sx._core, '_solve_complete_rational_roots', unnecessary_proof)
    expression = (x-1)*(x-2)*(x-3)
    X = sp.Symbol('x')
    assert foreign(sx.solve(expression, x, **options)) == sp.solve((X-1)*(X-2)*(X-3), X, **options)
