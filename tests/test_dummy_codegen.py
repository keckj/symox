"""Dummy identity survives lowering and extensible printer dispatch."""
import numpy as np
import pytest
import symox as s


@pytest.mark.parametrize('name', ['ccode', 'cxxcode', 'cudacode', 'torchcode', 'jaxcode', 'slangcode', 'syclcode'])
def test_standalone_codegen_uses_distinct_dummy_identifiers(name):
    first = s.Dummy('x', dummy_index=41)
    second = s.Dummy('x', dummy_index=42)
    result = getattr(s if name in ('ccode', 'cxxcode') else s.extra, name)(first + second**2)
    assert 'x_41' in result and 'x_42' in result


@pytest.mark.parametrize('name', ['torchcode', 'jaxcode', 'slangcode', 'syclcode'])
def test_exact_symbol_name_maps_distinguish_homonymous_dummies(name):
    first, second = s.Dummy('x'), s.Dummy('x')
    result = getattr(s.extra, name)(first + second, symbol_names={first: 'left', second: 'right'})
    assert 'left' in result and 'right' in result


@pytest.mark.parametrize('printer_name', ['PyTorchPrinter', 'JaxPrinter', 'SlangCodePrinter', 'SYCLCodePrinter'])
def test_symbol_subclass_hooks_reach_dummy_nodes_inside_expressions(printer_name):
    base = getattr(s.extra, printer_name)

    class Custom(base):
        def _print_Symbol(self, expr):
            return 'custom_' + str(expr.dummy_index)

    first, second = s.Dummy('x'), s.Dummy('x')
    result = Custom().doprint(first + second**2)
    assert 'custom_' + str(first.dummy_index) in result
    assert 'custom_' + str(second.dummy_index) in result


def test_native_numpy_evaluator_keeps_homonymous_dummy_inputs_distinct():
    first, second = s.Dummy('x'), s.Dummy('x')
    function = s.lambdify((first, second), first + 2*second, 'numpy')
    left, right = np.array([1., 2.]), np.array([3., 4.])
    np.testing.assert_allclose(function(left, right), left+2*right)


@pytest.mark.parametrize('name', ['torchcode_function', 'jaxcode_function', 'slangcode_function', 'syclcode_function'])
def test_function_export_maps_each_dummy_to_its_own_input(name):
    first, second = s.Dummy('x'), s.Dummy('x')
    source = getattr(s.extra, name)('evaluate', (first, second), first+2*second)
    assert '_sx_v_0' in source and '_sx_v_1' in source
    assert source
