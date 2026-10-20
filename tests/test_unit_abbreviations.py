"""StrPrinter's ``abbrev`` setting prints quantities like SymPy."""
import importlib

import pytest
import sympy as sp
import symox as sx

CASES = [
    lambda m, u, x: u.volt,
    lambda m, u, x: 2*u.volt/u.hertz,
    lambda m, u, x: x*u.volt,
    lambda m, u, x: u.radian/u.second,
    lambda m, u, x: u.volt*u.ampere,
    lambda m, u, x: 1/u.second,
    lambda m, u, x: u.gram*u.kilogram,
    lambda m, u, x: u.volt + x,
    lambda m, u, x: [u.volt, x],
    lambda m, u, x: m.sqrt(u.watt),
]


def _env(module):
    return module, importlib.import_module(module.__name__ + '.physics.units'), module.Symbol('x')


def _printer(module):
    return importlib.import_module(module.__name__ + '.printing.str').StrPrinter


@pytest.fixture
def restore_defaults():
    saved = {module: dict(_printer(module)._default_settings) for module in (sx, sp)}
    yield
    for module, settings in saved.items():
        _printer(module)._default_settings.clear()
        _printer(module)._default_settings.update(settings)


@pytest.mark.parametrize('index', range(len(CASES)))
def test_explicit_abbrev_setting_matches_sympy(index):
    def run(module):
        value = CASES[index](*_env(module))
        return module.sstr(value, abbrev=True), module.sstr(value, abbrev=False), str(value)
    assert run(sx) == run(sp)


@pytest.mark.parametrize('index', range(len(CASES)))
def test_default_abbrev_setting_matches_sympy(index, restore_defaults):
    # Lcapy sets StrPrinter._default_settings['abbrev'] to print V for volt.
    def run(module):
        value = CASES[index](*_env(module))
        _printer(module)._default_settings['abbrev'] = True
        return str(value), repr(value), module.sstr(value, abbrev=False)
    assert run(sx) == run(sp)


def test_str_printer_defaults_are_its_own():
    assert sx.printing.str.StrPrinter._default_settings == sp.printing.str.StrPrinter._default_settings
    assert 'abbrev' not in sx.printing.printer.Printer._default_settings


QUANTITY_QUERIES = [
    lambda m, u, x: m.Mul(7, u.ampere, evaluate=False),
    lambda m, u, x: x*u.volt,
    lambda m, u, x: m.sqrt(u.volt),
    lambda m, u, x: x,
    lambda m, u, x: u.volt/u.ampere,
    lambda m, u, x: 3*u.volt + x,
    lambda m, u, x: u.speed_of_light*x,
]


def _constant(module):
    return importlib.import_module(module.__name__ + '.physics.units.quantities').PhysicalConstant


@pytest.mark.parametrize('case', QUANTITY_QUERIES)
def test_quantity_class_queries_match_sympy(case):
    def query(module):
        m, u, x = _env(module)
        expr = case(m, u, x)
        return (expr.has(u.Quantity), sorted(map(str, expr.atoms(u.Quantity))),
                sorted(map(str, expr.atoms(_constant(module)))),
                sorted(isinstance(a, u.Quantity) for a in expr.args))
    assert query(sx) == query(sp)


def test_quantity_isinstance_matches_sympy():
    def query(module):
        m, u, x = _env(module)
        constant = _constant(module)
        return (isinstance(u.volt, u.Quantity), isinstance(u.volt, constant),
                isinstance(u.speed_of_light, constant), isinstance(x, u.Quantity),
                isinstance(3, u.Quantity), isinstance([1], u.Quantity),
                isinstance((x*u.volt).args[-1], u.Quantity) or isinstance((x*u.volt).args[0], u.Quantity))
    assert query(sx) == query(sp)
