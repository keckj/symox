"""Coverage accounting must retain omitted modules and unmeasurable entries."""
import json
import runpy
from types import SimpleNamespace
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
coverage=runpy.run_path(str(ROOT/'tools/sympy_coverage.py'))


def test_documented_families_are_in_inventory():
    inventory=coverage['load_inventory']('1.14.0')
    required={'sympy.physics.mechanics','sympy.physics.quantum.gate',
              'sympy.physics.units.quantities','sympy.stats',
              'sympy.crypto.crypto','sympy.holonomic','sympy.diffgeom',
              'sympy.liealgebras','sympy.categories'}
    assert required <= set(inventory['modules'])
    assert 'evaluate_pauli_product' in inventory['documented_exports']['sympy.physics.paulialgebra']
    with pytest.raises(ValueError,match='installed version'):
        coverage['load_inventory']('another-version')


def test_missing_native_module_is_zero():
    reference=SimpleNamespace(__name__='sympy.example',__all__=['one','two'])
    result=coverage['surface'](reference,None,reference.__all__)
    assert result['present']==0 and result['total']==2
    assert result['missing']==['one','two']
    assert coverage['percentage'](result)=='0.00%'
    assert coverage['percentage'](coverage['surface'](reference,None,[]))=='N/A'


def test_exports_preserve_documented_names_and_handle_missing_module_metadata():
    func=lambda:None
    func.__module__=None
    module=SimpleNamespace(__name__='example',helper=func,documented=42)
    assert coverage['exports'](module,['documented','absent'])==['documented']
    module.__all__=['helper']
    assert coverage['exports'](module,['documented'])==['documented','helper']


def test_missing_parent_and_internal_import_error_are_distinct(monkeypatch):
    import importlib
    def missing(name):
        raise ModuleNotFoundError('missing parent',name='symox.physics')
    monkeypatch.setattr(importlib,'import_module',missing)
    assert coverage['native_module']('symox.physics.units')[1:] == ('missing_module',None)
    def dependency(name):
        raise ModuleNotFoundError('missing dependency',name='unrelated_dependency')
    monkeypatch.setattr(importlib,'import_module',dependency)
    actual,status,error=coverage['native_module']('symox.physics.units')
    assert actual is None and status=='native_import_error'
    assert 'missing dependency' in error
