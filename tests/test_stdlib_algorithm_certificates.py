"""The source-equivalence evidence must become invalid after a code change."""
import sys
import pytest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from sympy_stdlib_certificates import certificates,verify_dependencies


def test_container_algorithm_source_certificates():
    report=certificates()
    assert len(report['records'])==5
    assert all(r['status']=='source_equivalent' for r in report['records'])


def test_extended_certificates_close_their_actual_standard_library_dependencies():
    report=certificates(extended=True)
    records={r['id']:r for r in report['records']}
    assert len(records)==21
    assert records['sympy.utilities.iterables.bracelets']['dependency_evidence']
    assert records['sympy.prod']['dependency_evidence']
    assert report['native_source_sha256']


def test_identical_wrapper_does_not_certify_a_changed_helper(monkeypatch):
    from sympy.utilities import iterables as reference
    from symox.utilities import iterables as native
    def broken(n,k,free=False):return iter(())
    monkeypatch.setitem(native.bracelets.__globals__,'necklaces',broken)
    with pytest.raises(ValueError,match='source differs'):
        verify_dependencies(reference.bracelets,native.bracelets)


def test_identical_body_does_not_certify_a_rebound_builtin(monkeypatch):
    from sympy.utilities import iterables as reference
    from symox.utilities import iterables as native
    monkeypatch.setitem(native.necklaces.__globals__,'range',lambda *args:())
    with pytest.raises(ValueError,match='unverified dependency range'):
        verify_dependencies(reference.necklaces,native.necklaces)
