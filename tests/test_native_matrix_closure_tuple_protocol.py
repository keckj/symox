"""Function closure tuples must be inspected without invoking subclass hooks."""
import importlib
import types

import pytest
import symox as sx

from symox.polys.matrices import dense
from symox.polys.matrices.domainmatrix import DomainMatrix


@pytest.mark.parametrize("route", ["lu", "solve", "hnf"])
@pytest.mark.parametrize("raises", [False, True])
def test_closure_tuple_subclass_keeps_complete_math_without_speculative_hooks(
    monkeypatch, route, raises
):
    observed = []

    class Closure(tuple):
        def __len__(self):
            observed.append("len")
            if raises:
                raise RuntimeError("foreign closure length")
            return tuple.__len__(self)

        def __iter__(self):
            observed.append("iter")
            if raises:
                raise RuntimeError("foreign closure iteration")
            return tuple.__iter__(self)

    name = "__floordiv__" if route == "hnf" else "__truediv__"
    original = sx._core.DomainElement.__dict__[name]
    clone = types.FunctionType(
        original.__code__, original.__globals__, original.__name__,
        original.__defaults__, Closure(original.__closure__),
    )
    assert type(clone.__closure__) is Closure
    monkeypatch.setattr(sx._core.DomainElement, name, clone)
    kernel_name = {"lu": "_ddm_ilu", "solve": "_ddm_ilu_solve", "hnf": "_ddm_hnf"}[route]
    kernel = getattr(sx._core, kernel_name)
    admissions = []

    def record(*args):
        result = kernel(*args)
        admissions.append(result)
        return result

    monkeypatch.setattr(sx._core, kernel_name, record)
    normalforms = importlib.import_module("symox.polys.matrices.normalforms")

    def qq(rows):
        return [[sx.QQ(value) for value in row] for row in rows]

    def run():
        if route == "lu":
            values = qq([[2, 1], [4, 3]])
            result = dense.ddm_ilu(values)
            assert result == [] and values == qq([[2, 1], [2, 1]])
            return result, [[str(value) for value in row] for row in values]
        if route == "solve":
            values = [[None], [None]]
            result = dense.ddm_ilu_solve(
                values, qq([[1, 0], [2, 1]]), qq([[2, 1], [0, 1]]), [], qq([[4], [10]])
            )
            assert result is None and values == qq([[1], [2]])
            return result, [[str(value) for value in row] for row in values]
        source = DomainMatrix.from_list([[12, 6, 4], [3, 9, 6], [2, 16, 14]], sx.ZZ)
        result = normalforms.hermite_normal_form(source)
        values = [[int(value) for value in row] for row in result.to_list()]
        assert values == [[10, 0, 2], [0, 15, 3], [0, 0, 2]]
        assert [[int(value) for value in row] for row in source.to_list()] == [
            [12, 6, 4], [3, 9, 6], [2, 16, 14]
        ]
        return str(result.domain), values

    observed.clear()
    actual = run()
    assert observed == [] and admissions == []
    if route == "lu":
        monkeypatch.setattr(dense, "_native_ilu", lambda *args: None)
    elif route == "solve":
        monkeypatch.setattr(dense, "_native_ilu_solve", lambda *args: False)
    else:
        monkeypatch.setattr(normalforms, "_native_hnf", lambda *args: None)
    expected = run()
    assert expected == actual and observed == [] and admissions == []
