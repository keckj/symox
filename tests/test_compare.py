"""Run every differential test module (tests/compare_*.py) under pytest."""

import importlib.util
import pathlib

import pytest

pytest.importorskip("sympy")

HERE = pathlib.Path(__file__).parent
MODULES = sorted(HERE.glob("compare_*.py"))


@pytest.mark.parametrize("path", MODULES, ids=[p.stem for p in MODULES])
def test_differential(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    bad = mod.run(verbose=False)
    assert not bad, "\n".join(f"{c}: symox={g} sympy={w}" for c, g, w in bad)
