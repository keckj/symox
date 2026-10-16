"""Regression checks for process isolation, reproducibility and serialization."""

import importlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys
import time

import pytest

import symox as sm

FUZZ_DIR = Path(__file__).resolve().parent / "fuzz"
sys.path.insert(0, str(FUZZ_DIR))
fuzz = importlib.import_module("fuzz_core")
triage = importlib.import_module("triage")


def isolated_worker(conn, state):
    while True:
        command = conn.recv()
        if command == "die":
            os._exit(19)
        if command == "hang":
            state.value = b"symox:test-hang"
            time.sleep(30)
        else:
            conn.send(("ok", command))


def test_worker_recovers_from_death_and_timeout():
    worker = fuzz.Worker(mp.get_context("spawn"), target=isolated_worker)
    try:
        assert worker.call("ready", 10) == ("ok", "ready")
        assert worker.call("die", 3)[0] == "died"
        assert worker.call("after death", 10) == ("ok", "after death")
        assert worker.call("hang", 0.1) == ("timeout", "symox:test-hang")
        assert worker.call("after timeout", 10) == ("ok", "after timeout")
    finally:
        worker.kill()
    assert not worker.proc.is_alive()
    assert worker.conn.closed


def test_integer_generator_never_nests_number_tuples():
    generator = fuzz.Gen(random.Random(42))
    atoms = [generator.atom() for _ in range(10000)]
    integers = [atom[1] for atom in atoms if atom[0] == "N"]
    assert any(abs(value) > 2**70 for value in integers)
    assert all(isinstance(value, int) for value in integers)


def test_infinity_comparisons_preserve_direction():
    fuzz.init_libs()
    assert fuzz.to_complex(fuzz._sympy.oo) == complex(math.inf)
    assert fuzz.to_complex(-fuzz._sympy.oo) == complex(-math.inf)
    assert not fuzz.close(complex(math.inf), complex(-math.inf), 1e-9)
    assert not fuzz.close(complex(math.inf), complex(0, math.inf), 1e-9)


@pytest.mark.parametrize("value", [1/3, -1/3.3, math.nextafter(1.0, 2.0), math.exp(0.5), 5e-324, 1e300])
def test_float_pickle_and_srepr_preserve_exact_binary_value(value):
    number = sm.Float(value)
    x = sm.Symbol("x", positive=True)
    for expression in (number, number*x + sm.Float(0.1)):
        restored = pickle.loads(pickle.dumps(expression))
        assert restored == expression
        assert hash(restored) == hash(expression)
        assert sm.sympify(sm.srepr(expression)) == expression
    assert float(pickle.loads(pickle.dumps(number))).hex() == value.hex()


@pytest.mark.parametrize("name", ["Eq", "Ne", "Lt", "Le", "Gt", "Ge"])
def test_relational_srepr_roundtrip(name):
    x, y = sm.symbols("x y")
    expression = getattr(sm, name)(x, y)
    assert sm.sympify(sm.srepr(expression)) == expression


def test_pickle_preserves_native_coefficients_and_nested_radicals():
    y, k = sm.symbols("y k")
    for expression in (sm.Float(1.25)*(y+k+sm.tanh(1)), sm.sqrt(sm.sqrt(12))):
        restored = pickle.loads(pickle.dumps(expression))
        assert restored == expression
        assert hash(restored) == hash(expression)


def cli(*args, timeout=20):
    return subprocess.run([sys.executable, str(FUZZ_DIR / "fuzz_core.py"), *map(str, args)],
                          text=True, capture_output=True, timeout=timeout)


@pytest.mark.parametrize("args", [
    ["--jobs", "0"], ["--timeout", "nan"], ["--time", "-1"],
    ["--seeds", "bad"], ["--mode", "missing"], ["--shrink-budget", "-1"],
    ["--replay", "core"], ["--replay", "core", "0", "missing"],
])
def test_invalid_cli_arguments_fail_before_spawning(args):
    assert cli(*args).returncode == 2


def test_wall_clock_limit_cancels_workers_and_shrinking(tmp_path):
    result = cli("--time", "0.1", "--timeout", "60", "--jobs", "3", "--out", tmp_path, timeout=12)
    assert result.returncode == 0, result.stderr
    assert json.loads((tmp_path / "stats.json").read_text())["interrupted"] is True


def test_finite_run_and_deduplication(tmp_path):
    for _ in range(2):
        result = cli("--seeds", "0:3", "--jobs", "3", "--timeout", "5", "--shrink-budget", "0", "--out", tmp_path)
        assert result.returncode == 0, result.stderr
        stats = json.loads((tmp_path / "stats.json").read_text())
        assert stats["cases"] == 3
        assert stats["deaths"] == stats["timeouts"] == 0
    assert stats["new"] == 0
    records = [json.loads(line) for line in (tmp_path / "findings.jsonl").read_text().splitlines()]
    assert len({record["sig"] for record in records}) == len(records)


def test_fail_on_is_explicit_and_opt_in(monkeypatch):
    monkeypatch.setattr(fuzz, "drive", lambda args: {"by_class": {"wrong": 1}})
    assert fuzz.main(["--seeds", "0:0"]) == 0
    assert fuzz.main(["--seeds", "0:0", "--fail-on", "wrong,crash"]) == 1


def test_triage_groups_checks_and_filters_classes():
    records = [dict(cls="wrong", min_src="x", check=check, min_msg="different value")
               for check in ("diff", "expand")]
    records.append(dict(cls="form", min_src="y", check="build", min_msg="different tree"))
    lines = list(triage.summarize(records, {"wrong"}))
    assert lines == ["[wrong] x   <diff,expand> x2", "      different value"]
