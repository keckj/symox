"""Replay selected Lie cases using the archived isolation/deadline protocol."""
import importlib.util
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location(
    "liegroup_audit", root / "audits/performance/2026-11-06-liegroup/scripts/corpus.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
selected = set(sys.argv[2:])
cases = [case for case in audit.cases() if case["id"] in selected]
assert {case["id"] for case in cases} == selected
audit.cases = lambda: cases
Path(sys.argv[1]).write_text(json.dumps(audit.replay("symox", 120), indent=2) + "\n")
