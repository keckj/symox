"""Check final receipts and consolidate their source/binary identities."""
import hashlib
import json

import validate

raw = validate.RAW
names = [
    "focused", "suite-run3", "cargo-test", "clippy", "ode", "ode-vs-wave45",
    "tsolve", "tsolve-compare", "unrad", "unrad-compare", "matcher",
    "matcher-report", "algfactor", "algfactor-compare", "algpoly",
    "algpoly-compare", "ratint2", "ratint2-compare", "polytools-probe",
    "pde", "pde-report", "liegroup", "liegroup-report", "has", "has-report",
    "fps-formal", "fps-extra", "fps-formal-compare", "fps-extra-compare",
    "limits2", "limits2-compare", "lie-boundary-baseline", "lie-boundary-final",
]
current = validate.identity()
expected = json.loads((raw / "build-input.json").read_text())
assert current["source_sha256"] == expected["source_sha256"]
commands = {}
for name in names:
    record = json.loads((raw / (name + ".command.json")).read_text())
    assert record["returncode"] == 0, (name, record)
    if name != "lie-boundary-baseline":
        assert record["extension_sha256"] == current["so_sha256"], name
    commands[name] = record
for phase in ["focused", "suite", "rust", "corpora", "recent", "boundaries"]:
    record = json.loads((raw / (phase + "-identity.json")).read_text())
    assert record["source_sha256"] == current["source_sha256"], phase
    assert record["so_sha256"] == current["so_sha256"], phase
suite = json.loads((raw / "suite-run3-summary.json").read_text())
assert suite["coverage"]["collection_matches"] and suite["coverage"]["execution_matches"]
assert not suite["failures"] and not suite["run_errors"]
before = {r["id"]: r for r in json.loads((raw / "lie-boundary-baseline.json").read_text())}
after = {r["id"]: r for r in json.loads((raw / "lie-boundary-final.json").read_text())}
assert before.keys() == after.keys() == {"dsolve:lg01", "dsolve:fsum"}
boundary = {}
for key in before:
    a, b = before[key]["result"], after[key]["result"]
    assert (a["ok"], a.get("error"), a.get("value")) == (b["ok"], b.get("error"), b.get("value"))
    boundary[key] = {"baseline": a, "final": b}
result = {
    "head": current["head"], "source_sha256": current["source_sha256"],
    "extension_sha256": current["so_sha256"], "suite_counts": suite["counts"],
    "suite_coverage": suite["coverage"], "commands": commands,
    "paired_lie_boundary_results": boundary,
    "artifact_sha256": {
        str(p.relative_to(raw)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(raw.rglob("*")) if p.is_file() and p.name != "final-validation.json"
    },
}
(raw / "final-validation.json").write_text(json.dumps(result, indent=2) + "\n")
print("Final receipts consistent:", len(commands), "successful commands;", suite["counts"])
