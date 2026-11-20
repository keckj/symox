"""Retain complete validation receipts and check their source/binary identity."""
from collections import Counter
import hashlib
import json
import shutil
import tarfile
import xml.etree.ElementTree as ET

import validate


def skipped_cases(xml):
    """Multiset of JUnit skip/xfail names and reasons, retaining duplicates."""
    cases = Counter()
    for case in xml.iter("testcase"):
        skipped = case.find("skipped")
        if skipped is None:
            continue
        reason = dict(skipped.attrib)
        if skipped.text and skipped.text.strip():
            reason["detail"] = skipped.text.strip()
        key = (case.get("classname", ""), case.get("name", ""),
               json.dumps(reason, sort_keys=True))
        cases[key] += 1
    return cases


def skip_counts(cases):
    counts = Counter()
    for (_, _, reason), count in cases.items():
        kind = "xfailed" if json.loads(reason).get("type") == "pytest.xfail" else "skipped"
        counts[kind] += count
    return {kind: counts[kind] for kind in ("skipped", "xfailed")}


def skip_rows(cases):
    return [{"test": [classname, name], "reason": json.loads(reason)}
            for classname, name, reason in sorted(cases.elements())]


raw = validate.RAW
current = validate.identity()
build = json.loads((raw / "build-input.json").read_text())
assert current["source_sha256"] == build["source_sha256"]
build_output = json.loads((raw / "build-output.json").read_text())
assert build_output["returncode"] == 0, build_output
assert build_output["source_unchanged"] is True, build_output
assert build_output["extension_sha256"] == current["so_sha256"], "Tested extension differs from the successful build output"
names = [
    "focused", "matcher47-probe", "polytools47-probe", "series47-boundary", "series47-cdir",
    "suite-final", "cargo-test", "clippy", "ode", "ode-vs-wave46",
    "tsolve", "tsolve-compare", "unrad", "unrad-compare", "matcher-sympy",
    "matcher", "matcher-report", "algfactor", "algfactor-compare", "algpoly",
    "algpoly-compare", "ratint2", "ratint2-compare", "pde", "pde-report",
    "liegroup", "liegroup-report", "has", "has-report", "fps-formal",
    "fps-extra", "fps-formal-compare", "fps-extra-compare", "limits2",
    "limits2-compare", "ode-paired-boundaries",
]
boundary = json.loads((raw / "lie-boundary-environment.json").read_text())
if boundary["cases"]:
    names.extend(["lie-boundary-baseline", "lie-boundary-final"])
    for label in ("baseline", "final"):
        rows = json.loads((raw / ("lie-boundary-" + label + ".json")).read_text())
        assert len(rows) == len(boundary["cases"])
        assert {row["id"] for row in rows} == set(boundary["cases"])
commands = {}
for name in names:
    record = json.loads((raw / (name + ".command.json")).read_text())
    assert record["returncode"] == 0, (name, record)
    expected_so = (json.loads((raw / "baseline.json").read_text()) if name == "lie-boundary-baseline" else current)
    expected_sha = expected_so.get("so_sha256", expected_so.get("extension_sha256"))
    assert record["extension_sha256"] == expected_sha, name
    commands[name] = record
for phase in ["focused", "lots", "suite", "rust", "corpora", "recent", "boundaries", "ode-boundaries"]:
    record = json.loads((raw / (phase + "-identity.json")).read_text())
    assert record["source_sha256"] == current["source_sha256"], phase
    assert record["so_sha256"] == current["so_sha256"], phase

ode_boundary = json.loads((raw / "ode-paired-boundaries/summary.json").read_text())
assert ode_boundary["complete"] and ode_boundary["rounds"] == 1
assert len(ode_boundary["keys"]) == 5 and len(ode_boundary["pairs"]) == 5
assert ode_boundary["extension_sha256"] == ode_boundary["extension_sha256_after"]
assert ode_boundary["extension_sha256"]["final"] == current["so_sha256"]
assert ode_boundary["extension_sha256"]["baseline"] == json.loads((raw / "baseline.json").read_text())["extension_sha256"]
for pair in ode_boundary["pairs"]:
    for label, call in pair["calls"].items():
        observation = call["observation"]
        assert call["returncode"] == 0 and observation["identity_matches"]
        assert observation["expected_extension_sha256"] == ode_boundary["extension_sha256"][label]

corpora_delta = json.loads((raw / "corpora-deltas.json").read_text())
assert corpora_delta["complete"]
assert set(corpora_delta["families"]) == {
    "tsolve", "unrad", "matcher", "algfactor", "algpoly", "ratint2",
    "pde", "liegroup", "fps-formal", "fps-extra", "ode",
}
limits_delta = json.loads((raw / "limits2-deltas.json").read_text())
assert limits_delta["status"] == "complete" and limits_delta["summary"]["cases"] == 1293
assert limits_delta["summary"]["source_line_deltas"] == 0
representation = json.loads((raw / "representation-review.json").read_text())
for family in ("pde", "liegroup"):
    rows = representation["families"][family]
    assert {row["key"] for row in rows} == {
        row["key"] for row in corpora_delta["families"][family]["representation_only_deltas"]
    }
    assert all(row["only_dummy_indices"] for row in rows)

# The lots phase also has post-command checks: the probe commands themselves
# can exit zero while reporting mismatches. Its pre-phase identity alone is
# not evidence that these post-command checks succeeded.
matcher = json.loads((raw / "matcher47-after.json").read_text())
polytools = json.loads((raw / "polytools47-after.json").read_text())
assert matcher["rows"] and all(row["public_agrees"] for row in matcher["rows"])
for row in matcher["rows"]:
    reference, native = row["sympy"], row["native"]
    if reference["ok"]:
        assert native.get("handled") and native.get("value") == reference["value"], row
    else:
        assert (native.get("type"), native.get("message")) == (reference["type"], reference["message"]), row
assert polytools and all(row["reference"] == row["actual"] for row in polytools)
cdir = json.loads((raw / "series47-cdir-after.json").read_text())
current_cases = cdir["current"]["cases"]
assert len(current_cases) == 8
assert current_cases.keys() == cdir["baseline"]["cases"].keys() == cdir["sympy"]["cases"].keys()
for name, row in current_cases.items():
    expected_engine = "sympy" if name in {"irrational_log_default", "irrational_log_positive"} else "baseline"
    expected = cdir[expected_engine]["cases"][name]["public_nseries"]
    assert row["public_nseries"] == expected, (name, row["public_nseries"], expected)

with tarfile.open(raw / "wave46-source.tar.gz") as archive:
    baseline_files = {member.name: hashlib.sha256(archive.extractfile(member).read()).hexdigest()
                      for member in archive if member.isfile()}
source_delta = {"added": [], "changed": [], "removed": []}
for path in sorted(current["files"].keys() | baseline_files.keys()):
    if path not in baseline_files:
        source_delta["added"].append(path)
    elif path not in current["files"]:
        source_delta["removed"].append(path)
    elif current["files"][path] != baseline_files[path]:
        source_delta["changed"].append(path)
assert not source_delta["removed"], source_delta
(raw / "source-delta.json").write_text(json.dumps(source_delta, indent=2)+"\n")
with tarfile.open(raw / "wave47-source.tar.gz", "w:gz") as archive:
    for path in current["files"]:
        archive.add(validate.ROOT / path, arcname=path, recursive=False)

evidence = validate.EVIDENCE / "suite-final"
suite = json.loads((evidence / "summary.json").read_text())
assert suite["coverage"]["collection_matches"] and suite["coverage"]["execution_matches"]
assert not suite["failures"] and not suite["run_errors"]
assert all(shard["returncode"] == 0 for shard in suite["shards"])
assert suite["identity"]["so_sha256"] == current["so_sha256"]
shutil.copyfile(evidence / "summary.json", raw / "suite-final-summary.json")
with tarfile.open(raw.parent.parent / "2026-11-09-wave46/raw/suite-receipts.tar.gz") as old:
    before = json.loads(old.extractfile("suite-run3/collection/execution.json").read())
    old_summary = json.loads(old.extractfile("suite-run3/summary.json").read())
    old_skips = Counter()
    for member in old:
        if member.isfile() and member.name.startswith("suite-run3/shard-") and member.name.endswith("/junit.xml"):
            old_skips.update(skipped_cases(ET.fromstring(old.extractfile(member).read())))
new_skips = Counter()
for path in sorted(evidence.glob("shard-*/junit.xml")):
    new_skips.update(skipped_cases(ET.parse(path).getroot()))
marker_delta = {
    "before": skip_counts(old_skips), "after": skip_counts(new_skips),
    "added": skip_rows(new_skips - old_skips),
    "removed": skip_rows(old_skips - new_skips),
}
assert marker_delta["before"] == {kind: old_summary["counts"].get(kind, 0) for kind in ("skipped", "xfailed")}
assert marker_delta["after"] == {kind: suite["counts"].get(kind, 0) for kind in ("skipped", "xfailed")}
(raw / "skip-xfail-delta.json").write_text(json.dumps(marker_delta, indent=2)+"\n")
after = json.loads((evidence / "collection/execution.json").read_text())
old_nodes, new_nodes = Counter(before["collected"]), Counter(after["collected"])
removed, added = sorted((old_nodes-new_nodes).elements()), sorted((new_nodes-old_nodes).elements())
renamed_old = "tests/test_native_matcher_parity.py::test_noncommutative_declines"
renamed_new = "tests/test_native_matcher_parity.py::test_commutative_wild_cannot_consume_noncommutative_prefix"
assert removed == [renamed_old], removed
assert renamed_new in added
new_files = {"tests/test_native_matcher_wave47.py", "tests/test_polytools_wave47.py", "tests/test_series_wave47.py"}
assert all(node == renamed_new or node.split("::")[0] in new_files for node in added)
delta = {"before": sum(old_nodes.values()), "after": sum(new_nodes.values()),
         "added": added, "removed": removed,
         "renamed_same_input": {"before": renamed_old, "after": renamed_new},
         "new_cases": len(added)-1}
(raw / "collection-delta.json").write_text(json.dumps(delta, indent=2)+"\n")
with tarfile.open(raw / "suite-receipts.tar.gz", "w:gz") as archive:
    for path in sorted(evidence.rglob("*")):
        if path.is_file() and (path.name in {"command.json", "execution.json", "junit.xml", "pytest.log", "result.json", "identity.json", "summary.json", "shards.json"}):
            archive.add(path, arcname="suite-final/"+str(path.relative_to(evidence)), recursive=False)
with tarfile.open(raw / "ode-receipts.tar.gz", "w:gz") as archive:
    for path in sorted((validate.EVIDENCE / "ode").rglob("*")):
        if path.is_file():
            archive.add(path, arcname="ode/"+str(path.relative_to(validate.EVIDENCE / "ode")), recursive=False)
result = {
    "head": current["head"], "source_sha256": current["source_sha256"],
    "extension_sha256": current["so_sha256"], "suite_counts": suite["counts"],
    "suite_coverage": suite["coverage"], "commands": commands,
    "collection_delta": {k:v for k,v in delta.items() if k not in {"added", "removed"}},
    "skip_xfail_delta": marker_delta,
    "lie_boundary_environment": boundary,
    "ode_boundary_observations": {
        "keys": ode_boundary["keys"], "pythonhashseed": ode_boundary["pythonhashseed"],
        "timeout_per_call_seconds": ode_boundary["timeout_per_call_seconds"],
        "memory_gib_per_worker": ode_boundary["memory_gib_per_worker"],
        "interpretation": ode_boundary["interpretation"],
    },
    "corpus_observation_deltas": {
        name: {"counts": family["counts"], "trace_counts": family["traces"]["counts"]}
        for name, family in corpora_delta["families"].items()
    },
    "limits_observation_deltas": limits_delta["summary"],
    "artifact_sha256": {str(p.relative_to(raw)): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted(raw.rglob("*")) if p.is_file() and p.name != "final-validation.json"},
}
(raw / "final-validation.json").write_text(json.dumps(result, indent=2)+"\n")
print("Consistent receipts:", len(commands), "commands completed with status 0;", suite["counts"], "new cases:", delta["new_cases"])
