#!/usr/bin/env python3
"""Compare recorded wave-46/47 corpora without replaying any computation.

Writes every changed payload, added/missing case and trace delta. Exit 2
means missing/invalid inputs or unexpected coverage, not a mathematical
verdict. Deltas alone do not make the script fail: reviewers classify them.
"""

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
ROOT = AUDIT.parents[2]
TIMES = {"seconds", "elapsed_seconds"}
TRACES = {"python_fallbacks", "trace", "native"}
SPECS = {
    "tsolve": ("tsolve-symox.json", "mapping", 115),
    "unrad": ("unrad-symox.json", "mapping", 90),
    "matcher": ("matcher-after.json", "id", 2311),
    "algfactor": ("algfactor-symox.jsonl", "case", 249),
    "algpoly": ("algpoly-symox.jsonl", "case", 139),
    "ratint2": ("ratint2-symox.jsonl", "case", 75),
    "pde": ("pde-symox.json", "id", 791),
    "liegroup": ("liegroup-symox.json", "id", 428),
    "fps-formal": ("fps-formal.json", "fps", 220),
    "fps-extra": ("fps-extra.json", "fps", 127),
    "ode": (None, "ode", 945),
}


def no_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key!r}")
        result[key] = value
    return result


def read_file(path, issues):
    if not path.is_file():
        issues.append({"kind": "missing_file", "path": str(path)})
        return None
    try:
        if path.suffix == ".jsonl":
            return [json.loads(line, object_pairs_hook=no_duplicate_keys)
                    for line in path.read_text().splitlines() if line.strip()]
        return json.loads(path.read_text(), object_pairs_hook=no_duplicate_keys)
    except (OSError, ValueError) as error:
        issues.append({"kind": "invalid_file", "path": str(path), "error": str(error)})
        return None


def load_rows(root, family, filename, kind, issues):
    files = sorted(root.glob("[0-9][0-9].json")) if kind == "ode" else [root / filename]
    if not files:
        issues.append({"kind": "missing_ode_parts", "path": str(root)})
    rows, locations = {}, {}
    for path in files:
        data = read_file(path, issues)
        if data is None:
            continue
        if kind in ("mapping", "fps"):
            if not isinstance(data, dict):
                issues.append({"kind": "invalid_schema", "path": str(path), "expected": "object"})
                continue
            entries = list(data.items())
        elif isinstance(data, list):
            entries = []
            for index, row in enumerate(data):
                try:
                    key = json.dumps([row["id"], row["method"]]) if kind == "ode" else str(row[kind])
                except (KeyError, TypeError) as error:
                    issues.append({"kind": "unmatched_row", "path": str(path), "index": index,
                                   "error": str(error), "row": row})
                    continue
                entries.append((key, row))
        else:
            issues.append({"kind": "invalid_schema", "path": str(path), "expected": "array"})
            continue
        for index, (key, row) in enumerate(entries):
            location = {"path": str(path), "index": index}
            if key in rows:
                issues.append({"kind": "duplicate_case", "family": family, "key": key,
                               "first_location": locations[key], "duplicate_location": location,
                               "first_row": rows[key], "duplicate_row": row})
                # Preserve the duplicate as an explicitly unmatched synthetic
                # row rather than silently overwriting either observation.
                key = f"DUPLICATE:{key}:{path.name}:{index}"
            rows[key], locations[key] = row, location
    return rows, [str(path) for path in files]


def strip_times(value):
    if isinstance(value, dict):
        return {key: strip_times(item) for key, item in value.items() if key not in TIMES}
    if isinstance(value, list):
        return [strip_times(item) for item in value]
    return value


def split_traces(row, family):
    # These corpus schemas keep execution traces at row level. Do not
    # reinterpret similarly named keys inside mathematical payloads.
    if not isinstance(row, dict):
        return row, {}
    keys = TRACES if family in ("pde", "liegroup") else TRACES - {"native"}
    # In matcher rows, `native` is a result implementation, not a trace.
    return ({key: value for key, value in row.items() if key not in keys},
            {key: value for key, value in row.items() if key in keys})


def compare(before, after):
    left, right = set(before), set(after)
    shared = sorted(left & right)
    changed = [key for key in shared if before[key] != after[key]]
    return {
        "counts": {"before": len(before), "after": len(after), "matched": len(shared),
                   "unchanged": len(shared) - len(changed), "changed": len(changed),
                   "missing_after": len(left - right), "added_after": len(right - left)},
        "missing_after": [{"key": key, "before": before[key]} for key in sorted(left - right)],
        "added_after": [{"key": key, "after": after[key]} for key in sorted(right - left)],
        "deltas": [{"key": key, "before": before[key], "after": after[key]} for key in changed],
    }


def corpus_normalizer(family):
    date = "2026-11-07" if family == "pde" else "2026-11-06"
    path = ROOT / f"audits/performance/{date}-{family}/scripts/corpus.py"
    spec = importlib.util.spec_from_file_location("normalizer_" + family, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module._normal, str(path)


def prepare(rows, family, normal, issues):
    main, traces, raw = {}, {}, {}
    for key, row in rows.items():
        if family.startswith("fps-"):
            if not isinstance(row, list) or len(row) < 3:
                issues.append({"kind": "invalid_fps_row", "family": family, "key": key, "row": row})
                main[key] = raw[key] = row
            else:
                main[key] = raw[key] = row[:2] + row[3:]
            traces[key] = {}
            continue
        payload, traces[key] = split_traces(strip_times(row), family)
        raw[key] = payload
        # Apply precisely the corpus helper to the rendered result value;
        # retain error class/message and all other result/row fields.
        if normal and isinstance(payload, dict) and isinstance(payload.get("result"), dict):
            result = payload["result"]
            if "value" in result:
                payload = dict(payload, result=dict(result, value=normal(result["value"])))
        main[key] = payload
    return main, traces, raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, default=ROOT / "audits/performance/2026-11-09-wave46/raw")
    parser.add_argument("--after", type=Path, default=AUDIT / "raw")
    parser.add_argument("--ode-before", type=Path, default=Path("/tmp/symox-wave47-baseline/ode"))
    parser.add_argument("--ode-after", type=Path, default=Path("/tmp/symox-wave47-evidence/ode"))
    parser.add_argument("--output", type=Path, default=AUDIT / "raw/corpora-deltas.json")
    parser.add_argument("--families", nargs="+", choices=list(SPECS), default=list(SPECS))
    args = parser.parse_args()
    report = {"complete": True, "ignored_timing_keys": sorted(TIMES),
              "fps_ignored_index": 2, "trace_keys_separated": sorted(TRACES), "families": {}}
    for family in args.families:
        filename, kind, expected = SPECS[family]
        issues = []
        before_root, after_root = ((args.ode_before, args.ode_after) if family == "ode"
                                   else (args.before, args.after))
        before, left_files = load_rows(before_root, family, filename, kind, issues)
        after, right_files = load_rows(after_root, family, filename, kind, issues)
        for label, rows in (("before", before), ("after", after)):
            if len(rows) != expected:
                issues.append({"kind": "unexpected_coverage", "side": label,
                               "expected": expected, "observed": len(rows)})
        normal, normal_source = corpus_normalizer(family) if family in ("pde", "liegroup") else (None, None)
        left, left_trace, left_raw = prepare(before, family, normal, issues)
        right, right_trace, right_raw = prepare(after, family, normal, issues)
        result = compare(left, right)
        result.update({"expected_cases": expected, "before_files": left_files,
                       "after_files": right_files, "issues": issues,
                       "normalizer": normal_source, "traces": compare(left_trace, right_trace)})
        if normal:
            raw_changes = compare(left_raw, right_raw)
            changed_keys = {row["key"] for row in result["deltas"]}
            result["representation_only_deltas"] = [
                row for row in raw_changes["deltas"] if row["key"] not in changed_keys
            ]
        if family == "matcher":
            result["implementations"] = {
                name: compare({key: row.get(name) for key, row in left.items()},
                              {key: row.get(name) for key, row in right.items()})
                for name in ("native", "public", "python_port", "pattern")
            }
        if family == "ode":
            result["verdicts"] = compare(
                {key: row.get("verdict") for key, row in left.items()},
                {key: row.get("verdict") for key, row in right.items()})
            result["before_verdict_counts"] = dict(Counter(row.get("verdict") for row in left.values()))
            result["after_verdict_counts"] = dict(Counter(row.get("verdict") for row in right.values()))
        report["families"][family] = result
        report["complete"] &= not issues and not result["missing_after"] and not result["added_after"]
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"complete": report["complete"], "families": {
        family: {"counts": result["counts"], "issues": len(result["issues"]),
                 "trace_changes": result["traces"]["counts"]["changed"]}
        for family, result in report["families"].items()}}, sort_keys=True))
    return 0 if report["complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
