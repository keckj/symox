"""Compare frozen wave-46 and wave-47 limits2 artifacts without replaying.

Case rows have the strict schema [value, seconds, source_line]. Only seconds,
at index 1, are excluded. Statistics and traces are compared separately.
Missing files, missing/extra cases, and malformed rows refuse the report.
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys


AUDIT = Path(__file__).resolve().parent
EXPECTED_COUNTS = {
    "limits": 554,
    "gruntz": 107,
    "series": 83,
    "nseries": 186,
    "order": 210,
    "formal": 85,
    "aseries": 9,
    "lseries": 5,
    "demidovich": 54,
}


class Refuse(ValueError):
    pass


def read(path, inputs):
    payload = path.read_bytes()
    inputs[str(path)] = hashlib.sha256(payload).hexdigest()
    return json.loads(payload)


def validate_cases(data, path, count):
    if not isinstance(data, dict) or len(data) != count:
        raise Refuse(f"{path}: expected {count} case rows")
    for key, row in data.items():
        if not isinstance(key, str) or not isinstance(row, list) or len(row) != 3:
            raise Refuse(f"{path}: malformed case {key!r}; expected [value, seconds, source_line]")
        value, seconds, source_line = row
        if (not isinstance(value, str)
                or type(seconds) not in (float, int) or not math.isfinite(seconds) or seconds < 0
                or type(source_line) is not int or source_line < 0):
            raise Refuse(f"{path}: invalid value/time/source_line types for {key!r}")


def validate_stats(data, path):
    if (not isinstance(data, dict)
            or not all(isinstance(key, str) and type(value) is int and value >= 0
                       for key, value in data.items())):
        raise Refuse(f"{path}: expected nonnegative integer statistics")


def validate_trace(data, path):
    if (not isinstance(data, list)
            or not all(isinstance(row, list) and len(row) == 2
                       and all(isinstance(item, str) for item in row) for row in data)):
        raise Refuse(f"{path}: expected [reason, expression] trace rows")


def counter_rows(counts):
    return [{"entry": list(entry), "count": count}
            for entry, count in sorted(counts.items())]


def compare(old_dir, new_dir, output):
    paths = {}
    for category in EXPECTED_COUNTS:
        for suffix in ("json", "stats.json", "trace.json"):
            paths[category, suffix] = (
                old_dir / f"wave46_{category}.{suffix}",
                new_dir / f"wave47_{category}.{suffix}",
            )
    missing = sorted(str(path) for pair in paths.values() for path in pair if not path.is_file())
    if missing:
        raise Refuse("missing required files: " + ", ".join(missing))

    inputs = {}
    categories = {}
    auxiliary = {}
    for category, count in EXPECTED_COUNTS.items():
        old_path, new_path = paths[category, "json"]
        old, new = read(old_path, inputs), read(new_path, inputs)
        validate_cases(old, old_path, count)
        validate_cases(new, new_path, count)
        if old.keys() != new.keys():
            raise Refuse(f"{category}: case keys differ; missing={sorted(old.keys() - new.keys())!r}, "
                         f"extra={sorted(new.keys() - old.keys())!r}")
        value_deltas = []
        source_line_deltas = []
        for key in sorted(old):
            if old[key][0] != new[key][0]:
                value_deltas.append({"case": key, "before": old[key][0], "after": new[key][0]})
            if old[key][2] != new[key][2]:
                source_line_deltas.append({"case": key, "before": old[key][2], "after": new[key][2]})
        categories[category] = {
            "cases": count,
            "same_case_keys": True,
            "value_deltas": value_deltas,
            "source_line_deltas": source_line_deltas,
        }

        old_path, new_path = paths[category, "stats.json"]
        old, new = read(old_path, inputs), read(new_path, inputs)
        validate_stats(old, old_path)
        validate_stats(new, new_path)
        auxiliary[f"{category}.stats.json"] = {
            "identical": old == new,
            "deltas": [{"counter": key, "before": old.get(key), "after": new.get(key)}
                       for key in sorted(old.keys() | new.keys()) if old.get(key) != new.get(key)],
        }

        old_path, new_path = paths[category, "trace.json"]
        old, new = read(old_path, inputs), read(new_path, inputs)
        validate_trace(old, old_path)
        validate_trace(new, new_path)
        before_counts, after_counts = Counter(map(tuple, old)), Counter(map(tuple, new))
        auxiliary[f"{category}.trace.json"] = {
            "identical": old == new,
            "before_entries": len(old),
            "after_entries": len(new),
            "reordered_only": old != new and before_counts == after_counts,
            "removed": counter_rows(before_counts - after_counts),
            "added": counter_rows(after_counts - before_counts),
            "reasons_before": dict(sorted(Counter(row[0] for row in old).items())),
            "reasons_after": dict(sorted(Counter(row[0] for row in new).items())),
        }

    summary = {
        "cases": sum(item["cases"] for item in categories.values()),
        "value_deltas": sum(len(item["value_deltas"]) for item in categories.values()),
        "source_line_deltas": sum(len(item["source_line_deltas"]) for item in categories.values()),
        "statistics_files_changed": sum(not item["identical"] for name, item in auxiliary.items()
                                        if name.endswith("stats.json")),
        "trace_files_changed": sum(not item["identical"] for name, item in auxiliary.items()
                                   if name.endswith("trace.json")),
    }
    return {
        "status": "complete",
        "old_directory": str(old_dir),
        "new_directory": str(new_dir),
        "output": str(output),
        "case_row_schema": ["value", "seconds", "source_line"],
        "ignored_case_indices": [1],
        "summary": summary,
        "categories": categories,
        "auxiliary_file_deltas": auxiliary,
        "input_sha256": inputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-dir", type=Path,
                        default=AUDIT.parent / "2026-11-09-wave46/raw/limits2")
    parser.add_argument("--new-dir", type=Path, default=AUDIT / "raw/limits2")
    parser.add_argument("--output", type=Path, default=AUDIT / "raw/limits2-deltas.json")
    args = parser.parse_args()
    try:
        report = compare(args.old_dir.resolve(), args.new_dir.resolve(), args.output.resolve())
    except (Refuse, OSError, ValueError) as error:
        # Replace any stale successful report with an explicit refusal.
        report = {"status": "refused", "error": str(error)}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report), file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "complete", **report["summary"], "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
