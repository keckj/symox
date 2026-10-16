"""Group saved fuzz findings by minimal case and operation."""

import argparse
from collections import defaultdict
import json


ORDER = {"died": 0, "timeout": 1, "crash": 2, "wrong": 3, "slow": 4, "weak": 5, "form": 6}


def summarize(records, classes=None):
    groups = defaultdict(list)
    for record in records:
        if classes is None or record["cls"] in classes:
            groups[(record["cls"], record["min_src"])].append(record)
    for (cls, src), rows in sorted(groups.items(), key=lambda pair: (ORDER.get(pair[0][0], 9), len(pair[0][1]), pair[0][1])):
        checks = ",".join(sorted({row["check"] for row in rows}))
        yield f"[{cls}] {src}   <{checks}> x{len(rows)}"
        yield f"      {rows[0]['min_msg'][:300]}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="findings.jsonl produced by fuzz_core.py")
    parser.add_argument("classes", nargs="?", help="comma-separated classes to include")
    args = parser.parse_args(argv)
    classes = set(args.classes.split(",")) if args.classes else None
    if classes is not None and not classes <= ORDER.keys():
        parser.error("unknown finding class")
    try:
        with open(args.path) as stream:
            records = [json.loads(line) for line in stream if line.strip()]
        for line in summarize(records, classes):
            print(line)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
