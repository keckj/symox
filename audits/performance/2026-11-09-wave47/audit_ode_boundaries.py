#!/usr/bin/env python3
"""Paired resource observations using the unmodified archived ODE Worker.

Preparation is cheap; invoking this script starts the selected heavy calls.
The controller acquires heavy.lock and runs fresh isolated workers for the
same five calls on baseline/final. Resource outcomes remain inconclusive.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time


AUDIT = Path(__file__).resolve().parent
ROOT = AUDIT.parents[2]
RUNNER = ROOT / "audits/sympy-1.14.0/corpus-2026-10-26/replay-ode.py"
CASES = RUNNER.with_name("ode-cases.json.gz")
KEYS = (
    "2nd_nonlinear_autonomous_conserved_03:default",
    "2nd_nonlinear_autonomous_conserved_04:default",
    "lie_group_01:default",
    "lie_group_03:default",
    "lie_group_03:hint",
)
TIMEOUT = 60
SUPERVISOR_TIMEOUT = 130  # archived startup bound 45 + call bound 60 + cleanup


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archived_runner():
    spec = importlib.util.spec_from_file_location("archived_ode_worker", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def item_for(key):
    case_id, method = key.rsplit(":", 1)
    cases = json.loads(gzip.decompress(CASES.read_bytes()))["cases"]
    selected = [case for case in cases if case["id"] == case_id]
    if len(selected) != 1 or method not in ("default", "hint"):
        raise ValueError(f"unmatched ODE key: {key}")
    return dict(selected[0], method=method)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def one_observation(args):
    """Internal helper; its environment selects the baseline/final import."""
    archive = archived_runner()
    record = {"key": args.one, "runtime": str(args.runtime),
              "expected_extension_sha256": args.expected_sha,
              "environment": {key: os.environ.get(key) for key in (
                  "PYTHONPATH", "PYTHONHASHSEED", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
              "worker_identities": [], "worker_stop_events": []}
    with args.out.with_suffix(".worker.stderr").open("w") as log:
        worker = archive.Worker("symox", TIMEOUT, log)
        original_stop = worker.stop

        def observed_stop():
            # Delegate to exactly Worker.stop; only record its effects.
            process = worker.process
            event = ({"pid": process.pid, "returncode_before_stop": process.poll()}
                     if process is not None else None)
            original_stop()
            if event is not None:
                event["returncode_after_stop"] = process.returncode
                record["worker_stop_events"].append(event)

        worker.stop = observed_stop
        started = time.monotonic()
        try:
            record["result"] = worker.run(item_for(args.one))
        except Exception as error:
            record["audit_error"] = {"type": type(error).__name__, "message": str(error)}
        finally:
            worker.stop()
            record["worker_identities"] = worker.identities
            record["wall_seconds"] = time.monotonic() - started
    record["identity_matches"] = bool(worker.identities) and all(
        identity.get("native_sha256") == args.expected_sha
        and Path(identity["module"]).resolve() == (args.runtime / "symox/__init__.py").resolve()
        for identity in worker.identities
    )
    write_json(args.out, record)
    return 0 if "audit_error" not in record and record["identity_matches"] else 2


def run_observation(out, key, label, runtime, expected_sha, round_index):
    stem = f"{key.replace(':', '-')}-round{round_index + 1}-{label}"
    path = out / (stem + ".json")
    env = dict(os.environ, PYTHONPATH=str(runtime), PYTHONHASHSEED="0",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    command = [sys.executable, str(Path(__file__).resolve()), "--one", key,
               "--runtime", str(runtime), "--expected-sha", expected_sha, "--out", str(path)]
    receipt = {"command": command, "key": key, "label": label, "round": round_index + 1,
               "stdout": stem + ".helper.stdout", "stderr": stem + ".helper.stderr"}
    started = time.monotonic()
    with (out / receipt["stdout"]).open("w") as stdout, (out / receipt["stderr"]).open("w") as stderr:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr,
                                   start_new_session=True)
        try:
            process.wait(timeout=SUPERVISOR_TIMEOUT)
        except subprocess.TimeoutExpired:
            receipt["audit_error"] = "SUPERVISOR_DEADLINE; entire helper process group killed"
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
    receipt["returncode"] = process.returncode
    receipt["wall_seconds"] = time.monotonic() - started
    receipt["observation_file"] = str(path)
    if path.is_file():
        receipt["observation"] = json.loads(path.read_text())
    else:
        receipt["audit_error"] = receipt.get("audit_error", "helper produced no observation")
    write_json(out / (stem + ".command.json"), receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, default=Path("/tmp/symox-wave47-baseline/python"))
    parser.add_argument("--final", type=Path, default=ROOT / "python")
    parser.add_argument("--lock", type=Path, default=Path("/home/keckj/symox-ci/heavy.lock"))
    parser.add_argument("--lock-held", action="store_true",
                        help="supervisor already holds heavy.lock; verify it is busy without blocking")
    parser.add_argument("--rounds", type=int, choices=(1, 2), default=1)
    parser.add_argument("--jobs", type=int, choices=(1, 2), default=1)
    parser.add_argument("--keys", nargs="+", choices=KEYS, default=list(KEYS))
    parser.add_argument("--one", choices=KEYS, help=argparse.SUPPRESS)
    parser.add_argument("--runtime", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--expected-sha", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.one:
        return one_observation(args)
    if len(set(args.keys)) != len(args.keys):
        parser.error("duplicate selected keys")
    # Resolve/validate every selected key before waiting for the heavy lock.
    items = {key: item_for(key) for key in args.keys}
    runtimes = {"baseline": args.baseline.resolve(), "final": args.final.resolve()}
    hashes = {label: sha(path / "symox/_core.abi3.so") for label, path in runtimes.items()}
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    archive = archived_runner()
    summary = {"complete": False, "keys": args.keys, "items": items,
               "rounds": args.rounds, "parallel_keys": args.jobs,
               "timeout_per_call_seconds": TIMEOUT, "memory_gib_per_worker": 8,
               "worker_startup_bound_seconds": 45, "supervisor_bound_seconds": SUPERVISOR_TIMEOUT,
               "pythonhashseed": "0", "fresh_worker_per_observation": True,
               "archived_runner": str(RUNNER), "archived_runner_sha256": sha(RUNNER),
               "worker_payload_sha256": hashlib.sha256(archive.WORKER.encode()).hexdigest(),
               "cases_sha256": sha(CASES), "audit_script_sha256": sha(Path(__file__)),
               "runtime_paths": {key: str(value) for key, value in runtimes.items()},
               "extension_sha256": hashes, "host": platform.node(), "python": sys.version,
               "lock": str(args.lock), "pairs": [],
               "interpretation": "Paired output/error/resource observations only. No regression verdict is predeclared. Fresh-worker controls do not reproduce the original reused worker history."}
    write_json(args.out / "summary.json", summary)
    waiting = time.monotonic()
    with args.lock.open("r") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            if not args.lock_held:
                summary["audit_error"] = "heavy.lock is busy; no observation ran"
                write_json(args.out / "summary.json", summary)
                return 3
            summary["lock_mode"] = "supervisor-held; verified busy by nonblocking flock"
        else:
            if args.lock_held:
                summary["audit_error"] = "--lock-held requires the supervising process to hold heavy.lock; no observation ran"
                write_json(args.out / "summary.json", summary)
                return 3
            summary["lock_mode"] = "audit acquired exclusive lock"
        summary["lock_wait_seconds"] = time.monotonic() - waiting
        # A queued audit must not silently observe another installed binary.
        if any(sha(runtimes[label] / "symox/_core.abi3.so") != expected for label, expected in hashes.items()):
            summary["audit_error"] = "extension identity changed while waiting for lock"
            write_json(args.out / "summary.json", summary)
            return 2

        def paired(key):
            pairs = []
            for index in range(args.rounds):
                order = ("baseline", "final") if index % 2 == 0 else ("final", "baseline")
                calls = {label: run_observation(args.out, key, label, runtimes[label], hashes[label], index)
                         for label in order}
                pair = {"key": key, "round": index + 1, "order": order, "calls": calls}
                left, right = (calls[label].get("observation", {}).get("result") for label in ("baseline", "final"))
                pair["paired_observation_verdict"] = archive.verdict(left, right) if left and right else "audit_failure_review_required"
                write_json(args.out / f"{key.replace(':', '-')}-round{index + 1}-pair.json", pair)
                pairs.append(pair)
            return pairs

        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            for pairs in executor.map(paired, args.keys):
                summary["pairs"].extend(pairs)
                write_json(args.out / "summary.json", summary)
        summary["extension_sha256_after"] = {
            label: sha(path / "symox/_core.abi3.so") for label, path in runtimes.items()
        }
        observations = [call for pair in summary["pairs"] for call in pair["calls"].values()]
        summary["complete"] = (len(summary["pairs"]) == len(args.keys) * args.rounds
                               and summary["extension_sha256_after"] == hashes
                               and all(call["returncode"] == 0 and "audit_error" not in call
                                       and call.get("observation", {}).get("identity_matches") for call in observations))
        write_json(args.out / "summary.json", summary)
    print(json.dumps({"complete": summary["complete"], "pairs": len(summary["pairs"]),
                      "verdicts": [pair["paired_observation_verdict"] for pair in summary["pairs"]]}))
    return 0 if summary["complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
