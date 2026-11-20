"""Build a frozen source tree; retain failed attempts and copy only success."""
import fcntl
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import validate

tag = sys.argv[1]
assert re.fullmatch(r"[a-z0-9_-]+", tag)
raw = validate.RAW
raw.mkdir(exist_ok=True)
with open("/home/keckj/symox-ci/heavy.lock", "rb") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    before = validate.identity()
    before["previous_so_sha256"] = before.pop("so_sha256")
    (raw / ("build-input-" + tag + ".json")).write_text(json.dumps(before, indent=2) + "\n")
    command = [validate.CARGO, "build", "--offline", "--release"]
    started = time.time()
    with (raw / ("build-" + tag + ".log")).open("w") as output:
        process = subprocess.run(command, cwd=validate.ROOT, env=validate.ENV,
                                 stdout=output, stderr=subprocess.STDOUT)
    unchanged = validate.identity()["source_sha256"] == before["source_sha256"]
    record = {
        "command": command, "returncode": process.returncode,
        "source_unchanged": unchanged, "elapsed_seconds": round(time.time() - started, 3),
        "environment": {key: validate.ENV[key] for key in ["CARGO_TARGET_DIR", "PYO3_PYTHON", "PATH"]},
    }
    if process.returncode == 0 and unchanged:
        built = Path(validate.TARGET) / "release/lib_core.so"
        shutil.copyfile(built, validate.ROOT / "python/symox/_core.abi3.so")
        record["extension_sha256"] = hashlib.sha256(built.read_bytes()).hexdigest()
        (raw / "build-input.json").write_text(json.dumps(before, indent=2) + "\n")
        (raw / "build-output.json").write_text(json.dumps(record, indent=2) + "\n")
    (raw / ("build-" + tag + ".command.json")).write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)
    if process.returncode:
        print((raw / ("build-" + tag + ".log")).read_text()[-12000:])
        raise SystemExit(process.returncode)
    if not unchanged:
        raise SystemExit("Sources changed during build; extension was not copied.")
