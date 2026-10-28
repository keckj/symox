#!/usr/bin/env bash
set -euo pipefail
: "${SYMOX_FROZEN_SNAPSHOT:?Set the immutable snapshot path supplied after root freeze}"
: "${SYMOX_NATIVE_SHA:?Set the complete frozen native binary SHA256}"
task_python="${SYMOX_AUDIT_PYTHON:-/home/keckj/symox-ci/venv/bin/python}"
task_archive="$(cd -- "$(dirname -- "$0")" && pwd)"
task_snapshot="$(realpath -- "$SYMOX_FROZEN_SNAPSHOT")"
case "$task_snapshot" in
  /home/keckj/Documents/symox|/home/keckj/Documents/symox/*) exit 23;;
esac
"$task_python" - "$task_snapshot" "$SYMOX_NATIVE_SHA" "$task_archive" <<'PY'
import hashlib,sys,pathlib,importlib.util,json
snapshot=pathlib.Path(sys.argv[1]);sha=sys.argv[2];archive=pathlib.Path(sys.argv[3])
assert len(sha)==64
assert hashlib.sha256((snapshot/'python/symox/_core.abi3.so').read_bytes()).hexdigest()==sha
assert hashlib.sha256((archive/'replay-integration.py').read_bytes()).hexdigest()=='9279d73a4f3d2980e6d4fcbf2e24fbbe007f28ca668024d3f5ebc02609d84236'
assert hashlib.sha256((archive/'replay-ode.py').read_bytes()).hexdigest()=='badcb0048e01375ac92b53cf450cbec1aa332b505982c9d49bc587e77481882e'
source=snapshot/'tools/integration_corpus.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='6f6262bdb0a886e0c79a5ae11eeea34742616a66b46fde46f74ceb066f6c1d94'
spec=importlib.util.spec_from_file_location('cases',source);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
cases=module.build_cases();assert len(cases)==1210 and len({c['id'] for c in cases})==1210
import sympy
assert sympy.__version__=='1.14.0'
(archive/'preflight.json').write_text(json.dumps({'snapshot':str(snapshot),'native_sha256':sha,'reference_version':sympy.__version__,'reference_module':sympy.__file__,'python':sys.executable,'expected_cases':1210,'corpus_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()},indent=2)+'\n')
PY
PYTHONPATH="$task_snapshot/python" "$task_python" -u "$task_archive/replay-integration.py" \
  --snapshot "$task_snapshot" --native-sha "$SYMOX_NATIVE_SHA" \
  --jobs 16 --timeout 60 --out "$task_archive/integration-r1-cb5c0c98" \
  2>&1 | tee "$task_archive/integration-r1-cb5c0c98.log"
