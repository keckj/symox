import json, os
from symox import _core

def pytest_sessionstart(session):
    if hasattr(_core, "_limit_trace"):
        _core._limit_trace(True)

def pytest_sessionfinish(session, exitstatus):
    out = os.environ.get("TRACE_OUT")
    if out and hasattr(_core, "_limit_trace"):
        json.dump(_core._limit_trace(False), open(out, "w"), indent=0)
