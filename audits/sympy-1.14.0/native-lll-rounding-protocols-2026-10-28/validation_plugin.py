"""Record collection and execution without changing test selection or outcomes."""
import json
import os
from pathlib import Path

_items = []
_files = set()
_finished = []


def pytest_collectreport(report):
    name = report.nodeid.split('::', 1)[0]
    if name.startswith('tests/') and name.endswith('.py'):
        _files.add(name)


def pytest_collection_finish(session):
    _items.extend(item.nodeid for item in session.items)
    output = Path(os.environ['SYMOX_VALIDATION_RECORD'])
    output.write_text(json.dumps({'files': sorted(_files), 'collected': _items}, indent=2))


def pytest_runtest_logfinish(nodeid, location):
    _finished.append(nodeid)


def pytest_sessionfinish(session, exitstatus):
    output = Path(os.environ['SYMOX_VALIDATION_RECORD'])
    output.write_text(json.dumps({'files': sorted(_files), 'collected': _items,
                                 'finished': _finished, 'exit_code': int(exitstatus)}, indent=2))
