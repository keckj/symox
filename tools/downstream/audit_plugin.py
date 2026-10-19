"""Record actual upstream pytest outcomes without changing their assertions."""
import json
import os
from pathlib import Path

REPORTS = []
COLLECTION = []
COLLECTION_SKIPS = []
os.environ.setdefault('SYMOX_DOWNSTREAM_PYTEST_ROOT_PID', str(os.getpid()))


def pytest_runtest_logreport(report):
    REPORTS.append(dict(nodeid=report.nodeid, phase=report.when, outcome=report.outcome,
                        duration=report.duration, wasxfail=getattr(report, 'wasxfail', None),
                        detail=str(report.longrepr)[:16000] if report.failed else ''))


def pytest_collectreport(report):
    if report.failed:
        COLLECTION.append(dict(nodeid=report.nodeid, detail=str(report.longrepr)[:16000]))
    elif report.skipped:
        COLLECTION_SKIPS.append(dict(nodeid=report.nodeid, detail=str(report.longrepr)[:2000]))


def pytest_sessionfinish(session, exitstatus):
    path = os.environ.get('SYMOX_DOWNSTREAM_TEST_REPORT')
    if path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if os.environ['SYMOX_DOWNSTREAM_PYTEST_ROOT_PID'] != str(os.getpid()):
            destination = destination.with_name(destination.stem+'.'+str(os.getpid())+destination.suffix)
        destination.write_text(json.dumps(dict(exitstatus=int(exitstatus), collected=session.testscollected,
                                               reports=REPORTS, collection_errors=COLLECTION,
                                               collection_skips=COLLECTION_SKIPS), indent=2))
