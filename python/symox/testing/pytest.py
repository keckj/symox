"""Compatible assertion helpers backed by pytest and native warning classes."""
import contextlib
import re
import sys
import warnings
import pytest as _pytest
from ..utilities.exceptions import SymPyDeprecationWarning
from ..utilities.exceptions import ignore_warnings  # noqa: F401  (SymPy re-exports it here)
ExceptionInfo=_pytest.ExceptionInfo
Failed=_pytest.fail.Exception
Skipped=_pytest.skip.Exception
XFail=_pytest.xfail.Exception
class XPass(AssertionError):pass


def _did_not_raise(error):
    # SymPy's own helper reports a missing exception as a bare DID NOT RAISE.
    if str(error).startswith('DID NOT RAISE'):raise Failed('DID NOT RAISE') from None
    raise error


class RaisesContext:
    def __init__(self,expectedException,match=None):self.context=_pytest.raises(expectedException,match=match)
    def __enter__(self):return self.context.__enter__()
    def __exit__(self,*exc):
        try:return self.context.__exit__(*exc)
        except Failed as error:_did_not_raise(error)


def raises(expectedException,code=None,**kwargs):
    if code is None:return RaisesContext(expectedException,**kwargs)
    if isinstance(code,str):
        raise TypeError('\'raises(xxx, "code")\' has been phased out; '
                        'change \'raises(xxx, "expression")\' '
                        'to \'raises(xxx, lambda: expression)\', '
                        '\'raises(xxx, "statement")\' '
                        'to \'with raises(xxx): statement\'')
    if not callable(code):raise TypeError('raises() expects a callable for the 2nd argument.')
    try:return _pytest.raises(expectedException,code,**kwargs)
    except Failed as error:_did_not_raise(error)

def skip(str):return _pytest.skip(str)
def SKIP(reason):return _pytest.mark.skip(reason=reason)
def XFAIL(func):return _pytest.mark.xfail(func)
def slow(func):func._slow=True;return func
def tooslow(func):func._tooslow=True;return func
def nocache_fail(func):func._nocache_fail=True;return func
def skip_under_pyodide(message=None):
    reason=message if isinstance(message,str) else 'unavailable on Pyodide'
    decorator=_pytest.mark.skipif(sys.platform=='emscripten',reason=reason)
    return decorator(message) if callable(message) else decorator

@contextlib.contextmanager
def warns(warningcls,*,match='',test_stacklevel=True):
    """SymPy's ``warns``: the expected category must be emitted and match
    ``match``; any other warning is an error, and with ``test_stacklevel``
    the warning must point at the calling ``test_*`` file (as in SymPy)."""
    import inspect, os
    with warnings.catch_warnings(record=True) as warnrec:
        warnings.simplefilter('error')
        warnings.filterwarnings('always',category=warningcls)
        yield warnrec
    if not any(issubclass(w.category,warningcls) for w in warnrec):
        msg=('Failed: DID NOT WARN.'
             ' No warnings of type %s was emitted.'
             ' The list of emitted warnings is: %s.')%(warningcls,[w.message for w in warnrec])
        raise Failed(msg)
    for w in warnrec:
        assert issubclass(w.category,warningcls)
        if not re.compile(match,re.IGNORECASE).match(str(w.message)):
            raise Failed(f"Failed: WRONG MESSAGE. A warning with of the correct category ({warningcls.__name__}) "
                         f"was issued, but it did not match the given match regex ({match!r})")
    if test_stacklevel:
        for f in inspect.stack():
            thisfile=f.filename
            file=os.path.split(thisfile)[1]
            if file.startswith('test_'):
                break
            elif file=='doctest.py':
                return
        else:
            raise RuntimeError('Could not find the file for the given warning to test the stacklevel')
        for w in warnrec:
            if w.filename!=thisfile:
                msg=f"""\
Failed: Warning has the wrong stacklevel. The warning stacklevel needs to be
set so that the line of code shown in the warning message is user code that
calls the deprecated code (the current stacklevel is showing code from
{w.filename} (line {w.lineno}), expected {thisfile})""".replace('\n',' ')
                raise Failed(msg)



@contextlib.contextmanager
def warns_deprecated_sympy():
    """Shorthand for ``warns(SymPyDeprecationWarning)``."""
    with warns(SymPyDeprecationWarning):
        yield

__all__=['ExceptionInfo','Failed','RaisesContext','SKIP','Skipped','XFAIL','XFail','XPass','nocache_fail','raises',
 'skip','skip_under_pyodide','slow','tooslow','warns','warns_deprecated_sympy']
