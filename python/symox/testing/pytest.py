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
    """SymPy's ``warns``: the expected category must be emitted and match ``match``.

    The stack-level check of the reference (which needs a ``test_*`` file
    on the stack) is not performed."""
    with warnings.catch_warnings(record=True) as warnrec:
        warnings.simplefilter('always')
        yield warnrec
    if not any(issubclass(w.category,warningcls) for w in warnrec):
        raise Failed('Failed: DID NOT WARN. No warnings of type %s was emitted.'
                     ' The list of emitted warnings is: %s.'%(warningcls,[w.message for w in warnrec]))
    for w in warnrec:
        if issubclass(w.category,warningcls) and not re.compile(match,re.IGNORECASE).match(str(w.message)):
            raise Failed(f"Failed: WRONG MESSAGE. A warning with of the correct category ({warningcls.__name__}) "
                         f"was issued, but it did not match the given match regex ({match!r})")
    for w in warnrec:
        if not issubclass(w.category,warningcls):
            warnings.warn_explicit(w.message,w.category,w.filename,w.lineno)
def warns_deprecated_sympy(*,match='',test_stacklevel=True):return warns(SymPyDeprecationWarning,match=match,test_stacklevel=test_stacklevel)

__all__=['ExceptionInfo','Failed','RaisesContext','SKIP','Skipped','XFAIL','XFail','XPass','nocache_fail','raises',
 'skip','skip_under_pyodide','slow','tooslow','warns','warns_deprecated_sympy']
