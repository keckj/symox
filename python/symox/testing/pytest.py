"""Compatible assertion helpers backed by pytest and native warning classes."""
import sys
import pytest as _pytest
from ..utilities.exceptions import SymPyDeprecationWarning
ExceptionInfo=_pytest.ExceptionInfo
Failed=_pytest.fail.Exception
Skipped=_pytest.skip.Exception
XFail=_pytest.xfail.Exception
class XPass(AssertionError):pass


class RaisesContext:
    def __init__(self,expectedException,match=None):self.context=_pytest.raises(expectedException,match=match)
    def __enter__(self):return self.context.__enter__()
    def __exit__(self,*exc):return self.context.__exit__(*exc)


def raises(expectedException,code=None,**kwargs):
    if code is None:return RaisesContext(expectedException,**kwargs)
    if not callable(code):raise TypeError('raises requires a callable or a context manager')
    return _pytest.raises(expectedException,code,**kwargs)

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

def warns(warningcls,*,match='',test_stacklevel=True):return _pytest.warns(warningcls,match=match or None)
def warns_deprecated_sympy(*,match='',test_stacklevel=True):return warns(SymPyDeprecationWarning,match=match,test_stacklevel=test_stacklevel)

__all__=['ExceptionInfo','Failed','RaisesContext','SKIP','Skipped','XFAIL','XFail','XPass','nocache_fail','raises',
 'skip','skip_under_pyodide','slow','tooslow','warns','warns_deprecated_sympy']
