"""Project test and doctest runners using pytest and the standard library.

These adapters run Symox tests and Python doctests; they do not embed SymPy's
custom collection, subprocess-management or assertion-rewriting engine.
"""
import base64
from collections import namedtuple
from contextlib import contextmanager
import doctest as _doctest
import fnmatch
import importlib
import os
from pathlib import Path
import pickle
import random
import subprocess as _subprocess
import sys
import warnings


class DependencyError(ImportError):pass
class TimeOutError(TimeoutError):pass
from .pytest import Skipped

class SymPyTestResults(_doctest.TestResults):__slots__=()
class SymPyOutputChecker(_doctest.OutputChecker):pass
class SymPyDocTestFinder(_doctest.DocTestFinder):pass
class SymPyDocTestRunner(_doctest.DocTestRunner):
    def __init__(self,*args,**kwargs):kwargs.setdefault('checker',SymPyOutputChecker());super().__init__(*args,**kwargs)


def get_sympy_dir():
    """Return the local Symox source directory for compatible test clients."""
    return str(Path(__file__).resolve().parents[1])
def convert_to_native_paths(paths):return [os.path.normpath(path) for path in paths]
def setup_pprint(disable_line_wrap=True):
    from ..interactive import init_printing
    init_printing(no_global=True,wrap_line=not disable_line_wrap)

def raise_on_deprecated():
    from ..utilities.exceptions import SymPyDeprecationWarning
    warnings.filterwarnings('error',category=SymPyDeprecationWarning)


def split_list(l,split,density=None):
    if split is None:return list(l)
    first,count=map(int,split.split('/') if isinstance(split,str) else split)
    if not 1<=first<=count:raise ValueError('test split must be a/b with 1 <= a <= b')
    size=len(l)//count;start=size*(first-1);end=size*first if first<count else len(l)
    return list(l[start:end])


def _paths(paths,default):return [str(p) for p in paths] if paths else [str(default)]

def test(*paths,subprocess=True,rerun=0,**kwargs):
    import pytest
    root=Path(__file__).resolve().parents[3]
    args=_paths(paths,root/'tests')
    for key,flag in [('verbose','-v'),('pdb','--pdb'),('failfast','-x')]:
        if kwargs.pop(key,False):args.append(flag)
    if kwargs.get('kw'):args+=['-k',str(kwargs.pop('kw'))]
    if kwargs.get('tb'):args+=['--tb',str(kwargs.pop('tb'))]
    if kwargs.pop('doctest',False):args.append('--doctest-modules')
    kwargs.pop('colors',None);kwargs.pop('force_colors',None);kwargs.pop('seed',None);kwargs.pop('slow',None);kwargs.pop('enhance_asserts',None)
    if kwargs:raise TypeError('unsupported pytest runner options: '+', '.join(kwargs))
    for attempt in range(int(rerun)+1):
        if subprocess:success=_subprocess.run([sys.executable,'-m','pytest',*args]).returncode==0
        else:success=pytest.main(args)==0
        if success:return True
    return False


def sympytestfile(filename,module_relative=True,name=None,package=None,globs=None,verbose=None,report=True,optionflags=0,extraglobs=None,raise_on_error=False,parser=None,encoding=None):
    kwargs=dict(module_relative=module_relative,name=name,package=package,globs=globs,verbose=verbose,report=report,optionflags=optionflags,extraglobs=extraglobs,raise_on_error=raise_on_error,encoding=encoding)
    if parser is not None:kwargs['parser']=parser
    result=_doctest.testfile(filename,**kwargs)
    return SymPyTestResults(*result)


def doctest(*paths,subprocess=True,rerun=0,**kwargs):
    paths=_paths(paths,Path(get_sympy_dir()))
    verbose=kwargs.pop('verbose',False);flags=kwargs.pop('optionflags',_doctest.ELLIPSIS|_doctest.NORMALIZE_WHITESPACE)
    if kwargs:raise TypeError('unsupported doctest options: '+', '.join(kwargs))
    failures=0
    for path in paths:
        p=Path(path)
        files=list(p.rglob('*.py')) if p.is_dir() else [p]
        for filename in files:
            if filename.suffix=='.py':
                try:
                    relative=filename.resolve().relative_to(Path(get_sympy_dir()).parent)
                    parts=list(relative.with_suffix('').parts)
                    if parts[-1]=='__init__':parts.pop()
                    module=importlib.import_module('.'.join(parts))
                except (ValueError,ImportError):
                    from importlib import util as importlib_util
                    spec=importlib_util.spec_from_file_location('_symox_doctest_'+str(abs(hash(str(filename)))),filename);module=importlib_util.module_from_spec(spec);spec.loader.exec_module(module)
                failures+=_doctest.testmod(module,verbose=verbose,optionflags=flags).failed
            else:failures+=sympytestfile(str(filename.resolve()),module_relative=False,verbose=verbose,optionflags=flags).failed
    return failures==0


def run_all_tests(test_args=(),test_kwargs=None,doctest_args=(),doctest_kwargs=None,examples_args=(),examples_kwargs=None):
    tests=test(*test_args,**(test_kwargs or {}));docs=doctest(*doctest_args,**(doctest_kwargs or {}))
    return tests and docs


def run_in_subprocess_with_hash_randomization(function,function_args=(),function_kwargs=None,command=sys.executable,module='symox.testing.runtests',force=False):
    if os.environ.get('PYTHONHASHSEED') is not None and not force:return None
    name=function if isinstance(function,str) else function.__name__
    payload=base64.b64encode(pickle.dumps((tuple(function_args),dict(function_kwargs or {})))).decode('ascii')
    code='import base64,pickle,importlib; args,kwargs=pickle.loads(base64.b64decode('+repr(payload)+')); result=getattr(importlib.import_module('+repr(module)+'),'+repr(name)+')(*args,**kwargs);raise SystemExit(0 if result else 1)'
    env=dict(os.environ,PYTHONHASHSEED=str(random.SystemRandom().randrange(2**32)))
    return _subprocess.run([command,'-c',code],env=env).returncode==0


class Reporter:
    def __init__(self):self.events=[];self.passed=self.failed=self.skipped=0
    def write(self,text,**kwargs):print(text,end='')
class PyTestReporter(Reporter):
    def __init__(self,verbose=False,tb='short',colors=True,force_colors=False,split=None):super().__init__();self.verbose,self.tb=verbose,tb
    def start(self,seed=None,msg='test process starts'):self.events.append(('start',seed,msg))
    def entering_filename(self,filename,n):self.events.append(('file',filename,n))
    def entering_test(self,function):self.events.append(('test',function))
    def leaving_filename(self):pass
    def root_dir(self,dir):self.directory=dir
    def test_pass(self,char='.'):self.passed+=1;self.write(char)
    def test_fail(self,exc_info):self.failed+=1;self.events.append(('failure',exc_info))
    def test_exception(self,exc_info):self.test_fail(exc_info)
    def test_skip(self,v=None):self.skipped+=1;self.events.append(('skip',v))
    def test_xfail(self):self.events.append(('xfail',))
    def test_xpass(self,v):self.events.append(('xpass',v))
    def import_error(self,filename,exc_info):self.events.append(('import_error',filename,exc_info));self.failed+=1
    def doctest_fail(self,name,error_msg):self.events.append(('doctest_failure',name,error_msg));self.failed+=1
    def finish(self):return self.failed==0


class SymPyTests:
    def __init__(self,reporter,kw='',post_mortem=False,seed=None,fast_threshold=None,slow_threshold=None):self.reporter,self.kw,self.files=reporter,kw,[]
    def get_test_files(self,dir,pat='test_*.py'):
        self.files=[str(path) for path in Path(dir).rglob(pat)];return self.files
    def matches(self,x):return self.kw in str(x)
    def test(self,sort=False,timeout=False,slow=False,enhance_asserts=False,fail_on_timeout=False):return test(*sorted(self.files) if sort else self.files,kw=self.kw,subprocess=True)
    def test_file(self,filename,**kwargs):return test(filename,kw=self.kw)
class SymPyDocTests:
    def __init__(self,reporter,normal):self.reporter,self.normal,self.files=reporter,normal,[]
    def get_test_files(self,dir,pat='*.py',init_only=True):
        self.files=[str(path) for path in Path(dir).rglob(pat) if not init_only or path.name=='__init__.py'];return self.files
    def test(self):return doctest(*self.files)
    def test_file(self,filename):return doctest(filename)

__all__=['DependencyError','PyTestReporter','Reporter','Skipped','SymPyDocTestFinder','SymPyDocTestRunner',
 'SymPyDocTests','SymPyOutputChecker','SymPyTestResults','SymPyTests','TimeOutError','convert_to_native_paths',
 'doctest','get_sympy_dir','raise_on_deprecated','run_all_tests','run_in_subprocess_with_hash_randomization',
 'setup_pprint','split_list','sympytestfile','test']

# Imported runner functions are API entry points, not pytest test cases.
test.__test__ = False
doctest.__test__ = False
