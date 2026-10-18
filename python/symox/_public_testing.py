"""Opt-in test runners; importing symox never imports pytest."""
from pathlib import Path


def test(*paths,**kwargs):
    import pytest
    arguments = [str(path) for path in paths]
    if not arguments:
        location = Path(__file__).resolve().parents[2]/'tests'
        if not location.is_dir():raise ValueError('test paths must be supplied outside a source checkout')
        arguments.append(str(location))
    if kwargs.pop('verbose',False):arguments.append('-v')
    if 'kw' in kwargs:arguments.extend(['-k',kwargs.pop('kw')])
    if kwargs.pop('tb',None):arguments.append('--tb=short')
    if kwargs:raise TypeError('unknown test options: '+', '.join(kwargs))
    return pytest.main(arguments) == 0


def doctest(*paths,**kwargs):
    import doctest as runner,importlib
    modules = []
    if not paths:
        import symox
        modules = [symox]
    failed = 0
    for path in paths:
        if isinstance(path,str) and Path(path).is_file():
            result = runner.testfile(path,module_relative=False,**kwargs)
            failed += result.failed
        else:modules.append(importlib.import_module(path) if isinstance(path,str) else path)
    for module in modules:failed += runner.testmod(module,**kwargs).failed
    return failed == 0
