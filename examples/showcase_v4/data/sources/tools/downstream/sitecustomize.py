"""Process-local, experimental SymPy import routing for downstream audits.

Never install this module on a user's normal PYTHONPATH. Reference SymPy
metadata stays installed for the control run, but its Python code is forbidden
in the native run. The facade's compatibility version is explicitly recorded.
"""
import atexit
import importlib
import importlib.abc
import importlib.util
import json
import os
from pathlib import Path
import sys
import types

BACKEND = os.environ.get('SYMOX_DOWNSTREAM_BACKEND')
ROUTED = {}
BLOCKED = []


def provenance():
    native = sys.modules.get('symox')
    reference = sys.modules.get('sympy')
    return dict(backend=BACKEND, pid=os.getpid(), target_api_version='1.14.0',
                native_version=getattr(native, '__version__', None),
                facade_version=getattr(reference, '__version__', None),
                routed_modules=dict(ROUTED), blocked_reference_execution=BLOCKED,
                loaded_symbolic_modules={name: getattr(module, '__file__', None)
                    for name, module in list(sys.modules.items())
                    if name == 'sympy' or name.startswith('sympy.')})


def save():
    destination = os.environ.get('SYMOX_DOWNSTREAM_PROVENANCE')
    if destination:
        original = Path(destination)
        path = original.with_name(original.stem+'.'+str(os.getpid())+original.suffix)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(provenance(), indent=2, sort_keys=True))


if BACKEND in ('sympy', 'symox'):
    atexit.register(save)
    os.environ.setdefault('MPLBACKEND', 'Agg')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

if BACKEND == 'symox':
    preloaded = [name for name in sys.modules if name == 'sympy' or name.startswith('sympy.')]
    if preloaded:
        BLOCKED.extend('preloaded:'+name for name in preloaded)
        save()
        sys.stderr.write('Reference SymPy was loaded before the native audit hook. Aborting.\n')
        os._exit(86)
    def audit(event, args):
        if event == 'exec':
            filename = getattr(args[0], 'co_filename', '')
            if '/site-packages/sympy/' in filename.replace('\\', '/'):
                BLOCKED.append(filename)
                raise ImportError('Reference SymPy code is forbidden in the native audit: '+filename)
    sys.addaudithook(audit)

    class NativeLoader(importlib.abc.Loader):
        def __init__(self, fullname, native):
            self.fullname, self.native = fullname, native

        def create_module(self, spec):
            if self.fullname in ('sympy', 'sympy.release'):
                facade = types.ModuleType(self.fullname)
                facade.__file__ = self.native.__file__
                facade.__path__ = list(getattr(self.native, '__path__', []))
                # Version gates describe the emulated namespace, not proof of
                # API parity. Keep the actual Symox version visible separately.
                facade.__version__ = '1.14.0'
                facade.__symox_version__ = importlib.import_module('symox').__version__
                facade.__getattr__ = lambda name: getattr(self.native, name)
                facade.__all__ = getattr(self.native, '__all__',
                    [name for name in vars(self.native) if not name.startswith('_')])
                # Like SymPy's package namespace, the facade's __dict__ holds
                # the public names (``from sympy import __dict__`` callers).
                for name, value in vars(self.native).items():
                    if not name.startswith('_'):
                        facade.__dict__.setdefault(name, value)
                return facade
            return self.native

        def exec_module(self, module):
            if self.fullname == 'sympy':
                # The native root has already initialized some child packages
                # and restored homonymous public callables (e.g. simplify).
                # Mirror those completed imports before routing later imports;
                # otherwise importlib rebinds the facade's callable to a module.
                for name, child in list(sys.modules.items()):
                    if name.startswith('symox.') and name != 'symox.release' and child is not None:
                        alias = 'sympy'+name[5:]
                        sys.modules.setdefault(alias, child)
                        ROUTED[alias] = name
            ROUTED[self.fullname] = self.native.__name__
            save()

    class NativeFinder(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname != 'sympy' and not fullname.startswith('sympy.'):
                return None
            native_name = 'symox'+fullname[5:]
            # A missing native path must fail; never fall through to SymPy.
            try:
                native = importlib.import_module(native_name)
            except ModuleNotFoundError as error:
                raise ModuleNotFoundError(f'{fullname} requires missing native module {native_name}',
                                          name=fullname) from error
            return importlib.util.spec_from_loader(fullname, NativeLoader(fullname, native),
                                                   is_package=hasattr(native, '__path__'))
    sys.meta_path.insert(0, NativeFinder())
