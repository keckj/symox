"""Ordered quantum operators, states, tensor products and representations."""
from ._core import *
from ._core import __all__

# Load the submodules named like their functions now, as the reference does, so
# a later ``import symox.physics.quantum.qapply`` keeps the function bound here.
import importlib as _importlib  # noqa: E402
for _name in ('qapply','represent'):
    _importlib.import_module(f'{__name__}.{_name}')
from ._core import qapply, represent  # noqa: E402
