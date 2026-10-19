"""Optional extensions must not change the compatible public import surface."""
import subprocess
import sys


def test_extra_source_exports_do_not_import_execution_libraries():
    subprocess.run([sys.executable, '-c', '''
import sys
import symox as s
from symox import extra
from symox.printing.pytorch import TorchPrinter, torch_code
x = s.Symbol('x')
for name in extra.__all__:
    assert not hasattr(s, name), name
    assert callable(getattr(extra, name)), name
assert 'torch.sin' in extra.torchcode(s.sin(x))
assert 'jnp.sin' in extra.jaxcode(s.sin(x))
assert 'sycl::sin' in extra.syclcode(s.sin(x))
assert 'sin' in extra.slangcode(s.sin(x))
assert 'sin' in extra.cudacode(s.sin(x))
assert TorchPrinter().doprint(x) == torch_code(x) == 'x'
assert not {'torch', 'jax', 'slangpy', 'sympy'}.intersection(sys.modules)
'''], check=True, capture_output=True, text=True)
