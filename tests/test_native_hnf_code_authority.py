"""An exact ModuleType and typed code table do not establish semantic origin."""
import json
import os
import subprocess
import sys

from test_native_hnf_function_bootstrap import SCRIPT


INJECTION = r'''
package = importlib.import_module('symox.polys.matrices')
proof = importlib.import_module('symox.polys.matrices._hnf_modular_native')
foreign = types.ModuleType(proof.__name__)
foreign.__dict__.update(proof.__dict__)
foreign._CODES = tuple((name, original.__code__ if name == target else code)
                       for name, code in proof._CODES)
assert type(foreign) is types.ModuleType
sys.modules[proof.__name__] = foreign
package._hnf_modular_native = foreign
'''


def test_foreign_code_table_cannot_authorize_changed_mathematical_helper():
    marker = "sys.modules.pop('symox.polys.matrices._hnf_native')"
    assert SCRIPT.count(marker) == 1
    source = SCRIPT.replace(marker, INJECTION + '\n' + marker)
    # Check this observation before the expected before-fix mismatch too.
    no_reference = "assert not any(name == 'sympy' or name.startswith('sympy.') for name in sys.modules)"
    source = source.replace(no_reference, '')
    source = source.replace('assert actual == expected,', no_reference + '\nassert actual == expected,')
    result = subprocess.run([sys.executable, '-c', source, '_gcdex', 'False'],
                            text=True, capture_output=True, env=os.environ.copy())
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads(result.stdout)
    assert data['events'] == ['gcdex']
    assert data['native_admissions'] == []
    assert data['reference_imported'] is False
