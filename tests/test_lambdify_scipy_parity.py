"""lambdify with and without SciPy, compared with SymPy 1.14 in the same
configuration.

SymPy's default modules are ['numpy', 'scipy'] when SciPy is importable and
['numpy'] otherwise; the NumPy printer then binds gamma, erf, ... to the
math module and prints the functions it does not support under their own
names (its generated function raises NameError). Each configuration runs in
a fresh interpreter, SciPy being made unimportable for the 'absent' one, so
the namespaces that SymPy and symox build once are those of the
configuration.

Compared: the generated source, signature and docstring, and the result
type, value or exception (type and message) of the calls. Documented
extensions (docs/algorithmic-parity.md, lambdaprint and lambdarr): float64
arrays and float arguments are evaluated natively where SymPy's generated
function raises (math.gamma on an array, a function printed without a
definition); for those calls only a value is required.
"""
import json
import os
import subprocess
import sys

import pytest

pytest.importorskip('numpy')
pytest.importorskip('sympy')

PROBE = r'''
import inspect, json, sys
if sys.argv[1] == 'absent':
    sys.modules['scipy'] = None
import numpy as np


def outcome(call):
    try:
        value = call()
    except Exception as error:
        return ['raise', type(error).__name__, str(error)]
    kind = type(value).__module__ + '.' + type(value).__name__
    if isinstance(value, np.ndarray):
        return ['value', kind, str(value.dtype), repr(value.tolist())]
    return ['value', kind, repr(value)]


def run(S):
    x, y = S.symbols('x y')
    cases = [
        ('gamma', S.gamma(x)), ('gamma+erf', S.gamma(x) + S.erf(x)), ('erfc', S.erfc(x)),
        ('loggamma', S.loggamma(x)), ('beta', S.beta(x, 2)), ('besselj', S.besselj(1, x)),
        ('zeta', S.zeta(x)), ('Si', S.Si(x)), ('airyai', S.airyai(x)),
        ('polygamma', S.polygamma(0, x)), ('sinc', S.sinc(x)), ('sinc2', S.sinc(2*x + 1)),
        ('sin', S.sin(x)), ('exp*gamma', S.exp(x)*S.gamma(x)), ('gamma**2', S.gamma(x)**2),
        ('besselj+sin', S.besselj(0, x) + S.sin(x)),
    ]
    args = [('float', lambda: 0.5), ('complex', lambda: 2.5 + 1j),
            ('array', lambda: np.array([0.5, 1.5])), ('carray', lambda: np.array([2.5 + 1j, 0.5 - 0.25j]))]
    out = {}
    for name, expr in cases:
        for modules in (None, 'numpy', 'scipy'):
            key = '%s|%s' % (name, modules)
            try:
                f = S.lambdify(x, expr, modules)
            except Exception as error:
                out[key] = {'lambdify': [type(error).__name__, str(error)]}
                continue
            record = {}
            try:
                record['source'] = inspect.getsource(f)
            except Exception as error:
                record['source'] = 'raise ' + type(error).__name__
            record['signature'] = str(inspect.signature(f))
            head, _, imports = (f.__doc__ or '').partition('Imported modules:\n\n')
            # SymPy lists the imports of one module in set order
            record['doc'] = [head, sorted(imports.split('\n'))]
            for label, arg in args:
                # twice: the first call may cache a dispatch decision
                record[label] = [outcome(lambda: f(arg())), outcome(lambda: f(arg()))]
            out[key] = record
    return out


if sys.argv[2] == 'sympy':
    import sympy as S
else:
    import symox as S
print(json.dumps(run(S)))
'''

# The calls where symox evaluates natively what SymPy's generated function
# rejects (the documented extensions): float arguments and float64 arrays.
EXTENSIONS = ('float', 'array')

# Remaining gaps recorded in docs/algorithmic-parity.md ("lambdify without
# SciPy"), by configuration: the forms SymPy rewrites (beta to a gamma ratio
# without scipy.special.beta) and SciPy's indexed forms (sici(x)[0],
# airy(x)[0]) are not reproduced, so there is no generated source and
# complex calls are not SymPy's; with SciPy the docstring of sinc lists
# 'from numpy import pi' where SymPy lists 'from scipy.constants import pi'.
GAPS = {
    'absent': {'beta|None': 'all', 'beta|numpy': 'all'},
    'installed': {'beta|numpy': 'all', 'Si|None': 'all', 'Si|scipy': 'all',
                  'airyai|None': 'all', 'airyai|scipy': 'all',
                  'sinc|scipy': 'doc', 'sinc2|scipy': 'doc'},
}


def _probe(configuration, library):
    done = subprocess.run([sys.executable, '-c', PROBE, configuration, library],
                          capture_output=True, text=True, env=dict(os.environ), timeout=300)
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


def _close(want, got):
    """Same outcome kind, type (and dtype), exception message; values equal
    to 1e-13 relative (native special-function kernels and SciPy's differ
    in the last digits; gamma is checked bit for bit in test_native_lambdarr)."""
    if want[0] != got[0] or want[1] != got[1]:
        return False
    if want[0] == 'raise':
        return want[2] == got[2]
    if want[1] == 'numpy.ndarray' and want[2] != got[2]:
        return False
    return _same_values(want[-1], got[-1])


def _same_values(want, got):
    import numpy as np
    a, b = (np.array(eval(v, {'np': np, 'inf': np.inf, 'nan': np.nan})) for v in (want, got))
    return a.shape == b.shape and bool(np.allclose(a, b, rtol=1e-13, atol=0, equal_nan=True))


@pytest.fixture(scope='module', params=['absent', 'installed'])
def configuration(request):
    if request.param == 'installed':
        pytest.importorskip('scipy')
    return request.param


@pytest.fixture(scope='module')
def probes(configuration):
    return configuration, _probe(configuration, 'sympy'), _probe(configuration, 'symox')


def test_lambdify_errors(probes):
    _, sympy, symox = probes
    assert sympy.keys() == symox.keys()
    for key in sympy:
        if 'lambdify' in sympy[key] or 'lambdify' in symox[key]:
            # 'scipy' without SciPy: SymPy's ImportError at lambdify time
            assert sympy[key] == symox[key], key


def test_sources_signatures_docstrings(probes):
    configuration, sympy, symox = probes
    gaps = GAPS[configuration]
    for key, expected in sympy.items():
        if 'lambdify' in expected or gaps.get(key) == 'all':
            continue
        for field in ('source', 'signature', 'doc'):
            if gaps.get(key) == field:
                continue
            assert symox[key][field] == expected[field], (key, field)


def test_call_results_and_errors(probes):
    configuration, sympy, symox = probes
    gaps = GAPS[configuration]
    for key, expected in sympy.items():
        if 'lambdify' in expected:
            continue
        for label in ('float', 'complex', 'array', 'carray'):
            for want, got in zip(expected[label], symox[key][label]):
                if want[0] == 'raise' and label in EXTENSIONS:
                    # documented extension: a native value instead of the error
                    assert got[0] == 'value' or got == want, (key, label, want, got)
                elif gaps.get(key) == 'all':
                    # recorded gap: a value is still SymPy's value
                    if want[0] == got[0] == 'value':
                        assert _same_values(want[-1], got[-1]), (key, label)
                else:
                    assert _close(want, got), (key, label, want, got)


def test_recorded_gaps_are_still_gaps(probes):
    # a gap that closes must be removed from GAPS (and from the parity notes)
    configuration, sympy, symox = probes
    for key, field in GAPS[configuration].items():
        fields = ('source', 'signature', 'doc') if field == 'all' else (field,)
        assert any(symox[key][f] != sympy[key][f] for f in fields), key
