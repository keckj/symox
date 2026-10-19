# Featured-project integration audit

Reference: **SymPy 1.14.0**. **The candidate does not pass the featured-project replacement gate.**

These results measure the stated upstream test selections, not universal compatibility or mathematical parity.
Each native run requires a clean reference run on the same unmodified upstream selection and environment.
A passing reference case that becomes skipped, expected to fail, or absent prevents native success.

The experimental process-local importer routes `sympy.*` to native `symox.*` and blocks execution of installed reference SymPy Python files.
Its facade advertises the **target namespace version** 1.14.0 for dependency version gates; the actual Symox version and wheel digest are recorded separately.
This audit shim is not a supported production replacement installer.

Sources come from the [official featured list](https://www.sympy.org/en/index.html). The yt unit-system component, yt unit tests and Devito symbolic tests have separate rows.
Commands, source revisions, wheel digests, dependency freezes, test reports and process provenance are in [the audit artifacts](../audits/downstream/summary.json).
Recorded assertions are checked against [the pinned upstream sources](../audits/downstream/source-integrity.json).
Local regression checks, their exact scope and the production fix commit are recorded in [the validation evidence](../audits/downstream/validation.json).

| Project | Version | SymPy 1.14 reference | Native Symox | Outcome |
|---|---|---|---|---|
| [Cadabra](../audits/downstream/results/cadabra.json) | 2.6.0.dev | 1 pass | 0 pass, 1 error | Native incompatibility |
| [ChemPy](../audits/downstream/results/chempy.json) | 0.10.2 | 573 pass, 28 skip, 5 xpass | 555 pass, 18 fail, 28 skip, 5 xpass | Native incompatibility |
| [Devito](../audits/downstream/results/devito.json) | 4.8.23 | 3719 pass, 4 fail, 5 skip, 4 xfail, 1 xpass | not run/completed | Reference fails; native not run |
| [Devito (symbolic algorithms)](../audits/downstream/results/devito_symbolic.json) | 4.8.23 | 136 pass | startup failure (exit 4) | Native incompatibility |
| [EinsteinPy](../audits/downstream/results/einsteinpy.json) | 0.4.0 | 287 pass, 9 xfail | 0 pass, 14 collection error | Native incompatibility |
| [galgebra](../audits/downstream/results/galgebra.json) | 0.6.0 | 115 pass, 2 fail, 2 skip | not run/completed | Reference fails; native not run |
| [Haarpy](../audits/downstream/results/haarpy.json) | 0.1.1 | 1107 pass | 729 pass, 378 fail | Native incompatibility |
| [LaTeX Expression](../audits/downstream/results/latexexpr.json) | 0.4.1 | 2 pass | 2 pass | No symbolic integration exercised |
| [Lcapy](../audits/downstream/results/lcapy.json) | 1.26 | 317 pass | 0 pass, 37 collection error | Native incompatibility |
| [NRPyLaTeX](../audits/downstream/results/nrpylatex.json) | 1.4.0 | 32 pass | 4 pass, 28 fail | Native incompatibility |
| [OctSymPy](../audits/downstream/results/octsympy.json) | 3.2.2+ | 131 pass | 0 pass, 1 error | Native incompatibility |
| [Optlang](../audits/downstream/results/optlang.json) | 1.9.1 | 346 pass, 49 skip | 0 pass, 17 collection error | Native incompatibility |
| [PyDy](../audits/downstream/results/pydy.json) | 0.9.4 | 76 pass | 17 pass, 13 fail, 46 error | Native incompatibility |
| [pyneqsys](../audits/downstream/results/pyneqsys.json) | 0.5.7 | 19 pass, 9 skip | 19 pass, 9 skip | Selected tests pass |
| [pyodesys](../audits/downstream/results/pyodesys.json) | 0.14.7 | 45 pass, 238 skip | 45 pass, 238 skip | Selected tests pass |
| [QMCPACK](../audits/downstream/results/qmcpack.json) | 4.4.9.dev | 1 pass | 0 pass, 1 error | Native incompatibility |
| [Quantum Programming in Python](../audits/downstream/results/quantum_programming.json) | — | not run/completed | not run/completed | Publication, not a package |
| [SageMath](../audits/downstream/results/sagemath.json) | 10.9 | 319 pass | 127 pass, 192 fail | Native incompatibility |
| [Scikit-fdiff](../audits/downstream/results/scikit_fdiff.json) | 0.7.0 | not run/completed | not run/completed | Requires an older SymPy |
| [SfePy](../audits/downstream/results/sfepy.json) | 2026.3 | 221 pass | 220 pass, 1 fail | Native incompatibility |
| [Spyder](../audits/downstream/results/spyder.json) | 6.1.7 | 1 pass | 1 pass | Selected tests pass |
| [Symbolic statistical modeling](../audits/downstream/results/symbolic_statistics.json) | — | not run/completed | not run/completed | Publication, not a package |
| [symjit](../audits/downstream/results/symjit.json) | 2.27.0 | 52 pass, 1 fail, 6 skip, 9 xfail | not run/completed | Reference fails; native not run |
| [PyTorch TorchInductor](../audits/downstream/results/torchinductor.json) | 2.14.1+cu130 | 236 pass | 0 pass, 1 collection error | Native incompatibility |
| [yt / unyt](../audits/downstream/results/yt.json) | 3.1.0 | 819 pass, 34 skip, 2 xfail | 0 pass, 12 collection error | Native incompatibility |
| [yt (main package)](../audits/downstream/results/yt_main.json) | 4.4.3 | 1 pass | 0 pass, 1 collection error | Native incompatibility |

## Scope and first blockers

### Cadabra

One unmodified upstream tests/sympy_cdb.cdb scalar bridge script, including all assertions (12 printed success markers plus the final unmarked assertion); not the entire tensor/C++ suite

```text
Traceback (most recent call last):
  File "<string>", line 15, in <module>
  File "<string>", line 9, in test01
TypeError: cannot convert object of type ExNode to a symox expression
Script ended execution due to an uncaught exception.
3x + (x)**2 + a (x)**2 + b x

```

### ChemPy

Upstream test tree selected in the pinned manifest.

```text
("sympy")
    def test_binary_rev():
        f = binary_rev(t, kf, kb, prod, major, minor, backend=sympy)
        dfdt = f.diff(t)
        num_dfdt = dfdt.subs(subsd)
        ans = kf * (minor - f) * (major - f) - kb * f
        # symbolic susbsitution fails:
>       assert abs(float(num_dfdt) - float(ans.subs(subsd))) < 2e-14
E       assert 7.105427357600714e-14 < 2e-14
E        +  where 7.105427357600714e-14 = abs((2.8721203367240307e-27 - 7.105427357601002e-14))
E        +    where 2.8721203367240307e-27 = float(sqrt(1093)*(-79 + sqrt(1093))*(sqrt(1093)*(-79 - sqrt(1093)) + (-79 + sqrt(1093))*(sqrt(1093) + 79)*exp(-2*sqrt(1093))...)*(-79 + sqrt(1093))*(sqrt(1093) + 79)*exp(-2*sqrt(1093))/(6*((-79 + sqrt(1093))*exp(-2*sqrt(1093)) + sqrt(1093) + 79)))
E        +    and   7.105427357601002e-14 = float(-7*(sqrt(1093)*(-79 - sqrt(1093)) + (-79 + sqrt(1093))*(sqrt(1093) + 79)*exp(-2*sqrt(1093)) + 79*sqrt(1093) + 6241)/(6...3) + 79)*exp(-2*sqrt(1093)) + 79*sqrt(1093) + 6241)/(6*((-79 + sqrt(1093))*exp(-2*sqrt(1093)) + sqrt(1093) + 79)) + 13))
E        +      where -7*(sqrt(1093)*(-79 - sqrt(1093)) + (-79 + sqrt(1093))*(sqrt(1093) + 79)*exp(-2*sqrt(1093)) + 79*sqrt(1093) + 6241)/(6...3) + 79)*exp(-2*sqrt(1093)) + 79*sqrt(1093) + 6241)/(6*((-79 + sqrt(1093))*exp(-2*sqrt(1093)) + sqrt(1093) + 79)) + 13) = subs({t: 2, kf: 3, kb: 7, major: 11, ...})
E        +        where subs = -kb*(sqrt(-4*kf*(kf*major*minor + kf*major*prod + kf*minor*prod + kf*prod**2) + (-kb - kf*major - kf*minor - 2*kf*prod...t(-4*kf*(kf*major*minor + kf*major*prod + kf*minor*prod + kf*prod**2) + (-kb - kf*major - kf*minor - 2*kf*prod)**2))))).subs

chempy/kinetics/tests/test_integrated.py:57: AssertionError
```

### Devito

Full published tests with four pytest workers, including MPI subprocess tests; original upstream maxfail policy retained

```text
item = <Function test_adjoint_F[(4, basic)-2]>

    @pytest.hookimpl(tryfirst=True, hookwrapper=True)
    def pytest_runtest_call(item):
        inside_pytest_marker = os.environ.get('DEVITO_PYTEST_FLAG', 0)
        with suppress(ValueError):
            inside_pytest_marker = int(inside_pytest_marker)

        if inside_pytest_marker:
            outcome = yield

        elif item.get_closest_marker("parallel"):
            # Spawn parallel processes to run test

            outcome = parallel(item, item.funcargs['mode'])
            if outcome:
                pytest.skip(f"{item} success in parallel")
            else:
>               pytest.fail(f"{item} failed in parallel")
E               Failed: <Function test_adjoint_F[(4, basic)-2]> failed in parallel

conftest.py:280: Failed
```

The four MPI cases that failed in the full parallel reference run were rerun without pytest workers: 4 passed. The full-suite reference remains unqualified; this is not evidence of a Symox regression.

### Devito (symbolic algorithms)

Complete upstream tests/test_symbolics.py, audited separately from the full suite with its four unstable parallel MPI cases

```text
ImportError while loading conftest '/tmp/symox-featured/sources/devito/devito-4.8.23/conftest.py'.
conftest.py:10: in <module>
    from devito import Eq, Revolver, configuration  # noqa
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
devito/__init__.py:11: in <module>
    from devito.types import NODE, CELL, Buffer  # noqa
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
devito/types/__init__.py:4: in <module>
    from .basic import *  # noqa
    ^^^^^^^^^^^^^^^^^^^^
devito/types/basic.py:11: in <module>
    from sympy.core.assumptions import _assume_rules
E   ImportError: cannot import name '_assume_rules' from 'symox.core.assumptions' (/tmp/symox-featured/envs/devito_symbolic/lib/python3.12/site-packages/symox/core/assumptions.py)

```

### EinsteinPy

Upstream test tree selected in the pinned manifest.

```text
tests/test_symbolic/test_christoffel.py:5: in <module>
    from einsteinpy.symbolic import ChristoffelSymbols, MetricTensor
src/einsteinpy/symbolic/__init__.py:2: in <module>
    from .constants import SymbolicConstant, get_constant
src/einsteinpy/symbolic/constants.py:39: in <module>
    c = SymbolicConstant("c", "Speed Of Light")
src/einsteinpy/symbolic/constants.py:27: in __new__
    instance._descriptive_name = descriptive_name
E   AttributeError: 'symox._core.Expr' object has no attribute '_descriptive_name'
```

### galgebra

Upstream test tree selected in the pinned manifest.

```text
sinstance(mat_rep, mv.Mv):     # Versor input
            if not mat_rep.is_versor:
                raise ValueError(mat_rep, 'is not a versor in Versor input for Lt!\n')
            V = mat_rep
            Vg_invol = V.g_invol()
            Vinv = V.inv()
            outermorphism = ga.lt(lambda x: Vg_invol * x * Vinv)
            self.lt_dict = simplify(outermorphism.lt_dict)
        # ## GSG code ends ###

        # ## GSG code starts ###
        elif isinstance(mat_rep, str):    # (One-letter) string input
            Amat = Symbolic_Matrix(mat_rep, coords=self.Ga.coords, f=f, mode=mode)
            if mode == 'g':
                self.__init__(Amat, ga=self.Ga)
            elif mode in ['s', 'a']:
                self.__init__(self.Ga.g_inv * Amat, ga=self.Ga)
        # ## GSG code ends ###

        elif callable(mat_rep):    # Linear multivector function input
            # Function is tested for linearity before use.
            F = mat_rep
            a = mv.Mv('a', 'vector', ga=self.Ga)
            b = mv.Mv('b', 'vector', ga=self.Ga)
            if F(a + b) != F(a) + F(b):
                raise ValueError('{} is not linear'.format(F))
            self.lt_dict = {}
            for base in self.Ga.basis:
                out = F(mv.Mv(base, ga=self.Ga))
                if not out.is_vector() and not out.is_zero():
                    raise ValueError('{} must return vectors'.format(F))
                self.lt_dict[base] = out.obj
        else:
>           raise TypeError("Unsupported argument type {}".format(type(mat_rep)))
E           TypeError: Unsupported argument type <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>

galgebra/lt.py:369: TypeError
```

### Haarpy

Upstream test tree selected in the pinned manifest.

```text
half_degree = 1

    @pytest.mark.parametrize("half_degree", range(1, 3))
    def test_weingarten_circular_orthogonal_hyperoctahedral_symbolic(half_degree):
        """Symbolic validation of COE Weingarten function against results shown in
        `Matsumoto. Weingarten calculus for matrix ensembles associated with compact symmetric spaces:
        <https://arxiv.org/abs/1301.5401>`_
        """
        if half_degree == 1:
>           for permutation in SymmetricGroup(2 * half_degree).generate():
                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
E           AttributeError: 'PermutationGroup' object has no attribute 'generate'. Did you mean: '_enumerate'?

haarpy/tests/test_circular_ensembles.py:55: AttributeError
```

### LaTeX Expression

Upstream test tree selected in the pinned manifest.

### Lcapy

Upstream test tree selected in the pinned manifest.

```text
eatured/git/lcapy/lcapy/tests/test_approximate.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/home/keckj/Documents/symox/tools/downstream/sitecustomize.py:93: in find_spec
    native = importlib.import_module(native_name)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
/home/keckj/.local/share/uv/python/cpython-3.12.9-linux-x86_64-gnu/lib/python3.12/importlib/__init__.py:90: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
<frozen importlib._bootstrap>:1387: in _gcd_import
    ???
<frozen importlib._bootstrap>:1360: in _find_and_load
    ???
<frozen importlib._bootstrap>:1324: in _find_and_load_unlocked
    ???
E   ModuleNotFoundError: No module named 'symox.core.logic'

The above exception was the direct cause of the following exception:
/home/keckj/.local/share/uv/python/cpython-3.12.9-linux-x86_64-gnu/lib/python3.12/importlib/__init__.py:90: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lcapy/tests/test_approximate.py:1: in <module>
    from lcapy import *
lcapy/__init__.py:18: in <module>
    from .functions import *
lcapy/functions.py:23: in <module>
    from .extrafunctions import Degrees as Degrees1
lcapy/extrafunctions.py:11: in <module>
    from sympy.core.logic import fuzzy_not
/home/keckj/Documents/symox/tools/downstream/sitecustomize.py:95: in find_spec
    raise ModuleNotFoundError(f'{fullname} requires missing native module {native_name}',
E   ModuleNotFoundError: sympy.core.logic requires missing native module symox.core.logic
```

### NRPyLaTeX

Upstream test tree selected in the pinned manifest.

```text
'Constant')):
                            return self._namespace[symbol]
                    else: self._define_tensor(IndexedSymbol(function))
                    return function
                self.scanner.reset(prev_state); self.scanner.lex()
            index = self._indexing_3()
            indexing.extend(index)
            symbol.extend(len(index) * ['U'])
            if self.accept('UNDERSCORE'):
                index, order = self._indexing_4()
                covariant = ';' in index
                for i in (',', ';'):
                    if i in index: index.remove(i)
                indexing.extend(index)
                symbol.extend((len(index) - order) * ['D'])
                if order > 0:
                    sentence = self.scanner.sentence
                    suffix = '_cd' if covariant else '_d'
                    symbol.append(suffix + order * 'D')
                    function = Function('Tensor')(Symbol(''.join(symbol)), *indexing)
                    old_latex = sentence[position:self.scanner.prev_state[0]]
                    new_latex = IndexedSymbol(function).latex_format(function)
                    self.scanner.sentence = sentence.replace(old_latex, new_latex)
                    self.scanner.marker = position
                    self.scanner.reset(sym_state)
                    return self._operator()
        symbol = ''.join(symbol)
        if symbol in self._namespace:
>           if isinstance(self._namespace[symbol], Function('Constant')):
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
E           TypeError: isinstance() arg 2 must be a type, a tuple of types, or a union

nrpylatex/core/parser.py:885: TypeError
```

### OctSymPy

Unmodified Octave BIST blocks in @sym/sym.m, diff.m, int.m and solve.m; not the complete Octave package

```text
Traceback (most recent call last):
  File "<stdin>", line 9, in <module>
ImportError: cannot import name 'Boolean' from 'symox.logic.boolalg' (/tmp/symox-featured/envs/sage-octave/lib/python3.12/site-packages/symox/logic/boolalg.py)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<stdin>", line 11, in octoutput_drv
NameError: name 'ET' is not defined

error: Python exception: ImportError: cannot import name 'Boolean' from 'symox.logic.boolalg' (/tmp/symox-featured/envs/sage-octave/lib/python3.12/site-packages/symox/logic/boolalg.py)
    occurred in python_header import block.
    Try "sympref reset" and repeat your command?
    (consider filing an issue at https://github.com/gnu-octave/symbolic/issues)
error: called from
    pycall_sympy__ at line 192 column 5
    valid_sym_assumptions at line 38 column 5
    assumptions at line 82 column 5
    syms at line 140 column 3
    /home/keckj/Documents/symox/tools/downstream/octave_bist.m at line 5 column 1

```

### Optlang

Upstream test tree selected in the pinned manifest.

```text
/home/keckj/.local/share/uv/python/cpython-3.12.9-linux-x86_64-gnu/lib/python3.12/importlib/__init__.py:90: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
<frozen importlib._bootstrap>:1387: in _gcd_import
    ???
<frozen importlib._bootstrap>:1360: in _find_and_load
    ???
<frozen importlib._bootstrap>:1310: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:488: in _call_with_frames_removed
    ???
<frozen importlib._bootstrap>:1387: in _gcd_import
    ???
<frozen importlib._bootstrap>:1360: in _find_and_load
    ???
<frozen importlib._bootstrap>:1310: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:488: in _call_with_frames_removed
    ???
<frozen importlib._bootstrap>:1387: in _gcd_import
    ???
<frozen importlib._bootstrap>:1360: in _find_and_load
    ???
<frozen importlib._bootstrap>:1331: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:935: in _load_unlocked
    ???
<frozen importlib._bootstrap_external>:999: in exec_module
    ???
<frozen importlib._bootstrap>:488: in _call_with_frames_removed
    ???
src/optlang/__init__.py:19: in <module>
    from optlang.util import list_available_solvers
src/optlang/util.py:27: in <module>
    from optlang import symbolics
src/optlang/symbolics.py:112: in <module>
    Expr = sympy.core.expr.Expr
           ^^^^^^^^^^
/home/keckj/Documents/symox/tools/downstream/sitecustomize.py:76: in <lambda>
    facade.__getattr__ = lambda name: getattr(self.native, name)
                                      ^^^^^^^^^^^^^^^^^^^^^^^^^^
E   AttributeError: module 'symox' has no attribute 'core'
```

### PyDy

Upstream test tree selected in the pinned manifest.

```text
self = <pydy.codegen.tests.test_c_code.TestCMatrixGenerator object at 0x7fed10e59280>

    def setup_method(self):

        self.prefix = 'boogly_bee'

>       sys = multi_mass_spring_damper(6, True, True)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pydy/codegen/tests/test_c_code.py:22:
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
pydy/models.py:113: in multi_mass_spring_damper
    return System(kane)
           ^^^^^^^^^^^^
pydy/system.py:104: in __init__
    self._constants_symbols = self._Kane_constant_symbols()
                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
pydy/system.py:653: in _Kane_constant_symbols
    from_eoms, from_sym_lists = self._Kane_inlist_insyms()
                                ^^^^^^^^^^^^^^^^^^^^^^^^^^
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

self = <pydy.system.System object at 0x7fed10e590a0>

    def _Kane_inlist_insyms(self):
        """TODO temporary."""

>       uaux = self.eom_method._uaux[:]
               ^^^^^^^^^^^^^^^^^^^^^
E       AttributeError: 'KanesMethod' object has no attribute '_uaux'

pydy/system.py:607: AttributeError
```

### pyneqsys

Complete published Python test tree; optional solver backend skips retained identically

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

### pyodesys

Complete published Python test tree; optional pycvodes/pygslodeiv2/pyodeint backends unavailable and explicitly skipped in both runs

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

### QMCPACK

Execute the real src/QMCWaveFunctions/tests/gen_matrix_ops.py reference generator and compare output; not the QMCPACK C++ or Monte Carlo suite

```text
1x1 matrix
Traceback (most recent call last):
  File "/tmp/symox-featured/git/qmcpack/src/QMCWaveFunctions/tests/gen_matrix_ops.py", line 240, in <module>
    test_1x1()
  File "/tmp/symox-featured/git/qmcpack/src/QMCWaveFunctions/tests/gen_matrix_ops.py", line 143, in test_1x1
    m1exp = exp(-m1)
            ^^^^^^^^
  File "/tmp/symox-featured/envs/qmcpack/lib/python3.12/site-packages/symox/__init__.py", line 124, in f
    return _core._apply(name, *args) if evaluate else _core._apply(name, *args, evaluate=False)
           ^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: cannot convert object of type Matrix to a symox expression

```

### Quantum Programming in Python

The featured link is a student thesis PDF, not an identifiable installable package/test distribution

### SageMath

All 319 upstream doctests in sage/interfaces/sympy.py; not the full SageMath test suite

```text
._sympysage_nan
   3 of   4 in sage.interfaces.sympy._sympysage_ninfty
   2 of   5 in sage.interfaces.sympy._sympysage_order
   3 of   4 in sage.interfaces.sympy._sympysage_pi
   4 of   9 in sage.interfaces.sympy._sympysage_piecewise
   3 of   4 in sage.interfaces.sympy._sympysage_pinfty
   5 of   9 in sage.interfaces.sympy._sympysage_polygamma
  10 of  14 in sage.interfaces.sympy._sympysage_polynomial
   4 of   8 in sage.interfaces.sympy._sympysage_polynomial_ring
   2 of   5 in sage.interfaces.sympy._sympysage_pow
   2 of   4 in sage.interfaces.sympy._sympysage_rational
   1 of   3 in sage.interfaces.sympy._sympysage_rational_field
   2 of   4 in sage.interfaces.sympy._sympysage_re
   5 of   9 in sage.interfaces.sympy._sympysage_real_interval
  12 of  15 in sage.interfaces.sympy._sympysage_relational
   2 of   6 in sage.interfaces.sympy._sympysage_rf
   2 of   4 in sage.interfaces.sympy._sympysage_symbol
   2 of   3 in sage.interfaces.sympy._sympysage_true
   3 of   4 in sage.interfaces.sympy._sympysage_uinfty
   2 of   6 in sage.interfaces.sympy._sympysage_ynm
   1 of   3 in sage.interfaces.sympy.check_all
   1 of   3 in sage.interfaces.sympy.check_expression
   2 of   5 in sage.interfaces.sympy.sympy_init
    [319 tests, 192 failures, 2.98s wall]
----------------------------------------------------------------------
python3 -m sage.doctest --warn-long 5.0 --random-seed=0 envs/sage-octave/lib/python3.12/site-packages/sage/interfaces/sympy.py  # 192 doctests failed
----------------------------------------------------------------------
Total time for all tests: 3.1 seconds
    cpu time: 1.7 seconds
    cumulative wall time: 3.0 seconds
Features detected for doctesting:

```

### Scikit-fdiff

Upstream test tree selected in the pinned manifest.

Excluded before native testing: sympy (<1.6).

### SfePy

Full unmodified test tree shipped in the compiled wheel, matched against the pinned source archive; slow gradient tests have a 600-second per-test budget in both backends

```text
y:684: in __init__
    sol = sm.solve(eqs)
          ^^^^^^^^^^^^^
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

f = array([[1.0*u0 - 1.0*u1]], dtype=object), symbols = (), flags = {}
Tuple = <class 'symox._compat.Tuple'>
reduce_inequalities = <function reduce_inequalities at 0x7fac01ace8e0>
sympify = <function make_sympify.<locals>.sympify at 0x7fac0094fec0>
sequences = (<class 'list'>, <class 'tuple'>, <class 'set'>, <class 'frozenset'>, <class 'symox._compat.Tuple'>)
equations = (array([[1.0*u0 - 1.0*u1]], dtype=object),)
prepared = [array([[1.0*u0 - 1.0*u1]], dtype=object)], relational = False
equation = array([[1.0*u0 - 1.0*u1]], dtype=object)
value = array([[1.0*u0 - 1.0*u1]], dtype=object)

    @wraps(native)
    def solve(f, *symbols, **flags):
        from . import Tuple, reduce_inequalities, sympify
        sequences=(list,tuple,set,frozenset,Tuple)
        if hasattr(f,'__next__'):
            f=list(f)
        equations=f if isinstance(f,sequences) else (f,)
        prepared=[];relational=False
        for equation in equations:
            value=sympify(equation) if isinstance(equation,str) else equation
            prepared.append(value)
            if getattr(value,'is_Relational',False) and value.rel_op!='==':
                relational=True
        if relational:
            targets=symbols[0] if len(symbols)==1 else symbols
            return reduce_inequalities(prepared,targets)
>       return native(f,*symbols,**flags)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^
E       TypeError: cannot convert object of type ndarray to a symox expression

envs/sfepy/lib/python3.12/site-packages/symox/_solve_relations.py:25: TypeError
```

### Spyder

Upstream IPython console SymPy integration tests; not the full GUI suite

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

### Symbolic statistical modeling

The featured link is the 2012 Symbolic Statistics with SymPy paper, not a separate package/test distribution

### symjit

Upstream test tree selected in the pinned manifest.

```text
self = <test_special_values.SpecialValues testMethod=test_unary_functions>

    def test_unary_functions(self):
        for name, (sym, ref) in special_unary().items():
            if name == "gamma":
                continue  # see test_gamma_at_zero
>           self.check_unary(name, sym, ref)

python/tests/test_special_values.py:104:
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
python/tests/test_special_values.py:83: in check_unary
    self.assertFalse(bad, f"{name} {options} {kind}: " + describe(bad, [points], got, want))
E   AssertionError: [16, 17] is not false : asinh {} vectorized: f(1.7976931348623157e+308) = inf, expected 710.4758600739439; f(-1.7976931348623157e+308) = -inf, expected -710.4758600739439
```

### PyTorch TorchInductor

Unmodified test/test_sympy_utils.py at the installed Torch commit; symbolic helpers used by TorchInductor, not the full PyTorch suite

```text
ImportError while importing test module '/tmp/symox-featured/git/torch/test/test_sympy_utils.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/usr/lib/python3.10/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
git/torch/test/test_sympy_utils.py:14: in <module>
    from sympy.core.relational import is_ge, is_gt, is_le, is_lt
E   ImportError: cannot import name 'is_ge' from 'symox.core.relational' (/tmp/symox-featured/envs/torchinductor/lib/python3.10/site-packages/symox/core/relational.py)
```

### yt / unyt

Full unyt unit-system test suite, the SymPy-consuming component identified by the featured list

```text
tlib._bootstrap>:488: in _call_with_frames_removed
    ???
<frozen importlib._bootstrap>:1387: in _gcd_import
    ???
<frozen importlib._bootstrap>:1360: in _find_and_load
    ???
<frozen importlib._bootstrap>:1310: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:488: in _call_with_frames_removed
    ???
<frozen importlib._bootstrap>:1387: in _gcd_import
    ???
<frozen importlib._bootstrap>:1360: in _find_and_load
    ???
<frozen importlib._bootstrap>:1331: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:935: in _load_unlocked
    ???
<frozen importlib._bootstrap_external>:999: in exec_module
    ???
<frozen importlib._bootstrap>:488: in _call_with_frames_removed
    ???
unyt/__init__.py:38: in <module>
    from unyt import physical_constants, unit_symbols
unyt/physical_constants.py:22: in <module>
    _add_constants(globals(), registry=_default_unit_registry)
unyt/unit_systems.py:91: in add_constants
    namespace[name] = quan.in_base(unit_system=registry.unit_system)
                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
unyt/array.py:1123: in in_base
    to_units = self.units.get_base_equivalent(unit_system)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
unyt/unit_object.py:636: in get_base_equivalent
    new_units = unit_system[self.dimensions]
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
unyt/unit_systems.py:256: in __getitem__
    units = _get_system_unit_string(key, self.units_map)
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
unyt/unit_systems.py:137: in _get_system_unit_string
    dim = list(factor.free_symbols)[0]
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
E   IndexError: list index out of range
```

### yt (main package)

Upstream yt/units/tests in addition to the separately audited full unyt suite; not the entire yt analysis suite

```text
<frozen importlib._bootstrap_external>:999: in exec_module
    ???
<frozen importlib._bootstrap>:488: in _call_with_frames_removed
    ???
envs/yt_main/lib/python3.12/site-packages/yt/__init__.py:12: in <module>
    import yt.units as units
envs/yt_main/lib/python3.12/site-packages/yt/units/__init__.py:1: in <module>
    from unyt.array import (
envs/yt_main/lib/python3.12/site-packages/unyt/__init__.py:38: in <module>
    from unyt import physical_constants, unit_symbols
envs/yt_main/lib/python3.12/site-packages/unyt/physical_constants.py:22: in <module>
    _add_constants(globals(), registry=_default_unit_registry)
envs/yt_main/lib/python3.12/site-packages/unyt/unit_systems.py:91: in add_constants
    namespace[name] = quan.in_base(unit_system=registry.unit_system)
                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
envs/yt_main/lib/python3.12/site-packages/unyt/array.py:1123: in in_base
    to_units = self.units.get_base_equivalent(unit_system)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
envs/yt_main/lib/python3.12/site-packages/unyt/unit_object.py:636: in get_base_equivalent
    new_units = unit_system[self.dimensions]
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
envs/yt_main/lib/python3.12/site-packages/unyt/unit_systems.py:256: in __getitem__
    units = _get_system_unit_string(key, self.units_map)
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
envs/yt_main/lib/python3.12/site-packages/unyt/unit_systems.py:137: in _get_system_unit_string
    dim = list(factor.free_symbols)[0]
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
E   IndexError: list index out of range
```

## Reproduction

See [the harness instructions](../tools/downstream/README.md) for runtime setup, dependency pins and native adapters.

```bash
python tools/downstream/audit.py --wheel /path/to/symox.whl --cache /tmp/symox-featured --out audits/downstream
# A CI gate for a specific integration:
python tools/downstream/audit.py --wheel /path/to/symox.whl --projects pyneqsys --require-native-pass
```

The ordinary audit exits successfully when evidence collection finishes, including negative results. `--require-native-pass` makes an unsuccessful integration fail the command.
Native runtimes use `tools/downstream/external.py`. Publications and incompatible dependency requirements are never counted as tested native packages.
