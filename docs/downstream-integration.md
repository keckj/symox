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
| [ChemPy](../audits/downstream/results/chempy.json) | 0.10.2 | 573 pass, 28 skip, 5 xpass | 556 pass, 17 fail, 28 skip, 5 xpass | Native incompatibility |
| [Devito](../audits/downstream/results/devito.json) | 4.8.23 | 3719 pass, 4 fail, 5 skip, 4 xfail, 1 xpass | not run/completed | Reference fails; native not run |
| [Devito (symbolic algorithms)](../audits/downstream/results/devito_symbolic.json) | 4.8.23 | 136 pass | startup failure (exit 4) | Native incompatibility |
| [EinsteinPy](../audits/downstream/results/einsteinpy.json) | 0.4.0 | 287 pass, 9 xfail | 0 pass, 14 collection error | Native incompatibility |
| [galgebra](../audits/downstream/results/galgebra.json) | 0.6.0 | 115 pass, 2 fail, 2 skip | not run/completed | Reference fails; native not run |
| [Haarpy](../audits/downstream/results/haarpy.json) | 0.1.1 | 1107 pass | 1107 pass | Selected tests pass |
| [LaTeX Expression](../audits/downstream/results/latexexpr.json) | 0.4.1 | 2 pass | 2 pass | No symbolic integration exercised |
| [Lcapy](../audits/downstream/results/lcapy.json) | 1.26 | 317 pass | 178 pass, 90 fail | Native incompatibility |
| [NRPyLaTeX](../audits/downstream/results/nrpylatex.json) | 1.4.0 | 32 pass | 4 pass, 28 fail | Native incompatibility |
| [OctSymPy](../audits/downstream/results/octsympy.json) | 3.2.2+ | 131 pass | 0 pass, 1 error | Native incompatibility |
| [Optlang](../audits/downstream/results/optlang.json) | 1.9.1 | 346 pass, 49 skip | 49 pass, 312 fail, 34 skip | Native incompatibility |
| [PyDy](../audits/downstream/results/pydy.json) | 0.9.4 | 76 pass | 17 pass, 13 fail, 46 error | Native incompatibility |
| [pyneqsys](../audits/downstream/results/pyneqsys.json) | 0.5.7 | 19 pass, 9 skip | 19 pass, 9 skip | Selected tests pass |
| [pyodesys](../audits/downstream/results/pyodesys.json) | 0.14.7 | 45 pass, 238 skip | 45 pass, 238 skip | Selected tests pass |
| [QMCPACK](../audits/downstream/results/qmcpack.json) | 4.4.9.dev | 1 pass | 0 pass, 1 error | Native incompatibility |
| [Quantum Programming in Python](../audits/downstream/results/quantum_programming.json) | — | not run/completed | not run/completed | Publication, not a package |
| [SageMath](../audits/downstream/results/sagemath.json) | 10.9 | 319 pass | 137 pass, 182 fail | Native incompatibility |
| [Scikit-fdiff](../audits/downstream/results/scikit_fdiff.json) | 0.7.0 | not run/completed | not run/completed | Requires an older SymPy |
| [SfePy](../audits/downstream/results/sfepy.json) | 2026.3 | 221 pass | 221 pass | Selected tests pass |
| [Spyder](../audits/downstream/results/spyder.json) | 6.1.7 | 1 pass | 1 pass | Selected tests pass |
| [Symbolic statistical modeling](../audits/downstream/results/symbolic_statistics.json) | — | not run/completed | not run/completed | Publication, not a package |
| [symjit](../audits/downstream/results/symjit.json) | 2.27.0 | 52 pass, 1 fail, 6 skip, 9 xfail | not run/completed | Reference fails; native not run |
| [PyTorch TorchInductor](../audits/downstream/results/torchinductor.json) | 2.14.1+cu130 | 236 pass | 0 pass, 1 collection error | Native incompatibility |
| [yt / unyt](../audits/downstream/results/yt.json) | 3.1.0 | 819 pass, 34 skip, 2 xfail | 819 pass, 34 skip, 2 xfail | Selected tests pass |
| [yt (main package)](../audits/downstream/results/yt_main.json) | 4.4.3 | 1 pass | 1 pass | Selected tests pass |

## Scope and first blockers

### Cadabra

One unmodified upstream tests/sympy_cdb.cdb scalar bridge script, including all assertions (12 printed success markers plus the final unmarked assertion); not the entire tensor/C++ suite

```text
Traceback (most recent call last):
  File "<string>", line 15, in <module>
  File "<string>", line 9, in test01
RuntimeError: SympyBridge::convert: do not know (yet) how to handle <class 'symox._core.Expr'>
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

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

### LaTeX Expression

Upstream test tree selected in the pinned manifest.

### Lcapy

Upstream test tree selected in the pinned manifest.

```text
mp = new._decompose_timedomain_expr(value)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lcapy/superposition.py:453: in _decompose_timedomain_expr
    sval = expr1.laplace()
           ^^^^^^^^^^^^^^^
lcapy/texpr.py:129: in laplace
    return self.LT(evaluate=evaluate, zero_initial_conditions=zero_initial_conditions,
lcapy/texpr.py:115: in LT
    result = laplace_transform(
lcapy/laplace.py:373: in laplace_transform
    return laplace_transformer.transform(expr, t, s,
lcapy/transformer.py:82: in transform
    return self.doit(expr, var, conjvar, evaluate, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lcapy/transformer.py:196: in doit
    ret = self.term(sterm, var, conjvar, **kwargs)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lcapy/laplace.py:363: in term
    return self.integrate_0(expr, t, s) * const
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lcapy/laplace.py:85: in integrate_0
    return self.integrate(expr, t, s, 0, sym.oo)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lcapy/laplace.py:65: in integrate
    self.error('Integration failed')
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

self = <lcapy.laplace.LaplaceTransformer object at 0x7f007874e3c0>
message = 'Integration failed'

    def error(self, message=''):
        if message == '':
            raise ValueError('Could not compute %s for %s' %
                             (self.name, self.expr))
>       raise ValueError('Could not compute %s for %s: %s' %
                         (self.name, self.expr, message))
E       ValueError: Could not compute Laplace transform for DiracDelta(t): Integration failed

lcapy/transformer.py:34: ValueError
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
self = <optlang.tests.test_container.ContainerTestCase testMethod=test_change_object_name>

    def test_change_object_name(self):
        var = Variable('blub')
>       self.model.add(var)

src/optlang/tests/test_container.py:190:
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

self = <optlang.interface.Model object at 0x7f3f135cfe30>, stuff = _blub
sloppy = False

    def add(self, stuff, sloppy=False):
        """Add variables and constraints.

        Parameters
        ----------
        stuff : iterable, Variable, Constraint
            Either an iterable containing variables and constraints or a single variable or constraint.

        sloppy : bool
            Check constraints for variables that are not part of the model yet.

        Returns
        -------
        None
        """
        if self._pending_modifications.toggle == 'remove':
            self.update()
            self._pending_modifications.toggle = 'add'
        if isinstance(stuff, Iterable):
            for elem in stuff:
                self.add(elem, sloppy=sloppy)
        elif isinstance(stuff, Variable):
            if stuff.__module__ != self.__module__:
>               raise TypeError("Cannot add Variable %s of interface type %s to model of type %s." % (
                    stuff, stuff.__module__, self.__module__))
E               TypeError: Cannot add Variable _blub of interface type symox._core to model of type optlang.interface.

src/optlang/interface.py:1399: TypeError
```

### PyDy

Upstream test tree selected in the pinned manifest.

```text
self = <pydy.codegen.tests.test_c_code.TestCMatrixGenerator object at 0x7fbea7764fb0>

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

self = <pydy.system.System object at 0x7fbeb834b020>

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
Matrix([[1.00000000000000]])

2x2 matrix
Traceback (most recent call last):
  File "/tmp/symox-featured/git/qmcpack/src/QMCWaveFunctions/tests/gen_matrix_ops.py", line 241, in <module>
    test_2x2()
  File "/tmp/symox-featured/git/qmcpack/src/QMCWaveFunctions/tests/gen_matrix_ops.py", line 152, in test_2x2
    m2exp = exp(m2)
            ^^^^^^^
  File "/tmp/symox-featured/envs/qmcpack/lib/python3.12/site-packages/symox/__init__.py", line 98, in __call__
    return self._factory(*args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/tmp/symox-featured/envs/qmcpack/lib/python3.12/site-packages/symox/__init__.py", line 151, in f
    return args[0].exp().as_immutable()
           ^^^^^^^^^^^^^
NotImplementedError: Cannot construct an exact analytic matrix function: numerical eigenvalues of float matrices are not implemented

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
   3 of   9 in sage.interfaces.sympy._sympysage_piecewise
   3 of   4 in sage.interfaces.sympy._sympysage_pinfty
   5 of   9 in sage.interfaces.sympy._sympysage_polygamma
   7 of  14 in sage.interfaces.sympy._sympysage_polynomial
   4 of   8 in sage.interfaces.sympy._sympysage_polynomial_ring
   2 of   5 in sage.interfaces.sympy._sympysage_pow
   2 of   4 in sage.interfaces.sympy._sympysage_rational
   1 of   3 in sage.interfaces.sympy._sympysage_rational_field
   2 of   4 in sage.interfaces.sympy._sympysage_re
   4 of   9 in sage.interfaces.sympy._sympysage_real_interval
  12 of  15 in sage.interfaces.sympy._sympysage_relational
   2 of   6 in sage.interfaces.sympy._sympysage_rf
   2 of   4 in sage.interfaces.sympy._sympysage_symbol
   2 of   3 in sage.interfaces.sympy._sympysage_true
   3 of   4 in sage.interfaces.sympy._sympysage_uinfty
   2 of   6 in sage.interfaces.sympy._sympysage_ynm
   1 of   3 in sage.interfaces.sympy.check_all
   1 of   3 in sage.interfaces.sympy.check_expression
   2 of   5 in sage.interfaces.sympy.sympy_init
    [319 tests, 182 failures, 1.62s wall]
----------------------------------------------------------------------
python3 -m sage.doctest --warn-long 5.0 --random-seed=0 envs/sage-octave/lib/python3.12/site-packages/sage/interfaces/sympy.py  # 182 doctests failed
----------------------------------------------------------------------
Total time for all tests: 1.7 seconds
    cpu time: 1.6 seconds
    cumulative wall time: 1.6 seconds
Features detected for doctesting:

```

### Scikit-fdiff

Upstream test tree selected in the pinned manifest.

Excluded before native testing: sympy (<1.6).

### SfePy

Full unmodified test tree shipped in the compiled wheel, matched against the pinned source archive; slow gradient tests have a 600-second per-test budget in both backends

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

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
j/Documents/symox/.venv/lib/python3.10/site-packages/torch/_inductor/standalone_compile.py:18: in <module>
    from torch._dynamo.aot_compile_types import BundledAOTAutogradSerializableCallable
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/_dynamo/__init__.py:13: in <module>
    from . import (
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/_dynamo/aot_compile.py:17: in <module>
    from torch._dynamo.convert_frame import GraphRuntimeEnv
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/_dynamo/convert_frame.py:62: in <module>
    from torch._dynamo.symbolic_convert import TensorifyState
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/_dynamo/symbolic_convert.py:54: in <module>
    from torch._dynamo.exc import (
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/_dynamo/exc.py:72: in <module>
    from .utils import counters
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/_dynamo/utils.py:68: in <module>
    import torch.fx.experimental.symbolic_shapes
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/fx/experimental/symbolic_shapes.py:80: in <module>
    from torch.utils._sympy.functions import (
/home/keckj/Documents/symox/.venv/lib/python3.10/site-packages/torch/utils/_sympy/functions.py:16: in <module>
    from sympy.core.operations import LatticeOp, ShortCircuit
/home/keckj/Documents/symox/tools/downstream/sitecustomize.py:105: in find_spec
    raise ModuleNotFoundError(f'{fullname} requires missing native module {native_name}',
E   ModuleNotFoundError: sympy.core.operations requires missing native module symox.core.operations
```

### yt / unyt

Full unyt unit-system test suite, the SymPy-consuming component identified by the featured list

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

### yt (main package)

Upstream yt/units/tests in addition to the separately audited full unyt suite; not the entire yt analysis suite

All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.

## Reproduction

See [the harness instructions](../tools/downstream/README.md) for runtime setup, dependency pins and native adapters.

```bash
python tools/downstream/audit.py --wheel /path/to/symox.whl --cache /tmp/symox-featured --out audits/downstream
# A CI gate for a specific integration:
python tools/downstream/audit.py --wheel /path/to/symox.whl --projects pyneqsys --require-native-pass
```

The ordinary audit exits successfully when evidence collection finishes, including negative results. `--require-native-pass` makes an unsuccessful integration fail the command.
Native runtimes use `tools/downstream/external.py`. Publications and incompatible dependency requirements are never counted as tested native packages.
