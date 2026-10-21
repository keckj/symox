"""Submodules initialized by SymPy's root import.

``import sympy`` loads these modules (SymPy 1.14.0), so callers reach them as
attributes (``sympy.core.symbol``, ``sympy.functions.elementary.integers``)
without importing them. The root import of symox loads the corresponding
symox modules that exist; a missing module stays an import error.
"""
import sys
from importlib import import_module
from importlib.util import LazyLoader, find_spec, module_from_spec

_MODULES = """
algebras algebras.quaternion assumptions assumptions.ask
assumptions.ask_generated assumptions.assume assumptions.cnf
assumptions.refine assumptions.relation assumptions.relation.binrel
assumptions.relation.equality calculus calculus.accumulationbounds
calculus.euler calculus.finite_diff calculus.singularities calculus.util
concrete concrete.expr_with_intlimits concrete.expr_with_limits
concrete.gosper concrete.products concrete.summations core
core._print_helpers core.add core.alphabets core.assumptions
core.assumptions_generated core.basic core.cache core.containers core.core
core.coreerrors core.decorators core.evalf core.expr core.exprtools
core.facts core.function core.intfunc core.kind core.logic core.mod core.mul
core.multidimensional core.numbers core.operations core.parameters
core.power core.random core.relational core.rules core.singleton
core.sorting core.symbol core.sympify core.traversal discrete
discrete.convolutions discrete.transforms external external.gmpy
external.importtools external.ntheory external.pythonmpq functions
functions.combinatorial functions.combinatorial.factorials
functions.combinatorial.numbers functions.elementary
functions.elementary._trigonometric_special functions.elementary.complexes
functions.elementary.exponential functions.elementary.hyperbolic
functions.elementary.integers functions.elementary.miscellaneous
functions.elementary.piecewise functions.elementary.trigonometric
functions.special functions.special.bessel functions.special.beta_functions
functions.special.bsplines functions.special.delta_functions
functions.special.elliptic_integrals functions.special.error_functions
functions.special.gamma_functions functions.special.hyper
functions.special.mathieu_functions functions.special.polynomials
functions.special.singularity_functions
functions.special.spherical_harmonics functions.special.tensor_functions
functions.special.zeta_functions geometry geometry.curve geometry.ellipse
geometry.entity geometry.exceptions geometry.line geometry.parabola
geometry.plane geometry.point geometry.polygon geometry.util integrals
integrals.deltafunctions integrals.integrals integrals.laplace
integrals.meijerint integrals.rationaltools integrals.singularityfunctions
integrals.transforms integrals.trigonometry interactive interactive.printing
interactive.session interactive.traversal logic logic.boolalg
logic.inference matrices matrices.decompositions matrices.dense
matrices.determinant matrices.eigen matrices.exceptions matrices.expressions
matrices.expressions._shape matrices.expressions.adjoint
matrices.expressions.blockmatrix matrices.expressions.companion
matrices.expressions.determinant matrices.expressions.diagonal
matrices.expressions.dotproduct matrices.expressions.funcmatrix
matrices.expressions.hadamard matrices.expressions.inverse
matrices.expressions.kronecker matrices.expressions.matadd
matrices.expressions.matexpr matrices.expressions.matmul
matrices.expressions.matpow matrices.expressions.permutation
matrices.expressions.sets matrices.expressions.slice
matrices.expressions.special matrices.expressions.trace
matrices.expressions.transpose matrices.graph matrices.immutable
matrices.inverse matrices.kind matrices.matrixbase matrices.reductions
matrices.repmatrix matrices.solvers matrices.sparse matrices.sparsetools
matrices.subspaces matrices.utilities multipledispatch
multipledispatch.conflict multipledispatch.core multipledispatch.dispatcher
multipledispatch.utils ntheory ntheory.continued_fraction ntheory.digits
ntheory.ecm ntheory.egyptian_fraction ntheory.factor_ ntheory.generate
ntheory.modular ntheory.multinomial ntheory.partitions_ ntheory.primetest
ntheory.qs ntheory.residue_ntheory parsing parsing.sympy_parser plotting
plotting.backends plotting.backends.base_backend
plotting.backends.matplotlibbackend
plotting.backends.matplotlibbackend.matplotlib plotting.backends.textbackend
plotting.backends.textbackend.text plotting.intervalmath
plotting.intervalmath.interval_arithmetic
plotting.intervalmath.interval_membership plotting.intervalmath.lib_interval
plotting.plot plotting.plot_implicit plotting.plotgrid plotting.pygletplot
plotting.series plotting.textplot plotting.utils polys polys.agca
polys.agca.extensions polys.agca.homomorphisms polys.agca.ideals
polys.agca.modules polys.appellseqs polys.compatibility polys.constructor
polys.densearith polys.densebasic polys.densetools polys.domainmatrix
polys.domains polys.domains.algebraicfield polys.domains.characteristiczero
polys.domains.complexfield polys.domains.compositedomain
polys.domains.domain polys.domains.domainelement
polys.domains.expressiondomain polys.domains.expressionrawdomain
polys.domains.field polys.domains.finitefield polys.domains.fractionfield
polys.domains.gaussiandomains polys.domains.gmpyfinitefield
polys.domains.gmpyintegerring polys.domains.gmpyrationalfield
polys.domains.groundtypes polys.domains.integerring
polys.domains.modularinteger polys.domains.polynomialring
polys.domains.pythonfinitefield polys.domains.pythonintegerring
polys.domains.pythonrational polys.domains.pythonrationalfield
polys.domains.rationalfield polys.domains.realfield polys.domains.ring
polys.domains.simpledomain polys.euclidtools polys.factortools
polys.fglmtools polys.fields polys.galoistools polys.groebnertools
polys.heuristicgcd polys.matrices polys.matrices._typing polys.matrices.ddm
polys.matrices.dense polys.matrices.dfm polys.matrices.domainmatrix
polys.matrices.domainscalar polys.matrices.eigen polys.matrices.exceptions
polys.matrices.linsolve polys.matrices.lll polys.matrices.normalforms
polys.matrices.rref polys.matrices.sdm polys.monomials polys.numberfields
polys.numberfields.basis polys.numberfields.exceptions
polys.numberfields.galois_resolvents polys.numberfields.galoisgroups
polys.numberfields.minpoly polys.numberfields.modules
polys.numberfields.primes polys.numberfields.subfield
polys.numberfields.utilities polys.orderings polys.orthopolys polys.partfrac
polys.polyclasses polys.polyconfig polys.polyerrors polys.polyfuncs
polys.polyoptions polys.polyquinticconst polys.polyroots polys.polytools
polys.polyutils polys.puiseux polys.rationaltools polys.ring_series
polys.rings polys.rootisolation polys.rootoftools polys.solvers
polys.specialpolys polys.sqfreetools printing printing.codeprinter
printing.conventions printing.defaults printing.dot printing.glsl
printing.gtk printing.jscode printing.julia printing.lambdarepr
printing.latex printing.maple printing.mathematica printing.mathml
printing.numpy printing.octave printing.precedence printing.pretty
printing.pretty.pretty printing.pretty.pretty_symbology
printing.pretty.stringpict printing.preview printing.printer printing.pycode
printing.python printing.rcode printing.repr printing.smtlib printing.str
printing.tableform printing.tree release series series.approximants
series.formal series.fourier series.gruntz series.limits series.limitseq
series.order series.residues series.sequences series.series
series.series_class sets sets.conditionset sets.contains sets.fancysets
sets.handlers sets.handlers.comparison sets.handlers.intersection
sets.handlers.union sets.ordinals sets.powerset sets.sets simplify
simplify.combsimp simplify.cse_main simplify.cse_opts simplify.epathtools
simplify.fu simplify.gammasimp simplify.hyperexpand simplify.powsimp
simplify.radsimp simplify.ratsimp simplify.simplify simplify.sqrtdenest
simplify.trigsimp solvers solvers.bivariate solvers.decompogen
solvers.deutils solvers.diophantine solvers.diophantine.diophantine
solvers.inequalities solvers.ode solvers.ode.ode solvers.ode.hypergeometric
solvers.ode.lie_group solvers.ode.nonhomogeneous
solvers.ode.riccati solvers.ode.single solvers.ode.subscheck
solvers.ode.systems solvers.pde solvers.polysys solvers.recurr
solvers.simplex solvers.solvers solvers.solveset strategies
strategies.branch strategies.branch.core strategies.branch.tools
strategies.branch.traverse strategies.core strategies.rl strategies.tools
strategies.traverse strategies.tree strategies.util tensor tensor.array
tensor.array.array_comprehension tensor.array.arrayop
tensor.array.dense_ndim_array tensor.array.mutable_ndim_array
tensor.array.ndim_array tensor.array.sparse_ndim_array tensor.functions
tensor.index_methods tensor.indexed utilities utilities.decorator
utilities.enumerative utilities.exceptions utilities.iterables
utilities.lambdify utilities.magic utilities.mathml utilities.memoization
utilities.misc utilities.source utilities.timeutils
""".split()


# Modules whose execution needs mpmath. The root import of symox does not
# load mpmath, so these are registered with a deferred execution: the module
# object exists as an attribute of its package and runs on first access.
_DEFERRED = frozenset(['core.evalf'])


def _register_deferred(qualified):
    if qualified in sys.modules:
        return
    spec = find_spec(qualified)
    if spec is None:
        raise ModuleNotFoundError(qualified, name=qualified)
    spec.loader = LazyLoader(spec.loader)
    module = module_from_spec(spec)
    sys.modules[qualified] = module
    spec.loader.exec_module(module)
    parent, _, child = qualified.rpartition('.')
    setattr(sys.modules[parent], child, module)


# Submodules that SymPy exposes as package attributes although the package
# also re-exports a function of the same name (sympy.core.evalf).
_MODULE_ATTRIBUTES = frozenset(['core.evalf'])


def _restore_exports(package):
    # Importing a submodule binds it on its package. A package whose public
    # name of the same spelling is a lazy re-export (symox._api_exports)
    # keeps that object, as SymPy's packages rebind theirs after importing
    # the submodule (sympy.matrices.expressions.trace is the function).
    prefix = package + '.'
    for qualified, module in list(sys.modules.items()):
        if module is None or not qualified.startswith(prefix):
            continue
        parent_name, _, child = qualified.rpartition('.')
        parent = sys.modules.get(parent_name)
        if parent is None or qualified[len(prefix):] in _MODULE_ATTRIBUTES:
            continue
        namespace = vars(parent)
        if child in namespace.get('_api_export_sources', ()) and namespace.get(child) is module:
            del namespace[child]


def install(package):
    for name in _MODULES:
        qualified = package + '.' + name
        try:
            if name in _DEFERRED:
                _register_deferred(qualified)
            else:
                import_module(qualified)
        except ModuleNotFoundError as error:
            # Only the absence of the module itself (or of a parent package)
            # is tolerated; a failing import inside an existing module is not.
            missing = error.name or ''
            if not (qualified == missing or qualified.startswith(missing + '.')):
                raise
    _restore_exports(package)
