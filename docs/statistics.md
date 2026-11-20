# Symbolic probability and statistics

[README](../README.md)

`from symox.stats import ...` exposes the 157 public names in SymPy 1.14.0's
`stats.__all__`. The distributions and algorithms run locally on native SymOx
expressions, NumPy and mpmath. They do not call SymPy at runtime.

This is API presence, not full SymPy behavioral compatibility. General integrals
and infinite sums can remain unevaluated, and the probability-space objects have
their own Python implementation. The limits below apply even when the coverage
row reads 100%.

## Basic calculations

```python
from symox import Eq, Rational
from symox.stats import Die, Normal, P, E, variance, covariance, sample

X, Y = Die('X', 6), Die('Y', 6)
P(X + Y > 10)                 # 1/12
E(X, X > 3)                  # 5
variance(X)                   # 35/12
covariance(X, Y)               # 0

Z = Normal('Z', 0, 1)
P(Z > 0)                     # 1/2
E(Z, Z > 0)                  # sqrt(2)/sqrt(pi)
E(Z, Eq(Z, 2))               # 2
sample(X + Y, size=(2, 3), seed=42)  # NumPy array
```

Use `Eq(X, value)` for an equality event. Python's `X == value` tests object
identity/equality, as it does for ordinary symbolic expressions.

Independent scalar constructors create separate probability spaces. Repeating a
scalar constructor with the same name and parameters returns the same variable.
A joint model keeps its components dependent. `random_symbols(expr)` returns a
set of the registered random variables.

`P` of `&`/`|` events and `E` of products factor over independent probability
spaces, and condition blocks independent of the query are dropped. A random
variable used as a condition means that it is observed, as in SymPy's
`Eq(Y, Y.symbol)`: with independent `X` and `Y`, `E(X + Y, Y)` is `Y`.
Overlapping `|` events of independent discrete variables use
inclusion-exclusion, where SymPy 1.14 adds the two probabilities.

## Scalar distributions

All public finite, discrete and continuous constructors are available, including
handmade laws. Examples include Bernoulli, Binomial, BetaBinomial,
Hypergeometric, the soliton laws, Poisson, Skellam, Normal, Beta, Gamma, StudentT,
Weibull, inverse Gaussian, noncentral Beta/Chi and the less common positive laws.
Constructor parameters follow the SymPy 1.14 names and defaults.

`density(X)` returns a distribution object with callable `pdf`/`pmf`, `set`,
`cdf`, `moment`, `mgf` and sampling methods. Finite distributions also expose
`dict`. A distribution's raw density formula is intended for its support;
`bounded_pdf(x)` explicitly applies the support restriction.

Finite probabilities, conditional probabilities and expectations use exact
weighted enumeration, limited to one million combined outcomes. Closed raw
moments are supplied for the usual laws; polynomial expectations combine those
moments without integration. Binomial, Normal and Poisson raw moments, and
polynomial expectations and variances over independent stock laws, run in
Rust; other cases keep the Python path. Other moments are represented by an exact native
`Integral` or `Sum`. Noncentral Beta uses its Poisson mixture. Heavy-tail moment
existence conditions matter: an undefined mean can be `nan`, and a divergent
positive moment can be `oo`.

`cdf`, `quantile`, `median`, characteristic functions, moment generating
functions and entropy have closed formulas where available. Finite CDFs are
ordered dictionaries. A symbolic quantile without an implemented inverse raises
`NotImplementedError`; arbitrary inverse-CDF solving is not guaranteed.
Densities and CDFs of the documented laws are written in the reference forms,
including the support `Piecewise` of the CDF; `cdf(X, **options)` forwards the
options to integration, as SymPy does. Some moments, `skewness`, `kurtosis` and
`entropy` come from simplified closed forms and can print differently from
SymPy's unsimplified integrals while being equal.

For a handmade law, supply a normalized density/mass function and its support:

```python
from symox import Symbol, Interval, Range, oo
from symox.stats import ContinuousRV, DiscreteRV, FiniteRV

x = Symbol('x')
k = Symbol('k', integer=True)
U = ContinuousRV(x, 2*x, Interval(0, 1))
D = DiscreteRV(k, 2**(-k), Range(1, oo))
F = FiniteRV('F', {0: Rational(1, 3), 1: Rational(2, 3)})
```

Handmade finite weights are checked by default. Continuous and infinite discrete
normalization is the caller's responsibility. Soliton laws accept a symbolic
support size, whose `density(...).dict` is a held `Density`; `E` and `P` of a
symbolic-size soliton raise `NotImplementedError`. `RobustSoliton` normalizes
over the support `1..k`, so its masses sum to 1; SymPy 1.14's normalizer does
not, and numeric values differ from the reference doctest.

## Conditioning and joint models

Finite conditioning is exact. Single-source continuous event conditioning uses
support intervals and normalized restricted densities; truncated Normal, Beta,
Gamma and Uniform polynomial moments have specialized formulas. Other cases
can remain exact integrals. `given` constructs a conditioned variable or a
transformed expression over that variable, and returns the expression unchanged
when the condition is independent of it. `given(Tuple(...), cond)` conditions
the joint space of dependent entries: finite spaces are filtered and
renormalized, and independent continuous variables under a linear `Eq` are
restricted to the hyperplane. A continuous inequality condition on such a tuple
raises SymPy's `NotImplementedError`.

Multivariate Normal and StudentT, asymmetric multivariate Laplace, Dirichlet,
Multinomial, NegativeMultinomial, Ewens, NormalGamma and generalized multivariate
log-Gamma models implement joint densities, component moments and joint draws.
`MultivariateBeta` aliases Dirichlet. `marginal_distribution(X, *indices)` returns
a native callable `Lambda`; omitted coordinates are integrated or summed when
no closed marginal is available. Dirichlet coordinates satisfy a simplex
constraint, and multinomial coordinates sum to the trial count; the joint
formulas must be used on those constrained supports.

```python
from symox import Matrix
from symox.stats import MultivariateNormal, Dirichlet

V = MultivariateNormal('V', [0, 0], [[1, Rational(1, 2)], [Rational(1, 2), 1]])
covariance(V[0], V[1])          # 1/2
E(V[0], Eq(V[1], 2))          # 1
variance(V[0], Eq(V[1], 2))    # 3/4

A = Dirichlet('A', 2, 3, 4)
E(A[0])                      # 2/9
covariance(A[0], A[1])        # -1/135
```

Point conditioning updates the full remaining Gaussian law using its covariance
matrix. NormalGamma has its conjugate conditional laws. Multinomial and Ewens
models condition by exact enumeration of their finite joint outcomes. Point
conditioning for other joint families (Dirichlet, for example) raises
`NotImplementedError`. Multivariate Normal, StudentT and Laplace laws with
`MatrixSymbol` parameters of symbolic dimension are refused at construction.
`Laplace(name, mu_list, sigma_matrix)` builds a multivariate Laplace law. Arbitrary
multivariate inequality conditioning and general independence proofs are not
complete. `independent` recognizes disjoint spaces and the supported Gaussian
and independent-component cases; a negative result is not a general proof of
dependence. Ewens partition enumeration is limited to sample size 30.

## Processes

Bernoulli, Poisson, Wiener and Gamma processes preserve dependence across
observation times. Wiener moments use covariance `Min(s, t)`. Poisson/Gamma
moments and path draws use independent increments. Wiener point conditions use
Gaussian conditioning; single-observation Poisson and Gamma bridges use
Binomial and scaled Beta laws. Conditions `Contains(t, interval)` on Poisson,
Wiener and Gamma processes follow the reference convention: the value at `t`
is the increment over the interval, and disjoint intervals are independent.
`BernoulliProcess.joint_distribution` returns a `JointDistributionHandmade`.
Higher non-Gaussian mixed moments require explicit ordered observation times.
More general process conditioning is not implemented.

Finite-state discrete and continuous Markov chains calculate transition powers
or matrix exponentials (`Q*D**t*Q**-1` or `Q*exp(t*D)*Q**-1` when
diagonalizable) and answer finite path queries. The conditioning event must
fix a state at the earliest observation time; later observations follow the
Markov property (SymPy 1.14 ignores intermediate ones). A transition `Matrix`
is checked like SymPy's: dimensions at construction, row sums (0 for a
generator) in `P`/`E` queries, with SymPy's `ValueError`. List inputs, a symox
extension, are checked at construction. Discrete chains accept integer, symbol
or `Str` states; a `MatrixSymbol` of explicit size is used entrywise, and one of
symbolic size has state space `Range(0, n)`. Queries at symbolic times return
the reference `Piecewise` on the time order. `TransitionMatrixOf`,
`GeneratorMatrixOf` and `StochasticStateSpaceOf` provide query-local overrides.

The finite chain computations run in Rust (`crates/symox-core/src/stats/markov.rs`)
for explicit transition matrices: `communication_classes` (with `decompose`,
`canonical_form`, `is_regular`, `is_ergodic` and `is_absorbing_chain`),
`stationary_distribution`, `fixed_row_vector` and `limiting_distribution` over
rational, Float and symbolic entries, `fundamental_matrix`,
`absorbing_probabilities`, and the path sums of `P` and `E` queries. Unique
stationary systems are solved as SymPy's `linsolve` does; a non-unique
solution is SymPy's parametric row (`1 - tau0, tau0`), computed on the Python
path. On a class where SymPy 1.14's period search does not terminate, the
kernel returns the true period. `ContinuousMarkovChain.transition_probabilities`
is the reference method.

A Markov chain constructor does not specify an initial probability law.
Unconditioned Markov probability/expectation queries remain held operators.
A sampled Markov trajectory starts at the first listed state unless an earliest
state is supplied in a conditional sampling query. `sample_stochastic_process`
returns an iterator sampled at integer times; increment processes start at zero.
Process queries concern observations from one process.

## Random matrices

Wishart and MatrixGamma implement positive-definite densities and Bartlett
sampling. MatrixNormal supports row/column covariance and joint Gaussian entry
moments. MatrixStudentT uses an inverse-Wishart mixture and StudentT entry
marginals. Numeric scale matrices must be symmetric positive definite. Singular
Wishart laws are outside the implemented density domain. Matrix-t means require
`nu > 1`, and covariances require `nu > 2`; higher mixed matrix-t moments are not
implemented.

`ExpectationMatrix`, `VarianceMatrix` and `CrossCovarianceMatrix` are matrix
expressions with the reference shapes, and their `expand` follows the reference
(linearity, constant factors, `A*Var(X)*A.T`). Products such as `A*X` of a
`MatrixSymbol` and a `RandomMatrixSymbol` keep their order. Applied to
matrix-distribution variables, where SymPy raises `ShapeError`, they evaluate
matrix statistics. For a column vector, variance is the full component
covariance matrix `E[(X-E[X]) * (X-E[X]).T]`. Row vectors likewise produce a
component covariance matrix. For general rectangular matrices, the additional
implemented convention aggregates column covariances using
`E[(X-E[X]).T * (X-E[X])]`. Transposes are ordinary transposes, not complex
conjugate transposes.

The classical Gaussian and circular ensembles support matrix sampling and
normalized joint eigenvalue/eigenphase laws for beta 1, 2 and 4. Symplectic draws
use a complex `2*n` by `2*n` representation of an `n` by `n` quaternion matrix;
Gaussian eigenvalues occur in pairs. Generic ensemble classes retain a symbolic
beta in the spectral formula and require a concrete classical beta for sampling.
`level_spacing_distribution` is the Wigner surmise, not an exact finite-size
spacing law. Circular matrices use Haar measure; requesting a density with
respect to ambient matrix-entry Lebesgue measure raises `NotImplementedError`.
`JointEigenDistribution` also represents an explicit random matrix's eigenvalue
pushforward using exact Dirac-delta expectations.

Some formulas intentionally correct reference behavior: Davis support starts at
its location, Kumaraswamy support is `[0, 1]`, Gompertz's MGF has the positive
exponent needed for its positive mean, NormalGamma's location marginal has the
correct StudentT scale, and Gaussian beta-ensemble spectral constants normalize
the implemented exponent. Tests include scalar reductions and spectral
normalization checks.

## Sampling and expression objects

`sample`, `sample_iter`, Monte Carlo `P`/`E`, and `sampling_density` support
explicit numeric parameters, NumPy sample shapes and reproducible integer or
NumPy-generator seeds. The `numpy` and `scipy` library selectors both use the
local NumPy implementation; SciPy is not required. Named laws use specialized
samplers where supplied, then numerical inverse-CDF sampling where needed.
Handmade joint sampling and general conditional bridge samplers are not supplied.
Conditional event sampling uses rejection with a finite work limit. Very rare
events and difficult tail integrals can fail explicitly. Empirical discrete
densities are frequency tables; continuous estimates use a Gaussian kernel.

Random scalars are Python probability-space objects coercible to independent
native `Dummy` expressions. They are not subclasses of SymPy's type hierarchy.
Inside native arithmetic they become random atoms that print as their names
(`X**2 + X`) and keep `pspace` and `distribution`; a random atom sorts as a
plain symbol, so terms can print in another order than SymPy's. `subs` on a
distribution parameter rebuilds the variable. Probability metadata is kept in
the Python process; serialization of a bare native expression does not
reconstruct its probability space.

Held `Probability`, `Expectation`, `Variance`, `Covariance`, `Moment` and
`CentralMoment` objects evaluate through `doit()`, and their `rewrite`
(Integral, Sum, Probability, Expectation), `expand`, `evaluate_integral` and
`evaluate_sum` follow the reference. Nested expectations evaluate with a deep
`doit`. `Variance(X).doit()` and `Covariance(X, Y).doit()` evaluate where SymPy
keeps them held. A narrow `Expr.doit` adapter also evaluates registered held
operators inside compound native expressions. The full SymPy assumptions and
probability-space class protocols are not implemented. Constructors and
queries can also be imported through `stats.crv_types`, `frv_types`,
`drv_types`, `joint_rv_types`, `matrix_distributions`,
`stochastic_process_types`, `random_matrix_models`, `symbolic_probability`,
`rv_interface` and the implemented helpers in `stats.rv`. `stats.crv` and
`stats.frv` expose the continuous and finite domain and probability-space
classes; `stats.drv` and `stats.joint_rv` are absent.
