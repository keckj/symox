"""Record scoped radical-control and complex-branch invariants, not global parity."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path

from sympy_parity import native_fingerprint


def review():
    import symox as sx
    import sympy as sp

    if sp.__version__ != '1.14.0':
        raise ValueError('Reference must be SymPy 1.14.0')
    root = Path(__file__).resolve().parents[1]
    files = [
        'crates/symox-core/src/assumptions_exact.rs',
        'crates/symox-core/src/assumptions.rs',
        'crates/symox-core/src/evalf_mp.rs',
        'crates/symox-core/src/poly/roots.rs',
        'crates/symox-core/src/solve.rs',
        'crates/symox-core/src/solve/recast.rs',
        'crates/symox-core/src/sets/solveset.rs',
        'crates/symox-py/src/solve.rs',
        'crates/symox-py/src/poly.rs',
        'crates/symox-py/src/lib.rs',
        'crates/symox-py/src/logic.rs',
        'python/symox/_public_special.py',
        'python/symox/_solve_relations.py',
        'python/symox/__init__.py',
        'python/symox/_solver_checks.py',
        'crates/symox-core/src/inspection.rs',
    ]
    hashes = {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files}
    common = dict(
        native_sources_sha256=hashes,
        proof_kind='Human invariant argument and independent numerical/exact regressions; not machine-checked.',
        universal_function_parity_certified=False,
    )
    records = [dict(
        common,
        id='sympy.N',
        status='principal_sqrt_and_cardano_pair_invariants_general_accuracy_pending',
        signature_reference=str(inspect.signature(sp.N)),
        signature_native=str(inspect.signature(sx.N)),
        evidence='tests/test_complex_power_branch_accuracy.py; tests/test_cubic_conjugate_accuracy.py',
        scope='Finite complex square roots and half-integer powers evaluated by the native MP path; rational affine Cardano pairs and reciprocal rational-power products satisfying exact norm identities. This is not a certificate of every evalf function, input or requested precision.',
        invariant_argument=(
            'For z=x+I*y with y!=0, set L=sqrt((abs(z)+abs(x))/2)>0 and '
            'S=y/(2*L). If x>=0 the principal square root is L+I*S; if x<0 '
            'it is abs(S)+I*sign(y)*L. The real part is positive and the '
            'imaginary part has the sign of y. Computing L uses a sum of '
            'nonnegative quantities; S uses the exact identity 2*Re(w)*Im(w)=y '
            'rather than a subtraction of nearly equal norms. For y=0 the '
            'principal real/positive-imaginary branch is selected explicitly. '
            'For c>0, abs(z)^2=c^3 and the principal u=z^(1/3), c/u=conjugate(u). '
            'Thus t=u+c/u is real in [sqrt(c),2*sqrt(c)] and is the largest '
            'root of t^3-3*c*t-2*Re(z). This polynomial is increasing above '
            'sqrt(c); evaluating it at a rational threshold certifies the '
            'sign of a rational affine expression in t. Native recognition '
            'requires exact rational norm equality and rational real part. '
            'For real rational exponents, a finite nonzero product T has '
            'squared modulus equal to the product of the powered squared '
            'moduli of its bases. When this is exactly rational q and b/a=q>0, '
            'a*T+b/T=a*(T+conjugate(T)) is real. Reciprocal factor vectors '
            'are matched after common base ordering; rational content may '
            'be extracted from additive bases only at small integer powers, '
            'where distributing the power is branch-independent. '
            'Only an exact structural real proof permits discarding a '
            'computed imaginary roundoff component. Genuine small imaginary '
            'inputs are never removed by an epsilon test. Cancellation in '
            'real and imaginary sums is tracked separately.'),
        precision='Independent cloned mpmath contexts resolve tiny components through scale 1000 before componentwise relative comparisons at 60 decimal digits. Quartic checks resolve exact real/imaginary structure in the reference before numerical evaluation and survive native pickle reconstruction. All three real cubic roots are compared with independent reference indexed roots through 200 digits. These tests do not prove correctly rounded output for all inputs; the internal 4096-bit extra-refinement budget, maxn and strict options remain obligations.',
        complexity='The stable square-root formula uses a constant number of MP arithmetic operations and square roots, without subtractive norm cancellation. Reciprocal-pair lookup uses a hash table, O(F) signature storage, and O(sum(f_i*log(f_i))) factor ordering for F stored factors across terms. Exact expression arithmetic, norm computations and coefficient bit sizes still count. No comparative whole-evalf bit-complexity certificate.',
    )]
    x, y, z = sx.symbols('x y z')
    assert (x+y-z).evalf(50,subs={x:1e16,y:1,z:1e16}) == sx.Float(1,50)
    substitution_scope=(
        'Sequential symbol substitutions in ordinary scalar arithmetic and '
        'functions, including nested powers, complex components and scalar '
        'matrix entries. Algebraic keys, binding expressions, Boolean/set '
        'control flow and the Python special-function adapter retain existing '
        'paths; their full substitution precision is not certified.')
    substitution_invariant=(
        'For each symbol mapping, unchanged leaves retain their exact stored '
        'values and changed arithmetic nodes are rebuilt with held native '
        'constructors. No addition, multiplication or power of replacement '
        'values occurs during this traversal. Traversing later mappings also '
        'transforms earlier inserted values, preserving sequential semantics. '
        'The resulting tree reaches the existing arbitrary-precision evaluator, '
        'which refines cancelling real and imaginary sums separately. The '
        'Python adapter delegates ordinary expressions directly to this path '
        'rather than calling eager .subs() first. Shared nodes are memoized '
        'with source references retained so pointer keys cannot be reused.')
    records[0]['evidence'] += '; tests/test_evalf_substitution_precision.py'
    records[0]['scope'] += ' Numerical substitution scope: '+substitution_scope
    records[0]['invariant_argument'] += ' '+substitution_invariant
    records.append(dict(common, id='sympy.Expr.evalf',
        status='held_scalar_substitution_and_refinement_scope_review',
        signature_reference=str(inspect.signature(sp.Expr.evalf)),
        signature_native=str(inspect.signature(sx.Expr.evalf)),
        evidence='tests/test_evalf_substitution_precision.py; tests/test_evalf_chop.py; tests/test_precision.py',
        scope=substitution_scope, invariant_argument=substitution_invariant,
        precision='Exact T+1-T=1 identities with identical finite binary literals through magnitude 10**1000, 15/50/120 requested decimal digits and nonlinear outer expressions. Separate-component relative checks resolve complex cancellation and a genuine 10**-100 branch-cut component using an independent 400-digit mpmath context. The reference receives a larger maxn only for test-oracle refinement; native maxn/strict are still unsupported. This is not a correctly rounded result or error-bound certificate for every numeric function or exhausted refinement budget.',
        complexity='The scalar eligibility scan memoizes expression DAG nodes. Each mapping visits and rebuilds changed DAG nodes once with a fresh memo, O(sum_i(V_i+E_i)) expected traversal work across intermediate DAGs, plus symbol comparisons and output storage. Floating exponents remain compact; there is no rationalization or exponent-sized allocation in substitution. The cost of subsequent multiprecision refinement and function evaluation remains separate.'))
    for name in ('solve', 'roots'):
        records.append(dict(
            common,
            id='sympy.'+name,
            status='radical_options_and_complete_QQ_fallback_scope_review',
            signature_reference=str(inspect.signature(getattr(sp, name))),
            signature_native=str(inspect.signature(getattr(sx, name))),
            evidence='tests/test_polynomial_root_controls.py'+('; tests/test_solve_output_conventions.py; tests/test_solve_factor_components.py; tests/test_solve_parameter_constraints.py; tests/test_solve_undetermined_dispatch.py; tests/test_solve_expression_targets.py; tests/test_solve_relational_dispatch.py' if name=='solve' else ''),
            scope=('Univariate rational-coefficient polynomial radical controls, including separately enabled quartics and exact indexed-root solve fallback; roots strict/filter/predicate behavior. '+('For solve, polynomial components, generic parameter constraints, rational coefficient identities, native sequences/output conventions, one commutative nonnumeric expression target with a rational recast and exact reconstruction certificate, and relational dispatch retaining all conditions through the native inequality engine. The relational family is bounded by that engine, including its remaining multivariate restrictions. ' if name=='solve' else '')+'Not every coefficient domain, parameter stratum, mixed target list, solver flag or quintic formula.'),
            invariant_argument=(
                'Cubic and quartic formula dispatch consults lexical thread-local '
                'options. A guard restores the prior options on exit or unwind; '
                'Python callbacks run after leaving that scope. The quartic '
                'resolvent may use its internal cubic regardless of whether '
                'external cubic formulas were disabled. Disabling radical '
                'formulas in solve preserves complete QQ polynomial solutions '
                'through the existing factor/indexed-root fallback, including '
                'multiplicity accounting and exact polynomial pole rejection. '
                'roots returns its supported radical subset; strict=True checks '
                'the retained multiplicity against the polynomial degree after '
                'filtering/predicates and raises UnsolvableFactorError if '
                'incomplete. Zero-root filter/predicate behavior follows the '
                'pinned reference. Formal Poly generator data are read directly '
                'rather than reconstructed by substitution.'+(
                ' A supplied finite polynomial product with nonzero constant '
                'coefficient vanishes iff at least one factor vanishes. Each '
                'positive-integral factor is therefore solved with every other '
                'equation retained, and the branches are combined without '
                'dropping an unsupported branch as an empty solution. Free '
                'factors independent of the selected unknowns are discarded '
                'only if their nonzeroness is proved. This preserves complete '
                'factor components where pivot division formerly lost one. '
                'When no proved-nonzero linear pivot already parametrizes '
                'the graph and a leading coefficient contains another unknown, '
                'native factorization can recover expanded product components '
                'before division. A proved-nonzero pivot is preferred over '
                'a singular one of equal degree. General coefficient-domain '
                'factorization and parameter degeneracies remain obligations. '
                'Identity assignments represent free '
                'parameters and are omitted from mapping outputs, while '
                'ordered tuple outputs retain their columns. Zero rows or '
                'parameter-only equations prove inconsistency only when their '
                'nonzeroness is proved; undecided residuals retain generic '
                'candidates and original-equation checking, not a false empty '
                'set. This does not prove all parameter strata. Bare scalar '
                'solve may first collect independent rational-function '
                'coefficients and accept a single nonempty solution mapping. '
                'Independent Symbol monomials form the certified basis; '
                'potentially dependent functional bases are excluded from '
                'this automatic attempt. Failed/multiple coefficient attempts '
                'fall back to ordinary algebraic solving. The collector uses '
                'the native system solver, outside the scalar radical scope. '
                'Native Tuple inputs are unfolded as equation/variable/initial '
                'point sequences, preserving container distinctions. A single '
                'non-Symbol scalar/system target can be represented by a fresh '
                'hygienic Dummy. Undefined-function derivative/binder '
                'coefficients remain separate exact algebraic atoms; they '
                'are not differentiated as replaced constants. Integer '
                'exponential templates use exp(k*u+delta)=exp(delta)*exp(u)**k. '
                'A recast is used only if rational in the new atom, if the '
                'atom actually occurs, and if restoring ALL original atoms '
                'reconstructs the original exactly or by ordinary ring '
                'expansion. Fractional lifting and forced power/log identities '
                'are not reconstruction proofs. Original expressions/metadata '
                'are restored in results and labels. This follows reference '
                'algebraic-atom conventions, not analytic function-range or '
                'differential-equation semantics; for example exp(x) as a '
                'requested atom may have the formal solution value zero. '
                'For exp(4*x+2)=5, the atom exp(2*x+1) has precisely the two '
                'values ±sqrt(5), by quadratic factorization and exact '
                'preimages; the reference can miss this affine-argument case. '
                'For relational solve inputs, the public dispatcher sends the '
                'entire original condition list to reduce_inequalities. It '
                'does not discard equations or inequalities or turn this result '
                'into an algebraic candidate list. Public aliases share this '
                'dispatcher; nonrelational inputs retain the native solver.' if name=='solve' else '')),
            precision='Exact polynomial coefficients and indexed roots; independent full root-count assertions, domain assumptions and numeric branch regressions. An undecidable residual is not a proof of absence. General exact verification and certified complex isolating boxes remain open.',
            complexity='Option lookup/restoration has constant overhead. Filtering is linear in the produced distinct roots plus predicate cost. '+('Product splitting can explore the product of factor counts across equations before duplicate/intersection elimination; this may exceed final output size. Conditional factor recovery adds native polynomial factorization cost. Coefficient dispatch adds structural rational classification, cancellation, sparse coefficient collection, canonical unknown sorting and a system solve; degree 10**6 is tested without dense degree allocation. Atom recasting adds DAG inspection, binder containment scans, exact substitution/reconstruction and optional expansion/cancellation; their expression growth and worst-case costs remain obligations. ' if name=='solve' else '')+'This does not certify factorization, isolation, algebraic expression growth or general solve bit complexity relative to SymPy.',
        ))
    records.append(dict(common,
        id='sympy.solve_univariate_inequality',
        status='three_real_cubic_roots_exact_sign_interval_scope_review',
        scope='x**3-3*x+1>0 over the reals with relational=False. One independently proved improvement where the pinned reference can refuse a result; not a certificate of all polynomial inequalities.',
        evidence='tests/test_cubic_conjugate_accuracy.py::test_cubic_inequality_contains_the_entire_exact_positive_set',
        invariant_argument='Let z=(-1+I*sqrt(3))/2=exp(2*pi*I/3) on its principal branch. The native endpoints z**(k/3)+z**(-k/3), k=4,7,1, are exactly 2*cos(8*pi/9), 2*cos(4*pi/9), 2*cos(2*pi/9). Cosine monotonicity puts them respectively in (-2,-1),(0,1),(1,2). Exact minimal polynomials and Sturm counts prove these are all three distinct real roots of t**3-3*t+1. A positive leading coefficient and simple roots give precisely (r1,r2) union (r3,infinity), with open endpoints. This proves the entire positive set, not a sample of points.',
        precision='Exact radical endpoints, principal-power conjugacy, rational isolating intervals and Sturm root counts. No machine-epsilon inference of missing branches.',
        complexity='A fixed cubic family and exact interval proof only. General root comparison, isolation and inequality bit complexity remain unproved.',
    ))
    return dict(
        schema=1, reference_version=sp.__version__,
        native_source_sha256=native_fingerprint(), records=records,
        limitations='Scoped invariant reviews only. Every universal function certificate remains false.',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    result = review()
    Path(args.out).write_text(json.dumps(result, indent=2)+'\n')
    print(f"{len(result['records'])} scoped radical/branch invariant reviews recorded")
