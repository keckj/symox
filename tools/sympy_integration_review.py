"""Scoped gamma-moment and Piecewise branch-protocol invariant records."""
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
    x, a = sx.symbols('x a')
    probe = sx.integrate(x**a*sx.exp(-x), (x, 0, sx.oo))
    assert probe.is_Piecewise
    assert probe.args[0].expr == sx.gamma(a+1)
    assert probe.args[0].cond == (sx.re(a)>-1)
    assert probe.func(*probe.args) == probe
    root = Path(__file__).resolve().parents[1]
    sources = [
        'crates/symox-core/src/integrate.rs', 'crates/symox-py/src/lib.rs',
        'python/symox/_matrixexpr_methods.py', 'python/symox/_public_calculus.py',
        'python/symox/_expression_replace.py', 'python/symox/_compat.py',
        'python/symox/series/fourier.py', 'python/symox/__init__.py',
    ]
    hashes = {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in sources}
    common = dict(native_sources_sha256=hashes,
        proof_kind='Human invariant argument with independent exact regressions; not machine-checked.',
        universal_function_parity_certified=False)
    records = []
    for identifier, reference, native in [
        ('sympy.integrate', sp.integrate, sx.integrate),
        ('sympy.Expr.integrate', sp.Expr.integrate, sx.Expr.integrate),
    ]:
        records.append(dict(common, id=identifier,
            status='gamma_moment_convergence_and_conditions_scope_review',
            signature_reference=str(inspect.signature(reference)),
            signature_native=str(inspect.signature(native)),
            evidence='tests/test_gamma_integral_conditions.py; tests/test_gamma_signed_powers.py; tests/test_integrate_poles.py',
            scope='On (0,infinity), c*x**a*exp(delta-beta*x**q), finite complex parameters, provably finite real q of known nonzero sign, and provably Re(beta)>0. Includes rational, symbolic signed-real and unbounded integer exponents. Whole-line extension requires proved even q and known integer parity of a, with convergence checked first. Sums use sufficient individual conditions and one original-integral fallback. Not all sign-unknown/complex q, oscillatory boundary regions or Risch/MeijerG strategy flags.',
            invariant_argument=(
                'For x>0, x**a=exp(a*log(x)) on the principal branch, and '
                'for q>0 absolute convergence at zero is Re(a)>-1 and '
                'exponential decay gives convergence at infinity. For q<0 '
                'the decaying endpoint is zero and convergence at infinity '
                'requires Re(a)<-1. Set '
                's=(a+1)/q and t=beta*x**q. The ray has argument strictly '
                'between -pi/2 and pi/2. Analytic contour rotation to the '
                'positive ray has vanishing small/large connecting arcs when '
                'Re(s)>0. A negative q reverses the integration direction, '
                'giving exp(delta)*Gamma(s)/(abs(q)*beta**s). beta**s '
                'uses the principal logarithm in this right half-plane. '
                'Unknown convergence remains a Piecewise with the ORIGINAL '
                'unevaluated integral as fallback; a provably false condition '
                'does not return a finite analytic continuation. For a sum, '
                'the conjunction of individual conditions licenses summing '
                'the values, while the original sum remains the fallback; '
                'individually divergent terms are never integrated separately '
                'into infinity-infinity. Bounds reversal negates the value '
                'without changing the condition. conds=none suppresses '
                'undecided displayed conditions; conds=separate returns the '
                'value and condition instead of a Piecewise. Multiple-limit '
                'separate rejects as in the reference.'),
            precision='Exact symbolic inputs retain native gamma values, powers and real-part inequalities, without a binary64 convergence threshold. Complex-rate comparisons normalize the known exact principal logarithm of 1+I against reference polar representation. Inexact input arithmetic retains its declared precision; this does not certify all numerical gamma evaluation or integration error bounds.',
            complexity='Phase coefficients are grouped sparsely by expression-valued monomial degree; degrees 10**6 and 2**100 are tested without dense degree allocation or reference imports. No i64 degree conversion bounds this moment reduction. Initial expansion, coefficient arithmetic, gamma evaluation and expression growth remain costs; no general asymptotic comparison against SymPy is certified.'))
    records.append(dict(common, id='sympy.Piecewise',
        status='constructor_class_and_paired_branch_protocol_scope_review',
        signature_reference=str(inspect.signature(sp.Piecewise)),
        signature_native=str(inspect.signature(sx.Piecewise)),
        evidence='tests/test_piecewise_branch_protocol.py; tests/test_piecewise_constructor_class.py; tests/test_expr_replace_options.py',
        scope='Finite supplied expression/condition branches: tuple/list/native-Tuple construction, grouped public args with expr/cond accessors, reconstruction, replacement, pickle and consumers covered by regressions. Public native-backed constructor class, instance predicates, fromiter, unbound methods and synchronized aliases. Held evaluation still drops false branches and branches after a true condition, and collapses a single true branch. Not nominal Function/Expr inheritance, every condition-conversion or option convention, or all inherited Basic/ExprCondPair protocols.',
        invariant_argument='Each internal branch (value,condition) yields one tuple-compatible public view referencing those exact two native expressions, with expr/cond/args/func access. The public Piecewise constructor class delegates construction to the native engine and tests instances by their native kind; all already loaded aliases of the native factory are synchronized to that class. The public head reconstructs these branches directly. Held evaluation applies an ordered linear pass removing false branches and stopping at a true branch; one remaining true branch returns its expression. Structural, replacement, folding, Fourier and legacy flat reconstruction adapters consume whole pairs rather than regrouping a nested vector as if it were flat. Internal Rust branch storage and exact expressions are retained.',
        precision='Argument presentation does not approximate branch values or conditions. Reference interchange is explicit; native reflection/replacement subprocesses forbid SymPy imports.',
        complexity='Public args creates O(k) pair views/references for k branches without recursively copying their DAGs; each pair field access is O(1). Instance checks inspect the native kind without reconstructing the expression. Held branch filtering uses O(k) additional references. Constructor condition simplification and downstream algorithms are separate costs.'))
    return dict(schema=1, reference_version=sp.__version__,
        native_source_sha256=native_fingerprint(Path(sx.__file__).resolve().parent), records=records,
        limitations='Scoped families and interface invariants only; global integration and Basic protocol parity remain unproved.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    result = review()
    Path(args.out).write_text(json.dumps(result, indent=2)+'\n')
    print(f"{len(result['records'])} scoped integration/branch protocol reviews recorded")
