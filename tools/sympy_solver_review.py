"""Scoped proofs for complete elementary trigonometric solution lattices."""
import hashlib
import inspect
import json
from pathlib import Path

from sympy_parity import native_fingerprint


def review():
    import sympy
    import symox
    if sympy.__version__!='1.14.0':raise ValueError('Reference must be SymPy 1.14.0')
    assert str(inspect.signature(sympy.solveset))==str(inspect.signature(symox.solveset))
    root=Path(__file__).resolve().parents[1]
    files=['crates/symox-core/src/sets/solveset.rs','crates/symox-core/src/special/lambert.rs']
    hashes={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files}
    return dict(schema=1,reference_version=sympy.__version__,native_source_sha256=native_fingerprint(),
        evidence='tests/test_trigonometric_solution_lattices.py',
        limitations='Elementary-family result invariants only. Not a proof of all solve/solveset inputs, parameter strata, all domains or comparative bit complexity.',
        records=[dict(id='sympy.solveset',
            status='trigonometric_lattice_invariant_review_complexity_comparison_pending',
            signature=str(inspect.signature(symox.solveset)),
            scope='Single sin/cos/tan/cot/sec/csc/exp of a*x+b equal to a finite numeric target, a provably nonzero, over Complexes; over Reals use exact range exclusion only for a provably real affine argument, or restrict the complete complex families. Does not cover arbitrary symbolic parameter degeneracies or every simplification of an intersection.',
            invariant_argument=(
                'Putting w=exp(I*z), sin(z)=k and cos(z)=k become quadratics '
                'w^2-2*I*k*w-1=0 and w^2-2*k*w+1=0. Their constant terms '
                'are nonzero, so both roots have logarithmic preimages, yielding '
                'the two exact inverse-function centers modulo 2*pi. Repeated '
                'centers may merge but no branch is dropped. Tangent inversion '
                'gives one center modulo pi except for the omitted values ±I. '
                'Substitution z=a*x+b is bijective when a!=0 and transforms '
                'each entire integer lattice, not a bounded sample. Cotangent '
                'at zero is the cosine-zero lattice; at nonzero k it reduces '
                'to tan(z)=1/k. Secant/cosecant at nonzero k reduce to '
                'cos(z)=1/k or sin(z)=1/k; those preimages cannot be poles. '
                'At zero they have no finite solutions. Exact nonzero proofs '
                'use assumptions or a certified nonzero real/imaginary part. '
                'Exp(z)=k has the complete lattice log(k)+2*pi*I*n for '
                'nonzero finite k and no solution for k=0. '
                'Changing a center by an integer multiple of its period '
                'preserves the full lattice. If normalization cannot be '
                'certified its original center is retained. Real sine/cosine '
                'targets outside [-1,1] are excluded only for a real affine '
                'argument, by rational endpoint proof. To restrict c+n*p to '
                'the reals, solve Im(c)+n*Im(p)=0: a real step with a real '
                'offset retains all indices, a real step with a certified '
                'nonreal offset retains none, and a certified nonzero '
                'imaginary step permits just n=-Im(c)/Im(p), provided it is '
                'an integer. Unresolved real restrictions remain exact set '
                'intersections. No finite sample is used to infer emptiness.'),
            precision='Exact symbolic ImageSets over unbounded integer indices; range/normalization decisions use outward rational enclosures and arbitrary-size integers. No finite branch-count sample or machine epsilon defines the returned families.',
            complexity='At most two symbolic lattice families per elementary inverse, with expression arithmetic and optional adaptive interval normalization. There is no group/solution-element enumeration or finite period-index cutoff. Symbolic preprocessing, coefficient bit growth, interval refinement and comparative reference complexity remain separate obligations.',
            improvements_over_reference=(
                'For sec/csc(z)=0 reference 1.14 returns {nan}; exact finite '
                'solution set is empty. For tan/cot(z)=±I reference can return '
                'infinite pseudo-roots; the omitted-value set is empty. These '
                'cases have independent analytic tests, not reference equality.'),
            native_sources_sha256=hashes,
            native_wrapper_sha256=hashlib.sha256(inspect.getsource(symox.solveset).encode()).hexdigest(),
            reference_wrapper_sha256=hashlib.sha256(inspect.getsource(sympy.solveset).encode()).hexdigest(),
            proof_kind='Human invariant argument plus exact whole-lattice coset comparisons; not machine-checked.',
            universal_function_parity_certified=False)])


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True);args=parser.parse_args()
    Path(args.out).write_text(json.dumps(review(),indent=2)+'\n')
    print('1 scoped trigonometric-lattice invariant review recorded')
