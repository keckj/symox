"""Scoped finite-group invariants; never universal API parity certificates."""
import hashlib
import inspect
import json
from pathlib import Path

from sympy_parity import native_fingerprint


CHAIN_ARGUMENT=(
    'Every stored transversal sends its base point to its dictionary key and '
    'fixes the earlier base points. Appending generators closes each orbit by '
    'BFS while retaining valid earlier representatives. The queued rectangles '
    'cover every orbit-point/generator pair exactly once as those two lists '
    'grow. For a pair (beta,g), U_beta*g*U_(beta*g)^-1 fixes the current base '
    'point. Sifting divides out successive transversal representatives. A '
    'nonidentity residue either expands a deficient orbit or supplies a moved '
    'point extending the base. At termination every Schreier generator lies '
    'in the next stabilizer, so Schreier lemma proves each stabilizer level '
    'and the final trivial stabilizer. Orbit-stabilizer then proves that the '
    'product of basic orbit lengths is exactly the group order. Sifting to '
    'identity is consequently equivalent to group membership.')

CHAIN_COMPLEXITY=(
    'For degree n and m initial generators, at most O(n^2) successful '
    'extensions occur: each adds a new point to a deficient basic orbit, or '
    'extends a base of length at most n. A level has at most m+O(n^2) '
    'generators and n orbit points. Persistent queues therefore process at '
    'most O(n^2*(m+n^2)) Schreier pairs. A pair and its sift take O(n^2) '
    'image operations. A conservative expected bound is O(m*n^4+n^6) image '
    'operations and O(m*n+n^3) stored image/index words, using amortized '
    'Python hashing. Reference incremental Schreier-Sims uses the same '
    'deficient-orbit extension criterion and restarts pair scans/rebuilds '
    'transversals after extensions. Native order/membership never enumerate '
    'all group elements. Tight comparative bounds, Python hash adversaries '
    'and integer bit costs are separate outstanding review obligations.')

REVIEWS={
    'order':(
        CHAIN_ARGUMENT,
        CHAIN_COMPLEXITY+' After construction the cached exact order is returned in O(1). '
        'For one generator cycle lengths give its order by exact LCM, in O(n) '
        'image operations plus integer GCD/LCM costs. Proven named S_n/A_n '
        'orders use exact factorial arithmetic.',
        'All validated finite permutation generators at a common finite degree; '
        'groups built by public constructors and proven named-family factories.'),
    'contains':(
        CHAIN_ARGUMENT+' A cyclic-group fast path checks a uniform shift on '
        'every disjoint cycle and combines the shift congruences by generalized '
        'CRT; consistency is equivalent to a single common generator power.',
        CHAIN_COMPLEXITY+' Subsequent chain membership takes O(n^2) image '
        'operations. Cyclic membership takes O(n) image operations plus exact '
        'CRT arithmetic. Named symmetric/alternating membership reduces to '
        'degree validity and exact parity.',
        'Validated permutation arguments; strict=True at the group degree, '
        'or strict=False with a valid identity-tail embedding.'),
    'schreier_sims':(
        CHAIN_ARGUMENT,CHAIN_COMPLEXITY,
        'Validated finite permutation groups; computes and caches a BSGS. '
        'Choice/order of an equally valid base, strong generators and '
        'transversals need not equal reference printing.'),
    'schreier_sims_incremental':(
        CHAIN_ARGUMENT+' Optional SLP tracking starts with singleton words '
        'for original generators. Persistent transversal paths append their '
        'edge generator; inverse paths reverse and invert those edges. The '
        'Schreier word and subsequent sift divisions therefore evaluate to '
        'the stored residue. A new strong generator references only earlier '
        'strong generators or their inverses. Adjacent inverse cancellation '
        'preserves that product. Reversing words on export gives public rmul '
        'order. Repeated base points add trivial stabilizer levels.',
        CHAIN_COMPLEXITY+' That bound assumes distinct supplied base points. '
        'For k supplied points including repetitions put r=k+n: orbit storage '
        'has at most r*n points and a sift has at most r levels. Optional '
        'word tracking materializes O(n*r) signed indices per Schreier pair '
        'and stores O(n^3*r) witness indices as a conservative bound for '
        'O(n^2) added generators. Words reference earlier strong generators '
        'without recursively expanding them into original generators; no '
        'group-element enumeration is introduced. These bounds include output '
        'storage and do not certify comparative reference SLP complexity.',
        'slp_dict=False or True; valid supplied base points (including '
        'repetitions) and validated common-degree generators. Strong '
        'generators and witness words may differ while certifying the same '
        'BSGS. Empty/identity-only lists have exact trivial-group outputs; '
        'reference 1.14 errors on some empty filtered lists.'),
    'generate_schreier_sims':(
        'The completed chain gives a unique transversal decomposition for each '
        'group element. Mixed-radix iteration visits each combination exactly '
        'once. Composing representatives from deepest stabilizer to the root '
        'therefore emits every element once and emits no outside element.',
        CHAIN_COMPLEXITY+' After the chain, cached partial products and radix '
        'sizes at least two give O(n) amortized image work per emitted array, '
        'plus public Permutation construction when af=False. Enumeration is '
        'lazy and does not preallocate storage proportional to the group order.',
        'af=False or af=True, validated finite groups. Enumeration order may '
        'differ because equivalent valid BSGS choices differ.'),
    'orbit':(
        'BFS starts at the supplied point, tuple or subset and repeatedly '
        'applies generators. The reached set is generator-closed. In a finite '
        'permutation group positive generator words already contain inverse '
        'actions, so this is exactly the group orbit. Unions repeat the same '
        'closure with multiple seeds.',
        'Point/union actions take O(n*m) image applications and O(n) point '
        'storage. Tuple/subset actions depend on the output orbit size and '
        'the number of supplied entries; they do not enumerate unrelated '
        'group elements. Hashing and output construction have their usual '
        'Python costs, as in the reference BFS implementation.',
        'Valid points, tuples or subsets with action="tuples", "sets" or '
        '"union". No certificate for arbitrary user-defined iterators.'),
    'orbit_transversal':(
        'BFS retains the generator word reaching each first-discovered point; '
        'multiplying the parent representative by that generator sends the '
        'starting point to the new point. One representative is retained '
        'per orbit point, hence the output is exactly an orbit transversal.',
        'At most n points and n*m point actions. At most n representative '
        'products of n images, giving O(n*m+n^2) image work and O(n^2) image '
        'storage; matches the reference representative-BFS bound.',
        'Valid point and pairs=False/True; representative choices may differ.'),
    'stabilizer':(
        'For a valid orbit transversal U, the products U_beta*g*U_(beta*g)^-1 '
        'fix the starting point. Schreier lemma states that those products '
        'generate its entire stabilizer. Identity products are omitted '
        'without changing the generated subgroup.',
        'Orbit-transversal work plus at most n*m products/inverses of n-image '
        'arrays, O(m*n^2) image work and output-dependent generator storage. '
        'Uses the same Schreier-generator construction as the reference; '
        'does not enumerate group elements.',
        'Valid point and validated finite generators; returned generating '
        'lists may differ while defining the same exact subgroup.'),
}


def review():
    import sympy
    from sympy.combinatorics import PermutationGroup as Reference
    from symox.combinatorics import PermutationGroup as Native
    from symox.combinatorics import _schreier
    if sympy.__version__!='1.14.0':raise ValueError('Reference version must be 1.14.0')
    chain_sha=hashlib.sha256(inspect.getsource(_schreier).encode()).hexdigest()
    reference_chain_sha=hashlib.sha256(inspect.getsource(Reference.schreier_sims_incremental).encode()).hexdigest()
    records=[]
    for name,(argument,complexity,scope) in REVIEWS.items():
        native,reference=getattr(Native,name),getattr(Reference,name)
        assert inspect.signature(native)==inspect.signature(reference),(name,inspect.signature(native),inspect.signature(reference))
        records.append(dict(id='sympy.combinatorics.PermutationGroup.'+name,
            status='finite_group_invariant_review_complexity_comparison_pending',
            signature=str(inspect.signature(native)),invariant_argument=argument,
            precision='Exact finite image arrays and unbounded Python integers; no numerical approximation or element-count cutoff.',
            complexity=complexity,scope=scope,
            native_method_sha256=hashlib.sha256(inspect.getsource(native).encode()).hexdigest(),
            reference_method_sha256=hashlib.sha256(inspect.getsource(reference).encode()).hexdigest(),
            native_chain_sha256=chain_sha,reference_chain_sha256=reference_chain_sha,
            proof_kind='human invariant argument and differential/independent closure tests; not machine-checked',
            universal_function_parity_certified=False))
    return dict(schema=1,reference_version=sympy.__version__,native_source_sha256=native_fingerprint(),
        evidence='tests/test_schreier_group_parity.py; tests/test_schreier_slp_parity.py',
        limitations='Finite-group result invariants do not certify all constructor/conversion protocols, inherited Basic methods or every invalid input. Comparative complexity certification remains pending.',
        records=records)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    result=review();Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(len(result['records']),'scoped finite-group invariant reviews recorded')
