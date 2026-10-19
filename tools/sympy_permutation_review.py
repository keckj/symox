"""Record scoped invariant/precision/complexity reviews of permutation algorithms.

This development tool imports the pinned reference to check signatures and
fingerprint its source. The invariant arguments are human reviews, not
machine-checked proofs or certificates for arbitrary Python object protocols.
"""
import hashlib
import inspect
import json
from pathlib import Path

from sympy_parity import native_fingerprint


REVIEWS={
    'inversion_vector':(
        'While visiting images from right to left, the Fenwick tree counts exactly the values already visited. Its prefix below the current value counts exactly the smaller images to its right. Dropping the last, always-zero count matches the reference representation.',
        'O(n log n) index operations and O(n) auxiliary words; reference O(n^2).'),
    'inversions':(
        'Each inverted pair has exactly one left endpoint. Summing the inversion vector therefore counts each inverted pair once and no other pair.',
        'O(n log n) index operations and O(n) auxiliary words, matching the reference large-input merge-count bound.'),
    'from_inversion_vector':(
        'At each position, selecting the specified zero-based order statistic among unused images creates exactly the requested count of smaller images remaining. Removal preserves the unused-image invariant. The last image is forced.',
        'O(n log n) index operations and O(n) auxiliary words; reference repeated list removal is O(n^2).'),
    'rank':(
        'Lexicographic rank is the factoradic value sum(v_i*(n-i-1)!). Horner evaluation of the inversion vector is this exact sum. Cached results preserve subsequent constant-time lookup.',
        'O(n log n) index operations plus O(n) exact growing-integer operations; reference O(n^2) index operations plus O(n) growing-integer operations. Cached lookup O(1) in both.'),
    'unrank_lex':(
        'Successive exact divisions extract factoradic digits of rank modulo n!. Choosing and removing each digit-th unused image reconstructs precisely the permutation with that lexicographic rank.',
        'O(n log n) order-statistic operations plus exact factorial/modulo and O(n) integer divisions; reference O(n^2) image adjustment plus factorial-sized arithmetic.'),
    'next_lex':(
        'The rightmost ascent precedes a decreasing suffix. Swapping its left entry with the smallest larger suffix entry and reversing the suffix gives the least greater lexicographic permutation. A decreasing permutation has no successor.',
        'O(n) time and output/working storage, matching the reference bound.'),
    'rank_nonlex':(
        'Each step records the last image s, moves the largest image to the last position, and updates the inverse by the same swap. The recurrence r_n=s+n*r_(n-1) is exactly the reference Myrvold-Ruskey recurrence; iterative reconstruction removes recursive stack limits.',
        'O(n) array operations plus O(n) exact growing-integer operations, matching the reference arithmetic work; no recursive Python stack.'),
    'unrank_nonlex':(
        'The digit r mod n identifies the image to exchange with the last slot, then r//n is the subproblem rank. Iterating the same swaps is the reference recurrence after reduction modulo n!.',
        'O(n) swaps plus factorial/modulo and O(n) exact divisions, matching the reference arithmetic bound; no recursive Python stack.'),
    'next_nonlex':(
        'The rank/unrank bijection in Myrvold-Ruskey order makes unrank(rank+1) the unique successor; equality with n! identifies the last permutation.',
        'Same ranking, factorial and unranking arithmetic bounds as the reference composition.'),
    'rank_trotterjohnson':(
        'For each largest image j, the tree counts its position among images <=j. In the insertion ordering the block direction alternates with the previous rank parity: rank_(j+1)=(j+1)*rank_j + (j-position if rank_j is even else position). This is the reference recurrence.',
        'O(n log n) index operations plus O(n) exact growing-integer operations; reference O(n^2) searching plus the same rank arithmetic.'),
    'unrank_trotterjohnson':(
        'Exact quotient/remainder division reverses the insertion-rank recurrence. Placing images from largest to smallest at the recorded order statistics reconstructs the unique permutation for every rank in [0,n!-1].',
        'O(n log n) index operations plus exact factorial/modulo and O(n) integer divisions; reference O(n^2) insertion and factorial-sized arithmetic.'),
    'next_trotterjohnson':(
        'An image moves left when its lower-image subpermutation is even and right when odd. Tree counts update those parities exactly. Exchanging the largest mobile image with its adjacent smaller neighbor is the Steinhaus-Johnson-Trotter successor; absence of a mobile image proves termination.',
        'O(n log n) index operations and O(n) words; reference repeatedly removes/searches subarrays and recomputes parity with an O(n^2) worst-case bound.'),
    'josephus':(
        'The current survivor rank records the front of the circular queue. Advancing by m-1 modulo its length and removing that order statistic is the same queue rotation/removal. The retained survivors are emitted from the next front. The deque branch performs the identical rotations directly.',
        'O(n) initialization/output plus O(k*min(m,log n)) index operations for k removals, with exact large-step modulo costs. The reference performs O(k*m) single deque rotations.'),
    'get_adjacency_distance':(
        'Each adjacent ordered pair occurs at most once in a valid permutation. Counting shared pairs equals the reference matrix inner product; subtracting from n-1 gives the same value, including -1 for degree zero.',
        'Expected O(n) time, worst-case O(n^2) under hash collisions, O(n) words; reference O(n^2) matrix work/storage.'),
    'get_positional_distance':(
        'The implementation sums the identical exact absolute differences of corresponding images after the same equal-degree precondition.',
        'O(n) arithmetic operations and O(1) auxiliary words; reference O(n) plus copied arrays.'),
    'get_precedence_distance':(
        'Map the first sequence through positions in the second. A pair disagrees in precedence exactly when it is an inversion in this relative sequence. Every pair is counted once.',
        'O(n log n) index operations and O(n) words; reference O(n^2) matrices and pair comparisons.'),
}


def review():
    import sympy
    from sympy.combinatorics import Permutation as Reference
    from symox.combinatorics import Permutation as Native
    if sympy.__version__!='1.14.0':raise ValueError('Reference version must be 1.14.0')
    records=[]
    for name,(argument,complexity) in REVIEWS.items():
        native,reference=getattr(Native,name),getattr(Reference,name)
        assert inspect.signature(native)==inspect.signature(reference),(name,inspect.signature(native),inspect.signature(reference))
        records.append(dict(id='sympy.combinatorics.Permutation.'+name,
            status='closed_integer_domain_review',
            signature=str(inspect.signature(native)),
            invariant_argument=argument,
            precision='Exact Python integers throughout; no floating-point conversion. Ranks may have factorial-sized bit lengths.',
            complexity=complexity,
            scope='Validated finite permutation image arrays; integer ranks/sizes; valid inversion vectors or a valid mutable inverse array when supplied; Josephus survivor threshold as documented. Resource availability is required.',
            native_method_sha256=hashlib.sha256(inspect.getsource(native).encode()).hexdigest(),
            reference_method_sha256=hashlib.sha256(inspect.getsource(reference).encode()).hexdigest(),
            proof_kind='human invariant argument plus executable differential regressions; not machine-checked',
            universal_function_parity_certified=False))
    return dict(schema=1,reference_version=sympy.__version__,native_source_sha256=native_fingerprint(),
        evidence='tests/test_permutation_algorithm_parity.py; tests/test_permutation_constructor_parity.py; tests/test_matrix_permutations.py',
        limitations='These closed-domain reviews do not certify arbitrary Python conversion/subclass protocols, every invalid argument shape, inherited Basic methods, all constructors or global SymPy parity. Complexity counts exact-integer operations separately from their bit costs.',
        improvements=[
            'next_lex on the empty permutation returns None instead of IndexError.',
            'next_trotterjohnson on the empty permutation returns None instead of the same permutation.',
            'Trotter-Johnson unranking extends periodically outside the documented valid rank interval, preserving a bijective image array instead of returning duplicate images.',
            'Iterative Myrvold-Ruskey ranking/unranking avoids the reference Python recursion-depth limit.'],
        records=records)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    report=review()
    Path(args.out).write_text(json.dumps(report,indent=2)+'\n')
    print(len(report['records']),'scoped permutation reviews recorded')
