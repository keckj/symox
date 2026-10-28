FFLU pre-existing mathematical failures and exact correction derivation

This is correctness evidence, not a speed comparison or a production port. The scripts read unmodified SymPy1.14.0 and the immutable Symox R11 snapshot at /home/keckj/.cache/symox-validation/native-normalforms-risch-lll-r11-20261028/after. Its native RELEASE identity is532c4b08d546d9d580306ad7e40317e2fbbbe80aa8d77d4749b35682d3ee9c90 and source identity1801af86581c822b537ed5899d0e7d03ee046db1987b39f13f422fe6fe44e1f9. No engine source, earlier delivery or main file was edited. No timers or benchmarks are present.

The DomainMatrix.fflu documentation specifies P*A=L*inv(D)*U. DDM.fflu returns P,L,D,U in the input domain; SDM/DomainMatrix sparse return sparse factors. The private DDM._fflu returns packed LU and a row permutation. P is on the left: a three-cycle witness has permutation[1,2,0] and verifies the documented identity.

Confirmed first witness, both engines, ZZ and QQ, all four private/public entries:

A=[[0,0,1,2],[0,0,2,3]]. Both return P=I, L=[[1,0],[2,3]], D=diag(1,3), U=[[0,0,1,2],[0,0,0,3]]. All advertised shapes, domains, ranks and triangularity happen to pass; nevertheless (L*inv(D)*U)[1,3]=7 whereas (P*A)[1,3]=3. Matching names, ranks and individual factors does not establish a valid decomposition.

Confirmed second witness:

A=[[0,1,2,3],[0,2,3,4],[0,3,4,6]]. Both engines raise ZeroDivisionError despite exact rank3. Symox says 'division by zero'; SymPy says 'integer division or modulo by zero' over ZZ and 'Zero divisor -1/0' over QQ. Those raw messages are preserved in witness-sympy.json and witness-symox.json; the messages already differ before a new port. Inputs, stored rows and stored cells remain unchanged.

Cause: Phase1 visits only min(m,n) columns, so later pivots are not eliminated. After skipped columns, it divides by packedLU[rank-1][rank-1], which need not be the preceding pivot. Phase2's existing formulas are mathematically sound when supplied a correct packed Phase1 result.

Coverage and counts, without hidden exclusions or deduplication:

There are203 fixture IDs:201 integer-valued fixtures each run over ZZ and QQ, plus2 fractional fixtures accepted only over QQ. Thus201*2+2=404 domain cases, not406. The203 fixture IDs represent164 numerically distinct matrices when exact entries and shape are canonicalized; intentional permutation/empty-case duplicates remain in the evidence. The201 integer fixtures comprise3 initial witnesses,48 row/column-permuted witnesses,140 shape/rank fixtures and10 delayed-pivot fixtures. All30 shape pairs m=0..4,n=0..5 are present, including0xn,mx0 and0x0. Every rank0..min(m,n) occurs for each shape, with two positional variants.

Each domain case runs private _fflu, DDM.fflu and DomainMatrix dense/sparse:1616 actual calls per engine. Both engines have1408 returned certified results,84 returned invalid decompositions and124 exceptions. The404-certificate mathematical prototype is a separate dataset and must not be counted as engine success. Likewise the original-Phase2 injection has404certificates for each engine and must not be counted as a complete engine FFLU call. No invalid or exceptional legacy case can be a speed baseline.

Every returned factorization is checked using fractions.Fraction alone: exact P*A=L*inv(D)*U, P permutation, correct shapes/domains/formats, L lower and nonsingular, D diagonal and nonsingular, U upper/echelon, and independent rank(A)=rank(U)=number of pivots. All sources/native identities are recorded unchanged. summary.json also records every invalid result/error with its independently certified input rank. All1056 public legacy results per engine that already satisfy the complete certificate retain exactly the same canonical factor entries in the mathematical correction.

Mathematical correction and justification:

Process every column until all rows have pivots. Select the first nonzero row at/below the current rank and swap full packed rows plus permutation entries. Retain the previous actual nonzero pivot p_prev, initially1. For each lower trailing cell use q=(p*cell-pivot_row_cell*multiplier)/p_prev; retain the multiplier in its packed pivot column. Update the previous pivot only after that pivot's elimination.

Let C be the already selected pivot columns. Before the next pivot, each active value is the determinant of the ordered minor using selected rows plus its row and columns C plus its column. Sylvester's determinant identity gives numerator=p_prev*the next minor. This proves exact divisibility over ZZ and the recurrence over QQ, independent of whether columns were skipped. Each chosen pivot is nonzero, so a consumed denominator is nonzero. Full row swaps also move stored lower multipliers, preserving the final row permutation convention.

With pivots p_i and p_-1=1, Phase2 extracts L's diagonal p_i and its stored lower multipliers. D_ii=p_(i-1)*p_i; unused L diagonals are1 and D retains the existing adjacent-diagonal-product convention. Eliminated residual rows are zero, so those trailing nonzero diagonal values do not change the product. The scaled elimination recurrence telescopes to P*A=L*inv(D)*U. Increasing pivot columns and exact zero residual rows certify rank completeness. Empty shapes use the same algebraic identities without a spurious inverse.

The prototype is only an executable mathematical derivation in proof code. It certifies all404domain cases, including1582 consumed exact divisions, and compares each pivot and each updated quotient against an independently evaluated determinant minor. Injecting the resulting packed LU into the unchanged original public Phase2 gives404 complete certificates on SymPy and404 on Symox, with every factor equal to the independently constructed prototype. The first corrected witness has L=[[1,0],[2,-1]], D=diag(1,-1), U=[[0,0,1,2],[0,0,0,-1]]. The second has L=[[1,0,0],[2,-1,0],[3,-2,-1]], D=diag(1,-1,1), U=[[0,1,2,3],[0,0,-1,-2],[0,0,0,-1]].

Future implementation contract:

The correction is an intentional mathematical divergence from erroneous SymPy1.14 behavior. It needs regression tests and parity/TODO documentation. Rust must port the admitted stock ZZ/QQ mathematical loops, with actual allocation/conversion/callback barriers retained. The original generic Python fallback must also receive the mathematical correction: unavailable native metadata must never restore either known bad result. Complete other-domain/custom-protocol behavior, errors outside these documented corrections, exact arithmetic, aliases and partial mutations must be preserved. Phase2 must be explicitly included before claiming public FFLU native; keeping its mathematics in Python would be a separate remaining port. No new build or port has been validated by this preparation.

Reproduce using the project's existing Python venv:

    PYTHONPATH=<R11>/after/python python reproduce.py --library sympy --fixtures fixtures-witness.json --out new-witness-sympy.json
    PYTHONPATH=<R11>/after/python python reproduce.py --library symox --fixtures fixtures-witness.json --out new-witness-symox.json
    python derive-correction.py
    PYTHONPATH=<R11>/after/python python verify-phase2.py --library sympy
    PYTHONPATH=<R11>/after/python python verify-phase2.py --library symox

Raw initial witness,201-neighbour and2-rational results, source identities, fixtures, prototype trace and validation scripts are included. checksums.sha256 binds all files. The certificate is finite evidence plus the stated algebraic derivation, not a claim of global Symox parity or an implemented Rust FFLU engine.
