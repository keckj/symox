Historical initial preparation, superseded by PUBLIC-REPORT.md and executed evidence below. No claims in this initial plan replace the later results.

FFLU inspection and proposed contract — preparation only

Snapshot read: native-normalforms-risch-lll-r11-20261028/after.
Source1801af86581c822b537ed5899d0e7d03ee046db1987b39f13f422fe6fe44e1f9; RELEASE532c4b08d546d9d580306ad7e40317e2fbbbe80aa8d77d4749b35682d3ee9c90. The R11 full gate is not yet green. No engine code, build, test or timing was started for this task. R1 HNF delivery, main and all performance preparations remain untouched.

Observed real routes

- DDM._fflu, ddm.py:968: Python Gauss–Bareiss elimination into an actual self.copy(), returning packed LU and permutation indices. Actual prefix reads shape/domain, calls self.copy(), constructs list(range(rows)), sets rank=0.
- DDM.fflu, ddm.py:1019: calls actual self._fflu(), then Python loops construct P, L, D, zero the packed multipliers in U and return P,L,D,U.
- DomainMatrix.fflu, domainmatrix.py:3401: actual self.from_rep capture, self.rep.fflu(), four from_rep calls. No domain gate or options here.
- SDM.fflu, sdm.py:1103: actual to_dfm_or_ddm().fflu(), then four to_sdm conversions. GROUND_TYPES='python' and DFM_dummy._supports_domain=False mean the stock route uses DDM.
- A separate _dfm.py FFLU exists for optional FLINT, but is not the active stock backend. It is outside the proposed port.
- SymPy1.14 installed source contains both identical Python phases. FFLU is a real API; this is not ordinary LU and not fraction-free RREF.

Scope to settle before implementation

The intended complete public claim is only stock homogeneous native ZZ/QQ on DDM and DomainMatrix dense/sparse, all shapes/ranks without artificial size/bit bounds. Other domains, custom arithmetic/types, aliases or unsafe protocols retain complete fallback. Stock native int cells are not silently conflated with native domain elements. Phase1-only work must be named Phase1 only; public FFLU is not ported until Phase2 arithmetic/public routes are covered.

The best split for preserving observations is a Rust elimination kernel plus separate Rust construction/publication kernels for P, L/U and D. Keep actual Python copy/zeros/from_rep/to_sdm allocation and conversion barriers in their existing order. P must be populated before the real L=zeros call; L/U must be published before D=zeros. Each barrier can execute callbacks/change state, so reauthenticate each native phase independently. Do not compute future phase values from stale pre-callback snapshots.

Private _fflu must retain its two-item result and packed-multiplier convention; public fflu retains four factors. Keep source cells/rows unchanged where the actual copy does so, exact output/domain/representation signatures, permutation object aliases, L-to-U pivot/multiplier aliases, first-D-to-L diagonal alias, distinct K.zero/K.one publication conventions and error/interrupt partial state. Raw malformed storage, copy aliases and overridden _fflu/zeros require the established generic behavior, not eager new errors. Audit/foreign module/foreign key cases learned from R1 must be admission controls. Unknown native element finalizers, arithmetic hooks and slot overrides require full fallback.

Mathematical concerns from source, NOT yet executed

Phase1 currently iterates only min(rows,cols) columns. It also chooses denominator LU[rank-1][rank-1] rather than the preceding actual pivot after skipped columns. Manual derivation predicts:

1. A=[[0,0,1,2],[0,0,2,3]]: Phase1 visits only columns0/1, so performs no elimination. Phase2 appears to return L=[[1,0],[2,3]], D=diag(1,3), U=[[0,0,1,2],[0,0,0,3]], P=I. L*inv(D)*U has final entry7 instead of A's3.
2. A=[[0,1,2,3],[0,2,3,4],[0,3,4,6]]: after the first nonzero pivot, the next elimination numerator is -1 while LU[0][0]=0, predicting exquo(-1,0).

Both must first be reproduced on exact immutable SymPy1.14 and Symox R11, with source/binary identities and Fraction certificates. They are not performance baselines. The parent explicitly requires a mathematically complete factorization, not preservation of a false PLDU identity. Any confirmed inherited bug needs an explicit documented intentional reference correction and regression contract.

A candidate mathematically correct algorithm visits every column until rank==rows, chooses the first nonzero pivot below rank, swaps full rows/perm, keeps the actual previous nonzero pivot as the denominator, and uses exact Bareiss quotients for all lower trailing entries. Phase2 then extracts retained multipliers by actual pivot column and completes unused L diagonal entries with one; D uses consecutive L diagonal products. This candidate needs algebraic proof and exact all-shape/rank testing before coding. The fallback contract remains an explicit open decision: a stock-domain correction cannot silently return the old wrong answer when native metadata is unavailable. A corrected generic fallback or a separately documented capability gap must be resolved before any complete public claim.

Required independent certificates

- P is an m-by-m permutation matrix and P*A reorders input rows according to the private permutation; verify the documented orientation rather than assuming A*P.
- L is m-by-m lower triangular with nonzero diagonal; D is diagonal/nonzero (including rank-deficient/all-zero cases using trailing ones); U is m-by-n row-echelon/upper triangular with strictly increasing pivot columns.
- Use fractions.Fraction alone to prove every entry of P*A == L*inv(D)*U; no float residuals. Empty shapes are interpreted algebraically without calling an inverse unnecessarily.
- Independent exact Gaussian rank(A)==rank(U)==number of nonzero U pivot rows; factor identity and invertible L/D certify completeness.
- For each ZZ quotient, record/prove numerator=denominator*quotient and denominator!=0 when the division is actually consumed. A zero denominator that is never consumed because a trailing range is empty must not trigger a speculative exception in an uncorrected legacy control.
- Compare all factor entries, shapes, domains, dense/sparse format and object aliases against SymPy where the reference factors satisfy the certificate; independently certify/document intentional corrected outputs where it does not.

Planned validation, not yet launched

After the parent's quiet-window release: exact preimplementation reproductions first. Then private stage from accepted corrected R11, standalone Rust/binding/helper and minimal ddm.py dispatch; no changes to dense.py, sdm.py or earlier deliveries unless required and explicitly scoped. Test DDM._fflu/DDM.fflu/SDM.fflu/DomainMatrix dense and sparse; dimensions0/1/2/3/4/8 including0xn/mx0, tall/wide, all-zero, signed/huge integers/rationals, delayed/empty pivot columns, row swaps at later ranks, dependent rows and every rank. Add native first-call assertions, NoReference workers, stock-signature comparisons, callback/namespace/audit/dictionary-key/error/mutation/alias controls, and unsupported GF/polynomial/inexact domain fallbacks. Then targeted plus ordered matrix/domain/protocol cohorts on local and Demeter, Rust/check/clippy, cumulative base/final manifests and safe patch replay. Parent runs the complete integrated release gate. No benchmark until baseline correctness/completeness, all release gates and an explicit quiet-window signal.
