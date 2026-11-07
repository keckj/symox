# algpoly: composite domains over QQ<a> in Rust (wave 44)

Lot `feat-algpoly-20261106`, based on main 583315a8 (wave 43). Scope:
`QQ<a>[t]` and `QQ<a>(t)` (elements, arithmetic, gcd, division, resultant,
sqf/factor, Poly over them), `construct_domain(extension=True)`,
unification, the 10 `ratint_logpart` cases left by ratint2, and the `EX`
element normalization.

## Harnesses

- `diffcorpus.py sympy|symox OUT.jsonl [cases...]` (139 cases): `R = A[t]` and
  `F = A.frac_field(t)` elements with `A = QQ.algebraic_field(sqrt(2))`,
  `Poly` over `R` and `F` (arithmetic, div, gcd, lcm, cofactors, cancel,
  resultant with and without `includePRS`, discriminant, sqf, factor,
  content, primitive, monic, clear_denoms, multivariate gcd), the public
  functions with `extension=`, `construct_domain(extension=True)`,
  unification with `QQ<a>`, `QQ<b>`, `QQ_I`, `ZZ_I[y]`, `RR`, `EX` and
  `QQ<b>(y)`, `EX` elements, and the 10 `ratint_logpart` cases over
  `QQ<sqrt(2)>[p]` and `QQ<sqrt(2)>(t)`. Each record keeps the printed
  result (or the exception) and a `value` that `compare.py` uses to tell a
  representation difference (DIFF) from a wrong value (WRONG). The timings
  are single first runs on a loaded machine and are only indicative.
  A case killed by the 60 s alarm is recorded as `CaseTimeout`.
- `compare.py ref.jsonl X.jsonl`: SAME / DIFF / DECLINED (NotImplementedError
  or missing) / WRONG (value differs from SymPy's, or an answer where SymPy
  raises).
- `upstream.py sympy|symox OUT.jsonl`: a transcription of the SymPy 1.14
  `test_domains.py` assertions involving algebraic fields
  (`test_Domain_unify` with `ALG`, and the `test_Domain_convert`
  composite-domain loop for every pair involving `QQ<sqrt(2)>` or
  `QQ<sqrt(3)>`), 728 records. The upstream module itself cannot be
  imported against symox (it needs `symox.external.gmpy`, rings, fields).
- `upstream_polytools.py sympy|symox OUT.jsonl`: the `test_polytools.py`
  functions that mention extensions, run statement by statement (a
  failing statement does not stop the function), 46 extension statements.

`raw/sympy.jsonl` is SymPy 1.14 (`~/symox-ci/venv`), `raw/before.jsonl` is the wave-43
`.so` (main 583315a8), and `raw/after.jsonl` is this branch. `raw/compare-*.txt` hold the
per-case tables. `raw/ratint2-after.jsonl` is the ratint2 corpus rerun on this
branch.

## Results

| corpus | before | after |
|---|---|---|
| algpoly corpus (139) | SAME 34, DIFF 10, DECLINED 95, WRONG 0 | SAME 134, DIFF 5, DECLINED 0, WRONG 0 |
| `test_domains.py` transcription (728) | 281 identical | 728 identical |
| `test_polytools.py` extension statements (46; SymPy 44) | 39 | 39 |
| ratint2 corpus (75) | 62 identical | 73 identical |

By group (algpoly corpus):

| group | cases | before SAME / DIFF / DECLINED | after SAME / DIFF | SymPy ms | after ms |
|---|---:|---|---|---:|---:|
| `R` elements | 14 | 13 / 1 / 0 | 14 / 0 | 54.1 | 21.4 |
| `F` elements | 19 | 1 / 2 / 16 | 19 / 0 | 90.1 | 1.7 |
| Poly over `R` | 34 | 0 / 0 / 34 | 33 / 1 | 575.2 | 106.3 |
| Poly over `F` | 25 | 0 / 1 / 24 | 22 / 3 | 722.8 | 80.0 |
| unify | 12 | 0 / 1 / 11 | 12 / 0 | 132.4 | 5.9 |
| construct_domain | 3 | 3 / 0 / 0 | 3 / 0 | 1.1 | 2.5 |
| public functions | 12 | 11 / 1 / 0 | 11 / 1 | 62.7 | 5.5 |
| `EX` elements | 10 | 6 / 4 / 0 | 10 / 0 | 8.4 | 0.9 |
| `ratint_logpart` | 10 | 0 / 0 / 10 | 10 / 0 | 1446.7 | 111.5 |

The timing columns cover the cases that succeed in both runs.

## The 5 differences left (no WRONG)

- `PR monic`, `PF gcd`, `PF lcm`, `PF sqf_list`: the values are equal and only the
  representative of the numerator/denominator pair differs. symox Poly
  coefficients over `QQ<a>(t)` are stored as expressions. Every result
  goes through `From<ExpressionPoly> for Polynomial`
  (`gaussian_ground::normalize_poly`), which converts each
  coefficient again from its expression into the canonical pair. SymPy keeps
  whatever pair its history produced (`dup_monic` cancels `c/lc`;
  `dmp_ff_prs_gcd` and `dup_lcm` use unnormalized subresultants, and its
  `sqf_list` even keeps coefficients with hundreds of digits).
  Example: SymPy `sqrt(2)*x*(t + sqrt(2))/(2*t)`, symox
  `x*(sqrt(2)*t/2 + 1)/t`.
- `factor(x**2 - 2*t**2, x, extension=sqrt(2))`: SymPy raises
  `CoercionFailed` (its `QQ<sqrt(2)>` ground cannot hold `t`). symox
  raises `ValueError: coefficient is not in the selected algebraic number
  field`, as it did before the lot (`AlgebraicPoly` path).

The test_polytools statements that fail both before and after involve no
composite algebraic ground: `CRootOf` extensions (`rootof` primitive elements,
issue 26808), `factor_list` with `sqrt(2)*pi`, `real_roots`/`all_roots`
with `extension=True` on a multivariate `Poly`, and `cancel(greedy=False)`.
SymPy fails 2 of the 46 statements itself.

Further gaps found during the lot. They also occur with the wave-43 `.so`,
so they are not regressions:

- the native `factor_list` of products with `sqrt(2)` fails its exact
  reconstruction certificate (`poly/checked_factors.rs`, `verify`). The
  failing cases are `(x - sqrt(2))*(x + sqrt(2))`, `(x - sqrt(2)*t)**2`
  and `(x - sqrt(2)*pi)*(x + sqrt(2)*pi)`;
- `cancel((x**2 - 2)/(x + sqrt(2)), greedy=False)` is returned unchanged,
  where SymPy gives `x - sqrt(2)`.
