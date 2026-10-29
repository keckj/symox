# hyperexpand shift operators and plans (wave 21, hyperx)

Harness: `bench.py` (11 public calls), `run.sh <label> <pythonpath|sympy> raw.jsonl`
(3 fresh processes per case; each records the first call and the median of 5
repeated calls), `summarize.py raw.jsonl` (median over processes).

- base: commit 82d674aa Python tree with its release `_core.abi3.so`
  (md5 fac9cd5a812d6af0de59eb4feb7f966f).
- after: this branch, release build from the same tree.
- sympy: SymPy 1.14 from `/home/keckj/Documents/symox/.venv`.

Machine was shared (load average about 7.5 on 8 cores), so the small ratios
are noisy. Times in seconds.

| case | base first | after first | ratio | base repeat | after repeat | ratio | sympy repeat | after/sympy |
|---|---|---|---|---|---|---|---|---|
| shift_a | 0.9384 | 0.0257 | 36.6x | 1.3543 | 0.0035 | 390.1x | 0.0144 | 0.24 |
| unshift_a | 3.4797 | 0.0361 | 96.4x | 4.2955 | 0.0135 | 317.5x | 0.0626 | 0.22 |
| unshift_b_build | 0.0320 | 0.0047 | 6.8x | 0.0122 | 0.0040 | 3.0x | 0.0049 | 0.82 |
| unshift_errors | 0.0251 | 0.0005 | 46.6x | 0.0007 | 0.0002 | 3.3x | 0.0007 | 0.30 |
| meijer_unshift | 0.0699 | 0.0405 | 1.7x | 0.0420 | 0.0146 | 2.9x | 0.0410 | 0.35 |
| plan_mixed | 0.0308 | 0.0027 | 11.2x | 0.0087 | 0.0022 | 3.9x | 0.0099 | 0.22 |
| plan_errors | 0.0005 | 0.0002 | 2.8x | 0.0003 | 0.0001 | 4.1x | 0.0000 | 1.55 |
| plan_apply | 38.7721 | 0.0575 | 674.1x | 39.4263 | 0.0329 | 1199.4x | 0.0760 | 0.43 |
| reduce_order | 0.0021 | 0.0025 | 0.9x | 0.0018 | 0.0017 | 1.0x | 0.0005 | 3.32 |
| meijerg_floats | 0.0210 | 0.0284 | 0.7x | 0.0012 | 0.0034 | 0.4x | 0.9504 | 0.00 |
| hyperexpand_plan | 0.0093 | 0.0089 | 1.0x | 0.0055 | 0.0055 | 1.0x | 0.1167 | 0.05 |

result differs: unshift_b_build base='Poly((864 - 1848*I)*(5082 + 2695*I)/(4161600*z)*_E' after='Poly((864 - 1848*I)*(5082 + 2695*I)/(4161600*z)*_x' sympy='Poly((864 - 1848*I)*(5082 + 2695*I)/(4161600*z)*_x'
result differs: unshift_errors base="['Cannot decrement a unit upper parameter', 'Singu" after="['Cannot decrement unit upper index.', 'Cannot inc" sympy="['Cannot decrement unit upper index.', 'Cannot inc"
result differs: plan_mixed base="['<symox.simplify.hyperexpand.ShiftB object at 0x7" after="['<Decrement lower 7/2.>', '<Decrement upper index" sympy="['<Decrement lower 7/2.>', '<Decrement upper index"
result differs: plan_errors base='no raise' after='Non-suitable parameters.' sympy='Non-suitable parameters.'
result differs: meijerg_floats base="TypeError: '3.00000000000000' object cannot be int" after='-2.36327180120735*z**(3/2)*hyper((-0.5, 3/2), (5/2' sympy='-2.36327180120735*z**(3/2)*hyper((-0.5, 3/2), (5/2'

Notes:
- `shift_a`, `unshift_a`, `plan_apply` include `verify_numerically`; the
  base cost was the evalf precision escalation on `f - g` (now SymPy's
  separate `.n()` plus `comp`).
- `plan_errors`, `unshift_errors`, `plan_mixed`: base returned no error or
  a different error text, `plan_mixed` had no SymPy `str`. `meijerg_floats`:
  base raised `TypeError`; after returns SymPy's result (3.4 ms repeated
  against SymPy 0.95 s; base 1.2 ms was the error path).
- Slowdowns: `meijerg_floats` against base (error path, see above);
  `reduce_order` is unchanged Python and still 3.3x slower than SymPy
  (1.7 ms against 0.5 ms), and `plan_errors` is 0.1 ms against SymPy
  0.05 ms (one native call plus error conversion).
