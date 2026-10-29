# RR[y] gcd scale, numeric signs, log forms (wave 22, rrnorm)

`bench.py` times 12 public calls (best of 3 in one process) and prints the
first 90 characters of each result. "before" (`before.txt`) is the main
checkout build (base c7289185); "after" (`after.txt`) is branch
`feat-rrnorm-20261029d`. Both are release builds, run under
`systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0`.

| call | before | after | result change |
|---|---|---|---|
| Poly.gcd RR[y] | 0.0002s | 0.0003s | 1.0*x*z+1.0*y -> 1.5*x*z+1.5*y (SymPy) |
| Poly.lcm RR[y] | 0.0002s | 0.0004s | 1.0 scale -> 1.5 scale (SymPy) |
| Poly.cofactors RR[y] | 0.0002s | 0.0003s | (g, 1.5, 3.0) -> (1.5*g, 1.0, 2.0) (SymPy) |
| gcd expr RR[y] | 0.0001s | 0.0002s | 1.0*x*z+1.0*y -> 1.5*x*z+1.5*y (SymPy) |
| Poly.sqf_list RR[y] | 0.0003s | 0.0006s | coeff 1.5 -> 1.0 (SymPy) |
| Poly.integrate RR[y] | 0.0001s | 0.0002s | domain RR[y] -> RR(y) (SymPy) |
| polylog(3,-1/2).is_extended_negative | 0.0000s | 0.0089s | None -> True (SymPy; costs one evalf(2)) |
| hyperexpand lerchphi | 0.1337s | 0.3284s | unchanged (slower: numerical signs of the polylog terms; SymPy 1.14 first call 0.42s, cached 0.08s) |
| solve exp(exp(exp(x)))+1 | 0.0165s | 0.0087s | log(log(pi)+I*pi/2) -> log(log(I*pi)) (SymPy) |
| solve exp(exp(x)+1)-4 | 0.0157s | 0.0085s | log(-1+log(4)) -> log(-1+2*log(2)) (SymPy) |
| simplify(log(-1+2*log(2))) | 0.0121s | 0.0029s | log(-1+log(4)) -> log(-1+2*log(2)) (SymPy) |
| integrate(atan(p*x)/x**2,(x,1,oo)) | 0.0013s | 0.0013s | unchanged (form gap, see docs) |

The new numeric sign answers cost one low-precision evalf per uncached query,
as in SymPy's `Expr._eval_is_extended_positive_negative`.

The hyperexpand call is the one measured slowdown: main at the time of the
rerun took 0.0953s, after takes 0.3284s, below SymPy's first call. A sum of
radicals alone does not take the numerical sign (the evaluation on large
solver roots made `nonlinsolve([x**3 + 1/y - 1, x - y], [x, y])` take 70 s
instead of 25 s before that gate).
