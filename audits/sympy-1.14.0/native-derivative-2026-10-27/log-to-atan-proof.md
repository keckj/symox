# Scope and branch argument

The selected inputs are exact stock univariate QQ/ZZ polynomials with the same native symbol generator and degree at most 128. Everything else keeps the full original algorithm. Zero inputs also decline, preserving the previous division errors. This is a dispatch bound, not a new mathematical limit.

The vectors use ascending coefficients, but all arithmetic is in the existing unbounded QQ field. The initial rotation `(f,g)=(-g,f)` when `degree(f)<degree(g)` is precisely the existing operation, rather than a new branch adjustment.

Over QQ[x], polynomial quotient and remainder are unique. The kernel repeats the reference half-extended-Euclid updates, starting with Bezout coefficients `(1,0)`, in the same order. Scaling the final coefficient and gcd by the inverse leading coefficient gives the identical monic `h` and representative `s`. It then computes `t=(h-s*a)/b` by exact polynomial division, as reference `dup_gcdex` does. No alternate solution of the Bezout identity is selected.

Consequently, each nonterminal step constructs exactly the same polynomial `u=(f*s+g*t)/h`, then recurses on the identical `(s,t)`. The maximal degree strictly decreases. The terminal step uses the identical quotient `p=f/g`. Failed invariants refuse the route rather than return a partial sum; interruptions propagate as errors.

The output is the unchanged collection of `2*atan(u)` terms and the terminal `2*atan(p)`. No atan addition formula, argument replacement or logarithmic branch normalization is introduced. With the stock evaluated constructors, the terms therefore agree pointwise wherever the reference expression is defined, for complex as well as real substitutions. This argument fixes the constants and branch cuts, whereas checking only derivatives would not. Differential tests independently compare the complete reference expression and values on real sides and complex points, as well as certify the exact rational derivative.

# Protocol and evaluation scope

The public adapter checks the stock atan factory, expression arithmetic, canonical Poly presentation adapter, recursion target and evaluated context. The binding checks exact native Poly storage, domain, generator and all methods/properties used by the recurrence, including the presentation adapter's native delegate. Real overrides, copied functools metadata, custom objects/domains and held evaluation retain the complete original route. Reloading the rationaltools module resets the known recursion function while retaining the original callback guards.

The Python additions are selection/protocol checks only. Polynomial division, extended gcd, exact quotients and construction of the entire atan sum execute in Rust. No runtime SymPy is used.

# Remaining reference differences and validation limits

The inherited zero-polynomial error wording is `polynomial division by zero`; SymPy1.14 uses `polynomial division`. The new route declines those inputs and does not change that existing gap. EX, algebraic, inexact, multivariate and custom-protocol computations retain their full previous implementation.

The separate ratint_ratpart evaluate(False) probe found identical recursive expression trees or errors for all eight cases between immutable releases576a and3ef. Seven successful results also satisfy the independent exact rational certificate A'+B=f/g. SymPy raises RecursionError on seven of these held-context probes and PolynomialError on the remaining input constructed in the context; these are recorded explicitly, not counted as matching successful results. There is no observed new regression justifying a separate guard patch in that probe.

DEBUG timings are not performance evidence. The release benchmark runner is prepared, but no timings were taken during concurrent validation load.
