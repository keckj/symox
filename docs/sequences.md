# Symbolic sequences

[README](../README.md) · [Series](series.md)

`sequence`, `SeqFormula`, `SeqPer`, `SeqAdd`, `SeqMul` and `EmptySequence` are
available at the package root and from `symox.series.sequences`. Entries and
linear algebra use native arithmetic; SymPy is only imported by explicit
`to_sympy()` conversion.

```python
from symox import sequence, SeqPer, symbols, oo

n, t = symbols('n t')
squares = sequence(n*n, (n, 0, oo))
assert squares[:5] == [0, 1, 4, 9, 16]
assert squares.coeff(3) == 9
assert squares.find_linear_recurrence(10) == [3, -3, 1]

periodic = SeqPer((1, 2), (n, 0, oo))
assert (squares+periodic)[:4] == [1, 3, 5, 11]
assert squares.coeff_mul(2)[:4] == [0, 2, 8, 18]
```

## Indices, periods and operations

Limits use `(index, start, stop)`, inclusive at both finite endpoints. A
two-element `(start, stop)` form infers an index from a single free symbol;
constants and constant periods use a fresh dummy index. Default limits are
`(index, 0, oo)`. One endpoint can be infinite, and endpoints must be integers,
infinities or symbolic bounds. Known reversed intervals produce `EmptySequence`.

`seq[position]` indexes positions starting at zero. `seq.coeff(point)` evaluates
at the actual sequence index. A sequence with start `-oo` iterates backwards
from its stop; for other sequences it iterates forwards from the start.
Negative positions are rejected. Finite slices return lists; infinite and
symbolic lengths require an explicit finite slice stop. Iterators also support
infinite sequences through ordinary lazy iteration.

`SeqPer` accepts nested periods and index-dependent entries such as `(n, n*n)`.
Its phase starts at the finite endpoint used for iteration. Addition and
multiplication intersect operand intervals. Formula operands combine formulas;
periodic operands align their actual origins before combining periods. Mixed
operands retain a `SeqAdd` or `SeqMul` object. Use `coeff_mul` for scalar scaling.
Empty sequences are neutral for sequence addition and absorbing for multiplication.

The objects expose `formula` or `periodical`, `period`, `gen`, `variables`,
`start`, `stop`, `length`, `interval`, `free_symbols`, `args`, `func`, copying,
hashing, pickling and structural sort keys. `subs` preserves the bound index
and renames it when needed to avoid capturing replacement symbols. Bounds
remain free. `xreplace` applies literal structural replacements, including the
index. `simplify` and `expand` transform generating expressions.

## Fitted linear recurrences

`seq.find_linear_recurrence(n, d=None, gfvar=None)` searches for the shortest
constant-coefficient relation on the first `n` terms, with order at most
`min(d, n//2)` when `d` is supplied. A result `[c1, c2, ...]` means
`a[k] = c1*a[k-1] + c2*a[k-2] + ...` on the fitted sample. The native exact
linear solver proposes coefficients, which are checked against every sample.
Symbolic coefficients are treated generically over the supported scalar field.

With `gfvar=t`, the result is `(coefficients, rational_generating_function)`
using those coefficients and the initial terms. No fitted relation returns
`[]`, or `([], None)` when a generating variable is requested.

A relation fitted to a finite prefix is not a proof for the infinite sequence.
In particular, sufficiently high orders can interpolate arbitrary finite data.
The rational generating function describes the continuation defined by the
fitted recurrence. Independent mathematical information is needed to identify
it with the original infinite formula.

## Limits

Slices materialize at most 100,000 terms. Periods have at most 4,096 entries;
larger combined periods retain a termwise operation without expanding them.
Recurrence searches use at most 256 terms and order 64. Periodic coefficients
need a concrete integer offset from their origin. Undecidable symbolic bounds
raise `NotImplementedError` when evaluating a coefficient.

This layer does not implement the complete SymPy `Basic` method surface for
sequence objects, general asymptotic sequence limits, plotting or code generation.
