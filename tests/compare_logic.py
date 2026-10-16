"""Differential tests for boolean logic (And, Or, ..., to_cnf, simplify_logic,
satisfiable, truth_table) against SymPy.

Every case is evaluated with SymPy and symox. Results are compared by their
printed form; when it differs the results are checked for logical
equivalence (and, for ``simplify_logic``, for equal size). ``SAT:`` cases
check that the symox model satisfies the formula (or that both report
unsatisfiability); ``ALLSAT:`` cases compare the sets of all models.
"""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
x & y
y & x
x | y & z
(x | y) & z
~x
~~x
~(x & y)
~(x | y)
x ^ y
x ^ y ^ z
x ^ x
x ^ True
Xor(x, y, True)
Xor(x, True, True)
Xor(x, y, x)
Xor(x < 1, x >= 1)
Xor(x < 1, 1 > x, y)
Implies(x, y)
x >> y
x << y
Implies(True, x)
Implies(x, x)
Implies(False, x)
Implies(x, True)
Implies(x, False)
Implies(x < 1, x >= 1)
Implies(x < 1, 1 > x)
Implies(x & y, z)
Equivalent(x, y)
Equivalent(x, x)
Equivalent(x, True)
Equivalent(x, False)
Equivalent(x, y, True)
Equivalent(x, y, False)
Equivalent(x, ~x)
Equivalent(x & y, z)
Equivalent(x < 1, x >= 1)
Equivalent(x < 1, 1 > x)
ITE(x, y, z)
ITE(True, y, z)
ITE(False, y, z)
ITE(x, True, False)
ITE(x, False, True)
ITE(x, y, y)
x & ~x
x | ~x
x & True
x & False
x | False
x | True
And(x, And(y, z))
Or(x, Or(y, Or(z, x)))
And()
Or()
And(x)
And(x < 1, y > 2)
Or(x < 1, x > 2)
~(x < 1)
~(x <= 1)
~(x > 1)
~(x >= 1)
~Eq(x, 1)
~Ne(x, 1)
And(x < 1, x >= 1)
Or(x < 1, x >= 1)
And(x < 1, 1 > x)
And(Eq(x, 1), y)
And(Ne(x, 1), y)
x & (x < 1)
Not(x & y | z)
~x & ~y
(x | y) & (~x | z) & a
Nand(x, y)
Nor(x, y)
Xnor(x, y)
Nand(x, False)
Nor(x, True)
And(x > 1, x < -1)
Or(x > 0, x < 0)
And(x < 1, x < 2)
(x > -oo) & (x < 1)
(x < 1) & (x > -oo)
And(y, x < 1, -oo < x)
And(x, y, x)
And(~x, x, y)
Or(x, y, x | z)
And(a, b) | And(a, c)
Not(True)
Not(False)
Not(Not(x) & y)
And(Eq(x + y, 1), Or(a, b))
Or(Eq(x + y, 1), And(a, b))
And(x + y > 1, a)
Xor(Eq(x + y, 1), Or(a, b))
Equivalent(Eq(x + y, 1), Or(a, b))
And(Implies(x, y), z)
Or(Implies(x, y), z)
And(Equivalent(x, y), z)
Xor(x & y, z)
Xor(x | y, z)
Or(Xor(x, y), z)
And(Xor(x, y), z)
Not(Xor(x, y))
Not(Implies(x, y))
Not(Equivalent(x, y))
And(ITE(x, y, z), a)
~ITE(x, y, z)
Implies(~x, y | z)
And(sin(x) > 0, a)
(x < 1) & (x < y)
Piecewise((1, x < 1), (2, (x > 2) & (x < 3)), (0, True))
Piecewise((1, (x > 2) | (x < -3)), (0, True))
Piecewise((1, ~(x > 2)), (0, True))
to_cnf((a & b) | c)
to_cnf(a >> b)
to_cnf(Equivalent(a, b))
to_cnf(Equivalent(a, b, c))
to_cnf(Xor(a, b, c))
to_cnf(~(a | b) | c)
to_cnf((a & b) | (c & d))
to_cnf((a & b) | (c & d) | (x & y))
to_cnf(ITE(a, b, c))
to_cnf(~ITE(a, b, c))
to_cnf(Implies(a & b, c | d))
to_cnf(a & (b | c))
to_cnf(a)
to_cnf(~(a & b & c))
to_cnf((a | b) & (a | b | c), True)
to_cnf((a & b) | (a & ~b), True)
to_dnf((a | b) & c)
to_dnf(Implies(a, b) & c)
to_dnf(Xor(a, b))
to_dnf((a | b) & (c | d))
to_dnf(~(a | b))
to_dnf(Equivalent(a, b))
to_dnf((a | b) & (a | c), True)
to_nnf(~(a & b))
to_nnf(Equivalent(a, b, c))
to_nnf(Xor(a, b))
to_nnf(~Xor(a, b))
to_nnf(Xor(a, b, c))
to_nnf(Implies(a, b))
to_nnf(ITE(a, b, c))
to_nnf(~ITE(a, b, c))
to_nnf(~Implies(a, b))
to_nnf(~Equivalent(a, b))
to_nnf((a | ~a) & b)
to_nnf((a | ~a) & b, False)
to_nnf(~(a & ~b) | (c >> d))
eliminate_implications(Implies(a, b))
eliminate_implications(Equivalent(a, b))
eliminate_implications(Equivalent(a, b, c))
is_cnf(a & (b | c))
is_cnf(a | (b & c))
is_cnf((a | b) & ~c)
is_cnf(~(a & b))
is_dnf(a & (b | c))
is_dnf((a & b) | c)
is_dnf(a ^ b)
is_nnf(~a & b)
is_nnf(~(a & b))
is_nnf(a & ~a)
is_nnf(a & ~a, False)
is_literal(a)
is_literal(~a)
is_literal(a & b)
is_literal(x < 1)
is_literal(x + y < 1)
distribute_and_over_or((a & b) | c)
distribute_or_over_and((a | b) & c)
simplify_logic((~x & ~y & ~z) | (~x & ~y & z))
simplify_logic((a & b) | (a & ~b))
simplify_logic(a | (a & b))
simplify_logic(a & (a | b))
simplify_logic((x & y) | (~x & z) | (y & z))
simplify_logic(Xor(a, b))
simplify_logic(Xor(a, b), 'cnf')
simplify_logic(Xor(a, b, c))
simplify_logic(Equivalent(a, b))
simplify_logic(Implies(a, b))
simplify_logic(ITE(a, b, c))
simplify_logic((x < 1) | (x >= 1))
simplify_logic((x < 1) & (x >= 1) | y)
simplify_logic((x < 1) & y | (x >= 1) & y)
simplify_logic(a & ~a)
simplify_logic(a | ~a)
simplify_logic(~(a & b) | c, 'dnf')
simplify_logic(~(a & b) | c, 'cnf')
simplify_logic((a | b) & (a | c))
simplify_logic((a | b) & (a | c), 'dnf')
simplify_logic((a & b) | (c & d), 'cnf')
simplify_logic(~((a | b) & ~c))
simplify_logic((a & ~b & c) | (a & b & c) | (~a & b & c) | (~a & ~b & ~c))
simplify_logic((a & b & c & d) | (a & b & c & ~d) | (a & b & ~c & d) | (~a & ~b))
simplify_logic(Or(And(a, b), And(~a, ~b), And(a, ~b, c)))
simplify_logic(Xor(a, b) | Xor(b, c))
simplify_logic(Equivalent(a, b) & Equivalent(b, c))
simplify_logic(Implies(a, b) & Implies(b, c) & Implies(c, a))
simplify_logic(a)
simplify_logic(x < 1)
simplify_logic(a & b, 'cnf')
simplify_logic(a | b, 'dnf')
SOPform([x, y, z], [[0, 0, 1], [0, 1, 1], [1, 1, 1]])
POSform([x, y, z], [[0, 0, 1], [0, 1, 1], [1, 1, 1]])
SOPform([a, b, c, d], [1, 3, 7, 11, 15], [0, 2, 5])
POSform([a, b, c, d], [1, 3, 7, 11, 15], [0, 2, 5])
SOPform([a, b, c, d], [0, 1, 2, 5, 6, 7, 8, 9, 10, 14])
POSform([a, b, c, d], [0, 1, 2, 5, 6, 7, 8, 9, 10, 14])
SOPform([x, y], [])
POSform([x, y], [[0, 0], [0, 1], [1, 0], [1, 1]])
SOPform([x, y], [[0, 0], [0, 1], [1, 0], [1, 1]])
SOPform([w, x, y, z], [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 1, 1]], [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]])
list(truth_table(x & y, [x, y]))
list(truth_table(x | y, [x, y], input=False))
list(truth_table(Implies(x, y) ^ z, [x, y, z]))
list(truth_table(x & y, [x]))
SAT:x | y
SAT:True
SAT:False
SAT:x
SAT:x & ~x
SAT:x ^ y
SAT:Equivalent(x, y, z) & ~x
SAT:(x | y) & (~x | y) & (x | ~y) & (~x | ~y)
SAT:(x | y | z) & (~x | ~y) & (~y | ~z) & (~x | ~z)
SAT:Implies(a, b) & Implies(b, c) & a & ~c
SAT:ITE(a, b, c) & ~b & ~c
SAT:Xor(a, b, c, d) & Equivalent(a, b) & ~c
SAT:(x < 1) & y
SAT:Nand(a, b) & Nor(c, d) & (a | c)
ALLSAT:x | y
ALLSAT:x & ~x
ALLSAT:Xor(a, b, c)
ALLSAT:Implies(a, b) & (b | c)
ALLSAT:Equivalent(a, b) | c
ALLSAT:ITE(a, b, c)
""".strip().splitlines()


def ns(mod):
    d = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
    if mod is sympy:
        from sympy.logic import boolalg
        from sympy.logic import inference
        for n in dir(boolalg):
            if not n.startswith("_"):
                d.setdefault(n, getattr(boolalg, n))
        d["satisfiable"] = inference.satisfiable
    d.update(zip("xyzwabcd", mod.symbols("x y z w a b c d")))
    return d


def to_sympy(v):
    if isinstance(v, (bool, int, list, tuple, dict)) or v is None:
        return v
    return sympy.sympify(str(v))


def equivalent(a, b):
    """Logical equivalence of two SymPy booleans."""
    try:
        return sympy.logic.inference.satisfiable(~sympy.Equivalent(a, b)) is False
    except Exception:  # noqa: BLE001
        return False


def size(e):
    return sympy.count_ops(e)


def check_model(formula, model):
    if model is False:
        return sympy.logic.inference.satisfiable(formula) is False
    if formula == sympy.true:
        return True
    rep = {}
    for k, v in model.items():
        rep[sympy.sympify(str(k))] = sympy.true if v else sympy.false
    val = formula.xreplace(rep)
    return val == sympy.true


def model_set(models):
    out = set()
    for m in models:
        if m is False:
            continue
        out.add(frozenset((str(k), bool(v)) for k, v in m.items()))
    return out


def run_case(c, sx, mx):
    if c.startswith("SAT:") or c.startswith("ALLSAT:"):
        allm = c.startswith("ALLSAT:")
        body = c.split(":", 1)[1]
        formula = eval(body, sx)
        try:
            got = eval(f"satisfiable({body}, all_models={allm})", mx)
            if allm:
                got = list(got)
        except Exception as e:  # noqa: BLE001
            return "error", f"{type(e).__name__}: {e}", ""
        want = eval(f"satisfiable({body}, all_models={allm})", sx)
        if allm:
            want = list(want)
            ok = model_set(got) == model_set(want)
            return ("same" if ok else "DIFFERENT"), str(got), str(want)
        if str(got) == str(want):
            return "same", str(got), str(want)
        ok = check_model(formula, got) and (got is False) == (want is False)
        return ("valid-model" if ok else "DIFFERENT"), str(got), str(want)
    want = eval(c, sx)
    try:
        got = eval(c, mx)
    except Exception as e:  # noqa: BLE001
        return "error", f"{type(e).__name__}: {e}", str(want)
    if str(got) == str(want):
        return "same", str(got), str(want)
    try:
        g = to_sympy(got)
    except Exception:  # noqa: BLE001
        return "unparseable", str(got), str(want)
    if isinstance(want, sympy.Basic) and isinstance(g, sympy.Basic):
        if g == want:
            return "same-tree", str(got), str(want)
        if equivalent(g, want):
            if c.startswith("simplify_logic") or c.startswith(("SOPform", "POSform")) or "True)" in c:
                # minimal forms: equivalent and no larger
                if size(g) <= size(want):
                    return "equivalent", str(got), str(want)
                return "DIFFERENT", str(got), str(want)
            return "equivalent", str(got), str(want)
    return "DIFFERENT", str(got), str(want)


def run(verbose=True):
    sx, mx = ns(sympy), ns(symox)
    bad, stats = [], {}
    for c in CASES:
        st, got, want = run_case(c, sx, mx)
        stats[st] = stats.get(st, 0) + 1
        if st not in ("same", "same-tree"):
            if verbose:
                print(f"[{st:11s}] {c:50s} symox: {got:40s} sympy: {want}")
            if st in ("DIFFERENT", "error", "unparseable"):
                bad.append((c, got, want))
    if verbose:
        print("\nsummary:", stats, f"(total {len(CASES)})")
    return bad


def test_logic_against_sympy():
    bad = run(verbose=False)
    assert not bad, "\n".join(f"{c}: symox={g} sympy={w}" for c, g, w in bad)


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
