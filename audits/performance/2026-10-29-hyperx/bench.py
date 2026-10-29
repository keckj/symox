"""hyperexpand shift operator / plan benchmark.

usage: python bench.py <module: symox|sympy> <case|all> [repeats]
Prints one JSON line per case: first call (fresh process when a single case is
given) and the median of the repeated calls, in seconds.
"""
import importlib, json, statistics, sys, time

mod = sys.argv[1]
M = importlib.import_module(mod)
H = importlib.import_module(mod + '.simplify.hyperexpand')
R = importlib.import_module(mod + '.core.random')
S, z = M.S, M.Symbol('z')
I, hyper, meijerg = M.I, M.hyper, M.meijerg
a1, a2, b1, b2, b3 = (S(n + 2) / 7 + S(n + 1) / 11 * I for n in range(5))
op = lambda f: z * f.diff(z)
tn = R.verify_numerically


def shift_a():
    h = hyper((a1, a2), (b1, b2, b3), z)
    return tn(H.ShiftA(a1).apply(h, op), hyper((a1 + 1, a2), (b1, b2, b3), z), z)


def unshift_a():
    h = hyper((a1, a2), (b1, b2, b3), z)
    return tn(H.UnShiftA((a1, a2), (b1, b2, b3), 1, z).apply(h, op), hyper((a1, a2 - 1), (b1, b2, b3), z), z)


def unshift_b_build():
    return H.UnShiftB((a1, a2), (b1, b2, b3), 2, z)._poly


def unshift_errors():
    out = []
    for f in (lambda: H.UnShiftA((1,), (0, -1, 1), 0, z), lambda: H.UnShiftB((0, 1), (1,), 0, z)):
        try:
            f()
            out.append('no raise')
        except ValueError as e:
            out.append(str(e))
    return out


def meijer_unshift():
    h = meijerg((a1,), (a2,), (b1,), (b2,), z)
    s = H.MeijerUnShiftC((a1,), (a2,), (b1,), (b2,), 0, z)
    return tn(s.apply(h, op), meijerg((a1,), (a2,), (b1,), (b2 - 1,), z), z)


def plan_mixed():
    hf = H.Hyper_Function
    return [str(o) for o in H.devise_plan(hf((1, 2, S.Half), (3, S(5) / 2)), hf((2, 3, S(5) / 2), (4, S(7) / 2)), z)]


def plan_errors():
    hf = H.Hyper_Function
    try:
        H.devise_plan(hf([2], [1]), hf([2], [2]), z)
        return 'no raise'
    except ValueError as e:
        return str(e)


def plan_apply():
    hf = H.Hyper_Function
    h = hyper([a1, a2], [b1], z)
    ops = H.devise_plan(hf((a1 + 1, a2 - 1), [b1]), hf((a1, a2), [b1]), z)
    return tn(H.apply_operators(h, ops, op), hyper((a1 + 1, a2 - 1), [b1], z), z)


def reduce_order():
    h = hyper((a1, 3), (b1, 1), z)
    nh, ops = H.reduce_order(H.Hyper_Function((a1, 3), (b1, 1)))
    return len(ops)


def meijerg_floats():
    try:
        return str(M.hyperexpand(meijerg(((3.0, 1), ()), ((S(3) / 2,), (0,)), z)))
    except Exception as e:
        return '%s: %s' % (type(e).__name__, e)


def hyperexpand_plan():
    return str(M.hyperexpand(hyper((1, 2, S.Half), (3, S(5) / 2), z)))[:60]


CASES = dict(shift_a=shift_a, unshift_a=unshift_a, unshift_b_build=unshift_b_build, unshift_errors=unshift_errors,
             meijer_unshift=meijer_unshift, plan_mixed=plan_mixed, plan_errors=plan_errors, plan_apply=plan_apply,
             reduce_order=reduce_order, meijerg_floats=meijerg_floats, hyperexpand_plan=hyperexpand_plan)

names = list(CASES) if sys.argv[2] == 'all' else [sys.argv[2]]
reps = int(sys.argv[3]) if len(sys.argv) > 3 else 5
for n in names:
    t = time.perf_counter(); r = CASES[n](); first = time.perf_counter() - t
    ts = []
    for _ in range(reps):
        t = time.perf_counter(); CASES[n](); ts.append(time.perf_counter() - t)
    print(json.dumps(dict(case=n, first=first, median=statistics.median(ts), result=str(r)[:120])), flush=True)
