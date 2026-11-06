"""classify.py TAG [RAW]: assertion sides of the statement harness whose
symox str differs from SymPy 1.14, grouped by cause.

The cause of a key is the first rule of CAUSES whose key pattern (and,
optionally, SymPy / symox value pattern) matches; unmatched keys are
reported as 'unclassified'.  The rules encode the analysis of the
2026-11-05 limits lot and are applied unchanged to every tag, so a before /
after pair shows which causes shrank."""
import collections, glob, json, os, re, sys

CAUSES = [
    # (cause, key regex, sympy value regex, symox value regex)
    ("sympy-better-or-equal: symox value correct, SymPy wrong/timeout", r"issue_5172:2|eval_special_slow:2", None, None),
    ("gruntz API: compare/sign/mrv/rewrite/mrv_leadterm not exposed", r"mrv_leadterm|demidovich.py:test_leadterm", None, r"^EXC"),
    ("gruntz API: pure gruntz (I*a coefficient -> c0*oo, errors)", r"gruntz_I|gruntz_eval_special_fail", None, None),
    ("Float in limit: nsimplify (_real_to_rational) at entry", r"issue_5164|issue_5383|limit_with_Float|floor:14|ceiling:14|issue_26991", None, None),
    ("Bessel at imaginary infinity", r"bessel_functions_at_infinity", None, None),
    ("heuristics() API missing", r"test_heuristic", None, None),
    ("Limit.free_symbols with non-symbol variable", r"issue_9205", None, None),
    ("held where SymPy evaluates (re/im, complex direction)", r"issue_21756", None, None),
    ("limit evaluates where SymPy holds (validated or proposed deviation)", r"issue_19770|issue_7535", None, None),
    ("value form: rational/trig/sign canonical form differs (same value)", r"test_limits.py:test_(rational|issue_14393|issue_16708|issue_17751|issue_6364|issue_5740|basic1:50|issue_8634|series_AccumBounds|issue_13715|issue_21785)|aseries_trig", None, None),
    ("test not evaluated (statement harness)", r"sympyissue_22986", None, None),
    ("formal power series (fps) gaps", r"test_formal", None, None),
    ("Order: API/attributes (as_expr_variables, leadterm, ...)", r"test_order.py", None, r"^EXC"),
    ("Order: noncommutative/conjugate/subs/multivariate forms", r"test_order.py", None, None),
    ("series of undefined functions / Subs / Derivative", r"issue_3978|series_of_Subs|issue_11313|issue_5223:(7|8|9|10|12|2)$", None, None),
    ("series/nseries: symbolic or log exponents (NotImplementedError)", r"test_(n?series|aseries)\.py", None, r"^EXC"),
    ("series/nseries: value form or order term differs", r"test_(n?series|aseries)\.py", None, None),
]


def classify(k, ref, got):
    for cause, kp, sp, gp in CAUSES:
        if re.search(kp, k) and (sp is None or re.search(sp, ref)) and (gp is None or re.search(gp, got)):
            return cause
    return "unclassified"


def main():
    tag = sys.argv[1]
    args = [a for a in sys.argv[2:] if not a.startswith('-')]
    raw = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'raw')
    groups = collections.defaultdict(list)
    same = total = 0
    for f in sorted(glob.glob(os.path.join(raw, 'tl', 'sympy_*.json'))):
        name = os.path.basename(f)[6:-5]
        ref = json.load(open(f))
        p = os.path.join(raw, 'tl', '%s_%s.json' % (tag, name))
        d = json.load(open(p)) if os.path.exists(p) else {}
        for k, v in sorted(ref.items()):
            total += 1
            got = d.get(k, ['-'])[0]
            if got == v[0]:
                same += 1
                continue
            groups[classify(k, v[0], got)].append(k)
    print('%s: %d/%d identical, %d differ' % (tag, same, total, total - same))
    for cause, ks in sorted(groups.items(), key=lambda t: -len(t[1])):
        print('  %4d  %s' % (len(ks), cause))
        if '-v' in sys.argv:
            for k in ks:
                print('          ' + k)


main()
