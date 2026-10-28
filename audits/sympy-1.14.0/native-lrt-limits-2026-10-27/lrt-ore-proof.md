# Frozen LRT/Ore selection correction

Python protocol selection only; Rust and mathematical algorithms unchanged.
Two source files: rationaltools.py and test_lrt_stock_forwarding.py.

The eight ordinary LRT failures were caused by `resultant is _LRT_RESULTANT`.
The initial bootstrap imports the native builtin before the public wrt/sort
adapter is installed. After reload, rationaltools imports that official Python
wrapper from symox.polys. Every other LRT guard remained true. The captured
builtin is retained; the wrapper is accepted only with the exact bootstrap
code/globals, default signature and closure (name='resultant', native=the captured
builtin). LRT always calls it with only includePRS=True, so the verified branch
forwards directly. User @wraps wrappers and changed delegates retain the full
old algorithm and its callbacks; __wrapped__ is not a certificate.

Both Atan and LRT now validate stored _op_priority definitions through builtin
MRO/dict descriptors. The definitions are compared to bootstrap captures without
calling properties. This check precedes the no-Ore fast return: the ordinary
arithmetic delegate also reads priority. Custom attribute hooks likewise decline.
The observed quartic Atan priority callback and the LRT callbacks are compared
to forced old routes, including matching KeyboardInterrupt messages/traces.

LRT supplies the same verified-Ore machinery with its additional native scalar
classes and arithmetic methods. This covers the observed One multiplication and
Symbol power callbacks while keeping Atan's existing operator capture separate.
The public Expr facade's as_poly override is checked on the methods actually
resolved by those scalar classes, rather than only the native storage base.
These two existing facade tests had been masked by the always-false resultant
guard after reload; the original 157 tests are unchanged.

LRT requires the original bootstrap capture flag. Normal reload preserves the
canonical captures; removing/reimporting the module after bootstrap keeps the
complete old route instead of blessing newly captured user operators. Wild and
host/custom inputs retain the explicit full fallback already in the base.

This is a bounded correction for the demonstrated routing and callback cases.
No universal claim about arbitrary mutation of library internals is made. No
new Python mathematics, runtime SymPy, Rust code or recompiled binary is added.
