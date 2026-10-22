"""Structural pattern matching adapted from SymPy 1.14; see LICENSE-SYMPY.

This is the reference ``matches`` protocol (sympy/core/basic.py,
operations.py, add.py, mul.py, power.py, symbol.py and function.py): the
same traversal order, the same handling of exact and wild parts, the same
rewriting of integer powers and coefficients, and the same ``old`` flag.
``expr.match(pattern)`` is ``matches(pattern, expr)``.
"""

__all__ = ['matches']


def _api():
    from . import (S, Add, Mul, Dummy, ordered, signsimp, collect)
    return S, Add, Mul, Dummy, ordered, signsimp, collect


def _is_wild(node):
    return getattr(node, 'is_Wild', False) is True


def _wild_function_counts():
    from ._public_core import _wild_functions
    return _wild_functions


def _is_wild_function(node):
    try:
        return node in _wild_function_counts()
    except TypeError:
        return False


def _has_wild(node):
    pending = [node]
    seen = set()
    while pending:
        item = pending.pop()
        if _is_wild(item):
            return True
        try:
            if item in seen:
                continue
            seen.add(item)
        except TypeError:
            pass
        pending.extend(getattr(item, 'args', ()))
    return False


def _is_expr(node):
    from .core.expr import Expr
    return isinstance(node, Expr)


def _same_class(pattern, expr):
    if getattr(pattern, 'is_Atom', False):
        return type(expr) is type(pattern) or isinstance(expr, type(pattern))
    try:
        return expr.func == pattern.func
    except AttributeError:
        return isinstance(expr, type(pattern))


def _nonzero(expr):
    # bool(Basic) in the reference: numbers are false at zero.
    return not (getattr(expr, 'is_Number', False) and expr == 0)


def _coeff_isneg(a):
    if getattr(a, 'is_MatMul', False):
        a = a.args[0]
    if a.is_Mul:
        a = a.args[0]
    return bool(a.is_Number and a.is_extended_negative)


def matches(pattern, expr, repl_dict=None, old=False):
    """``pattern.matches(expr, repl_dict, old)`` of the reference."""
    custom = _custom_matches(pattern)
    if custom is not None:
        return custom(expr, repl_dict, old=old)
    if _is_wild(pattern):
        if _is_wild_function(pattern):
            return _wild_function(pattern, expr, repl_dict)
        return _wild(pattern, expr, repl_dict)
    if not hasattr(pattern, 'is_Add'):
        return _native_nonscalar(pattern, expr, repl_dict)
    if pattern.is_Add:
        return _commutative(pattern, expr, repl_dict, old)
    if pattern.is_Mul:
        return _mul(pattern, expr, repl_dict, old)
    if pattern.is_Pow:
        return _pow(pattern, expr, repl_dict, old)
    if getattr(pattern, 'is_Atom', False):
        if pattern == expr:
            return {} if repl_dict is None else dict(repl_dict)
        return None
    return _basic(pattern, expr, repl_dict, old)


def _native_nonscalar(pattern, expr, repl_dict):
    # Matrix expression patterns keep the native structural matcher.
    from ._expr_methods import _native_match
    try:
        result = _native_match(expr, pattern)
    except TypeError:
        result = {} if expr == pattern else None
    if result is None:
        return None
    bindings = {} if repl_dict is None else dict(repl_dict)
    bindings.update(result)
    return bindings


def _custom_matches(pattern):
    # A user subclass defining its own matches keeps it.
    method = getattr(type(pattern), 'matches', None)
    if method is None:
        return None
    module = getattr(method, '__module__', '') or ''
    if module.startswith('symox') or module == 'builtins':
        return None
    return getattr(pattern, 'matches')


def _wild(pattern, expr, repl_dict):
    if any(expr.has(x) for x in pattern.exclude):
        return None
    if not all(f(expr) for f in pattern.properties):
        return None
    repl_dict = {} if repl_dict is None else dict(repl_dict)
    repl_dict[pattern] = expr
    return repl_dict


def _wild_function(pattern, expr, repl_dict):
    if not (getattr(expr, 'is_Function', False) or getattr(expr, 'is_UndefinedFunction', False)):
        return None
    from ._public_core import _wild_nargs
    nargs = _wild_nargs(pattern)
    if nargs is not None and len(expr.args) not in nargs:
        return None
    repl_dict = {} if repl_dict is None else dict(repl_dict)
    repl_dict[pattern] = expr
    return repl_dict


def _basic(pattern, expr, repl_dict, old):
    from . import sympify
    expr = sympify(expr)
    if not _same_class(pattern, expr):
        return None
    repl_dict = {} if repl_dict is None else dict(repl_dict)
    if pattern == expr:
        return repl_dict
    pargs, eargs = pattern.args, expr.args
    if len(pargs) != len(eargs):
        return None
    d = repl_dict
    for arg, other_arg in zip(pargs, eargs):
        if arg == other_arg:
            continue
        if getattr(arg, 'is_Relational', False):
            try:
                d = matches(arg.xreplace(d), other_arg, d, old=old)
            except TypeError:
                d = None
        else:
            d = matches(arg.xreplace(d), other_arg, d, old=old)
        if d is None:
            return None
    return d


def _pow(pattern, expr, repl_dict, old):
    from . import sympify
    S = _api()[0]
    expr = sympify(expr)
    if repl_dict is None:
        repl_dict = {}
    # pattern = 1: the exponent can match 0
    if expr == 1 and expr.is_Integer:
        d = matches(pattern.exp, S.Zero, repl_dict)
        if d is not None:
            return d
    if not _is_expr(expr):
        return None
    b, e = expr.as_base_exp()
    sb, se = pattern.as_base_exp()
    if sb.is_Symbol and se.is_Integer and _nonzero(expr):
        if e.is_rational:
            return matches(sb, b**(e/se), repl_dict)
        return matches(sb, expr**(1/se), repl_dict)
    d = dict(repl_dict)
    d = matches(pattern.base, b, d)
    if d is None:
        return None
    d = matches(pattern.exp.xreplace(d), e, d)
    if d is None:
        return _basic(pattern, expr, repl_dict, old)
    return d


def _new_rawargs(pattern, args):
    S, Add, Mul = _api()[:3]
    if not args:
        return S.Zero if pattern.is_Add else S.One
    if len(args) == 1:
        return args[0]
    return (Add if pattern.is_Add else Mul)(*args, evaluate=False)


def _make_args(pattern, expr):
    Add, Mul = _api()[1:3]
    return (Add if pattern.is_Add else Mul).make_args(expr)


def _add_simple(pattern, expr, repl_dict):
    # (w+3).matches(x+5) -> {w: x+2}
    coeff, terms = pattern.as_coeff_add()
    if len(terms) == 1:
        return matches(terms[0], expr - coeff, repl_dict)
    return None


def _mul_simple(pattern, expr, repl_dict):
    # (w*3).matches(x*5) -> {w: x*5/3}
    Mul = _api()[2]
    coeff, terms = pattern.as_coeff_Mul()
    terms = Mul.make_args(terms)
    if len(terms) == 1:
        return matches(terms[0], _mul_inverse(expr, coeff), repl_dict)
    return None


def _add_inverse(lhs, rhs):
    S, Add, Mul, Dummy, ordered, signsimp, collect = _api()
    inf = (S.Infinity, S.NegativeInfinity)
    if lhs.has(*inf) or rhs.has(*inf):
        oo = Dummy('oo')
        reps = {S.Infinity: oo, S.NegativeInfinity: -oo}
        ireps = {v: k for k, v in reps.items()}
        eq = lhs.xreplace(reps) - rhs.xreplace(reps)
        if eq.has(oo):
            eq = eq.replace(lambda x: x.is_Pow and x.base == oo, lambda x: x.base)
        rv = eq.xreplace(ireps)
    else:
        rv = lhs - rhs
    srv = signsimp(rv)
    return srv if srv.is_Number else rv


def _mul_inverse(lhs, rhs):
    S, Add, Mul, Dummy, ordered, signsimp, collect = _api()
    if lhs == rhs:
        return S.One

    def check(left, right):
        if left.is_Float and right.is_comparable:
            return left.__add__(0) == right.evalf().__add__(0)
        return False

    if check(lhs, rhs) or check(rhs, lhs):
        return S.One
    if any(i.is_Pow or i.is_Mul for i in (lhs, rhs)):
        # a literal I does not combine with a power of -1
        d = Dummy('I')
        _i = {S.ImaginaryUnit: d}
        i_ = {d: S.ImaginaryUnit}
        a = lhs.xreplace(_i).as_powers_dict()
        b = rhs.xreplace(_i).as_powers_dict()
        blen = len(b)
        for bi in tuple(b.keys()):
            if bi in a:
                a[bi] -= b.pop(bi)
                if not a[bi]:
                    a.pop(bi)
        if len(b) != blen:
            lhs = Mul(*[k**v for k, v in a.items()]).xreplace(i_)
            rhs = Mul(*[k**v for k, v in b.items()]).xreplace(i_)
    rv = lhs/rhs
    srv = signsimp(rv)
    return srv if srv.is_Number else rv


def _numeric_first(term):
    if term.is_Mul and term.args[0].is_Number:
        return term.args[0]
    return 0


def _commutative(pattern, expr, repl_dict, old):
    S, Add, Mul, Dummy, ordered, signsimp, collect = _api()
    if _is_expr(pattern) and not _is_expr(expr):
        return None
    if repl_dict is None:
        repl_dict = {}
    if pattern == expr:
        return repl_dict
    if pattern.is_Add:
        d = _add_simple(pattern, expr, repl_dict)
    else:
        d = _mul_simple(pattern, expr, repl_dict)
    if d is not None:
        return d

    # eliminate the exact part: (2+a+w1+w2).matches(expr) -> (w1+w2).matches(expr-a-2)
    wild_part, exact_part = [], []
    for p in pattern.args:
        (wild_part if _has_wild(p) and not expr.has(p) else exact_part).append(p)
    if not exact_part:
        wild_part = list(ordered(wild_part))
        if pattern.is_Add:
            wild_part = sorted(wild_part, key=_numeric_first)
    else:
        exact = _new_rawargs(pattern, exact_part)
        free = expr.free_symbols
        if free and (exact.free_symbols - free):
            # symbols of the exact part are absent from the expression
            return None
        if pattern.is_Add:
            newexpr = _add_inverse(expr, exact)
        else:
            newexpr = _mul_inverse(expr, exact)
        if not old and (expr.is_Add or expr.is_Mul):
            check = newexpr
            if _coeff_isneg(check):
                check = -check
            if check.count_ops() > expr.count_ops():
                return None
        newpattern = _new_rawargs(pattern, wild_part)
        return matches(newpattern, newexpr, repl_dict)

    identity = S.Zero if pattern.is_Add else S.One
    i = 0
    saw = set()
    while expr not in saw:
        saw.add(expr)
        args = tuple(ordered(_make_args(pattern, expr)))
        if pattern.is_Add and expr.is_Add:
            args = tuple(sorted(args, key=_numeric_first))
        expr_list = (identity,) + args
        for last_op in reversed(expr_list):
            for w in reversed(wild_part):
                d1 = matches(w, last_op, repl_dict)
                if d1 is not None:
                    d2 = matches(pattern.xreplace(d1), expr, d1)
                    if d2 is not None:
                        return d2

        if i == 0:
            if pattern.is_Mul:
                # make e**i look like Mul
                if expr.is_Pow and expr.exp.is_Integer:
                    if expr.exp > 0:
                        expr = Mul(*[expr.base, expr.base**(expr.exp - 1)], evaluate=False)
                    else:
                        expr = Mul(*[1/expr.base, expr.base**(expr.exp + 1)], evaluate=False)
                    i += 1
                    continue
            elif pattern.is_Add:
                # make i*e look like Add
                c, e = expr.as_coeff_Mul()
                if abs(c) > 1:
                    if c > 0:
                        expr = Add(*[e, (c - 1)*e], evaluate=False)
                    else:
                        expr = Add(*[-e, (c + 1)*e], evaluate=False)
                    i += 1
                    continue

                # try collection on non-Wild symbols
                was = expr
                did = set()
                for w in reversed(wild_part):
                    c = Mul(*[f for f in Mul.make_args(w) if not _has_wild(f)])
                    free = c.free_symbols - did
                    if free:
                        did.update(free)
                        expr = collect(expr, free)
                if expr != was:
                    continue

        break
    return None


def _mul(pattern, expr, repl_dict, old):
    from . import sympify
    Mul = _api()[2]
    expr = sympify(expr)
    if pattern.is_commutative and expr.is_commutative:
        return _commutative(pattern, expr, repl_dict, old)
    if pattern.is_commutative is not expr.is_commutative:
        return None

    # both are noncommutative
    c1, nc1 = pattern.args_cnc()
    c2, nc2 = expr.args_cnc()
    c1, c2 = [c or [1] for c in [c1, c2]]
    comm_mul_self = Mul(*c1)
    comm_mul_expr = Mul(*c2)
    repl_dict = matches(comm_mul_self, comm_mul_expr, repl_dict, old)
    if not repl_dict and c1 != c2:
        return None
    nc1 = _expand_pows(nc1)
    nc2 = _expand_pows(nc2)
    repl_dict = _noncomm(nc1, nc2, repl_dict)
    return repl_dict or None


def _expand_pows(arg_list):
    new_args = []
    for arg in arg_list:
        if arg.is_Pow and arg.exp > 0:
            new_args.extend([arg.base] * int(arg.exp))
        else:
            new_args.append(arg)
    return new_args


def _noncomm(nodes, targets, repl_dict=None):
    repl_dict = {} if repl_dict is None else dict(repl_dict)
    agenda = []
    state = (0, 0)
    node_ind, target_ind = state
    wildcard_dict = {}
    while target_ind < len(targets) and node_ind < len(nodes):
        node = nodes[node_ind]
        if _is_wild(node):
            _add_wildcard(wildcard_dict, state)
        states_matches = _new_states(wildcard_dict, state, nodes, targets)
        if states_matches:
            new_states, new_matches = states_matches
            agenda.extend(new_states)
            if new_matches:
                for match in new_matches:
                    repl_dict[match] = new_matches[match]
        if not agenda:
            return None
        state = agenda.pop()
        node_ind, target_ind = state
    return repl_dict


def _add_wildcard(dictionary, state):
    node_ind, target_ind = state
    if node_ind in dictionary:
        begin, end = dictionary[node_ind]
        dictionary[node_ind] = (begin, target_ind)
    else:
        dictionary[node_ind] = (target_ind, target_ind)


def _new_states(dictionary, state, nodes, targets):
    node_ind, target_ind = state
    node = nodes[node_ind]
    target = targets[target_ind]
    if target_ind >= len(targets) - 1 and node_ind < len(nodes) - 1:
        return None
    if _is_wild(node):
        match_attempt = _match_wilds(dictionary, node_ind, nodes, targets)
        if match_attempt:
            other_node_inds = [ind for ind in dictionary if nodes[ind] == node]
            for ind in other_node_inds:
                other_begin, other_end = dictionary[ind]
                curr_begin, curr_end = dictionary[node_ind]
                other_targets = targets[other_begin:other_end + 1]
                current_targets = targets[curr_begin:curr_end + 1]
                for curr, other in zip(current_targets, other_targets):
                    if curr != other:
                        return None
            new_state = [(node_ind, target_ind + 1)]
            if node_ind < len(nodes) - 1:
                new_state.append((node_ind + 1, target_ind + 1))
            return new_state, match_attempt
    else:
        if node_ind >= len(nodes) - 1 and target_ind < len(targets) - 1:
            return None
        match_attempt = matches(node, target)
        if match_attempt:
            return [(node_ind + 1, target_ind + 1)], match_attempt
        elif node == target:
            return [(node_ind + 1, target_ind + 1)], None
        else:
            return None
    return None


def _match_wilds(dictionary, wildcard_ind, nodes, targets):
    Mul = _api()[2]
    wildcard = nodes[wildcard_ind]
    begin, end = dictionary[wildcard_ind]
    terms = targets[begin:end + 1]
    mult = Mul(*terms) if len(terms) > 1 else terms[0]
    return matches(wildcard, mult)
