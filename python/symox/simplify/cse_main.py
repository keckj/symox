"""Compatible imports backed by existing native implementations."""

from collections import defaultdict
import symox as S

# FuncArgTracker, Unevaluated, match_common_args and opt_cse are ported from
# SymPy 1.14 ``sympy/simplify/cse_main.py`` (BSD-3-Clause, see LICENSE-SYMPY).
from symox import Basic, Mul, Add, Pow, MatrixExpr, MatMul, MatAdd, MatPow, Inverse
from symox.core.containers import OrderedSet
from symox.utilities.iterables import iterable, topological_sort, numbered_symbols, sift  # noqa: F401
from symox.matrices.expressions.matexpr import MatrixElement  # noqa: F401
from symox import S as _S

from symox.core.exprtools import factor_terms  # noqa: E402
from symox.simplify import cse_opts  # noqa: E402

# (preprocessor, postprocessor) pairs which are commonly useful (SymPy's
# ``basic_optimizations``).
basic_optimizations = [(cse_opts.sub_pre, cse_opts.sub_post),
                       (factor_terms, None)]


class FuncArgTracker:

    def __init__(self, funcs):
        # To minimize the number of symbolic comparisons, all function arguments
        # get assigned a value number.
        self.value_numbers = {}
        self.value_number_to_value = []

        # Both of these maps use integer indices for arguments / functions.
        self.arg_to_funcset = []
        self.func_to_argset = []

        for func_i, func in enumerate(funcs):
            func_argset = OrderedSet()

            for func_arg in func.args:
                arg_number = self.get_or_add_value_number(func_arg)
                func_argset.add(arg_number)
                self.arg_to_funcset[arg_number].add(func_i)

            self.func_to_argset.append(func_argset)

    def get_args_in_value_order(self, argset):
        return [self.value_number_to_value[argn] for argn in sorted(argset)]

    def get_or_add_value_number(self, value):
        nvalues = len(self.value_numbers)
        value_number = self.value_numbers.setdefault(value, nvalues)
        if value_number == nvalues:
            self.value_number_to_value.append(value)
            self.arg_to_funcset.append(OrderedSet())
        return value_number

    def stop_arg_tracking(self, func_i):
        for arg in self.func_to_argset[func_i]:
            self.arg_to_funcset[arg].remove(func_i)


    def get_common_arg_candidates(self, argset, min_func_i=0):
        count_map = defaultdict(lambda: 0)
        if not argset:
            return count_map

        funcsets = [self.arg_to_funcset[arg] for arg in argset]
        # As an optimization below, we handle the largest funcset separately from
        # the others.
        largest_funcset = max(funcsets, key=len)

        for funcset in funcsets:
            if largest_funcset is funcset:
                continue
            for func_i in funcset:
                if func_i >= min_func_i:
                    count_map[func_i] += 1

        # We pick the smaller of the two containers (count_map, largest_funcset)
        # to iterate over to reduce the number of iterations needed.
        (smaller_funcs_container,
         larger_funcs_container) = sorted(
                 [largest_funcset, count_map],
                 key=len)

        for func_i in smaller_funcs_container:
            # Not already in count_map? It can't possibly be in the output, so
            # skip it.
            if count_map[func_i] < 1:
                continue

            if func_i in larger_funcs_container:
                count_map[func_i] += 1

        return {k: v for k, v in count_map.items() if v >= 2}

    def get_subset_candidates(self, argset, restrict_to_funcset=None):
        iarg = iter(argset)

        indices = OrderedSet(
            fi for fi in self.arg_to_funcset[next(iarg)])

        if restrict_to_funcset is not None:
            indices &= restrict_to_funcset

        for arg in iarg:
            indices &= self.arg_to_funcset[arg]

        return indices

    def update_func_argset(self, func_i, new_argset):
        new_args = OrderedSet(new_argset)
        old_args = self.func_to_argset[func_i]

        for deleted_arg in old_args - new_args:
            self.arg_to_funcset[deleted_arg].remove(func_i)
        for added_arg in new_args - old_args:
            self.arg_to_funcset[added_arg].add(func_i)

        self.func_to_argset[func_i].clear()
        self.func_to_argset[func_i].update(new_args)


class Unevaluated:

    def __init__(self, func, args):
        self.func = func
        self.args = args

    def __str__(self):
        return "Uneval<{}>({})".format(
                self.func, ", ".join(str(a) for a in self.args))

    def _symox_(self):
        return self.as_unevaluated_basic()

    def as_unevaluated_basic(self):
        return self.func(*self.args, evaluate=False)

    @property
    def free_symbols(self):
        return set().union(*[a.free_symbols for a in self.args])

    __repr__ = __str__


def match_common_args(func_class, funcs, opt_subs):

    # Sort to ensure that whole-function subexpressions come before the items
    # that use them.
    funcs = sorted(funcs, key=lambda f: len(f.args))
    arg_tracker = FuncArgTracker(funcs)

    changed = OrderedSet()

    for i in range(len(funcs)):
        common_arg_candidates_counts = arg_tracker.get_common_arg_candidates(
                arg_tracker.func_to_argset[i], min_func_i=i + 1)

        # Sort the candidates in order of match size.
        # This makes us try combining smaller matches first.
        common_arg_candidates = OrderedSet(sorted(
                common_arg_candidates_counts.keys(),
                key=lambda k: (common_arg_candidates_counts[k], k)))

        while common_arg_candidates:
            j = common_arg_candidates.pop(last=False)

            com_args = arg_tracker.func_to_argset[i].intersection(
                    arg_tracker.func_to_argset[j])

            if len(com_args) <= 1:
                # This may happen if a set of common arguments was already
                # combined in a previous iteration.
                continue

            # For all sets, replace the common symbols by the function
            # over them, to allow recursive matches.

            diff_i = arg_tracker.func_to_argset[i].difference(com_args)
            if diff_i:
                # com_func needs to be unevaluated to allow for recursive matches.
                com_func = Unevaluated(
                        func_class, arg_tracker.get_args_in_value_order(com_args))
                com_func_number = arg_tracker.get_or_add_value_number(com_func)
                arg_tracker.update_func_argset(i, diff_i | OrderedSet([com_func_number]))
                changed.add(i)
            else:
                # Treat the whole expression as a CSE.
                #
                # The reason this needs to be done is somewhat subtle. Within
                # tree_cse(), to_eliminate only contains expressions that are
                # seen more than once. The problem is unevaluated expressions
                # do not compare equal to the evaluated equivalent. So
                # tree_cse() won't mark funcs[i] as a CSE if we use an
                # unevaluated version.
                com_func_number = arg_tracker.get_or_add_value_number(funcs[i])

            diff_j = arg_tracker.func_to_argset[j].difference(com_args)
            arg_tracker.update_func_argset(j, diff_j | OrderedSet([com_func_number]))
            changed.add(j)

            for k in arg_tracker.get_subset_candidates(
                    com_args, common_arg_candidates):
                diff_k = arg_tracker.func_to_argset[k].difference(com_args)
                arg_tracker.update_func_argset(k, diff_k | OrderedSet([com_func_number]))
                changed.add(k)

        if i in changed:
            opt_subs[funcs[i]] = Unevaluated(func_class,
                arg_tracker.get_args_in_value_order(arg_tracker.func_to_argset[i]))

        arg_tracker.stop_arg_tracking(i)


def opt_cse(exprs, order='canonical'):
    opt_subs = {}

    adds = OrderedSet()
    muls = OrderedSet()

    seen_subexp = set()
    collapsible_subexp = set()

    def _find_opts(expr):

        if not isinstance(expr, (Basic, Unevaluated)):
            return

        if getattr(expr, 'is_Atom', False) or getattr(expr, 'is_Order', False):
            return

        if iterable(expr):
            list(map(_find_opts, expr))
            return

        if expr in seen_subexp:
            return expr
        seen_subexp.add(expr)

        list(map(_find_opts, expr.args))

        if not isinstance(expr, MatrixExpr) and expr.could_extract_minus_sign():
            # XXX -expr does not always work rigorously for some expressions
            # containing UnevaluatedExpr.
            # https://github.com/sympy/sympy/issues/24818
            if isinstance(expr, Add):
                neg_expr = Add(*(-i for i in expr.args))
            else:
                neg_expr = -expr

            if not neg_expr.is_Atom:
                opt_subs[expr] = Unevaluated(Mul, (_S.NegativeOne, neg_expr))
                seen_subexp.add(neg_expr)
                expr = neg_expr

        if isinstance(expr, (Mul, MatMul)):
            if len(expr.args) == 1:
                collapsible_subexp.add(expr)
            else:
                muls.add(expr)

        elif isinstance(expr, (Add, MatAdd)):
            if len(expr.args) == 1:
                collapsible_subexp.add(expr)
            else:
                adds.add(expr)

        elif isinstance(expr, Inverse):
            # Do not want to treat `Inverse` as a `MatPow`
            pass

        elif isinstance(expr, (Pow, MatPow)):
            base, exp = expr.base, expr.exp
            if exp.could_extract_minus_sign():
                opt_subs[expr] = Unevaluated(Pow, (Pow(base, -exp), -1))

    for e in exprs:
        if isinstance(e, (Basic, Unevaluated)):
            _find_opts(e)

    # Handle collapsing of multinary operations with single arguments
    edges = [(s, s.args[0]) for s in collapsible_subexp
             if s.args[0] in collapsible_subexp]
    for e in reversed(topological_sort((collapsible_subexp, edges))):
        opt_subs[e] = opt_subs.get(e.args[0], e.args[0])

    # split muls into commutative
    commutative_muls = OrderedSet()
    for m in muls:
        c, nc = m.args_cnc(cset=False)
        if c:
            c_mul = m.func(*c)
            if nc:
                if c_mul == 1:
                    new_obj = m.func(*nc)
                else:
                    if isinstance(m, MatMul):
                        new_obj = m.func(c_mul, *nc, evaluate=False)
                    else:
                        new_obj = m.func(c_mul, m.func(*nc), evaluate=False)
                opt_subs[m] = new_obj
            if len(c) > 1:
                commutative_muls.add(c_mul)

    match_common_args(Add, adds, opt_subs)
    match_common_args(Mul, commutative_muls, opt_subs)

    return opt_subs


def reps_toposort(r):
    pending=[(symbol,S.sympify(expression)) for symbol,expression in r];keys={symbol for symbol,expression in pending};done=set();ordered=[]
    while pending:
        ready=next((i for i,(symbol,expression) in enumerate(pending) if not (expression.free_symbols&keys)-done),None)
        if ready is None:raise ValueError('cyclic common-subexpression replacements')
        symbol,expression=pending.pop(ready);ordered.append((symbol,expression));done.add(symbol)
    return ordered


def cse_separate(r,e):
    replacements=list(r);expressions=[]
    for expression in e:
        if getattr(expression,'class_name','')=='Equality' and expression.lhs.is_Symbol:replacements.append((expression.lhs,expression.rhs))
        else:expressions.append(expression)
    return [reps_toposort(replacements),expressions]


def cse_release_variables(r, e):
    # Verbatim port of the reference cse_release_variables.
    if not r:
        return r, e
    s, p = zip(*r)
    esyms = S.symbols('_:%d' % len(e))
    syms = list(esyms)
    s = list(s)
    in_use = set(s)
    p = list(p)
    e = [(e[i], syms[i]) for i in range(len(e))]
    e, syms = zip(*sorted(e,
        key=lambda x: -sum(p[s.index(i)].count_ops()
        for i in x[0].free_symbols & in_use)))
    syms = list(syms)
    p += e
    rv = []
    i = len(p) - 1
    while i >= 0:
        _p = p.pop()
        c = in_use & _p.free_symbols
        if c:
            rv.extend([(s, None) for s in sorted(c, key=str)])
        if i >= len(r):
            rv.append((syms.pop(), _p))
        else:
            rv.append((s[i], _p))
        in_use -= c
        i -= 1
    rv.reverse()
    return rv, esyms


def preprocess_for_cse(expr,optimizations):
    for before,after in optimizations:
        if before is not None:expr=before(expr)
    return expr

def postprocess_for_cse(expr,optimizations):
    for before,after in reversed(optimizations):
        if after is not None:expr=after(expr)
    return expr


def tree_cse(exprs, symbols, opt_subs=None, order='canonical', ignore=()):
    """Perform raw CSE on expression tree, taking opt_subs into account.

    Port of SymPy 1.14 ``tree_cse``.
    """
    from symox.polys.rootoftools import RootOf
    from symox.matrices.expressions import MatrixSymbol
    from symox.matrices.expressions.matexpr import MatrixElement
    from symox.core.sorting import ordered

    if opt_subs is None:
        opt_subs = {}

    ## Find repeated sub-expressions

    to_eliminate = set()

    seen_subexp = set()
    excluded_symbols = set()

    def _find_repeated(expr):
        if not isinstance(expr, (Basic, Unevaluated)):
            return

        if isinstance(expr, RootOf):
            return

        if isinstance(expr, Basic) and (
                getattr(expr, 'is_Atom', False) or
                getattr(expr, 'is_Order', False) or
                isinstance(expr, (MatrixSymbol, MatrixElement))):
            if expr.is_Symbol:
                excluded_symbols.add(expr.name)
            return

        if iterable(expr):
            args = expr

        else:
            if expr in seen_subexp:
                for ign in ignore:
                    if ign in expr.free_symbols:
                        break
                else:
                    to_eliminate.add(expr)
                    return

            seen_subexp.add(expr)

            if expr in opt_subs:
                expr = opt_subs[expr]

            args = expr.flat() if _is_concrete_matrix(expr) else expr.args

        list(map(_find_repeated, args))

    for e in exprs:
        if isinstance(e, Basic):
            _find_repeated(e)

    ## Rebuild tree

    # Remove symbols from the generator that conflict with names in the expressions.
    symbols = (_ for _ in symbols if _.name not in excluded_symbols)

    replacements = []

    subs = {}

    def _rebuild(expr):
        if not isinstance(expr, (Basic, Unevaluated)):
            return expr

        if _is_concrete_matrix(expr):
            # The reference ImmutableMatrix is a Basic with args
            # (rows, cols, Tuple(*entries)); rebuild entrywise.
            return type(expr)(expr.rows, expr.cols,
                              [_rebuild(arg) for arg in expr.flat()])

        if not expr.args:
            return expr

        if iterable(expr):
            new_args = [_rebuild(arg) for arg in expr.args]
            return expr.func(*new_args)

        if expr in subs:
            return subs[expr]

        orig_expr = expr
        if expr in opt_subs:
            expr = opt_subs[expr]

        if _is_concrete_matrix(expr):
            new_expr = _rebuild(expr)
            if orig_expr in to_eliminate:
                return _eliminate(orig_expr, new_expr)
            return new_expr

        # If enabled, parse Muls and Adds arguments by order to ensure
        # replacement order independent from hashes
        if order != 'none':
            if isinstance(expr, (Mul, MatMul)):
                c, nc = expr.args_cnc()
                if c == [1]:
                    args = nc
                else:
                    args = list(ordered(c)) + nc
            elif isinstance(expr, (Add, MatAdd)):
                args = list(ordered(expr.args))
            else:
                args = expr.args
        else:
            args = expr.args

        new_args = list(map(_rebuild, args))
        if isinstance(expr, Unevaluated) or new_args != args:
            new_expr = expr.func(*new_args)
        else:
            new_expr = expr

        if orig_expr in to_eliminate:
            return _eliminate(orig_expr, new_expr)
        else:
            return new_expr

    def _eliminate(orig_expr, new_expr):
        try:
            sym = next(symbols)
        except StopIteration:
            raise ValueError("Symbols iterator ran out of symbols.")

        if isinstance(orig_expr, MatrixExpr):
            sym = MatrixSymbol(sym.name, orig_expr.rows,
                orig_expr.cols)

        subs[orig_expr] = sym
        replacements.append((sym, new_expr))
        return sym

    reduced_exprs = []
    for e in exprs:
        if isinstance(e, Basic):
            reduced_e = _rebuild(e)
        else:
            reduced_e = e
        reduced_exprs.append(reduced_e)
    return replacements, reduced_exprs


def _is_concrete_matrix(expr):
    from symox import _core
    return isinstance(expr, _core.Matrix)


def _basic_optimizations():
    return basic_optimizations


def cse(exprs, symbols=None, optimizations=None, postprocess=None,
        order='canonical', ignore=(), list=True):
    """Perform common subexpression elimination on an expression.

    Port of SymPy 1.14 ``cse``: ``opt_cse`` finds the Add/Mul/Pow
    rewrites, ``tree_cse`` the repeated subtrees.
    """
    from symox.core.containers import Tuple
    from symox.matrices import (MatrixBase, Matrix, ImmutableMatrix,
                                SparseMatrix, ImmutableSparseMatrix)
    from symox.utilities.iterables import numbered_symbols
    from symox import Symbol, sympify

    if not list:
        return _cse_homogeneous(exprs,
            symbols=symbols, optimizations=optimizations,
            postprocess=postprocess, order=order, ignore=ignore)

    if isinstance(exprs, (int, float)):
        exprs = sympify(exprs)

    # Handle the case if just one expression was passed.
    if isinstance(exprs, (Basic, MatrixBase)) or getattr(exprs, '_symox_matrixexpr', False):
        exprs = [exprs]

    copy = exprs
    temp = []
    for e in exprs:
        if isinstance(e, (Matrix, ImmutableMatrix)):
            temp.append(Tuple(*e.flat()))
        elif isinstance(e, (SparseMatrix, ImmutableSparseMatrix)):
            temp.append(Tuple(*e.todok().items()))
        else:
            temp.append(e)
    exprs = temp
    del temp

    if optimizations is None:
        optimizations = []
    elif optimizations == 'basic':
        optimizations = _basic_optimizations()

    # Preprocess the expressions to give us better optimization opportunities.
    reduced_exprs = [preprocess_for_cse(e, optimizations) for e in exprs]

    if symbols is None:
        symbols = numbered_symbols(cls=Symbol)
    else:
        # In case we get passed an iterable with an __iter__ method instead of
        # an actual iterator.
        symbols = iter(symbols)

    # Find other optimization opportunities.
    opt_subs = opt_cse(reduced_exprs, order)

    # Main CSE algorithm.
    replacements, reduced_exprs = tree_cse(reduced_exprs, symbols, opt_subs,
                                           order, ignore)

    # Postprocess the expressions to return the expressions to canonical form.
    exprs = copy
    replacements = [(sym, postprocess_for_cse(subtree, optimizations))
                    for sym, subtree in replacements]
    reduced_exprs = [postprocess_for_cse(e, optimizations)
                     for e in reduced_exprs]

    # Get the matrices back
    for i, e in enumerate(exprs):
        if isinstance(e, (Matrix, ImmutableMatrix)):
            reduced_exprs[i] = Matrix(e.rows, e.cols, reduced_exprs[i])
            if isinstance(e, ImmutableMatrix):
                reduced_exprs[i] = reduced_exprs[i].as_immutable()
        elif isinstance(e, (SparseMatrix, ImmutableSparseMatrix)):
            m = SparseMatrix(e.rows, e.cols, {})
            for k, v in reduced_exprs[i]:
                m[tuple(int(j) for j in k)] = v
            if isinstance(e, ImmutableSparseMatrix):
                m = m.as_immutable()
            reduced_exprs[i] = m

    if postprocess is None:
        return replacements, reduced_exprs

    return postprocess(replacements, reduced_exprs)


def _cse_homogeneous(exprs, **kwargs):
    """Same as ``cse`` but the ``reduced_exprs`` are returned with the
    same type as ``exprs`` or a sympified version of the same."""
    from symox import sympify
    if isinstance(exprs, str):
        replacements, reduced_exprs = _cse_homogeneous(
            sympify(exprs), **kwargs)
        return replacements, repr(reduced_exprs)
    if isinstance(exprs, (list, tuple, set)):
        replacements, reduced_exprs = cse(exprs, **kwargs)
        return replacements, type(exprs)(reduced_exprs)
    if isinstance(exprs, dict):
        keys = list(exprs.keys())
        replacements, values = cse([exprs[k] for k in keys], **kwargs)
        reduced_exprs = dict(zip(keys, values))
        return replacements, reduced_exprs

    try:
        replacements, (reduced_exprs,) = cse(exprs, **kwargs)
    except TypeError:  # For example 'mpf' objects
        return [], exprs
    else:
        return replacements, reduced_exprs


__all__ = ['cse', 'numbered_symbols', 'basic_optimizations', 'FuncArgTracker','Unevaluated','cse_release_variables','cse_separate','match_common_args','opt_cse',
 'postprocess_for_cse','preprocess_for_cse','reps_toposort','tree_cse','cse','numbered_symbols']
