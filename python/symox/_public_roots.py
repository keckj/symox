"""Indexed roots of univariate polynomials and symmetric sums over roots.

``CRootOf`` follows the reference's ``sympy.polys.rootoftools`` step by step:
the polynomial is preprocessed to an integer basis, factored, and its roots
are located by certified isolation (continued fractions for the real roots,
rectangles with edge Sturm sequences for the others). The isolating
intervals live in caches keyed by the irreducible factor and are refined in
the same order as the reference, so indices, assumptions and rational
approximations agree with it. A root object is a native expression with the
reference's arguments ``(expr, Integer(index))``.
"""
import operator
from itertools import chain

from . import _core
from ._compat import _ExprType

_sums = {}


def _polynomial(f,x=None,rational=True):
    from . import Poly,Symbol,sympify
    if isinstance(f,Poly):
        if x is None:x = f.gens[0]
        f = f.as_expr()
    f = sympify(f)
    if x is None:
        symbols = sorted(f.free_symbols,key=str)
        if len(symbols) != 1:raise ValueError('a univariate polynomial is required')
        x = symbols[0]
    x = sympify(x)
    polynomial = Poly(f,x)
    if len(polynomial.gens) != 1 or polynomial.degree() < 1:raise ValueError('polynomial degree must be positive')
    if rational:
        _require_rational(polynomial)
    return polynomial


def _require_rational(polynomial):
    if any(c.is_Rational is not True for c in polynomial.all_coeffs()):
        raise NotImplementedError('algebraic root isolation currently requires rational coefficients')


class _pure_key_dict:
    """A dictionary keyed by ``PurePoly``; a univariate expression is
    converted to one."""

    def __init__(self):
        self._dict = {}

    def __getitem__(self, k):
        from . import Expr, PurePoly
        if not isinstance(k, PurePoly):
            if not (isinstance(k, Expr) and len(k.free_symbols) == 1):
                raise KeyError
            k = PurePoly(k, expand=False)
        return self._dict[k]

    def __setitem__(self, k, v):
        from . import Expr, PurePoly
        if not isinstance(k, PurePoly):
            if not (isinstance(k, Expr) and len(k.free_symbols) == 1):
                raise ValueError('expecting univariate expression')
            k = PurePoly(k, expand=False)
        self._dict[k] = v

    def __contains__(self, k):
        try:
            self[k]
            return True
        except KeyError:
            return False


_reals_cache = _pure_key_dict()
_complexes_cache = _pure_key_dict()


def _pure_factors(poly):
    from . import PurePoly
    _, factors = poly.factor_list()
    return [(PurePoly(f, expand=False), m) for f, m in factors]


def _imag_count_of_factor(f):
    """Return the number of imaginary roots for irreducible
    univariate polynomial ``f``.
    """
    from . import Dummy, I, Poly, oo
    terms = [(i, j) for (i,), j in f.terms()]
    if any(i % 2 for i, j in terms):
        return 0
    # update signs
    even = [(i, I**i*j) for i, j in terms]
    even = Poly.from_dict(dict(even), Dummy('x'))
    return int(even.count_roots(-oo, oo))


def _ordered(seq):
    from .core.sorting import ordered
    return ordered(seq)


def rootof(f, x, index=None, radicals=True, expand=True):
    """An indexed root of a univariate polynomial.

    Returns either a :obj:`ComplexRootOf` object or an explicit
    expression involving radicals.
    """
    return CRootOf(f, x, index=index, radicals=radicals, expand=expand)


class RootOf(metaclass=_ExprType):
    """Represents a root of a univariate polynomial."""
    __module__ = 'symox.polys.rootoftools'
    _matches = staticmethod(lambda expr:expr.class_name == 'ComplexRootOf')

    def __new__(cls, f, x, index=None, radicals=True, expand=True):
        return rootof(f, x, index=index, radicals=radicals, expand=expand)


class ComplexRootOf(metaclass=_ExprType):
    """Represents an indexed complex root of a polynomial."""
    __module__ = 'symox.polys.rootoftools'
    _matches = staticmethod(lambda expr:expr.class_name == 'ComplexRootOf')

    def __new__(cls, f, x, index=None, radicals=False, expand=True):
        """ Construct an indexed complex root of a polynomial."""
        from . import PurePoly, sympify
        from .polys.polyerrors import PolynomialError
        from .polys.polyroots import preprocess_roots
        x = sympify(x)

        if index is None and x.is_Integer:
            x, index = None, x
        elif index is not None:
            index = sympify(index)

        if index is not None and index.is_Integer:
            index = int(index)
        else:
            raise ValueError("expected an integer root index, got %s" % index)

        if x is None:
            poly = PurePoly(f, greedy=False, expand=expand)
        else:
            poly = PurePoly(f, x, greedy=False, expand=expand)

        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are allowed")

        if not poly.gen.is_Symbol:
            # PurePoly(sin(x) + 1) == PurePoly(x + 1) but the roots differ.
            raise PolynomialError("generator must be a Symbol")

        degree = poly.degree()

        if degree <= 0:
            raise PolynomialError("Cannot construct CRootOf object for %s" % f)

        if index < -degree or index >= degree:
            raise IndexError("root index out of [%d, %d] range, got %d" %
                             (-degree, degree - 1, index))
        elif index < 0:
            index += degree

        dom = poly.get_domain()

        if not dom.is_Exact:
            poly = poly.to_exact()

        roots = cls._roots_trivial(poly, radicals)

        if roots is not None:
            return roots[index]

        coeff, poly = preprocess_roots(poly)
        dom = poly.get_domain()

        if not dom.is_ZZ:
            raise NotImplementedError("CRootOf is not supported over %s" % dom)

        root = cls._indexed_root(poly, index, lazy=True)
        return coeff * cls._postprocess_root(root, radicals)

    @classmethod
    def _new(cls, poly, index):
        """Construct new ``CRootOf`` object from raw data. """
        from . import PurePoly
        poly = PurePoly(poly)
        # The native object carries its assumptions: locate the root.
        if poly not in _reals_cache:
            cls._indexed_root(poly, index)
        reals = _reals_cache[poly]
        if index < len(reals):
            kind = (True, bool(reals[index].a < 0), False)
        else:
            if poly not in _complexes_cache:
                cls._indexed_root(poly, index)
            ivl = _complexes_cache[poly][index - len(reals)]
            kind = (False, False, bool(ivl.ax*ivl.bx <= 0))
        return _core.complex_root_node(poly.as_expr(), index, kind)

    @classmethod
    def real_roots(cls, poly, radicals=True):
        """Get real roots of a polynomial. """
        return cls._get_roots("_real_roots", poly, radicals)

    @classmethod
    def all_roots(cls, poly, radicals=True):
        """Get real and complex roots of a polynomial. """
        return cls._get_roots("_all_roots", poly, radicals)

    @classmethod
    def _get_reals_sqf(cls, currentfactor, use_cache=True):
        """Get real root isolating intervals for a square-free factor."""
        from .polys.rootisolation import dup_isolate_real_roots_sqf
        if use_cache and currentfactor in _reals_cache:
            real_part = _reals_cache[currentfactor]
        else:
            _reals_cache[currentfactor] = real_part = \
                dup_isolate_real_roots_sqf(
                    currentfactor.rep.to_list(), currentfactor.rep.dom, blackbox=True)

        return real_part

    @classmethod
    def _get_complexes_sqf(cls, currentfactor, use_cache=True):
        """Get complex root isolating intervals for a square-free factor."""
        from .polys.rootisolation import dup_isolate_complex_roots_sqf
        if use_cache and currentfactor in _complexes_cache:
            complex_part = _complexes_cache[currentfactor]
        else:
            _complexes_cache[currentfactor] = complex_part = \
                dup_isolate_complex_roots_sqf(
                currentfactor.rep.to_list(), currentfactor.rep.dom, blackbox=True)
        return complex_part

    @classmethod
    def _get_reals(cls, factors, use_cache=True):
        """Compute real root isolating intervals for a list of factors. """
        reals = []

        for currentfactor, k in factors:
            try:
                if not use_cache:
                    raise KeyError
                r = _reals_cache[currentfactor]
                reals.extend([(i, currentfactor, k) for i in r])
            except KeyError:
                real_part = cls._get_reals_sqf(currentfactor, use_cache)
                new = [(root, currentfactor, k) for root in real_part]
                reals.extend(new)

        reals = cls._reals_sorted(reals)
        return reals

    @classmethod
    def _get_complexes(cls, factors, use_cache=True):
        """Compute complex root isolating intervals for a list of factors. """
        complexes = []

        for currentfactor, k in _ordered(factors):
            try:
                if not use_cache:
                    raise KeyError
                c = _complexes_cache[currentfactor]
                complexes.extend([(i, currentfactor, k) for i in c])
            except KeyError:
                complex_part = cls._get_complexes_sqf(currentfactor, use_cache)
                new = [(root, currentfactor, k) for root in complex_part]
                complexes.extend(new)

        complexes = cls._complexes_sorted(complexes)
        return complexes

    @classmethod
    def _reals_sorted(cls, reals):
        """Make real isolating intervals disjoint and sort roots. """
        cache = {}

        for i, (u, f, k) in enumerate(reals):
            for j, (v, g, m) in enumerate(reals[i + 1:]):
                u, v = u.refine_disjoint(v)
                reals[i + j + 1] = (v, g, m)

            reals[i] = (u, f, k)

        reals = sorted(reals, key=lambda r: r[0].a)

        for root, currentfactor, _ in reals:
            if currentfactor in cache:
                cache[currentfactor].append(root)
            else:
                cache[currentfactor] = [root]

        for currentfactor, root in cache.items():
            _reals_cache[currentfactor] = root

        return reals

    @classmethod
    def _refine_imaginary(cls, complexes):
        from .utilities import sift
        sifted = sift(complexes, lambda c: c[1])
        complexes = []
        for f in _ordered(sifted):
            nimag = _imag_count_of_factor(f)
            if nimag == 0:
                # refine until xbounds are neg or pos
                for u, f, k in sifted[f]:
                    while u.ax*u.bx <= 0:
                        u = u._inner_refine()
                    complexes.append((u, f, k))
            else:
                # refine until all but nimag xbounds are neg or pos
                potential_imag = list(range(len(sifted[f])))
                while True:
                    assert len(potential_imag) > 1
                    for i in list(potential_imag):
                        u, f, k = sifted[f][i]
                        if u.ax*u.bx > 0:
                            potential_imag.remove(i)
                        elif u.ax != u.bx:
                            u = u._inner_refine()
                            sifted[f][i] = u, f, k
                    if len(potential_imag) == nimag:
                        break
                complexes.extend(sifted[f])
        return complexes

    @classmethod
    def _refine_complexes(cls, complexes):
        """return complexes such that no bounding rectangles of non-conjugate
        roots would intersect. In addition, assure that neither ay nor by is
        0 to guarantee that non-real roots are distinct from real roots in
        terms of the y-bounds.
        """
        # get the intervals pairwise-disjoint.
        # If rectangles were drawn around the coordinates of the bounding
        # rectangles, no rectangles would intersect after this procedure.
        for i, (u, f, k) in enumerate(complexes):
            for j, (v, g, m) in enumerate(complexes[i + 1:]):
                u, v = u.refine_disjoint(v)
                complexes[i + j + 1] = (v, g, m)

            complexes[i] = (u, f, k)

        # refine until the x-bounds are unambiguously positive or negative
        # for non-imaginary roots
        complexes = cls._refine_imaginary(complexes)

        # make sure that all y bounds are off the real axis
        # and on the same side of the axis
        for i, (u, f, k) in enumerate(complexes):
            while u.ay*u.by <= 0:
                u = u.refine()
            complexes[i] = u, f, k
        return complexes

    @classmethod
    def _complexes_sorted(cls, complexes):
        """Make complex isolating intervals disjoint and sort roots. """
        complexes = cls._refine_complexes(complexes)
        # XXX don't sort until you are sure that it is compatible
        # with the indexing method but assert that the desired state
        # is not broken
        C, F = 0, 1  # location of ComplexInterval and factor
        fs = {i[F] for i in complexes}
        for i in range(1, len(complexes)):
            if complexes[i][F] != complexes[i - 1][F]:
                # if this fails the factors of a root were not
                # contiguous because a discontinuity should only
                # happen once
                fs.remove(complexes[i - 1][F])
        for i, cmplx in enumerate(complexes):
            # negative im part (conj=True) comes before
            # positive im part (conj=False)
            assert cmplx[C].conj is (i % 2 == 0)

        # update cache
        cache = {}
        # -- collate
        for root, currentfactor, _ in complexes:
            cache.setdefault(currentfactor, []).append(root)
        # -- store
        for currentfactor, root in cache.items():
            _complexes_cache[currentfactor] = root

        return complexes

    @classmethod
    def _reals_index(cls, reals, index):
        """
        Map initial real root index to an index in a factor where
        the root belongs.
        """
        i = 0

        for j, (_, currentfactor, k) in enumerate(reals):
            if index < i + k:
                poly, index = currentfactor, 0

                for _, currentfactor, _ in reals[:j]:
                    if currentfactor == poly:
                        index += 1

                return poly, index
            else:
                i += k

    @classmethod
    def _complexes_index(cls, complexes, index):
        """
        Map initial complex root index to an index in a factor where
        the root belongs.
        """
        i = 0
        for j, (_, currentfactor, k) in enumerate(complexes):
            if index < i + k:
                poly, index = currentfactor, 0

                for _, currentfactor, _ in complexes[:j]:
                    if currentfactor == poly:
                        index += 1

                index += len(_reals_cache[poly])

                return poly, index
            else:
                i += k

    @classmethod
    def _count_roots(cls, roots):
        """Count the number of real or complex roots with multiplicities."""
        return sum(k for _, _, k in roots)

    @classmethod
    def _indexed_root(cls, poly, index, lazy=False):
        """Get a root of a composite polynomial by index. """
        factors = _pure_factors(poly)

        # If the given poly is already irreducible, then the index does not
        # need to be adjusted, and we can postpone the heavy lifting of
        # computing and refining isolating intervals until that is needed.
        if lazy and len(factors) == 1 and factors[0][1] == 1:
            return factors[0][0], index

        reals = cls._get_reals(factors)
        reals_count = cls._count_roots(reals)

        if index < reals_count:
            return cls._reals_index(reals, index)
        else:
            complexes = cls._get_complexes(factors)
            return cls._complexes_index(complexes, index - reals_count)

    @classmethod
    def _real_roots(cls, poly):
        """Get real roots of a composite polynomial. """
        factors = _pure_factors(poly)

        reals = cls._get_reals(factors)
        reals_count = cls._count_roots(reals)

        roots = []

        for index in range(0, reals_count):
            roots.append(cls._reals_index(reals, index))

        return roots

    @classmethod
    def _all_roots(cls, poly, use_cache=True):
        """Get real and complex roots of a composite polynomial. """
        factors = _pure_factors(poly)

        reals = cls._get_reals(factors, use_cache=use_cache)
        reals_count = cls._count_roots(reals)

        roots = []

        for index in range(0, reals_count):
            roots.append(cls._reals_index(reals, index))

        complexes = cls._get_complexes(factors, use_cache=use_cache)
        complexes_count = cls._count_roots(complexes)

        for index in range(0, complexes_count):
            roots.append(cls._complexes_index(complexes, index))

        return roots

    @classmethod
    def _roots_trivial(cls, poly, radicals):
        """Compute roots in linear, quadratic and binomial cases. """
        from .polys.polyroots import roots_binomial, roots_linear, roots_quadratic
        if poly.degree() == 1:
            return roots_linear(poly)

        if not radicals:
            return None

        if poly.degree() == 2:
            return roots_quadratic(poly)
        elif poly.length() == 2 and poly.TC():
            return roots_binomial(poly)
        else:
            return None

    @classmethod
    def _preprocess_roots(cls, poly):
        """Take heroic measures to make ``poly`` compatible with ``CRootOf``."""
        from .polys.polyroots import preprocess_roots
        dom = poly.get_domain()

        if not dom.is_Exact:
            poly = poly.to_exact()

        coeff, poly = preprocess_roots(poly)
        dom = poly.get_domain()

        if not dom.is_ZZ:
            raise NotImplementedError(
                "sorted roots not supported over %s" % dom)

        return coeff, poly

    @classmethod
    def _postprocess_root(cls, root, radicals):
        """Return the root if it is trivial or a ``CRootOf`` object. """
        poly, index = root
        roots = cls._roots_trivial(poly, radicals)

        if roots is not None:
            return roots[index]
        else:
            return cls._new(poly, index)

    @classmethod
    def _get_roots(cls, method, poly, radicals):
        """Return postprocessed roots of specified kind. """
        from . import Dummy, symbols
        from .polys.polyerrors import PolynomialError
        from .utilities import numbered_symbols
        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are allowed")

        dom = poly.get_domain()

        # get rid of gen and it's free symbol
        d = Dummy()
        poly = poly.subs(poly.gen, d)
        x = symbols('x')
        # see what others are left and select x or a numbered x
        # that doesn't clash
        free_names = {str(i) for i in poly.free_symbols}
        for x in chain((symbols('x'),), numbered_symbols('x')):
            if x.name not in free_names:
                poly = poly.replace(d, x)
                break

        if dom.is_QQ or dom.is_ZZ:
            return cls._get_roots_qq(method, poly, radicals)
        elif dom.is_AlgebraicField or dom.is_ZZ_I or dom.is_QQ_I:
            return cls._get_roots_alg(method, poly, radicals)
        else:
            # XXX: not sure how to handle ZZ[x] which appears in some tests.
            # This should not be supported but it is not clear how to remove
            # it...
            return cls._get_roots_qq(method, poly, radicals)

    @classmethod
    def _get_roots_qq(cls, method, poly, radicals):
        """Return postprocessed roots of specified kind
        for polynomials with rational coefficients. """
        coeff, poly = cls._preprocess_roots(poly)
        roots = []

        for root in getattr(cls, method)(poly):
            roots.append(coeff*cls._postprocess_root(root, radicals))

        return roots

    @classmethod
    def _get_roots_alg(cls, method, poly, radicals):
        """Return postprocessed roots of specified kind
        for polynomials with algebraic coefficients. It assumes
        the domain is already an algebraic field. First it
        finds the roots using _get_roots_qq, then uses the
        square-free factors to filter roots and get the correct
        multiplicity.
        """
        # Existing QQ code can handle non-square-free polynomials
        # The lifted polynomial has the same roots as poly
        roots = cls._get_roots_qq(method, poly.lift(), radicals)
        subroots = {}
        for f, m in poly.sqf_list()[1]:
            if method == "_real_roots":
                roots_filt = f.which_real_roots(roots)
            elif method == "_all_roots":
                roots_filt = f.which_all_roots(roots)
            for r in roots_filt:
                subroots[r] = m

        roots_seen = set()
        roots_flat = []
        for r in roots:
            if r in subroots and r not in roots_seen:
                m = subroots[r]
                roots_flat.extend([r] * m)
                roots_seen.add(r)

        return roots_flat

    @classmethod
    def clear_cache(cls):
        """Reset cache for reals and complexes.

        The intervals used to approximate a root instance are updated
        as needed. When a request is made to see the intervals, the
        most current values are shown. `clear_cache` will reset all
        CRootOf instances back to their original state.
        """
        global _reals_cache, _complexes_cache
        _reals_cache = _pure_key_dict()
        _complexes_cache = _pure_key_dict()


CRootOf = ComplexRootOf


# Methods of the native root objects.

def _root_poly(self):
    from . import PurePoly
    return PurePoly(self.args[0], expand=False)


def _root_index(self):
    return int(self.args[1])


def _root_expr(self):
    return self.args[0]


def _ensure_reals_init(self):
    """Ensure that our poly has entries in the reals cache. """
    poly = self.poly
    if poly not in _reals_cache:
        ComplexRootOf._indexed_root(poly, self.index)


def _ensure_complexes_init(self):
    """Ensure that our poly has entries in the complexes cache. """
    poly = self.poly
    if poly not in _complexes_cache:
        ComplexRootOf._indexed_root(poly, self.index)


def _get_interval(self):
    """Internal function for retrieving isolation interval from cache. """
    self._ensure_reals_init()
    poly, index = self.poly, self.index
    if self.is_real:
        return _reals_cache[poly][index]
    else:
        reals_count = len(_reals_cache[poly])
        self._ensure_complexes_init()
        return _complexes_cache[poly][index - reals_count]


def _set_interval(self, interval):
    """Internal function for updating isolation interval in cache. """
    self._ensure_reals_init()
    poly, index = self.poly, self.index
    if self.is_real:
        _reals_cache[poly][index] = interval
    else:
        reals_count = len(_reals_cache[poly])
        self._ensure_complexes_init()
        _complexes_cache[poly][index - reals_count] = interval


def _reset(self):
    """
    Reset all intervals
    """
    ComplexRootOf._all_roots(self.poly, use_cache=False)


def eval_approx(self, n, return_mpmath=False):
    """Evaluate this complex root to the given precision.

    This uses secant method and root bounds are used to both
    generate an initial guess and to check that the root
    returned is valid. If ever the method converges outside the
    root bounds, the bounds will be made smaller and updated.
    """
    from mpmath import findroot, mpc, mpf, workprec
    from mpmath.libmp.libmpf import dps_to_prec
    from . import Dummy, Float, I
    from .utilities import lambdify
    prec = dps_to_prec(n)
    with workprec(prec):
        g = self.poly.gen
        if not g.is_Symbol:
            d = Dummy('x')
            if self.is_imaginary:
                d *= I
            func = lambdify(d, self.expr.subs(g, d))
        else:
            expr = self.expr
            if self.is_imaginary:
                expr = self.expr.subs(g, I*g)
            func = lambdify(g, expr)

        interval = self._get_interval()
        while True:
            if self.is_real:
                a = mpf(str(interval.a))
                b = mpf(str(interval.b))
                if a == b:
                    root = a
                    break
                x0 = mpf(str(interval.center))
                x1 = x0 + mpf(str(interval.dx))/4
            elif self.is_imaginary:
                a = mpf(str(interval.ay))
                b = mpf(str(interval.by))
                if a == b:
                    root = mpc(mpf('0'), a)
                    break
                x0 = mpf(str(interval.center[1]))
                x1 = x0 + mpf(str(interval.dy))/4
            else:
                ax = mpf(str(interval.ax))
                bx = mpf(str(interval.bx))
                ay = mpf(str(interval.ay))
                by = mpf(str(interval.by))
                if ax == bx and ay == by:
                    # the sign of the imaginary part will be assigned
                    # according to the desired index using the fact that
                    # roots are sorted with negative imag parts coming
                    # before positive (and all imag roots coming after real
                    # roots)
                    root = mpc(ax, ay)
                    break
                x0 = mpc(*map(str, interval.center))
                x1 = x0 + mpc(*map(str, (interval.dx, interval.dy)))/4
            try:
                # without a tolerance, this will return when (to within
                # the given precision) x_i == x_{i-1}
                root = findroot(func, (x0, x1))
                # If the (real or complex) root is not in the 'interval',
                # then keep refining the interval. This happens if findroot
                # accidentally finds a different root outside of this
                # interval because our initial estimate 'x0' was not close
                # enough. It is also possible that the secant method will
                # get trapped by a max/min in the interval; the root
                # verification by findroot will raise a ValueError in this
                # case and the interval will then be tightened -- and
                # eventually the root will be found.
                #
                # It is also possible that findroot will not have any
                # successful iterations to process (in which case it
                # will fail to initialize a variable that is tested
                # after the iterations and raise an UnboundLocalError).
                if self.is_real or self.is_imaginary:
                    if not bool(root.imag) == self.is_real and (
                            a <= root <= b):
                        if self.is_imaginary:
                            root = mpc(mpf('0'), root.real)
                        break
                elif (ax <= root.real <= bx and ay <= root.imag <= by):
                    break
            except (UnboundLocalError, ValueError):
                pass
            interval = interval.refine()

    # update the interval so we at least (for this precision or
    # less) don't have much work to do to recompute the root
    self._set_interval(interval)
    if return_mpmath:
        return root
    return (_float_new(root.real._mpf_, prec) +
        I*_float_new(root.imag._mpf_, prec))


def _float_new(mpf_, prec):
    """The reference ``Float._new``: zero is ``S.Zero``, anything else a
    Float holding exactly ``mpf_`` at precision ``prec``."""
    from . import Float, S
    from mpmath import make_mpf
    from mpmath.libmp import fzero
    if mpf_ == fzero:
        return S.Zero
    return Float(make_mpf(mpf_), precision=prec)


def eval_rational(self, dx=None, dy=None, n=15):
    """
    Return a Rational approximation of ``self`` that has real
    and imaginary component approximations that are within ``dx``
    and ``dy`` of the true values, respectively. Alternatively,
    ``n`` digits of precision can be specified.

    The interval is refined with bisection and is sure to
    converge. The root bounds are updated when the refinement
    is complete so recalculation at the same or lesser precision
    will not have to repeat the refinement and should be much
    faster.
    """
    from . import I, Rational, S
    dy = dy or dx
    if dx:
        rtol = None
        dx = dx if isinstance(dx, Rational) else Rational(str(dx))
        dy = dy if isinstance(dy, Rational) else Rational(str(dy))
    else:
        # assume n is the number of desired digits
        rtol = S(10)**-(n + 2)  # +2 for guard digits
    interval = self._get_interval()
    while True:
        if self.is_real:
            if rtol:
                dx = abs(interval.center*rtol)
            interval = interval.refine_size(dx=dx)
            c = interval.center
            real = Rational(c)
            imag = S.Zero
            if not rtol or interval.dx < abs(c*rtol):
                break
        elif self.is_imaginary:
            if rtol:
                dy = abs(interval.center[1]*rtol)
                dx = 1
            interval = interval.refine_size(dx=dx, dy=dy)
            c = interval.center[1]
            imag = Rational(c)
            real = S.Zero
            if not rtol or interval.dy < abs(c*rtol):
                break
        else:
            if rtol:
                dx = abs(interval.center[0]*rtol)
                dy = abs(interval.center[1]*rtol)
            interval = interval.refine_size(dx, dy)
            c = interval.center
            real, imag = map(Rational, c)
            if not rtol or (
                    interval.dx < abs(c[0]*rtol) and
                    interval.dy < abs(c[1]*rtol)):
                break

    # update the interval so we at least (for this precision or
    # less) don't have much work to do to recompute the root
    self._set_interval(interval)
    return real + I*imag


def _eval_evalf(self, prec, **kwargs):
    """Evaluate this complex root to the given precision."""
    # all kwargs are ignored
    # ``eval_rational(n=prec_to_dps(prec))._evalf(prec)``: each rational
    # part rounded once to ``prec`` bits.
    return root_float(self, prec)


def root_float(node, prec):
    """The value of a root object at ``prec`` bits, as the reference's
    ``_eval_evalf`` gives it: a rational approximation rounded once."""
    from mpmath.libmp.libmpf import prec_to_dps
    from . import Float, I, S
    value = node.eval_rational(n=prec_to_dps(prec))
    re, im = value.as_real_imag()
    re = Float(re, precision=prec) if re else S.Zero
    im = Float(im, precision=prec) if im else S.Zero
    return re + I*im


def _eval_is_eq(lhs, rhs):
    """The reference's ``_eval_is_eq`` for a root object and a number."""
    from . import Rational, sympify
    from .core.relational import is_le
    if getattr(rhs, 'class_name', None) == 'ComplexRootOf':
        return lhs == rhs
    if not rhs.is_number:
        return None
    if not rhs.is_finite:
        return False
    z = lhs.expr.subs(lhs.expr.free_symbols.pop(), rhs).is_zero
    if z is False:  # all roots will make z True but we don't know
        # whether this is the right root if z is True
        return False
    o = rhs.is_real, rhs.is_imaginary
    s = lhs.is_real, lhs.is_imaginary
    assert None not in s  # this is part of initial refinement
    if o != s and None not in o:
        return False
    re, im = rhs.as_real_imag()
    if lhs.is_real:
        if im:
            return False
        i = lhs._get_interval()
        a, b = [Rational(str(_)) for _ in (i.a, i.b)]
        return sympify(a <= rhs and rhs <= b)
    i = lhs._get_interval()
    r1, r2, i1, i2 = [Rational(str(j)) for j in (
        i.ax, i.bx, i.ay, i.by)]
    return is_le(r1, re) and is_le(re, r2) and is_le(i1, im) and is_le(im, i2)


class RootSum(metaclass=_ExprType):
    """Represents a sum of all roots of a univariate polynomial.

    The reference ``sympy.polys.rootoftools.RootSum``. The native node holds
    ``Lambda(gen, poly_expr)`` and ``fun``, so the generator of the
    polynomial is bound; ``args`` gives the reference ``(expr, fun, gen)``.
    """
    _matches = staticmethod(lambda expr:expr.class_name == 'RootSum')

    def __new__(cls, expr, func=None, x=None, auto=True, quadratic=False):
        """Construct a new ``RootSum`` instance of roots of a polynomial."""
        from . import Lambda, S, Add
        from .polys.polyerrors import MultivariatePolynomialError
        from .polys.polyroots import roots_linear, roots_quadratic
        coeff, poly = cls._transform(expr, x)

        if not poly.is_univariate:
            raise MultivariatePolynomialError(
                "only univariate polynomials are allowed")

        if func is None:
            func = Lambda(poly.gen, poly.gen)
        else:
            is_func = getattr(func, 'is_Function', False) or isinstance(func, Lambda)
            if is_func and _accepts_one(func):
                if not isinstance(func, Lambda):
                    func = Lambda(poly.gen, func(poly.gen))
            else:
                raise ValueError(
                    "expected a univariate function, got %s" % func)

        var, expr = func.variables[0], func.expr

        if coeff is not S.One:
            expr = expr.subs(var, coeff*var)

        deg = poly.degree()

        if not expr.has(var):
            return deg*expr

        if expr.is_Add:
            add_const, expr = expr.as_independent(var)
        else:
            add_const = S.Zero

        if expr.is_Mul:
            mul_const, expr = expr.as_independent(var)
        else:
            mul_const = S.One

        func = Lambda(var, expr)

        rational = cls._is_func_rational(poly, func)
        factors, terms = _pure_factors(poly), []

        for poly, k in factors:
            if poly.degree() == 1:
                term = func(roots_linear(poly)[0])
            elif quadratic and poly.degree() == 2:
                term = sum(map(func, roots_quadratic(poly)))
            else:
                if not rational or not auto:
                    term = cls._new(poly, func, auto)
                else:
                    term = cls._rational_case(poly, func)

            terms.append(k*term)

        return mul_const*Add(*terms) + deg*add_const

    @classmethod
    def _new(cls, poly, func, auto=True):
        """Construct new raw ``RootSum`` instance. """
        from . import Lambda
        node = _core._apply('RootSum', Lambda(poly.gen, poly.as_expr()), func, evaluate=False)
        _sums[node] = poly, auto
        return node

    @classmethod
    def new(cls, poly, func, auto=True):
        """Construct new ``RootSum`` instance. """
        if not func.expr.has(*func.variables):
            return func.expr

        rational = cls._is_func_rational(poly, func)

        if not rational or not auto:
            return cls._new(poly, func, auto)
        else:
            return cls._rational_case(poly, func)

    @classmethod
    def _transform(cls, expr, x):
        """Transform an expression to a polynomial. """
        from . import PurePoly
        from .polys.polyroots import preprocess_roots
        poly = PurePoly(expr, x, greedy=False)
        return preprocess_roots(poly)

    @classmethod
    def _is_func_rational(cls, poly, func):
        """Check if a lambda is a rational function. """
        var, expr = func.variables[0], func.expr
        return expr.is_rational_function(var)

    @classmethod
    def _rational_case(cls, poly, func):
        """Handle the rational function case. """
        from . import symbols, together, Poly, QQ, factor
        from .polys.polyerrors import GeneratorsNeeded
        from .polys.polyfuncs import symmetrize, viete
        roots = symbols('r:%d' % poly.degree())
        var, expr = func.variables[0], func.expr

        f = sum(expr.subs(var, r) for r in roots)
        p, q = together(f).as_numer_denom()

        domain = QQ[roots]

        p = p.expand()
        q = q.expand()

        try:
            p = Poly(p, domain=domain, expand=False)
        except GeneratorsNeeded:
            p, p_coeff = None, (p,)
        else:
            p_monom, p_coeff = zip(*p.terms())

        try:
            q = Poly(q, domain=domain, expand=False)
        except GeneratorsNeeded:
            q, q_coeff = None, (q,)
        else:
            q_monom, q_coeff = zip(*q.terms())

        coeffs, mapping = symmetrize(p_coeff + q_coeff, formal=True)
        formulas, values = viete(poly, roots), []

        for (sym, _), (_, val) in zip(mapping, formulas):
            values.append((sym, val))

        for i, (coeff, _) in enumerate(coeffs):
            coeffs[i] = coeff.subs(values)

        n = len(p_coeff)

        p_coeff = coeffs[:n]
        q_coeff = coeffs[n:]

        if p is not None:
            p = Poly(dict(zip(p_monom, p_coeff)), *p.gens).as_expr()
        else:
            (p,) = p_coeff

        if q is not None:
            q = Poly(dict(zip(q_monom, q_coeff)), *q.gens).as_expr()
        else:
            (q,) = q_coeff

        return factor(p/q)


def _accepts_one(func):
    nargs = getattr(func, 'nargs', None)
    if nargs is None:
        return True
    from . import Integer
    try:
        return 1 in nargs
    except TypeError:
        # the nominal FiniteSet of ``Lambda.nargs`` is not iterable here
        return Integer(1) in getattr(nargs, 'args', ())


def _pure_factors(poly):
    from . import PurePoly
    _, factors = poly.factor_list()
    return [(PurePoly(f, expand=False), m) for f, m in factors]


def _root_sum_poly(self):
    from . import PurePoly
    stored = _sums.get(self)
    if stored is not None:
        return stored[0]
    (gen,), expr = _root_sum_lambda(self)
    return PurePoly(expr, gen)


def _root_sum_lambda(self):
    lam = _core.Expr.args.__get__(self)[0]
    return tuple(lam.variables), lam.expr


def _root_sum_fun(self):
    return _core.Expr.args.__get__(self)[1]


def _root_sum_auto(self):
    stored = _sums.get(self)
    return True if stored is None else stored[1]


def _root_sum_expr(self):
    return self.poly.as_expr()


def _root_sum_args(self):
    return (self.expr, self.fun, self.poly.gen)


def _root_sum_is_commutative(self):
    return True


def root_sum_doit(node, **hints):
    from . import Add
    from .polys.polyroots import roots
    if not hints.get('roots', True):
        return node

    _roots = roots(node.poly, multiple=True)

    if len(_roots) < node.poly.degree():
        return node
    else:
        return Add(*[node.fun(r) for r in _roots])


def root_sum_evalf(node, prec):
    from . import Add
    from mpmath.libmp.libmpf import prec_to_dps
    from .polys.polyerrors import DomainError, PolynomialError
    try:
        _roots = node.poly.nroots(n=prec_to_dps(prec))
    except (DomainError, PolynomialError):
        return node
    else:
        return Add(*[node.fun(r) for r in _roots])


def root_sum_derivative(node, x):
    from . import Lambda
    var, expr = node.fun.args
    if getattr(var, 'is_Tuple', False) or isinstance(var, tuple):
        (var,) = tuple(var)
    func = Lambda(var, expr.diff(x))
    return RootSum.new(node.poly, func, node.auto)


def install():
    from ._node_types import node_class
    cls = node_class('ComplexRootOf')
    for name, function in (('poly', _root_poly), ('index', _root_index), ('expr', _root_expr),
                           ('func', lambda self: ComplexRootOf)):
        type.__setattr__(cls, name, property(function))
    root_sum = node_class('RootSum')
    for name, function in (('poly', _root_sum_poly), ('fun', _root_sum_fun), ('auto', _root_sum_auto),
                           ('expr', _root_sum_expr), ('args', _root_sum_args),
                           ('is_commutative', _root_sum_is_commutative)):
        type.__setattr__(root_sum, name, property(function))
    type.__setattr__(root_sum, '_eval_evalf', lambda self, prec: root_sum_evalf(self, prec))
    type.__setattr__(root_sum, '_eval_derivative', lambda self, x: root_sum_derivative(self, x))
    for method in (_ensure_reals_init, _ensure_complexes_init, _get_interval, _set_interval,
                   _reset, eval_approx, eval_rational, _eval_evalf, _eval_is_eq):
        type.__setattr__(cls, method.__name__, method)
