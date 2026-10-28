"""This module implements tools for integrating rational functions. """

from symox.core.function import Lambda
from symox.core.numbers import I
from symox.core.singleton import S
from symox.core.symbol import (Dummy, Symbol, symbols)
from symox.functions.elementary.exponential import log
from symox.functions.elementary.trigonometric import atan
from symox.polys.polyerrors import DomainError
from symox.polys.polyroots import roots
from symox.polys.polytools import cancel
from symox.polys.rootoftools import RootSum
from symox.polys import Poly, resultant, ZZ
from symox.polys import _wrt_option as _lrt_wrt


def _atan_function_states(functions):
    from types import FunctionType
    pending=list(functions);seen=set();states=[]
    while pending:
        function=pending.pop()
        if type(function) is not FunctionType or id(function) in seen:continue
        seen.add(id(function));closure=function.__closure__
        contents=tuple(cell.cell_contents for cell in closure or ())
        states.append((function,function.__code__,closure,contents,function.__defaults__,function.__kwdefaults__))
        pending.extend(value for value in contents if type(value) is FunctionType)
    return tuple(states)


# These guards only select a native route. The full general algorithm below
# remains responsible for custom constructors, Poly methods and held evaluation.
from symox import _core as _rational_core
from symox import _poly_constructor as _lrt_constructor
import symox as _lrt_public
import types as _lrt_types
from symox.core.sympify import sympify as _lrt_stock_sympify
from symox import _node_types as _atan_nodes
from symox._public_core import _evaluation as _atan_evaluation
_ATAN_HEAD = globals().get('_ATAN_HEAD', atan)
_ATAN_POLY = globals().get('_ATAN_POLY', Poly)
_ATAN_TYPE = globals().get('_ATAN_TYPE', type)
_ATAN_GETATTR = globals().get('_ATAN_GETATTR', getattr)
_ATAN_FACTORY = globals().get('_ATAN_FACTORY', atan.__dict__.get('_symox_factory'))
_ATAN_CALL = globals().get('_ATAN_CALL', type(atan).__dict__.get('__call__'))
_ATAN_AS_EXPR = globals().get('_ATAN_AS_EXPR', Poly.as_expr)
_ATAN_APPLY = globals().get('_ATAN_APPLY', _rational_core._apply)
_ATAN_OPERATORS = globals().get('_ATAN_OPERATORS', tuple(
    (cls, name, getattr(cls, name))
    for cls in (_rational_core.Expr, atan, _atan_nodes.node_class('atan'), _atan_nodes.node_class('Add'), _atan_nodes.node_class('Mul'), type(S.Zero), type(S.Pi))
    for name in ('__add__', '__radd__', '__mul__', '__rmul__')))
_LRT_RESULTANT = globals().get('_LRT_RESULTANT', _rational_core.resultant)
if '_LRT_WRT_CODE' not in globals():
    _lrt_wrt_factory = _lrt_wrt.__dict__.get('_wrap')
    _LRT_WRT_CODE = next((code for code in (
        _lrt_wrt_factory.__code__.co_consts if type(_lrt_wrt_factory) is _lrt_types.FunctionType else ())
        if type(code) is _lrt_types.CodeType and code.co_name == 'wrapper'), None)
    del _lrt_wrt_factory
_LRT_WRT_GLOBALS = globals().get('_LRT_WRT_GLOBALS', _lrt_wrt.__dict__)
_LRT_CONSTRUCT = _lrt_constructor._native_construct
_LRT_SYMPIFY = globals().get('_LRT_SYMPIFY', _lrt_stock_sympify)
_LRT_AS_POLY_CODE = next(code for code in _lrt_public._install_poly_methods.__code__.co_consts
                         if isinstance(code, _lrt_types.CodeType) and code.co_name == 'as_poly')
_ATAN_ALL = globals().get('_ATAN_ALL', all)
_ATAN_BOOTSTRAP_STATES = globals().get('_ATAN_BOOTSTRAP_STATES',
    _atan_function_states(original for _,_,original in _ATAN_OPERATORS))
# Symox imports this adapter during bootstrap, before later user overrides.
# A removed/reimported module cannot bless arbitrary newly captured operators.
import sys as _atan_sys
_ATAN_BOOTSTRAP_CAPTURE = globals().get('_ATAN_BOOTSTRAP_CAPTURE',
    bool(getattr(getattr(_atan_sys.modules.get('symox'),'__spec__',None),'_initializing',False)))
_ATAN_GENS = globals().get('_ATAN_GENS', Poly.gens)
_ATAN_SCALARS = globals().get('_ATAN_SCALARS', tuple(_atan_nodes.node_class(name) for name in
    ('Symbol','Dummy','Integer','Rational','Zero','One','NegativeOne','Half','Pi','Add','Mul','Pow','atan')))
_ATAN_OBJECT_GETATTRIBUTE = globals().get('_ATAN_OBJECT_GETATTRIBUTE', object.__getattribute__)
_ATAN_MRO_DESCRIPTOR = globals().get('_ATAN_MRO_DESCRIPTOR', type.__dict__['__mro__'])
_ATAN_DICT_DESCRIPTOR = globals().get('_ATAN_DICT_DESCRIPTOR', type.__dict__['__dict__'])
_ATAN_MISSING = globals().get('_ATAN_MISSING', object())


def _atan_raw_attribute(cls,name):
    # Inspect the stored definition; neither instance nor metaclass
    # descriptors may run during route selection.
    for base in _ATAN_MRO_DESCRIPTOR.__get__(cls):
        namespace=_ATAN_DICT_DESCRIPTOR.__get__(base)
        if name in namespace:return namespace[name]
    return _ATAN_MISSING


_ATAN_PRIORITIES = globals().get('_ATAN_PRIORITIES', tuple(
    (cls,_atan_raw_attribute(cls,'_op_priority')) for cls in _ATAN_SCALARS))
_LRT_OPERATORS = globals().get('_LRT_OPERATORS', _ATAN_OPERATORS + tuple(
    (cls,name,getattr(cls,name)) for cls in (_rational_core.Expr,)+_ATAN_SCALARS
    for name in ('__add__','__radd__','__sub__','__rsub__','__mul__','__rmul__',
                 '__truediv__','__rtruediv__','__pow__','__rpow__','__neg__')
    if all(cls is not previous or name!=operation for previous,operation,_ in _ATAN_OPERATORS)))
_LRT_BOOTSTRAP_STATES = globals().get('_LRT_BOOTSTRAP_STATES',
    _atan_function_states(original for _,_,original in _LRT_OPERATORS))




def _atan_states_match(states):
    from types import FunctionType
    if _ATAN_TYPE(states) is not tuple:return False
    try:
        for state in states:
            if _ATAN_TYPE(state) is not tuple or len(state)!=6:return False
            function,code,closure,contents,defaults,kwdefaults=state
            if (_ATAN_TYPE(function) is not FunctionType or _ATAN_TYPE(contents) is not tuple
                    or function.__code__ is not code or function.__closure__ is not closure
                    or function.__defaults__ is not defaults or function.__kwdefaults__ is not kwdefaults
                    or len(closure or ())!=len(contents)
                    or not _ATAN_ALL(cell.cell_contents is value for cell,value in zip(closure or (),contents))):return False
    except (AttributeError,ValueError,TypeError):return False
    return True


def _atan_forwarding_proof(proof,operation_code,name):
    import builtins
    from types import FunctionType
    operation,delegate,states,namespace,ore,mro,isinstance_,getattr_=proof
    if (_ATAN_TYPE(operation) is not FunctionType or operation.__code__ is not operation_code
            or operation.__globals__ is not namespace or _ATAN_TYPE(namespace) is not dict
            or namespace.get('_OreOperator') is not ore
            or _ATAN_TYPE(ore) is not _ATAN_TYPE or ore.__mro__ is not mro
            or namespace.get('isinstance',builtins.isinstance) is not isinstance_
            or namespace.get('getattr',builtins.getattr) is not getattr_):return False
    try:
        cells=dict(zip(operation_code.co_freevars,operation.__closure__ or ()))
        if (set(cells)!= {'original','name'} or cells['original'].cell_contents is not delegate
                or _ATAN_TYPE(cells['name'].cell_contents) is not str
                or cells['name'].cell_contents!=name):return False
    except (ValueError,TypeError):return False
    return _atan_states_match(states)


def _atan_vector_proof(proof,operation_code,name,vector_types):
    import builtins
    from types import FunctionType
    operation,delegate,states,namespace,basis,mro,isinstance_,operation_name=proof
    if (_ATAN_TYPE(operation) is not FunctionType or operation.__code__ is not operation_code
            or operation.__globals__ is not namespace or _ATAN_TYPE(namespace) is not dict
            or _ATAN_TYPE(basis) is not _ATAN_TYPE or basis.__mro__ is not mro
            or namespace.get('isinstance',builtins.isinstance) is not isinstance_
            or operation_name!=name or not _atan_states_match(states)):return False
    try:
        cells=dict(zip(operation_code.co_freevars,operation.__closure__ or ()))
        if (set(cells)!={'original','name','_BasisDependent'}
                or cells['original'].cell_contents is not delegate
                or cells['_BasisDependent'].cell_contents is not basis
                or _ATAN_TYPE(cells['name'].cell_contents) is not str
                or cells['name'].cell_contents!=name):return False
    except (ValueError,TypeError):return False
    # Exclude these basis classes once for every scalar MRO below, after
    # checking all operator layers. No state survives this selector call.
    vector_types[id(basis)]=basis
    return True


def _atan_forwarding_chain(current,original,name,ore_proofs,ore_code,vector_proofs,vector_code,vector_types):
    seen=set();forwarded=0
    while id(current) not in seen:
        seen.add(id(current))
        ore=next((proof for proof in ore_proofs if current is proof[0]),None)
        vector=next((proof for proof in vector_proofs if current is proof[0]),None)
        if ore is not None:
            if (_ATAN_GETATTR(current,'_symox_ore_proof',None) is not ore
                    or not _atan_forwarding_proof(ore,ore_code,name)):return False,forwarded
            delegate=ore[1];forwarded|=1
        elif vector is not None:
            if (_ATAN_GETATTR(current,'_symox_vector_proof',None) is not vector
                    or not _atan_vector_proof(vector,vector_code,name,vector_types)):return False,forwarded
            delegate=vector[1];forwarded|=2
        else:return current is original,forwarded
        if current is original:return True,forwarded
        current=delegate
    return False,forwarded


def _atan_stock_operators(f,g,operators=None,states=None):
    import sys
    from types import CodeType,FunctionType,GetSetDescriptorType
    if operators is None:operators=_ATAN_OPERATORS
    if states is None:states=_ATAN_BOOTSTRAP_STATES
    if _ATAN_TYPE(operators) is not tuple or not _atan_states_match(states):return False
    recurrence=sys.modules.get('symox.holonomic.recurrence')
    proofs=()
    if recurrence is not None:
        proofs=recurrence.__dict__.get('_SCALAR_DISPATCH_PROOFS')
        installer=recurrence.__dict__.get('_install_scalar_dispatch')
        if (_ATAN_TYPE(proofs) is not tuple or _ATAN_TYPE(installer) is not FunctionType
                or proofs is not _ATAN_GETATTR(installer,'_stock_proofs',None)):return False
        factory=next((code for code in installer.__code__.co_consts
                      if _ATAN_TYPE(code) is CodeType and code.co_name=='factory'),None)
        if factory is None:return False
        operation_code=next((code for code in factory.co_consts
                            if _ATAN_TYPE(code) is CodeType and code.co_name=='operation'),None)
        if operation_code is None:return False
    vector=sys.modules.get('symox.vector._adapters')
    vector_proofs=();vector_code=None
    if vector is not None:
        namespace=vector.__dict__;vector_proofs=namespace.get('_SCALAR_DISPATCH_PROOFS')
        installer=namespace.get('install');vector_code=namespace.get('_SCALAR_OPERATION_CODE')
        if (_ATAN_TYPE(vector_proofs) is not tuple or _ATAN_TYPE(installer) is not FunctionType
                or installer.__code__ is not namespace.get('_SCALAR_INSTALL_CODE')
                or vector_proofs is not _ATAN_GETATTR(installer,'_stock_proofs',None)
                or _ATAN_TYPE(vector_code) is not CodeType):return False
    if (not _ATAN_ALL(_ATAN_TYPE(p) is tuple and len(p)==8 for p in proofs)
            or not _ATAN_ALL(_ATAN_TYPE(p) is tuple and len(p)==8 for p in vector_proofs)):return False
    verified={};forwarded=0;vector_types={}
    for cls,name,original in operators:
        current=_ATAN_GETATTR(cls,name)
        key=(id(current),id(original),name)
        if key in verified:
            if not verified[key]:return False
            continue
        valid,has_forwarder=_atan_forwarding_chain(current,original,name,proofs,
                                                  operation_code if recurrence is not None else None,
                                                  vector_proofs,vector_code,vector_types)
        verified[key]=valid
        if not valid:return False
        forwarded|=has_forwarder
    # Inspect each scalar MRO and each base dictionary once. Priority and
    # attribute hooks apply to all delegates; holonomic markers apply only
    # when an actual Ore layer was authenticated. A vector-only chain never
    # reads that marker in its normal scalar path.
    has_ore=forwarded&1
    bases=tuple(vector_types.values())
    for cls,priority in _ATAN_PRIORITIES:
        mro=_ATAN_MRO_DESCRIPTOR.__get__(cls)
        value=_ATAN_MISSING;attribute=_ATAN_MISSING
        found_value=False;found_attribute=False
        for base in mro:
            namespace=_ATAN_DICT_DESCRIPTOR.__get__(base)
            if not found_value and '_op_priority' in namespace:
                value=namespace['_op_priority'];found_value=True
            if not found_attribute and '__getattribute__' in namespace:
                attribute=namespace['__getattribute__'];found_attribute=True
            if '__getattr__' in namespace:return False
            if has_ore and '_symox_holonomic_function' in namespace:return False
            if bases:
                for basis in bases:
                    if base is basis:return False
        if value is not priority or attribute is not _ATAN_OBJECT_GETATTRIBUTE:return False
        if has_ore and _ATAN_GETATTR(cls,'__getattribute__') is not _ATAN_OBJECT_GETATTRIBUTE:return False
    if not forwarded:return True
    if f is None and g is None:return True
    if (_ATAN_TYPE(f) is not _ATAN_POLY or _ATAN_TYPE(g) is not _ATAN_POLY
            or _ATAN_TYPE(_ATAN_GENS) is not GetSetDescriptorType or Poly.gens is not _ATAN_GENS):return False
    # Read the stock native descriptor for selection only. Host generators
    # retain the full algorithm, including their arithmetic/attribute hooks.
    for polynomial in (f,g):
        for gen in _ATAN_GENS.__get__(polynomial,Poly):
            cls=_ATAN_TYPE(gen)
            if cls is not _ATAN_SCALARS[0] and cls is not _ATAN_SCALARS[1]:return False
    return True


def _stock_lrt_as_poly():
    # The public Expr facade can override a scalar method without changing
    # the storage base's descriptor. Check the methods actually resolved by
    # the native scalar classes produced by LRT, without invoking descriptors.
    for cls in (_rational_core.Expr,)+_ATAN_SCALARS:
        method=_atan_raw_attribute(cls,'as_poly')
        if (type(method) is not _lrt_types.FunctionType
                or method.__code__ is not _LRT_AS_POLY_CODE
                or method.__globals__ is not _lrt_public.__dict__):return False
    return True


def _stock_lrt_resultant():
    if resultant is _LRT_RESULTANT:return True
    # The ordering-options adapter is installed after this module's first
    # bootstrap import. Reload then sees that official wrapper instead of
    # its original builtin. LRT supplies only includePRS: the verified
    # wrapper's first branch forwards directly, without invoking any hook.
    if (type(resultant) is not _lrt_types.FunctionType
            or resultant.__code__ is not _LRT_WRT_CODE
            or resultant.__globals__ is not _LRT_WRT_GLOBALS
            or resultant.__defaults__ is not None
            or resultant.__kwdefaults__ is not None):return False
    try:
        cells=dict(zip(_LRT_WRT_CODE.co_freevars,resultant.__closure__ or ()))
        return (set(cells)=={'name','native'}
                and type(cells['name'].cell_contents) is str
                and cells['name'].cell_contents=='resultant'
                and cells['native'].cell_contents is _LRT_RESULTANT)
    except (AttributeError,ValueError,TypeError):return False


def _stock_atan_protocol(f=None,g=None):
    return (_ATAN_BOOTSTRAP_CAPTURE and atan is _ATAN_HEAD and Poly is _ATAN_POLY and _atan_evaluation.get()
            and _ATAN_TYPE(atan).__dict__.get('__call__') is _ATAN_CALL
            and atan.__dict__.get('_symox_factory') is _ATAN_FACTORY
            and Poly.as_expr is _ATAN_AS_EXPR
            and _rational_core._apply is _ATAN_APPLY
            and _atan_stock_operators(f,g)
            and log_to_atan is _ATAN_RECURSION)


def ratint(f, x, **flags):
    """
    Performs indefinite integration of rational functions.

    Explanation
    ===========

    Given a field :math:`K` and a rational function :math:`f = p/q`,
    where :math:`p` and :math:`q` are polynomials in :math:`K[x]`,
    returns a function :math:`g` such that :math:`f = g'`.

    Examples
    ========

    >>> from symox.integrals.rationaltools import ratint
    >>> from symox.abc import x

    >>> ratint(36/(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2), x)
    (12*x + 6)/(x**2 - 1) + 4*log(x - 2) - 4*log(x + 1)

    References
    ==========

    .. [1] M. Bronstein, Symbolic Integration I: Transcendental
       Functions, Second Edition, Springer-Verlag, 2005, pp. 35-70

    See Also
    ========

    symox.integrals.integrals.Integral.doit
    symox.integrals.rationaltools.ratint_logpart
    symox.integrals.rationaltools.ratint_ratpart

    """
    if isinstance(f, tuple):
        p, q = f
    else:
        p, q = f.as_numer_denom()

    p, q = Poly(p, x, composite=False, field=True), Poly(q, x, composite=False, field=True)

    coeff, p, q = p.cancel(q)
    poly, p = p.div(q)

    result = poly.integrate(x).as_expr()

    if p.is_zero:
        return coeff*result

    g, h = ratint_ratpart(p, q, x)

    P, Q = h.as_numer_denom()

    P = Poly(P, x)
    Q = Poly(Q, x)

    q, r = P.div(Q)

    result += g + q.integrate(x).as_expr()

    if not r.is_zero:
        symbol = flags.get('symbol', 't')

        if not isinstance(symbol, Symbol):
            t = Dummy(symbol)
        else:
            t = symbol.as_dummy()

        L = ratint_logpart(r, Q, x, t)

        real = flags.get('real')

        if real is None:
            if isinstance(f, tuple):
                p, q = f
                atoms = p.atoms() | q.atoms()
            else:
                atoms = f.atoms()

            for elt in atoms - {x}:
                if not elt.is_extended_real:
                    real = False
                    break
            else:
                real = True

        eps = S.Zero

        if not real:
            for h, q in L:
                _, h = h.primitive()
                eps += RootSum(
                    q, Lambda(t, t*log(h.as_expr())), quadratic=True)
        else:
            for h, q in L:
                _, h = h.primitive()
                R = log_to_real(h, q, x, t)

                if R is not None:
                    eps += R
                else:
                    eps += RootSum(
                        q, Lambda(t, t*log(h.as_expr())), quadratic=True)

        result += eps

    return coeff*result


def ratint_ratpart(f, g, x):
    """
    Horowitz-Ostrogradsky algorithm.

    Explanation
    ===========

    Given a field K and polynomials f and g in K[x], such that f and g
    are coprime and deg(f) < deg(g), returns fractions A and B in K(x),
    such that f/g = A' + B and B has square-free denominator.

    Examples
    ========

        >>> from symox.integrals.rationaltools import ratint_ratpart
        >>> from symox.abc import x, y
        >>> from sympy import Poly
        >>> ratint_ratpart(Poly(1, x, domain='ZZ'),
        ... Poly(x + 1, x, domain='ZZ'), x)
        (0, 1/(x + 1))
        >>> ratint_ratpart(Poly(1, x, domain='EX'),
        ... Poly(x**2 + y**2, x, domain='EX'), x)
        (0, 1/(x**2 + y**2))
        >>> ratint_ratpart(Poly(36, x, domain='ZZ'),
        ... Poly(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2, x, domain='ZZ'), x)
        ((12*x + 6)/(x**2 - 1), 12/(x**2 - x - 2))

    See Also
    ========

    ratint, ratint_logpart
    """
    f = Poly(f, x)
    g = Poly(g, x)

    from symox._core import _ratint_ratpart
    native = _ratint_ratpart(f, g, x)
    if native is not None:
        return native

    from symox.solvers.solvers import solve

    u, v, _ = g.cofactors(g.diff())

    n = u.degree()
    m = v.degree()

    A_coeffs = [ Dummy('a' + str(n - i)) for i in range(0, n) ]
    B_coeffs = [ Dummy('b' + str(m - i)) for i in range(0, m) ]

    C_coeffs = A_coeffs + B_coeffs

    A = Poly(A_coeffs, x, domain=ZZ[C_coeffs])
    B = Poly(B_coeffs, x, domain=ZZ[C_coeffs])

    H = f - A.diff()*v + A*(u.diff()*v).quo(u) - B*u

    result = solve(H.coeffs(), C_coeffs)

    A = A.as_expr().subs(result)
    B = B.as_expr().subs(result)

    rat_part = cancel(A/u.as_expr(), x)
    log_part = cancel(B/v.as_expr(), x)

    return rat_part, log_part


def ratint_logpart(f, g, x, t=None):
    r"""
    Lazard-Rioboo-Trager algorithm.

    Explanation
    ===========

    Given a field K and polynomials f and g in K[x], such that f and g
    are coprime, deg(f) < deg(g) and g is square-free, returns a list
    of tuples (s_i, q_i) of polynomials, for i = 1..n, such that s_i
    in K[t, x] and q_i in K[t], and::

                           ___    ___
                 d  f   d  \  `   \  `
                 -- - = --  )      )   a log(s_i(a, x))
                 dx g   dx /__,   /__,
                          i=1..n a | q_i(a) = 0

    Examples
    ========

    >>> from symox.integrals.rationaltools import ratint_logpart
    >>> from symox.abc import x
    >>> from sympy import Poly
    >>> ratint_logpart(Poly(1, x, domain='ZZ'),
    ... Poly(x**2 + x + 1, x, domain='ZZ'), x)
    [(Poly(x + 3*_t/2 + 1/2, x, domain='QQ[_t]'),
    ...Poly(3*_t**2 + 1, _t, domain='ZZ'))]
    >>> ratint_logpart(Poly(12, x, domain='ZZ'),
    ... Poly(x**2 - x - 2, x, domain='ZZ'), x)
    [(Poly(x - 3*_t/8 - 1/2, x, domain='QQ[_t]'),
    ...Poly(-_t**2 + 16, _t, domain='ZZ'))]

    See Also
    ========

    ratint, ratint_ratpart
    """
    f, g = Poly(f, x), Poly(g, x)

    t = t or Dummy('t')
    if (_ATAN_BOOTSTRAP_CAPTURE and Poly is _rational_core.Poly and Poly is _ATAN_POLY
            and _stock_lrt_resultant()
            and Poly.__init__ is object.__init__
            and _lrt_public.Poly is _rational_core.Poly
            and _lrt_constructor.Poly is _rational_core.Poly
            and _lrt_constructor.construct is _LRT_CONSTRUCT
            and _lrt_public.sympify is _LRT_SYMPIFY
            and type(_LRT_RESULTANT) is _lrt_types.BuiltinFunctionType
            and _LRT_RESULTANT.__name__ == 'resultant'
            and _LRT_RESULTANT.__module__ == 'symox._core'
            and _stock_lrt_as_poly()
            and _atan_evaluation.get()
            and (_ATAN_TYPE(t) is _ATAN_SCALARS[0] or _ATAN_TYPE(t) is _ATAN_SCALARS[1])
            and _atan_stock_operators(f,g,_LRT_OPERATORS,_LRT_BOOTSTRAP_STATES)):
        native = _rational_core._ratint_logpart(f, g, x, t)
        if native is not None:
            return native
    a, b = g, f - g.diff()*Poly(t, x)

    res, R = resultant(a, b, includePRS=True)
    res = Poly(res, t, composite=False)

    assert res, "BUG: resultant(%s, %s) cannot be zero" % (a, b)

    R_map, H = {}, []

    for r in R:
        R_map[r.degree()] = r

    def _include_sign(c, sqf):
        if c.is_extended_real and (c < 0) == True:
            h, k = sqf[0]
            c_poly = c.as_poly(h.gens)
            sqf[0] = h*c_poly, k

    C, res_sqf = res.sqf_list()
    _include_sign(C, res_sqf)

    for q, i in res_sqf:
        _, q = q.primitive()

        if g.degree() == i:
            H.append((g, q))
        else:
            h = R_map[i]
            h_lc = Poly(h.LC(), t, field=True)

            c, h_lc_sqf = h_lc.sqf_list(all=True)
            _include_sign(c, h_lc_sqf)

            for a, j in h_lc_sqf:
                h = h.quo(Poly(a.gcd(q)**j, x))

            inv, coeffs = h_lc.invert(q), [S.One]

            for coeff in h.coeffs()[1:]:
                coeff = coeff.as_poly(inv.gens)
                T = (inv*coeff).rem(q)
                coeffs.append(T.as_expr())

            h = Poly(dict(list(zip(h.monoms(), coeffs))), x)

            H.append((h, q))

    return H


def log_to_atan(f, g):
    """
    Convert complex logarithms to real arctangents.

    Explanation
    ===========

    Given a real field K and polynomials f and g in K[x], with g != 0,
    returns a sum h of arctangents of polynomials in K[x], such that:

                   dh   d         f + I g
                   -- = -- I log( ------- )
                   dx   dx        f - I g

    Examples
    ========

        >>> from symox.integrals.rationaltools import log_to_atan
        >>> from symox.abc import x
        >>> from sympy import Poly, sqrt, S
        >>> log_to_atan(Poly(x, x, domain='ZZ'), Poly(1, x, domain='ZZ'))
        2*atan(x)
        >>> log_to_atan(Poly(x + S(1)/2, x, domain='QQ'),
        ... Poly(sqrt(3)/2, x, domain='EX'))
        2*atan(2*sqrt(3)*x/3 + sqrt(3)/3)

    See Also
    ========

    log_to_real
    """
    if _stock_atan_protocol(f,g):
        native = _rational_core._log_to_atan(f, g)
        if native is not None:
            return native

    if f.degree() < g.degree():
        f, g = -g, f

    f = f.to_field()
    g = g.to_field()

    if g.degree() == 0 and f.domain.is_EX:
        # Division by a nonzero constant is exact. SymPy's EX elements are
        # canonicalised by cancel, so its f.div(g) leaves a zero remainder;
        # symox's EX arithmetic does not always recognise (a/c)*c - a as 0
        # for nested radicals, so divide by the constant directly.
        c = g.LC()
        p = Poly([coeff/c for coeff in f.all_coeffs()], *f.gens, domain=f.domain)
        return 2*atan(p.as_expr())

    p, q = f.div(g)

    if q.is_zero:
        return 2*atan(p.as_expr())
    else:
        s, t, h = g.gcdex(-f)
        u = (f*s + g*t).quo(h)
        A = 2*atan(u.as_expr())

        return A + log_to_atan(s, t)


_ATAN_RECURSION = log_to_atan


def _get_real_roots(f, x):
    """get real roots of f if possible"""
    rs = roots(f, filter='R')

    try:
        num_roots = f.count_roots()
    except DomainError:
        return rs
    else:
        if len(rs) == num_roots:
            return rs
        else:
            return None


def log_to_real(h, q, x, t):
    r"""
    Convert complex logarithms to real functions.

    Explanation
    ===========

    Given real field K and polynomials h in K[t,x] and q in K[t],
    returns real function f such that:
                          ___
                  df   d  \  `
                  -- = --  )  a log(h(a, x))
                  dx   dx /__,
                         a | q(a) = 0

    Examples
    ========

        >>> from symox.integrals.rationaltools import log_to_real
        >>> from symox.abc import x, y
        >>> from sympy import Poly, S
        >>> log_to_real(Poly(x + 3*y/2 + S(1)/2, x, domain='QQ[y]'),
        ... Poly(3*y**2 + 1, y, domain='ZZ'), x, y)
        2*sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/3
        >>> log_to_real(Poly(x**2 - 1, x, domain='ZZ'),
        ... Poly(-2*y + 1, y, domain='ZZ'), x, y)
        log(x**2 - 1)/2

    See Also
    ========

    log_to_atan
    """
    from symox.simplify.radsimp import collect
    u, v = symbols('u,v', cls=Dummy)

    H = h.as_expr().xreplace({t: u + I*v}).expand()
    Q = q.as_expr().xreplace({t: u + I*v}).expand()

    H_map = collect(H, I, evaluate=False)
    Q_map = collect(Q, I, evaluate=False)

    a, b = H_map.get(S.One, S.Zero), H_map.get(I, S.Zero)
    c, d = Q_map.get(S.One, S.Zero), Q_map.get(I, S.Zero)

    R = Poly(resultant(c, d, v), u)

    R_u = _get_real_roots(R, u)

    if R_u is None:
        return None

    result = S.Zero

    for r_u in R_u.keys():
        C = Poly(c.xreplace({u: r_u}), v)
        if not C:
            # t was split into real and imaginary parts
            # and denom Q(u, v) = c + I*d. We just found
            # that c(r_u) is 0 so the roots are in d
            C = Poly(d.xreplace({u: r_u}), v)
            # we were going to reject roots from C that
            # did not set d to zero, but since we are now
            # using C = d and c is already 0, there is
            # nothing to check
            d = S.Zero

        R_v = _get_real_roots(C, v)

        if R_v is None:
            return None

        R_v_paired = [] # take one from each pair of conjugate roots
        for r_v in R_v:
            if r_v not in R_v_paired and -r_v not in R_v_paired:
                if r_v.is_negative or r_v.could_extract_minus_sign():
                    R_v_paired.append(-r_v)
                elif not r_v.is_zero:
                    R_v_paired.append(r_v)

        for r_v in R_v_paired:

            D = d.xreplace({u: r_u, v: r_v})

            if D.evalf(chop=True) != 0:
                continue

            A = Poly(a.xreplace({u: r_u, v: r_v}), x)
            B = Poly(b.xreplace({u: r_u, v: r_v}), x)

            AB = (A**2 + B**2).as_expr()

            result += r_u*log(AB) + r_v*log_to_atan(A, B)

    R_q = _get_real_roots(q, t)

    if R_q is None:
        return None

    for r in R_q.keys():
        result += r*log(h.as_expr().subs(t, r))

    return result
