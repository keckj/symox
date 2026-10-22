"""Public sympification options and container protocols, without SymPy."""
from . import _core
from ._node_types import native_module
import io
import tokenize
import re

# User converters take precedence over the built-in container conversions.
converter={}
_unusual_space = re.compile(r'[^\S \t\f]')


def _unicode_space(source):
    # Retain the Rust parser's support for Unicode whitespace, while leaving
    # whitespace inside quoted strings untouched.
    text=source.replace('\n','')
    # Python's tokenizer already handles ordinary spaces, tabs and form feeds.
    # Only inspect token positions when some other whitespace needs replacing;
    # quoted Unicode/ASCII whitespace must still remain untouched.
    if _unusual_space.search(text) is None:
        return text
    try:
        positions=[token.start[1] for token in tokenize.generate_tokens(io.StringIO(text).readline)
                   if token.type==tokenize.ERRORTOKEN and token.string.isspace()]
    except tokenize.TokenError:return text
    for index in reversed(positions):text=text[:index]+' '+text[index+1:]
    return text


def _reference_object(a):
    # Objects of the reference library keep their native conversion path.
    for cls in type(a).__mro__:
        module=native_module(cls)
        if isinstance(module,str) and (module=='sympy' or module.startswith('sympy.')):return True
    return False


def _forwarded_expression(a):
    try:marker=getattr(a,'_symox_forwarded_',None)
    except Exception:return None
    if isinstance(marker,tuple) and len(marker)==1 and isinstance(marker[0],_core.Expr):return marker[0]
    return None


def make_sympify(original):
    def sympify(a,locals=None,convert_xor=True,strict=False,rational=False,evaluate=None):
        from ._public_core import SympifyError,_evaluation,Dict
        from ._basic import Basic
        if isinstance(a, Basic):return a
        if isinstance(a,(_core.Expr,Dict)):return a
        for cls in type(a).__mro__:
            conversion=converter.get(cls)
            if conversion is not None:return conversion(a)
        from .core.sympify import CantSympify
        if isinstance(a,CantSympify):raise SympifyError(a)
        if a is None:
            if strict:raise SympifyError(a)
            return None
        if isinstance(a,str):
            if strict:raise SympifyError(a)
            from .parsing.sympy_parser import (parse_expr,standard_transformations,
                                                convert_xor as xor_transform,rationalize)
            transformations=standard_transformations
            if convert_xor:transformations+=(xor_transform,)
            if rational:transformations+=(rationalize,)
            state=_evaluation.get() if evaluate is None else evaluate
            try:
                return parse_expr(_unicode_space(a),local_dict=locals,
                                  transformations=transformations,evaluate=state)
            except (SyntaxError,ValueError,TypeError,NotImplementedError,tokenize.TokenError) as error:
                if isinstance(error,SympifyError):raise
                raise SympifyError(a,error) from error
        from . import Tuple,FiniteSet,Dict,ImmutableMatrix,Float,Integer,Range
        from decimal import Decimal
        if isinstance(a,Decimal):return Float(a)
        if isinstance(a,range):return Range(a.start,a.stop,a.step)
        if isinstance(a,tuple):return Tuple(*(sympify(value) for value in a))
        if isinstance(a,(set,frozenset)):return FiniteSet(*(sympify(value) for value in a))
        if isinstance(a,dict):return Dict(a)
        if isinstance(a,list):
            if strict:raise SympifyError(a)
            return type(a)(sympify(value,locals=locals,convert_xor=convert_xor,
                                  rational=rational,evaluate=evaluate) for value in a)
        if getattr(a,'_symox_sparse',False):
            # SymPy keeps the sparse representation: sympify(SparseMatrix)
            # is an ImmutableSparseMatrix, which printers dispatch on.
            return a if getattr(a,'_immutable',False) else a.as_immutable()
        if isinstance(a,_core.Matrix):
            return ImmutableMatrix(a)
        if getattr(a,'is_NDimArray',False) and getattr(a,'_mutable',False):
            # SymPy converts mutable arrays to their immutable counterpart.
            return a.as_immutable()
        if not isinstance(a,type) and not _reference_object(a) and _forwarded_expression(a) is not None:
            # SymPy returns proxies forwarding attribute access to one of its
            # expressions unchanged (they expose its ``__sympy__`` marker).
            return a
        mpmath_number=_from_mpmath(a)
        if mpmath_number is not None:return mpmath_number
        try:
            return original(a)
        except (ValueError,TypeError) as error:
            module=native_module(type(a)) or ''
            if (module=='symox' or module.startswith('symox.')) and (hasattr(a,'args') or callable(a)):
                # Respect explicit _symox_ conversions first (units and
                # algebraic coefficient wrappers rely on them). Structural
                # objects without that conversion remain native; never call
                # their reference-library interoperability hook implicitly.
                return a
            if getattr(a,'__sympy__',None) is True:
                # A rejected symbolic object must retain its validation error;
                # converting it to float would discard assumptions and exactness.
                raise SympifyError(a,error) from error
            hook=getattr(a,'_sympy_',None)
            if hook is not None:
                result=hook()
                # SymPy returns the hook result unchanged, whatever its type;
                # a hook returning its input keeps the native cycle error.
                if result is a:raise
                return result
            if hasattr(a,'__float__'):return Float(float(a))
            if hasattr(a,'__int__'):return Integer(int(a))
            if type(error) not in (ValueError,TypeError):raise
            raise SympifyError(a,error) from error
    return sympify


def _from_mpmath(value):
    """SymPy's ``sympify_mpmath``: ``Expr._from_mpmath(x, x.context.prec)``."""
    module=type(value).__module__ or ''
    if not module.startswith('mpmath'):return None
    from . import Float, I
    prec=getattr(getattr(value,'context',None),'prec',None)
    if prec is None:return None
    if hasattr(value,'_mpf_'):
        return Float(value._mpf_,precision=prec)
    if hasattr(value,'_mpc_'):
        real,imag=value._mpc_
        return Float(real,precision=prec)+Float(imag,precision=prec)*I
    return None
