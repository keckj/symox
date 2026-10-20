"""Pure polynomials over the native polynomial representation.

Generator names are presentation data when two pure polynomial operands are
unified. Coefficient arithmetic and polynomial algorithms stay in the native
Poly backends; results retain the public PurePoly type.
"""
from functools import wraps
import inspect
from collections.abc import Mapping

from ._core import Poly as _Poly


_ALIGNED = frozenset(('add','sub','mul','div','rem','quo','exquo','pdiv','prem',
    'pquo','pexquo','gcd','lcm','cofactors','gcdex','half_gcdex','resultant',
    'subresultants','compose','invert','unify','rat_clear_denoms','eq','ne'))


def _restore(cls,terms,gens,domain):return cls(terms,*gens,domain=domain)


def _wrap(value,cls):
    if isinstance(value,_Poly):
        return value if isinstance(value,cls) else _Poly.__new__(cls,value)
    if isinstance(value,tuple):return tuple(_wrap(v,cls) for v in value)
    if isinstance(value,list):return [_wrap(v,cls) for v in value]
    if isinstance(value,dict):return {_wrap(k,cls):_wrap(v,cls) for k,v in value.items()}
    return value


class _PurePolyType(type):
    def __instancecheck__(cls,value):
        if super().__instancecheck__(value):return True
        from . import _CharPolyData
        return isinstance(value,_CharPolyData)


class PurePoly(_Poly,metaclass=_PurePolyType):
    """A native polynomial whose generators are ignored during unification."""
    __slots__=()

    def __new__(cls,rep,*gens,**args):
        return _Poly.__new__(cls,rep,*gens,**args)

    @classmethod
    def _from_coeffs(cls,coeffs,gen,domain,expr=None):
        # Preserve the existing charpoly data representation for coefficient
        # domains which do not yet have a native Poly backend. Their algorithms
        # remain an explicit parity gap rather than a mislabeled typed domain.
        from . import _CharPolyData
        def fallback():
            expression=expr
            if expression is None:
                from . import _core
                # Match the previous native CharPoly.as_expr conversion even
                # inside Python's evaluate(False) context.
                terms=tuple(_core._node('Mul',(coefficient,
                    _core._node('Pow',(gen,len(coeffs)-i-1),True)),True)
                    for i,coefficient in enumerate(coeffs))
                expression=_core._node('Add',terms,True)
            return _CharPolyData._from_coeffs(coeffs,gen,domain,expression)
        if str(domain) in ('EX','EXRAW'):
            return fallback()
        try:
            result=_Poly.__new__(cls,list(coeffs),gen,domain=domain)
            result.domain
            return result
        except (NotImplementedError,ValueError,TypeError):
            return fallback()

    def _with_gens(self,*gens):return _wrap(_Poly._with_gens(self,*gens),type(self))

    @property
    def func(self):return type(self)

    @property
    def free_symbols(self):
        return set().union(*(coefficient.free_symbols for coefficient in self.as_dict().values()))

    def as_expr(self,*gens):
        if not gens:return _Poly.as_expr(self)
        from . import Add,Mul,sympify
        from .polys.polyerrors import GeneratorsError
        if len(gens)==1 and isinstance(gens[0],Mapping):
            replacements=gens[0]
            if any(g not in self.gens for g in replacements):
                raise GeneratorsError('replacement generator is not in the polynomial')
            gens=tuple(replacements.get(g,g) for g in self.gens)
        elif len(gens)!=len(self.gens):
            raise ValueError('invalid number of generators')
        values=tuple(sympify(g,strict=True) for g in gens)
        return Add(*(coefficient*Mul(*(g**e for g,e in zip(values,monomial)))
                     for monomial,coefficient in self.as_dict().items()))

    def _aligned(self,other):
        if not isinstance(other,_Poly):return other
        if len(self.gens)!=len(other.gens):
            from . import UnificationFailed
            raise UnificationFailed('pure polynomials require the same number of generators')
        return _Poly._with_gens(other,*self.gens)

    def __eq__(self,other):
        from .polys.polyerrors import PolynomialError,DomainError,CoercionFailed
        try:
            if not isinstance(other,_Poly):other=_Poly(other,*self.gens,domain=self.domain)
            if len(self.gens)!=len(other.gens):return False
            other=self._aligned(other)
            domain=self.domain.unify(other.domain,symbols=self.gens)
            left=_Poly.set_domain(self,domain)
            right=_Poly.set_domain(other,domain)
            return left.as_dict()==right.as_dict()
        except (PolynomialError,DomainError,CoercionFailed,ValueError,TypeError):return False

    def __ne__(self,other):return not self==other

    def __hash__(self):
        # Exact-domain unification must not give equal PurePoly objects
        # different hashes. Comparisons through inexact rounding, or with
        # ordinary Poly objects, still inherit nontransitive reference rules.
        return hash((PurePoly,len(self.gens),tuple(sorted(self.as_dict().items()))))

    def __repr__(self):return 'Pure'+_Poly.__repr__(self)
    __str__=__repr__

    def __reduce__(self):return _restore,(type(self),self.as_dict(),self.gens,self.domain)
    def __reduce_ex__(self,protocol):return self.__reduce__()

    def to_sympy(self):
        # Explicit interchange only; native mathematical operations never
        # import the reference library.
        import sympy
        kwargs=({'modulus':self.modulus} if self.modulus is not None else
                {'domain':str(self.domain)})
        return sympy.PurePoly({monomial:coefficient.to_sympy()
                              for monomial,coefficient in self.as_dict().items()},
                             *(g.to_sympy() for g in self.gens),**kwargs)

    _sympy_=to_sympy


def _method(name,function):
    @wraps(function)
    def method(self,*args,**kwargs):
        if name in _ALIGNED:
            args=tuple(self._aligned(arg) for arg in args)
            kwargs={key:self._aligned(value) for key,value in kwargs.items()}
        return _wrap(function(self,*args,**kwargs),type(self))
    return method


def _binary(name,reverse=False):
    function=getattr(_Poly,name)
    @wraps(function)
    def operation(self,other):
        if reverse and isinstance(other,_Poly) and not isinstance(other,PurePoly):
            return getattr(_Poly,{'__radd__':'__add__','__rsub__':'__sub__',
                                  '__rmul__':'__mul__'}[name])(other,self)
        return _wrap(function(self,self._aligned(other)),type(self))
    return operation


def install():
    """Install after all Python Poly method adapters have been attached."""
    for name in dir(_Poly):
        if name.startswith('_') or name in PurePoly.__dict__:continue
        descriptor=inspect.getattr_static(_Poly,name)
        if isinstance(descriptor,classmethod):
            function=getattr(_Poly,name)
            def factory(cls,*args,_function=function,**kwargs):
                return _wrap(_function(*args,**kwargs),cls)
            factory=wraps(descriptor.__func__)(factory)
            setattr(PurePoly,name,classmethod(factory))
        elif isinstance(descriptor,property):
            def getter(self,_descriptor=descriptor):return _wrap(_descriptor.__get__(self,type(self)),type(self))
            setattr(PurePoly,name,property(getter,doc=descriptor.__doc__))
        else:
            function=getattr(_Poly,name)
            if callable(function):setattr(PurePoly,name,_method(name,function))
    for name in ('__add__','__sub__','__mul__','__floordiv__','__mod__','__divmod__'):
        setattr(PurePoly,name,_binary(name))
    for name in ('__radd__','__rsub__','__rmul__'):setattr(PurePoly,name,_binary(name,True))
    for name in ('__neg__','__pos__','__abs__','__pow__'):
        setattr(PurePoly,name,_method(name,getattr(_Poly,name)))
