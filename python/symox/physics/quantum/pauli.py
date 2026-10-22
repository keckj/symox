"""Pauli operators and the sigma_z eigenstates (port of sympy.physics.quantum.pauli)."""
import symox as S
from ._core import Operator,Ket,Bra,ComplexSpace,_DUALS
from .._algebra import Atom,OperatorExpr

__all__=['SigmaX','SigmaY','SigmaZ','SigmaMinus','SigmaPlus','SigmaZKet','SigmaZBra','qsimplify_pauli']


def _check_format(options):
    fmt=options.get('format','sympy')
    if fmt!='sympy':raise NotImplementedError('Representation in format %s not implemented.'%fmt)


class SigmaOpBase(Operator):
    """Pauli sigma operator, base class."""
    _matrix=None
    def __init__(self,*args,**hints):
        super().__init__(*(args or (False,)))
    @property
    def name(self):return self.labels[0]
    @property
    def use_name(self):return bool(self.labels[0]) is not False
    @classmethod
    def default_args(cls):return (False,)
    def _atom_repr(self):return type(self).__name__+'()'
    def _same(self,other):return isinstance(other,SigmaOpBase) and other.name==self.name
    def _multiply_atom(self,other):
        # The reference Mul merges equal neighbours into a power, which
        # _eval_power reduces: sigma_x**2 = 1, sigma_minus**2 = 0.
        if other==self:return (S.Integer(1) if self._involution else S.Integer(0)),()
        return None
    def _commutator(self,other):
        rule=self._commutators.get(type(other).__name__) if isinstance(other,SigmaOpBase) else None
        if rule is None:return None
        # Operators acting on differently named spins commute.  The reference
        # omits this check in some sigma_minus/sigma_plus rules.
        if not self._same(other):return S.Integer(0)
        return rule(self)
    def _anticommutator(self,other):
        rule=self._anticommutators.get(type(other).__name__) if isinstance(other,SigmaOpBase) else None
        # Differently named spins commute, so {A,B} = 2*A*B; the reference
        # returns the same-spin value regardless of the names.
        if rule is None or not self._same(other):return None
        return rule(self)
    def _represent_default_basis(self,**options):
        _check_format(options);return S.Matrix(self._matrix)
    def _apply_to(self,state,**options):
        if not isinstance(state,SigmaZKet):return None
        return self._on_ket(state)


class SigmaX(SigmaOpBase):
    is_hermitian=True
    _involution=True
    _matrix=[[0,1],[1,0]]
    _commutators={'SigmaY':lambda s:2*S.I*SigmaZ(s.name),'SigmaZ':lambda s:-2*S.I*SigmaY(s.name)}
    _anticommutators={'SigmaY':lambda s:S.Integer(0),'SigmaZ':lambda s:S.Integer(0)}
    def _on_ket(self,ket):return SigmaZKet(1) if ket.n==0 else SigmaZKet(0)


class SigmaY(SigmaOpBase):
    is_hermitian=True
    _involution=True
    _matrix=[[0,-S.I],[S.I,0]]
    _commutators={'SigmaZ':lambda s:2*S.I*SigmaX(s.name),'SigmaX':lambda s:-2*S.I*SigmaZ(s.name)}
    _anticommutators={'SigmaX':lambda s:S.Integer(0),'SigmaZ':lambda s:S.Integer(0)}
    def _on_ket(self,ket):return S.I*SigmaZKet(1) if ket.n==0 else -S.I*SigmaZKet(0)


class SigmaZ(SigmaOpBase):
    is_hermitian=True
    _involution=True
    _matrix=[[1,0],[0,-1]]
    _commutators={'SigmaX':lambda s:2*S.I*SigmaY(s.name),'SigmaY':lambda s:-2*S.I*SigmaX(s.name)}
    _anticommutators={'SigmaX':lambda s:S.Integer(0),'SigmaY':lambda s:S.Integer(0)}
    def _on_ket(self,ket):return ket if ket.n==0 else -1*ket


class SigmaMinus(SigmaOpBase):
    _involution=False
    _matrix=[[0,0],[1,0]]
    _commutators={'SigmaX':lambda s:-SigmaZ(s.name),'SigmaY':lambda s:S.I*SigmaZ(s.name),'SigmaZ':lambda s:2*s}
    _anticommutators={'SigmaZ':lambda s:S.Integer(0),'SigmaX':lambda s:S.Integer(1),'SigmaY':lambda s:-S.I,'SigmaPlus':lambda s:S.Integer(1)}
    def _dual(self):return SigmaPlus(self.name)
    def _on_ket(self,ket):return SigmaZKet(1) if ket.n==0 else S.Integer(0)


class SigmaPlus(SigmaOpBase):
    _involution=False
    _matrix=[[0,1],[0,0]]
    _commutators={'SigmaX':lambda s:SigmaZ(s.name),'SigmaY':lambda s:S.I*SigmaZ(s.name),'SigmaZ':lambda s:-2*s,'SigmaMinus':lambda s:SigmaZ(s.name)}
    _anticommutators={'SigmaZ':lambda s:S.Integer(0),'SigmaX':lambda s:S.Integer(1),'SigmaY':lambda s:S.I,'SigmaMinus':lambda s:S.Integer(1)}
    def _dual(self):return SigmaMinus(self.name)
    def _on_ket(self,ket):return S.Integer(0) if ket.n==0 else SigmaZKet(0)


class _SigmaZState:
    is_orthogonal=True
    def __init__(self,n):
        if n not in (0,1):raise ValueError('n must be 0 or 1')
        super().__init__(n);self.hilbert_space=ComplexSpace(2)
    @property
    def n(self):return self.labels[0]


class SigmaZKet(_SigmaZState,Ket):
    @classmethod
    def dual_class(cls):return SigmaZBra
    def _represent_default_basis(self,**options):
        _check_format(options);return S.Matrix([[1],[0]]) if self.n==0 else S.Matrix([[0],[1]])


class SigmaZBra(_SigmaZState,Bra):
    @classmethod
    def dual_class(cls):return SigmaZKet
    def _represent_default_basis(self,**options):
        _check_format(options);return S.Matrix([[1,0]]) if self.n==0 else S.Matrix([[0,1]])


_DUALS[SigmaZKet]=SigmaZBra
_DUALS[SigmaZBra]=SigmaZKet


def _half(a,sign):return (S.Integer(1)+sign*SigmaZ(a.name))/2


_PRODUCTS={
    ('SigmaX','SigmaY'):lambda a:S.I*SigmaZ(a.name),
    ('SigmaX','SigmaZ'):lambda a:-S.I*SigmaY(a.name),
    ('SigmaX','SigmaMinus'):lambda a:_half(a,1),
    ('SigmaX','SigmaPlus'):lambda a:_half(a,-1),
    ('SigmaY','SigmaX'):lambda a:-S.I*SigmaZ(a.name),
    ('SigmaY','SigmaZ'):lambda a:S.I*SigmaX(a.name),
    ('SigmaY','SigmaMinus'):lambda a:-S.I*_half(a,1),
    ('SigmaY','SigmaPlus'):lambda a:S.I*_half(a,-1),
    ('SigmaZ','SigmaX'):lambda a:S.I*SigmaY(a.name),
    ('SigmaZ','SigmaY'):lambda a:-S.I*SigmaX(a.name),
    ('SigmaZ','SigmaMinus'):lambda a:-SigmaMinus(a.name),
    ('SigmaZ','SigmaPlus'):lambda a:SigmaPlus(a.name),
    ('SigmaMinus','SigmaX'):lambda a:_half(a,-1),
    ('SigmaMinus','SigmaY'):lambda a:-S.I*_half(a,-1),
    ('SigmaMinus','SigmaZ'):lambda a:SigmaMinus(a.name),
    ('SigmaMinus','SigmaPlus'):lambda a:_half(a,-1),
    ('SigmaPlus','SigmaX'):lambda a:_half(a,1),
    ('SigmaPlus','SigmaY'):lambda a:S.I*_half(a,1),
    ('SigmaPlus','SigmaZ'):lambda a:-SigmaPlus(a.name),
    ('SigmaPlus','SigmaMinus'):lambda a:_half(a,1),
}


def _qsimplify_pauli_product(a,b):
    """Product of two Pauli operators of the same spin (reference helper)."""
    if a==b:return a*b
    return _PRODUCTS[(type(a).__name__,type(b).__name__)](a)


def _split(value):
    """``args_cnc`` of a product: (scalar, single operator) or (1, compound)."""
    if not isinstance(value,OperatorExpr):return S.sympify(value),None
    if isinstance(value,Atom):return S.Integer(1),value
    if len(value._terms)==1:
        (word,coefficient),=value._terms.items()
        if len(word)==1:return coefficient,word[0]
    return S.Integer(1),value


def qsimplify_pauli(e):
    """Simplify products of Pauli operators acting on the same spin."""
    if isinstance(e,Atom) or not isinstance(e,OperatorExpr):return e
    result=S.Integer(0)
    for word,coefficient in e._terms.items():
        factors=[];pending=list(word);scalar=coefficient
        while pending:
            current=pending.pop(0)
            while pending and isinstance(current,SigmaOpBase) and isinstance(pending[0],SigmaOpBase) and current.name==pending[0].name:
                factor,current=_split(_qsimplify_pauli_product(current,pending.pop(0)))
                scalar=scalar*factor
                if current is None:current=S.Integer(1)
            factors.append(current)
        term=scalar
        for factor in factors:term=term*factor
        result=result+term
    return result
