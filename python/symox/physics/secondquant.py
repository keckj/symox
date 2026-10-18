"""Vacuum boson/fermion Fock algebra, normal ordering and Wick contractions."""
from itertools import combinations, combinations_with_replacement
import symox as S
from ._algebra import Atom, OperatorExpr, _terms
from .quantum._core import (Operator, StateBase, KetBase, BraBase, Dagger,
    InnerProduct, qapply, FockSpace)

KroneckerDelta = S.KroneckerDelta


class BosonicOperator(Operator):
    is_boson = True
    is_creator = False
    def __init__(self, state):self.state=S.sympify(state);super().__init__(self.state)
    def _atom_repr(self):return ('Bd' if self.is_creator else 'B')+'(%s)'%self.state
    def _commutator(self, other):
        if not isinstance(other,BosonicOperator):return None
        if self.is_creator == other.is_creator:return S.Integer(0)
        return (-1 if self.is_creator else 1)*KroneckerDelta(self.state,other.state)
    def _apply_to(self, ket, **options):
        if not isinstance(ket,FockStateBosonKet) or not self.state.is_Integer:return None
        i=int(self.state)
        if not 0<=i<len(ket.occupations):raise ValueError('mode index lies outside the Fock state')
        values=list(ket.occupations);n=values[i]
        if self.is_creator:values[i]+=1;factor=S.sqrt(n+1)
        else:
            if n==0:return S.Integer(0)
            values[i]-=1;factor=S.sqrt(n)
        return factor*type(ket)(values)
class AnnihilateBoson(BosonicOperator):
    def _dual(self):return CreateBoson(self.state)
class CreateBoson(BosonicOperator):
    is_creator=True
    def _dual(self):return AnnihilateBoson(self.state)
class _FermionOperator(Operator):
    is_boson=False
    is_creator=False
    def __init__(self,state):self.state=S.sympify(state);super().__init__(self.state)
    def _dual(self):return AnnihilateFermion(self.state) if self.is_creator else CreateFermion(self.state)
    def _atom_repr(self):return ('Fd' if self.is_creator else 'F')+'(%s)'%self.state
    def _multiply_atom(self,other):
        if isinstance(other,type(self)) and self.state==other.state:return S.Integer(0),()
        return None
    def _apply_to(self,ket,**options):
        if not isinstance(ket,FockStateFermionKet) or not self.state.is_Integer:return None
        # Fermionic states are occupied orbital lists, in increasing order.
        orbital=self.state;values=list(ket.occupations)
        if self.is_creator:
            if orbital in values:return S.Integer(0)
            position=sum(int(v)<int(orbital) for v in values);values.insert(position,orbital)
        else:
            if orbital not in values:return S.Integer(0)
            position=values.index(orbital);values.pop(position)
        return (-1)**position*type(ket)(values,fermi_level=ket.fermi_level)
class AnnihilateFermion(_FermionOperator):pass
class CreateFermion(_FermionOperator):is_creator=True
_LADDER=(BosonicOperator,_FermionOperator)


class FockState(StateBase):
    def __init__(self,occupations):
        self.occupations=tuple(map(S.sympify,occupations));StateBase.__init__(self,*self.occupations);self.hilbert_space=FockSpace()
    def _subs_parameters(self,mapping,**kwargs):return type(self)([v.subs(mapping,**kwargs) for v in self.occupations])
    def _inner_with(self,ket):
        if type(self) not in (FockStateBosonBra,FockStateFermionBra) or self.statistics!=getattr(ket,'statistics',None):return S.Integer(0)
        if len(self.occupations)!=len(ket.occupations):return S.Integer(0)
        return S.Mul(*(KroneckerDelta(a,b) for a,b in zip(self.occupations,ket.occupations)))
class FockStateKet(FockState,KetBase):pass
class FockStateBra(FockState,BraBase):pass
class FockStateBosonKet(FockStateKet):
    statistics='boson'
    def __init__(self,occupations):
        if any(S.sympify(n).is_integer is False or S.sympify(n).is_negative is True for n in occupations):raise ValueError('bosonic occupation numbers must be nonnegative integers')
        super().__init__(occupations)
    def _dual(self):return FockStateBosonBra(self.occupations)
class FockStateBosonBra(FockStateBra):
    statistics='boson'
    def _dual(self):return FockStateBosonKet(self.occupations)
class FockStateFermionKet(FockStateKet):
    statistics='fermion'
    def __init__(self,occupations,fermi_level=0):
        values=tuple(map(S.sympify,occupations))
        if len(set(values))!=len(values):raise ValueError('a fermionic orbital cannot be occupied twice')
        if any(v.is_Integer is not True for v in values):raise NotImplementedError('fermionic Fock states currently require explicit integer orbitals')
        if tuple(sorted(values,key=int))!=values:raise ValueError('fermionic occupied orbitals must be sorted')
        self.fermi_level=fermi_level
        if fermi_level!=0:raise NotImplementedError('filled Fermi-sea states are not supported; use the empty-vacuum convention')
        super().__init__(values)
    def _dual(self):return FockStateFermionBra(self.occupations,fermi_level=self.fermi_level)
class FockStateFermionBra(FockStateBra):
    statistics='fermion'
    def __init__(self,occupations,fermi_level=0):
        ket=FockStateFermionKet(occupations,fermi_level);self.fermi_level=ket.fermi_level;super().__init__(ket.occupations)
    def _dual(self):return FockStateFermionKet(self.occupations,fermi_level=self.fermi_level)

BBra=FockStateBosonBra;BKet=FockStateBosonKet;FBra=FockStateFermionBra;FKet=FockStateFermionKet
B=AnnihilateBoson;Bd=CreateBoson;F=AnnihilateFermion;Fd=CreateFermion

def apply_operators(expression):return qapply(expression)


class Commutator:
    def __new__(cls,a,b):
        if a==b:return S.Integer(0)
        result=getattr(a,'_commutator',lambda b:None)(b)
        if result is not None:return result
        obj=super().__new__(cls);obj.a,obj.b=a,b;return obj
    def doit(self,**hints):return wicks(self.a*self.b-self.b*self.a)
    def _symox_(self):return S.sympify(self.doit())


def contraction(a,b):
    if not isinstance(a,_LADDER) or not isinstance(b,_LADDER):raise TypeError('contractions require creation or annihilation operators')
    if a.is_boson!=b.is_boson or a.is_creator or not b.is_creator:return S.Integer(0)
    return KroneckerDelta(a.state,b.state)


def _normal_word(word, contractions):
    for i in range(len(word)-1):
        a,b=word[i:i+2]
        if not isinstance(a,_LADDER) or not isinstance(b,_LADDER):raise TypeError('normal ordering accepts ladder-operator polynomials')
        if a.is_boson!=b.is_boson:raise NotImplementedError('mixed bosonic and fermionic normal ordering is not supported')
        same=a.is_creator==b.is_creator
        if same and not a.is_boson and a.state==b.state:return S.Integer(0)
        # Canonical orbital order gives deterministic antisymmetric fermion words.
        reverse=same and str(a.state)>str(b.state)
        if (not a.is_creator and b.is_creator) or reverse:
            sign=1 if a.is_boson else -1
            result=sign*_normal_word(word[:i]+(b,a)+word[i+2:],contractions)
            if contractions and not same:
                factor=contraction(a,b)
                if factor!=0:result+=factor*_normal_word(word[:i]+word[i+2:],contractions)
            return result
    result=S.Integer(1)
    for atom in word:result=result*atom
    return result


class NO(Atom):
    def __new__(cls,expression):
        if not isinstance(expression,OperatorExpr):return S.sympify(expression)
        return super().__new__(cls)
    def __init__(self,expression):self.expression=expression;super().__init__(expression)
    def _key(self):return type(self),tuple((w,c) for w,c in self.expression._terms.items())
    def doit(self,**hints):return sum((c*_normal_word(w,False) for w,c in self.expression._terms.items()),S.Integer(0))
    def _atom_repr(self):return ':'+str(self.expression)+':'
    def _subs_parameters(self,mapping,**kwargs):return NO(self.expression.subs(mapping,**kwargs))


def _fully_contracted(word):
    if not word:return S.Integer(1)
    result=S.Integer(0);first=word[0]
    for j in range(1,len(word)):
        factor=contraction(first,word[j])
        if factor!=0:
            sign=1 if first.is_boson else (-1)**(j-1)
            result+=sign*factor*_fully_contracted(word[1:j]+word[j+1:])
    return result

def wicks(expression,**kwargs):
    if isinstance(expression,NO):return expression.doit()
    result=S.Integer(0)
    for word,coefficient in _terms(expression).items():
        if any(not isinstance(a,_LADDER) for a in word):raise TypeError('Wick expansion accepts ladder-operator polynomials')
        value=_fully_contracted(word) if kwargs.get('keep_only_fully_contracted',False) else _normal_word(word,True)
        result+=coefficient*value
    return evaluate_deltas(result) if kwargs.get('simplify_kronecker_deltas',False) else result


class BosonicBasis:
    def __init__(self,states):self._states=list(states)
    def __len__(self):return len(self._states)
    def __getitem__(self,index):return self._states[index]
    def index(self,state):return self._states.index(state)
class VarBosonicBasis(BosonicBasis):
    def __init__(self,n_max):
        if int(n_max)!=n_max or n_max<1:raise ValueError('basis dimension must be positive')
        super().__init__([BKet([n]) for n in range(n_max)])
class FixedBosonicBasis(BosonicBasis):
    def __init__(self,n_particles,n_levels):
        if n_particles<0 or n_levels<1:raise ValueError('invalid particle or level count')
        states=[]
        for occupied in combinations_with_replacement(range(n_levels),n_particles):states.append(BKet([occupied.count(i) for i in range(n_levels)]))
        super().__init__(states)

def matrix_rep(operator,basis):
    basis=list(basis);result=S.zeros(len(basis))
    for j,ket in enumerate(basis):
        applied=apply_operators(operator*ket)
        for i,bra in enumerate(map(Dagger,basis)):result[i,j]=apply_operators(bra*applied)
    return result


def evaluate_deltas(expression):
    if isinstance(expression,OperatorExpr):return OperatorExpr._make({w:evaluate_deltas(c) for w,c in expression._terms.items()})
    expression=S.expand(expression)
    if expression.is_Add:return S.Add(*(evaluate_deltas(t) for t in expression.args))
    deltas=[n for n in S.preorder_traversal(expression) if getattr(n,'class_name','')=='KroneckerDelta']
    for delta in deltas:
        a,b=delta.args[:2]
        if a==b:expression=expression.xreplace({delta:S.Integer(1)});continue
        # Only eliminate a delta when the eliminated index appears elsewhere.
        rest=S.simplify(expression/delta)
        if a.is_Symbol and rest.has(a):expression=delta*rest.subs(a,b)
        elif b.is_Symbol and rest.has(b):expression=delta*rest.subs(b,a)
    return expression

def _antisym_sort(values):
    values=tuple(map(S.sympify,values))
    if len(set(values))!=len(values):return 0,values
    order=sorted(range(len(values)),key=lambda i:str(values[i]));inversions=sum(order[i]>order[j] for i in range(len(order)) for j in range(i+1,len(order)))
    return (-1)**inversions,tuple(values[i] for i in order)

def AntiSymmetricTensor(symbol,upper,lower):
    su,upper=_antisym_sort(upper);sl,lower=_antisym_sort(lower)
    if not su or not sl:return S.Integer(0)
    return su*sl*S._core._apply('AntiSymmetricTensor',S.Symbol(symbol) if isinstance(symbol,str) else S.sympify(symbol),S.Tuple(*upper),S.Tuple(*lower),evaluate=False)


def substitute_dummies(expression,new_indices=False,pretty_indices=None):
    # Native Dummy identity, unlike symbol spelling, determines dummy status.
    symbols=sorted((s for s in expression.free_symbols if getattr(s,'class_name','')=='Dummy'),key=str)
    mapping={s:S.Dummy(f'd_{i}') if new_indices else S.Symbol(f'd_{i}') for i,s in enumerate(symbols)}
    return expression.xreplace(mapping)


class PermutationOperator(Atom):
    def __init__(self,i,j):self.i,self.j=S.sympify(i),S.sympify(j);super().__init__(self.i,self.j)
    def get_permuted(self,expression):return expression-expression.xreplace({self.i:self.j,self.j:self.i})
    def doit(self,**hints):return self

def simplify_index_permutations(expression,permutation_operators):
    # Factor f(i,j)-f(j,i) against each requested antisymmetrizer.
    for p in permutation_operators:
        if isinstance(expression,OperatorExpr):continue
        for term in S.Add.make_args(S.expand(expression)):
            swapped=term.xreplace({p.i:p.j,p.j:p.i})
            if S.expand(expression).has(-swapped):
                expression=p*term+(expression-term+swapped);break
    return expression

__all__=['Dagger','KroneckerDelta','BosonicOperator','AnnihilateBoson','CreateBoson','AnnihilateFermion','CreateFermion','FockState','FockStateBra','FockStateKet','FockStateBosonKet','FockStateBosonBra','FockStateFermionKet','FockStateFermionBra','BBra','BKet','FBra','FKet','F','Fd','B','Bd','apply_operators','InnerProduct','BosonicBasis','VarBosonicBasis','FixedBosonicBasis','Commutator','matrix_rep','contraction','wicks','NO','evaluate_deltas','AntiSymmetricTensor','substitute_dummies','PermutationOperator','simplify_index_permutations']
