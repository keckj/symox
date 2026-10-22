"""Vacuum boson/fermion Fock algebra, normal ordering and Wick contractions."""
from itertools import combinations, combinations_with_replacement
import symox as S
from ._algebra import Atom, OperatorExpr, _terms
from .quantum._core import (Operator, StateBase, KetBase, BraBase, Dagger,
    InnerProduct, qapply, FockSpace)

KroneckerDelta = S.KroneckerDelta


class SecondQuantizationError(Exception):pass
class AppliesOnlyToSymbolicIndex(SecondQuantizationError):pass
class ContractionAppliesOnlyToFermions(SecondQuantizationError):pass
class ViolationOfPauliPrinciple(SecondQuantizationError):pass
class SubstitutionOfAmbigousOperatorFailed(SecondQuantizationError):pass
class WicksTheoremDoesNotApply(SecondQuantizationError):pass


# Fermi-level restrictions of the dummies this module introduces. Native
# symbols carry no above_fermi/below_fermi assumption, so the quasi-particle
# dummies created by contractions and bracket removal are recorded here.
_FERMI_LEVEL = {}
_LEVEL_KEY = {'above_fermi':'above','below_fermi':'below'}


def _assumed(index, key):
    """Whether ``index`` carries the reference ``above_fermi``/``below_fermi`` assumption."""
    try:
        if getattr(index,'assumptions0',{}).get(key):return True
    except Exception:
        pass
    try:return _FERMI_LEVEL.get(index)==_LEVEL_KEY[key]
    except TypeError:return False


def _fermi_dummy(name, level, template=None):
    assumptions = {}
    if template is not None:
        try:assumptions = {k:v for k,v in template.assumptions0.items() if k not in _LEVEL_KEY}
        except Exception:assumptions = {}
    try:index = S.Dummy(name,**assumptions)
    except Exception:index = S.Dummy(name)
    _FERMI_LEVEL[index] = level
    return index


def _delta(a, b):
    """KroneckerDelta that vanishes across the Fermi level, like the reference."""
    a,b = S.sympify(a),S.sympify(b)
    if (_assumed(a,'below_fermi') and _assumed(b,'above_fermi')) or (_assumed(a,'above_fermi') and _assumed(b,'below_fermi')):
        return S.Integer(0)
    return KroneckerDelta(a,b)


class BosonicOperator(Operator):
    is_boson = True
    is_creator = False
    def __init__(self, state):self.state=S.sympify(state);super().__init__(self.state)
    @property
    def is_symbolic(self):return not self.state.is_Integer
    def _atom_repr(self):return ('CreateBoson' if self.is_creator else 'AnnihilateBoson')+'(%s)'%self.state
    def __str__(self):return ('b+' if self.is_creator else 'b')+'(%s)'%self.state
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
    def _atom_repr(self):return ('CreateFermion' if self.is_creator else 'AnnihilateFermion')+'(%s)'%self.state
    def __str__(self):return ('f+' if self.is_creator else 'f')+'(%s)'%self.state
    @property
    def is_symbolic(self):return not self.state.is_Integer
    @property
    def is_restricted(self):
        if _assumed(self.state,'below_fermi'):return -1
        if _assumed(self.state,'above_fermi'):return 1
        return 0
    @property
    def is_above_fermi(self):return not _assumed(self.state,'below_fermi')
    @property
    def is_below_fermi(self):return not _assumed(self.state,'above_fermi')
    @property
    def is_only_below_fermi(self):return self.is_below_fermi and not self.is_above_fermi
    @property
    def is_only_above_fermi(self):return self.is_above_fermi and not self.is_below_fermi
    # Quasi-particle character relative to the Fermi vacuum.
    @property
    def is_q_creator(self):
        if self.is_creator:return 1 if self.is_above_fermi else 0
        return -1 if self.is_below_fermi else 0
    @property
    def is_q_annihilator(self):
        if self.is_creator:return -1 if self.is_below_fermi else 0
        return 1 if self.is_above_fermi else 0
    @property
    def is_only_q_creator(self):return self.is_only_above_fermi if self.is_creator else self.is_only_below_fermi
    @property
    def is_only_q_annihilator(self):return self.is_only_below_fermi if self.is_creator else self.is_only_above_fermi
    def _sortkey(self):
        h=hash(self);label=str(self.state)
        if self.is_only_q_creator:return (1,label,h)
        if self.is_only_q_annihilator:return (4,label,h)
        if not self.is_creator:return (3,label,h)
        return (2,label,h)
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

def apply_operators(e):return qapply(e)


class Commutator:
    def __new__(cls,a,b):
        if a==b:return S.Integer(0)
        result=getattr(a,'_commutator',lambda b:None)(b)
        if result is not None:return result
        obj=super().__new__(cls);obj.a,obj.b=a,b;return obj
    def doit(self,**hints):return wicks(self.a*self.b-self.b*self.a)
    def _symox_(self):return S.sympify(self.doit())


def _boson_contraction(a,b):
    # Bosonic contractions are a symox extension: the reference raises.
    if a.is_creator or not b.is_creator:return S.Integer(0)
    return KroneckerDelta(a.state,b.state)


def contraction(a,b):
    """Contraction of two fermionic operators relative to the Fermi vacuum."""
    if isinstance(a,_FermionOperator) and isinstance(b,_FermionOperator):
        if not a.is_creator and b.is_creator:
            if _assumed(b.state,'below_fermi') or _assumed(a.state,'below_fermi'):return S.Integer(0)
            if _assumed(b.state,'above_fermi') or _assumed(a.state,'above_fermi'):return _delta(a.state,b.state)
            return _delta(a.state,b.state)*_delta(b.state,_fermi_dummy('a','above'))
        if a.is_creator and not b.is_creator:
            if _assumed(b.state,'above_fermi') or _assumed(a.state,'above_fermi'):return S.Integer(0)
            if _assumed(b.state,'below_fermi') or _assumed(a.state,'below_fermi'):return _delta(a.state,b.state)
            return _delta(a.state,b.state)*_delta(b.state,_fermi_dummy('i','below'))
        return S.Integer(0)
    if isinstance(a,BosonicOperator) and isinstance(b,BosonicOperator):return _boson_contraction(a,b)
    raise ContractionAppliesOnlyToFermions(*(isinstance(v,_FermionOperator) for v in (a,b)))


def _legacy_contraction(a,b):
    if a.is_boson!=b.is_boson or a.is_creator or not b.is_creator:return S.Integer(0)
    return KroneckerDelta(a.state,b.state)


def _normal_word(word, contractions):
    # Vacuum normal ordering, kept for bosonic strings (unsupported by the reference).
    for i in range(len(word)-1):
        a,b=word[i:i+2]
        if not isinstance(a,_LADDER) or not isinstance(b,_LADDER):raise TypeError('normal ordering accepts ladder-operator polynomials')
        if a.is_boson!=b.is_boson:raise NotImplementedError('mixed bosonic and fermionic normal ordering is not supported')
        same=a.is_creator==b.is_creator
        if same and not a.is_boson and a.state==b.state:return S.Integer(0)
        reverse=same and str(a.state)>str(b.state)
        if (not a.is_creator and b.is_creator) or reverse:
            sign=1 if a.is_boson else -1
            result=sign*_normal_word(word[:i]+(b,a)+word[i+2:],contractions)
            if contractions and not same:
                factor=_legacy_contraction(a,b)
                if factor!=0:result+=factor*_normal_word(word[:i]+word[i+2:],contractions)
            return result
    return _prod(word)


def _prod(factors):
    result=S.Integer(1)
    for atom in factors:result=result*atom
    return result


def _sqkey_operator(sq_operator):return sq_operator._sortkey()


def _sqkey_index(index):
    h=hash(index);label=str(index)
    if getattr(index,'is_Dummy',False):
        if _assumed(index,'above_fermi'):return (20,label,h)
        if _assumed(index,'below_fermi'):return (21,label,h)
        return (22,label,h)
    if _assumed(index,'above_fermi'):return (10,label,h)
    if _assumed(index,'below_fermi'):return (11,label,h)
    return (12,label,h)


def _sort_anticommuting_fermions(string1, key=_sqkey_operator):
    """Bubble sort of anticommuting fermions by ``key``, with the swap count."""
    verified=False;sign=0
    rng=list(range(len(string1)-1));rev=list(range(len(string1)-3,-1,-1))
    keys=list(map(key,string1));key_val=dict(zip(keys,string1))
    while not verified:
        verified=True
        for i in rng:
            left,right=keys[i],keys[i+1]
            if left==right:raise ViolationOfPauliPrinciple([left,right])
            if left>right:verified=False;keys[i:i+2]=[right,left];sign+=1
        if verified:break
        for i in rev:
            left,right=keys[i],keys[i+1]
            if left==right:raise ViolationOfPauliPrinciple([left,right])
            if left>right:verified=False;keys[i:i+2]=[right,left];sign+=1
    return [key_val[k] for k in keys],sign


class NO(Atom):
    """Normal-ordered product relative to the Fermi vacuum, as in the reference."""
    def __new__(cls,arg):
        if isinstance(arg,NO):return arg
        if not isinstance(arg,OperatorExpr):return S.sympify(arg)
        if isinstance(arg,Atom) and not isinstance(arg,NO):return arg
        result=S.Integer(0)
        for word,coefficient in arg._terms.items():result=result+coefficient*cls._from_word(word)
        if type(result) is OperatorExpr:result.__class__=_NormalOrderedSum
        return result
    def __init__(self,*args,**kwargs):pass
    @classmethod
    def _raw(cls,ops,legacy=False):
        obj=Atom.__new__(cls);obj.ops=tuple(ops);obj._legacy=legacy;Atom.__init__(obj,*obj.ops);return obj
    @classmethod
    def _from_word(cls,word):
        if not word:return S.Integer(1)
        seq=[];found=False
        for fac in word:
            if isinstance(fac,NO):seq.extend(fac.ops);found=True
            else:seq.append(fac)
        if found:return cls._from_word(tuple(seq))
        if len(seq)==1:return seq[0]
        if not all(isinstance(v,_FermionOperator) for v in seq):
            # Bosonic normal ordering is a symox extension of the reference.
            return cls._raw(seq,legacy=True)
        try:newseq,sign=_sort_anticommuting_fermions(seq)
        except ViolationOfPauliPrinciple:return S.Integer(0)
        result=cls._raw(newseq)
        return -result if sign%2 else result
    def _key(self):return type(self),self.ops
    @property
    def expression(self):return _prod(self.ops)
    def _parameter_symbols(self):return set().union(*(a._parameter_symbols() for a in self.ops))
    def _subs_parameters(self,mapping,**kwargs):return NO(_prod(a._subs_parameters(mapping,**kwargs) for a in self.ops))
    def _atom_repr(self):return 'NO(%s)'%'*'.join(a._atom_repr() for a in self.ops)
    def __repr__(self):return self._atom_repr()
    def __str__(self):return ':%s:'%'*'.join(a._atom_repr() for a in self.ops)
    @property
    def args(self):return (self.expression,)
    @property
    def has_q_creators(self):return self.ops[0].is_q_creator
    @property
    def has_q_annihilators(self):return self.ops[-1].is_q_annihilator
    def doit(self,**hints):
        if self._legacy:return _normal_word(self.ops,False)
        if hints.get('remove_brackets',True):return self._remove_brackets()
        return NO(_prod(a.doit(**hints) if hasattr(a,'doit') else a for a in self.ops))
    def _remove_brackets(self):
        subslist={}
        for i in self.iter_q_creators():
            op=self.ops[i]
            if op.is_q_annihilator:
                if getattr(op.state,'is_Dummy',False):
                    below=_fermi_dummy('i','below',op.state);above=_fermi_dummy('a','above',op.state)
                    cls=type(op)
                    subslist[i]=cls(below)*_delta(below,op.state)+cls(above)*_delta(above,op.state)
                else:raise SubstitutionOfAmbigousOperatorFailed(op)
        if not subslist:return self.expression
        result=NO(_prod(subslist.get(i,op) for i,op in enumerate(self.ops)))
        return _doit_terms(result)
    def _expand_operators(self):return NO(self._remove_brackets())
    def __getitem__(self,i):
        if isinstance(i,slice):return list(self.ops[i])
        return self.ops[i]
    def __len__(self):return len(self.ops)
    def iter_q_annihilators(self):
        for i in range(len(self.ops)-1,-1,-1):
            if self.ops[i].is_q_annihilator:yield i
            else:break
    def iter_q_creators(self):
        for i in range(len(self.ops)):
            if self.ops[i].is_q_creator:yield i
            else:break
    def get_subNO(self,i):return NO(_prod(self.ops[:i]+self.ops[i+1:]))


class _NormalOrderedSum(OperatorExpr):
    # A sum or multiple of normal-ordered products whose doit reaches the
    # products, as Mul.doit and Add.doit do in the reference.
    def doit(self,**hints):return _doit_terms(self,**hints)


def _doit_terms(expr,**hints):
    """Apply ``doit`` to the normal-ordered factors of every term."""
    if not isinstance(expr,OperatorExpr):return expr
    result=S.Integer(0)
    for word,coefficient in expr._terms.items():
        term=coefficient
        for atom in word:term=term*(atom.doit(**hints) if isinstance(atom,NO) else atom)
        result=result+term
    return result


def _fully_contracted(word):
    if not word:return S.Integer(1)
    result=S.Integer(0);first=word[0]
    for j in range(1,len(word)):
        factor=_legacy_contraction(first,word[j])
        if factor!=0:
            sign=1 if first.is_boson else (-1)**(j-1)
            result+=sign*factor*_fully_contracted(word[1:j]+word[j+1:])
    return result


def _get_contractions(string1, keep_only_fully_contracted=False):
    """All contractions of an operator string, as in the reference Wick expansion."""
    result=[] if keep_only_fully_contracted and string1 else [NO(_prod(string1))]
    for i in range(len(string1)-1):
        for j in range(i+1,len(string1)):
            c=contraction(string1[i],string1[j])
            if c!=0:
                coeff=-c if (j-i+1)%2 else c
                oplist=string1[i+1:j]+string1[j+1:]
                if oplist:result.append(coeff*NO(_prod(string1[:i])*_get_contractions(oplist,keep_only_fully_contracted)))
                else:result.append(coeff*NO(_prod(string1[:i])))
        if keep_only_fully_contracted:break
    total=S.Integer(0)
    for value in result:total=total+value
    return total


def _is_zero(e):
    if isinstance(e,OperatorExpr):return not e._terms
    try:return S.sympify(e)==0
    except Exception:return False


def wicks(e,**kw_args):
    """Wick expansion relative to the Fermi vacuum, as in the reference."""
    if _is_zero(e):return S.Integer(0)
    opts={'simplify_kronecker_deltas':False,'expand':True,'simplify_dummies':False,'keep_only_fully_contracted':False}
    opts.update(kw_args)
    if isinstance(e,(NO,_FermionOperator)):
        return S.Integer(0) if opts['keep_only_fully_contracted'] else e
    if not isinstance(e,OperatorExpr):return S.sympify(e)
    e=_doit_terms(e)
    if not isinstance(e,OperatorExpr):return e
    terms=list(e._terms.items())
    if len(terms)>1:
        result=S.Integer(0)
        for word,coefficient in terms:result=result+wicks(coefficient*_prod(word),**kw_args)
        return substitute_dummies(result) if opts['simplify_dummies'] else result
    word,coefficient=terms[0];n=len(word)
    if n==0:result=coefficient
    elif n==1:
        if opts['keep_only_fully_contracted']:return S.Integer(0)
        result=e
    elif isinstance(word[0],BosonicOperator):
        # Bosonic Wick expansion is a symox extension (the reference raises).
        if any(not isinstance(a,_LADDER) for a in word):raise TypeError('Wick expansion accepts ladder-operator polynomials')
        result=coefficient*(_fully_contracted(word) if opts['keep_only_fully_contracted'] else _normal_word(word,True))
    else:
        result=coefficient*_get_contractions(tuple(word),keep_only_fully_contracted=opts['keep_only_fully_contracted'])
    if opts['simplify_kronecker_deltas']:result=evaluate_deltas(result)
    return result


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

def matrix_rep(op,basis):
    basis=list(basis);result=S.zeros(len(basis))
    for j,ket in enumerate(basis):
        applied=apply_operators(op*ket)
        for i,bra in enumerate(map(Dagger,basis)):result[i,j]=apply_operators(bra*applied)
    return result


def _kd_above(d):
    a,b=d.args[:2];return not (_assumed(a,'below_fermi') or _assumed(b,'below_fermi'))


def _kd_below(d):
    a,b=d.args[:2];return not (_assumed(a,'above_fermi') or _assumed(b,'above_fermi'))


def _kd_preferred_position(d):
    if not _kd_above(d):return 0 if _assumed(d.args[0],'below_fermi') else 1
    if not _kd_below(d):return 0 if _assumed(d.args[0],'above_fermi') else 1
    return 0


def _kd_preferred_index(d):return d.args[1] if _kd_preferred_position(d) else d.args[0]
def _kd_killable_index(d):return d.args[0] if _kd_preferred_position(d) else d.args[1]


def _kd_equal_information(d):
    a,b=d.args[:2]
    if _assumed(a,'below_fermi') and _assumed(b,'below_fermi'):return True
    if _assumed(a,'above_fermi') and _assumed(b,'above_fermi'):return True
    return _kd_below(d) and _kd_above(d)


def _is_delta(node):return isinstance(node,S.Expr) and getattr(node,'class_name','')=='KroneckerDelta'


def _clean_deltas(e):
    """Zero the deltas that connect indices on opposite sides of the Fermi level."""
    if isinstance(e,OperatorExpr):
        result=S.Integer(0)
        for word,coefficient in e._terms.items():result=result+_clean_deltas(coefficient)*_prod(word)
        return result
    if not isinstance(e,S.Expr):return e
    replacements={n:S.Integer(0) for n in S.preorder_traversal(e) if _is_delta(n) and _delta(*n.args[:2])==0}
    return e.xreplace(replacements) if replacements else e


def _evaluate_mul(factors):
    """The reference Mul branch of evaluate_deltas on an explicit factor list."""
    indices={};deltas=[]
    for factor in factors:
        for s in factor.free_symbols:indices[s]=indices[s]+1 if s in indices else 0
        if _is_delta(factor):deltas.append(factor)
    e=_prod(factors)
    for d in deltas:
        killable,preferred=_kd_killable_index(d),_kd_preferred_index(d)
        if killable.is_Symbol and indices.get(killable):
            e=_clean_deltas(e.subs(killable,preferred))
            if len(deltas)>1:return evaluate_deltas(e)
        elif preferred.is_Symbol and indices.get(preferred) and _kd_equal_information(d):
            e=_clean_deltas(e.subs(preferred,killable))
            if len(deltas)>1:return evaluate_deltas(e)
    return e


def evaluate_deltas(e):
    """Eliminate KroneckerDelta factors over summed (repeated) indices, as in the reference."""
    if isinstance(e,OperatorExpr):
        result=S.Integer(0)
        for word,coefficient in e._terms.items():
            for term in S.Add.make_args(S.expand(coefficient)):
                result=result+_evaluate_mul(list(S.Mul.make_args(term))+list(word))
        return result
    e=S.sympify(e)
    if e.is_Add:return S.Add(*(evaluate_deltas(t) for t in e.args))
    if e.is_Mul:return _evaluate_mul(list(e.args))
    return e

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


def _has_index(expr,index):
    if isinstance(expr,OperatorExpr):return index in expr.free_symbols
    return S.sympify(expr).has(index)


class PermutationOperator(Atom):
    """Antisymmetrizing index permutation P(ij), as in the reference."""
    def __init__(self,i,j):
        i,j=sorted(map(S.sympify,(i,j)),key=S.default_sort_key)
        self.i,self.j=i,j;super().__init__(i,j)
    def get_permuted(self,expr):
        i,j=self.i,self.j
        if _has_index(expr,i) and _has_index(expr,j):
            tmp=S.Dummy()
            expr=expr.subs(i,tmp);expr=expr.subs(j,i);expr=expr.subs(tmp,j)
            return -expr
        return expr
    def doit(self,**hints):return self


def simplify_index_permutations(expr,permutation_operators):
    """Collect f(ij) - f(ji) pairs into P(ij)*f(ij), as in the reference."""
    def _get_indices(expr,ind):
        result=[]
        for arg in getattr(expr,'args',()):
            if arg in ind:result.append(arg)
            elif getattr(arg,'args',()):result.extend(_get_indices(arg,ind))
        return result
    def _choose_one_to_keep(a,b,ind):return min(a,b,key=lambda x:S.default_sort_key(S.Tuple(*_get_indices(x,ind))))
    if isinstance(expr,OperatorExpr):return expr
    expr=S.expand(expr)
    if not expr.is_Add:return expr
    # Terms are (applied permutation operators, scalar term) pairs; the
    # operators are kept outside the scalar so the terms stay hashable.
    terms={((),t) for t in expr.args}
    for P in permutation_operators:
        new_terms=set();on_hold=set()
        while terms:
            perms,term=terms.pop();permuted=(perms,P.get_permuted(term))
            if permuted in terms|on_hold:
                if permuted in terms:terms.remove(permuted)
                else:on_hold.remove(permuted)
                new_terms.add((perms+(P,),_choose_one_to_keep(term,permuted[1],(P.i,P.j))))
            else:
                permuted1=permuted;permuted=(perms,substitute_dummies(permuted[1]))
                if permuted1==permuted:on_hold.add((perms,term))
                elif permuted in terms|on_hold:
                    if permuted in terms:terms.remove(permuted)
                    else:on_hold.remove(permuted)
                    new_terms.add((perms+(P,),_choose_one_to_keep(term,permuted[1],(P.i,P.j))))
                else:new_terms.add((perms,term))
        terms=new_terms|on_hold
    result=S.Integer(0)
    for perms,term in terms:result=result+(_prod(perms)*term if perms else term)
    return result

__all__=['Dagger','KroneckerDelta','BosonicOperator','AnnihilateBoson','CreateBoson','AnnihilateFermion','CreateFermion','FockState','FockStateBra','FockStateKet','FockStateBosonKet','FockStateBosonBra','FockStateFermionKet','FockStateFermionBra','BBra','BKet','FBra','FKet','F','Fd','B','Bd','apply_operators','InnerProduct','BosonicBasis','VarBosonicBasis','FixedBosonicBasis','Commutator','matrix_rep','contraction','wicks','NO','evaluate_deltas','AntiSymmetricTensor','substitute_dummies','PermutationOperator','simplify_index_permutations','SecondQuantizationError','AppliesOnlyToSymbolicIndex','ContractionAppliesOnlyToFermions','ViolationOfPauliPrinciple','SubstitutionOfAmbigousOperatorFailed','WicksTheoremDoesNotApply']
