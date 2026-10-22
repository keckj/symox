"""Vacuum boson/fermion Fock algebra, normal ordering and Wick contractions."""
from itertools import combinations, combinations_with_replacement
import symox as S
from ._algebra import Atom, OperatorExpr, _terms
from .quantum._core import (Operator, StateBase, KetBase, BraBase, Dagger,
    qapply, FockSpace)

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
    assumptions[level+'_fermi'] = True
    try:return S.Dummy(name,**assumptions)
    except Exception:
        assumptions.pop(level+'_fermi')
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


class _SqAtoms:
    def atoms(self,*types):return set().union(*(S.sympify(v).atoms(*types) for v in self.labels))


class BosonicOperator(_SqAtoms,Operator):
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
    def apply_operator(self, state):
        if not self.is_symbolic and isinstance(state,FockStateKet):
            element=self.state
            if self.is_creator:return S.sqrt(state[element]+1)*state.up(element)
            return S.sqrt(state[element])*state.down(element)
        return self*state
    def _apply_to(self, ket, **options):
        # apply_operators leaves symbolic labels alone, like the reference _apply_Mul.
        if self.is_symbolic or not isinstance(ket,FockStateKet):return None
        return self.apply_operator(ket)
class AnnihilateBoson(BosonicOperator):
    def _dual(self):return CreateBoson(self.state)
class CreateBoson(BosonicOperator):
    is_creator=True
    def _dual(self):return AnnihilateBoson(self.state)
class _FermionOperator(_SqAtoms,Operator):
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
    def _ladder(self,ket):return ket.up(self.state) if self.is_creator else ket.down(self.state)
    def apply_operator(self,state):
        if isinstance(state,FockStateFermionKet):return self._ladder(state)
        if isinstance(state,OperatorExpr) and len(state._terms)==1:
            (word,coefficient),=state._terms.items()
            if word and isinstance(word[0],FockStateFermionKet):return coefficient*self._ladder(word[0])*_prod(word[1:])
        return self*state
    def _apply_to(self,ket,**options):
        if self.is_symbolic or not isinstance(ket,FockStateFermionKet):return None
        return self.apply_operator(ket)
class AnnihilateFermion(_FermionOperator):pass
class CreateFermion(_FermionOperator):is_creator=True
_LADDER=(BosonicOperator,_FermionOperator)


def _tuple_text(values):
    values=tuple(values)
    return '(%s,)'%values[0] if len(values)==1 else '(%s)'%', '.join(map(str,values))


class FockState(StateBase):
    lbracket=rbracket=''
    def __init__(self,occupations):
        if getattr(self,'_ready',False):return
        self.occupations=tuple(map(S.sympify,occupations));StateBase.__init__(self,*self.occupations);self.hilbert_space=FockSpace()
    def __getitem__(self,i):return self.occupations[int(i)]
    @property
    def args(self):return (S.Tuple(*self.occupations),)
    def __len__(self):return len(self.occupations)
    def _labels(self):return self.occupations
    def __repr__(self):return 'FockState(%s)'%_tuple_text(self.occupations)
    def __str__(self):return '%s%s%s'%(self.lbracket,_tuple_text(self._labels()),self.rbracket)
    # Products print the class-name form, like the reference StrPrinter.
    def _atom_repr(self):return '%s(%s)'%(type(self).__name__,_tuple_text(self.occupations))
    def _subs_parameters(self,mapping,**kwargs):return type(self)([v.subs(mapping,**kwargs) for v in self.occupations])
    def _inner_with(self,ket):
        if isinstance(ket,FockStateKet):return InnerProduct(self,ket)
        return StateBase._inner_with(self,ket)
class FockStateKet(FockState,KetBase):lbracket,rbracket='|','>'
class FockStateBra(FockState,BraBase):lbracket,rbracket='<','|'


class BosonState(FockState):
    def up(self,i):
        i=int(i);values=list(self.occupations);values[i]=values[i]+1
        return type(self)(values)
    def down(self,i):
        i=int(i);values=list(self.occupations)
        if values[i]==0:return S.Integer(0)
        values[i]=values[i]-1
        return type(self)(values)


def _only_below_fermi(i,fermi_level):
    i=S.sympify(i)
    if i.is_number:return bool(i<=fermi_level)
    return _assumed(i,'below_fermi')


def _only_above_fermi(i,fermi_level):
    i=S.sympify(i)
    if i.is_number:return bool(i>fermi_level)
    if _assumed(i,'above_fermi'):return True
    return not fermi_level


class FermionState(FockState):
    """Excitations of the Fermi sea: holes (orbitals at or below ``fermi_level``) and particles.

    The reference stores ``fermi_level`` on the class (the last construction
    wins); here every state keeps its own level.
    """
    def __new__(cls,occupations,fermi_level=0):
        occupations=list(map(S.sympify,occupations));sign=0
        if len(occupations)>1:
            try:occupations,sign=_sort_anticommuting_fermions(occupations,key=_sqkey_index)
            except ViolationOfPauliPrinciple:return S.Integer(0)
        if sum(_only_below_fermi(i,fermi_level) for i in occupations)>fermi_level:return S.Integer(0)
        obj=super().__new__(cls);obj.fermi_level=fermi_level
        FockState.__init__(obj,occupations);obj._ready=True
        return -obj if sign%2 else obj
    def __init__(self,*args,**kwargs):pass
    def up(self,i):
        i=S.sympify(i);present=i in self.occupations
        if _only_above_fermi(i,self.fermi_level):return S.Integer(0) if present else self._add_orbit(i)
        if _only_below_fermi(i,self.fermi_level):return self._remove_orbit(i) if present else S.Integer(0)
        if present:return KroneckerDelta(i,S.Dummy('i',below_fermi=True))*self._remove_orbit(i)
        return KroneckerDelta(i,S.Dummy('a',above_fermi=True))*self._add_orbit(i)
    def down(self,i):
        i=S.sympify(i);present=i in self.occupations
        if _only_above_fermi(i,self.fermi_level):return self._remove_orbit(i) if present else S.Integer(0)
        if _only_below_fermi(i,self.fermi_level):return S.Integer(0) if present else self._add_orbit(i)
        if present:return KroneckerDelta(i,S.Dummy('i',below_fermi=True))*self._add_orbit(i)
        return KroneckerDelta(i,S.Dummy('a',above_fermi=True))*self._remove_orbit(i)
    def _remove_orbit(self,i):
        values=list(self.occupations);pos=values.index(i);del values[pos]
        state=type(self)(values,self.fermi_level)
        return -state if pos%2 else state
    def _add_orbit(self,i):return type(self)((i,)+self.occupations,self.fermi_level)
    def _negate_holes(self,occupations):
        return tuple(-i if _only_below_fermi(i,self.fermi_level) and i.is_number else i for i in occupations)
    def _labels(self):return self._negate_holes(self.occupations)
    def __repr__(self):
        if self.fermi_level:return 'FockStateKet(%s, fermi_level=%s)'%(_tuple_text(self.occupations),self.fermi_level)
        return 'FockStateKet(%s)'%_tuple_text(self.occupations)
    def _subs_parameters(self,mapping,**kwargs):return type(self)([v.subs(mapping,**kwargs) for v in self.occupations],self.fermi_level)


class FockStateBosonKet(BosonState,FockStateKet):
    statistics='boson'
    def __init__(self,occupations):
        if any(S.sympify(n).is_integer is False or S.sympify(n).is_negative is True for n in occupations):raise ValueError('bosonic occupation numbers must be nonnegative integers')
        super().__init__(occupations)
    def _dual(self):return FockStateBosonBra(self.occupations)
class FockStateBosonBra(BosonState,FockStateBra):
    statistics='boson'
    def _dual(self):return FockStateBosonKet(self.occupations)
class FockStateFermionKet(FermionState,FockStateKet):
    statistics='fermion'
    def _dual(self):return FockStateFermionBra(self.occupations,self.fermi_level)
class FockStateFermionBra(FermionState,FockStateBra):
    statistics='fermion'
    def _dual(self):return FockStateFermionKet(self.occupations,self.fermi_level)

BBra=FockStateBosonBra;BKet=FockStateBosonKet;FBra=FockStateFermionBra;FKet=FockStateFermionKet
B=AnnihilateBoson;Bd=CreateBoson;F=AnnihilateFermion;Fd=CreateFermion


class InnerProduct:
    """Overlap of Fock states, evaluated at construction like the reference."""
    def __new__(cls,bra,ket):
        if not isinstance(bra,FockStateBra):raise TypeError('must be a bra')
        if not isinstance(ket,FockStateKet):raise TypeError('must be a ket')
        return cls.eval(bra,ket)
    @classmethod
    def eval(cls,bra,ket):
        result=S.Integer(1)
        for i,j in zip(bra.occupations,ket.occupations):
            result*=KroneckerDelta(i,j)
            if result==0:break
        return result


def apply_operators(e):return qapply(e)


def _sort_key(e):
    # Approximates the reference Basic.sort_key: products sort before operators.
    if isinstance(e,Atom):return ((5,0,type(e).__name__),(len(e.labels),tuple(_sort_key(v) for v in e.labels)))
    if isinstance(e,OperatorExpr):
        (word,coefficient),=e._terms.items()
        return ((3,0,'Mul'),(len(word),tuple(_sort_key(v) for v in word)))
    return ((0,0,''),S.default_sort_key(e))


class Commutator(Atom):
    """The reference commutator [A,B]; it evaluates at construction when it can."""
    def __new__(cls,a,b):
        result=cls.eval(a,b)
        if result is not None:return result
        obj=super().__new__(cls);Atom.__init__(obj,a,b);return obj
    def __init__(self,*args):pass
    @classmethod
    def eval(cls,a,b):
        if _is_zero(a) or _is_zero(b) or a==b:return S.Integer(0)
        if any(not isinstance(v,OperatorExpr) and S.sympify(v).is_commutative is not False for v in (a,b)):return S.Integer(0)
        terms_a,terms_b=_commutator_terms(a),_commutator_terms(b)
        if len(terms_a)>1:return _sum(cls(c*w,b) for c,w in terms_a)
        if len(terms_b)>1:return _sum(cls(a,c*w) for c,w in terms_b)
        (ca,wa),=terms_a;(cb,wb),=terms_b
        if wa==1 or wb==1:return S.Integer(0)
        if ca*cb!=1:return ca*cb*cls(wa,wb)
        a,b=wa,wb
        if isinstance(a,BosonicOperator) and isinstance(b,BosonicOperator):
            if isinstance(b,CreateBoson) and isinstance(a,AnnihilateBoson):return KroneckerDelta(a.state,b.state)
            if isinstance(a,CreateBoson) and isinstance(b,AnnihilateBoson):return -KroneckerDelta(a.state,b.state)
            return S.Integer(0)
        if isinstance(a,_FermionOperator) and isinstance(b,_FermionOperator):return wicks(a*b)-wicks(b*a)
        if _sort_key(a)>_sort_key(b):return -cls(b,a)
        return None
    def doit(self,**hints):
        a,b=self.labels
        if hints.get('wicks'):
            a=a.doit(**hints);b=b.doit(**hints)
            try:return wicks(a*b)-wicks(b*a)
            except (ContractionAppliesOnlyToFermions,WicksTheoremDoesNotApply):pass
        result=a*b-b*a
        return result.doit(**hints) if hasattr(result,'doit') else result
    def _parameter_symbols(self):return set().union(*(v.free_symbols for v in self.labels))
    def _subs_parameters(self,mapping,**kwargs):return Commutator(*(v.subs(mapping,**kwargs) for v in self.labels))
    def _atom_repr(self):return 'Commutator(%s, %s)'%self.labels
    def __repr__(self):return 'Commutator(%s,%s)'%self.labels
    def __str__(self):return '[%s,%s]'%self.labels


def _commutator_terms(e):
    """The expanded (commuting part, noncommuting part) terms of a commutator argument."""
    if isinstance(e,OperatorExpr):
        return [(c,_prod(w)) for w,c in e._terms.items()] or [(S.Integer(0),S.Integer(1))]
    result=[]
    for term in S.Add.make_args(S.expand(S.sympify(e))):
        c,nc=term.args_cnc();result.append((S.Mul(*c),S.Mul(*nc)))
    return result


def _sum(values):
    total=S.Integer(0)
    for value in values:total=total+value
    return total


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
    def subs(self,*args,**kwargs):
        # Operator replacements reach the factors, like the reference NO.subs.
        mapping=args[0] if len(args)==1 else {args[0]:args[1]}
        if isinstance(mapping,(list,tuple)):mapping=dict(mapping)
        if not any(isinstance(k,OperatorExpr) for k in mapping):return super().subs(mapping,**kwargs)
        scalars={k:v for k,v in mapping.items() if not isinstance(k,OperatorExpr)}
        result=NO(_prod(mapping.get(op,op) for op in self.ops))
        return result.subs(scalars,**kwargs) if scalars and isinstance(result,OperatorExpr) else result
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


def _expand_coefficients(e):
    if isinstance(e,OperatorExpr):
        terms={w:S.expand(c) for w,c in e._terms.items()}
        return e if all(terms[w]==c for w,c in e._terms.items()) else OperatorExpr(terms)
    return S.expand(e)


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
    if opts['expand']:result=_expand_coefficients(result)
    if opts['simplify_kronecker_deltas']:result=evaluate_deltas(result)
    return result


class BosonicBasis:pass


class VarBosonicBasis:
    """Single-mode bosonic basis with occupations 0, ..., n_max-1."""
    def __init__(self,n_max):
        self.n_max=n_max;self._build_states()
    def _build_states(self):
        self.basis=[FockStateBosonKet([i]) for i in range(self.n_max)];self.n_basis=len(self.basis)
    def index(self,state):return self.basis.index(state)
    def state(self,i):return self.basis[i]
    def __getitem__(self,i):return self.state(i)
    def __len__(self):return len(self.basis)
    def __repr__(self):return repr(self.basis)


class FixedBosonicBasis(BosonicBasis):
    """Basis of ``n_particles`` bosons over ``n_levels`` levels, in the reference order."""
    def __init__(self,n_particles,n_levels):
        self.n_particles=n_particles;self.n_levels=n_levels
        self._build_particle_locations();self._build_states()
    def _build_particle_locations(self):
        # The reference nests "for i0 in range(n_levels) for i1 in range(i0 + 1) ...".
        if self.n_particles==0:locations=[() for _ in range(self.n_levels)]
        else:
            locations=[(i,) for i in range(self.n_levels)]
            for _ in range(self.n_particles-1):locations=[t+(j,) for t in locations for j in range(t[-1]+1)]
        self.particle_locations=locations
    def _build_states(self):
        self.basis=[]
        for tuple_of_indices in self.particle_locations:
            occ_numbers=self.n_levels*[0]
            for level in tuple_of_indices:occ_numbers[level]+=1
            self.basis.append(FockStateBosonKet(occ_numbers))
        self.n_basis=len(self.basis)
    def index(self,state):return self.basis.index(state)
    def state(self,i):return self.basis[i]
    def __getitem__(self,i):return self.state(i)
    def __len__(self):return len(self.basis)
    def __repr__(self):return repr(self.basis)

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
            for term in S.Add.make_args(S.expand(_clean_deltas(coefficient))):
                result=result+_evaluate_mul(list(S.Mul.make_args(term))+list(word))
        return result
    e=_clean_deltas(S.sympify(e))
    if e.is_Add:return S.Add(*(evaluate_deltas(t) for t in e.args))
    if e.is_Mul:return _evaluate_mul(list(e.args))
    return e

class AntiSymmetricTensor(S.Function):
    """Tensor antisymmetric in its upper and in its lower indices, as in the reference."""
    # Host hooks of a symox class, so eval canonicalizes again after subs.
    _symox_user_hooks=True
    def __new__(cls,symbol,upper,lower):
        return super().__new__(cls,S.sympify(symbol),S.Tuple(*upper),S.Tuple(*lower))
    @classmethod
    def eval(cls,symbol,upper,lower):
        # Indices sort like the reference _sqkey_index; this reruns on subs.
        try:
            su,signu=_sort_anticommuting_fermions(list(upper),key=_sqkey_index)
            sl,signl=_sort_anticommuting_fermions(list(lower),key=_sqkey_index)
        except ViolationOfPauliPrinciple:
            return S.Integer(0)
        if tuple(su)==tuple(upper) and tuple(sl)==tuple(lower):return None
        result=cls(symbol,su,sl)
        return -result if (signu+signl)%2 else result
    @property
    def symbol(self):return self.args[0]
    @property
    def upper(self):return self.args[1]
    @property
    def lower(self):return self.args[2]


def _is_dummy(s):return getattr(s,'is_Dummy',False) is True


def _factor_dummies(factor):
    if isinstance(factor,OperatorExpr):return {s for s in factor.free_symbols if _is_dummy(s)}
    return set(S.sympify(factor).atoms(S.Dummy))


def _key_text(factor):
    """The reference KeyPrinter text: every Dummy prints as ``(name_index)``."""
    if isinstance(factor,_LADDER):
        return '%s(%s)'%(type(factor).__name__ if not isinstance(factor,BosonicOperator) else ('CreateBoson' if factor.is_creator else 'AnnihilateBoson'),_key_text(factor.state))
    if isinstance(factor,NO):return 'NO(%s)'%'*'.join(_key_text(op) for op in factor.ops)
    if isinstance(factor,OperatorExpr):
        return factor._subs_parameters({d:S.Symbol(_key_text(d)) for d in _factor_dummies(factor)})._atom_repr()
    factor=S.sympify(factor)
    if _is_dummy(factor):return '(%s_%i)'%(factor.name,factor.dummy_index)
    dummies=factor.atoms(S.Dummy)
    return str(factor.xreplace({d:S.Symbol(_key_text(d)) for d in dummies}) if dummies else factor)


def _term_factors(term):
    """Mul.make_args of one term: commuting factors, then the operator word."""
    if isinstance(term,OperatorExpr):
        if not term._terms:return [S.Integer(0)]
        if len(term._terms)>1:raise NotImplementedError('dummy ordering of a sum of operator words')
        (word,coefficient),=term._terms.items()
        return list(S.Mul.make_args(coefficient))+list(word)
    return list(S.Mul.make_args(S.sympify(term)))


def _level_mask(d):
    if _assumed(d,'below_fermi'):return '0'
    if _assumed(d,'above_fermi'):return '1'
    return '2'


def _has_dups(values):
    values=list(values);return len(set(values))!=len(values)


def _get_ordered_dummies(mul,verbose=False):
    """Canonical order of the dummies of one term, ported from the reference."""
    args=_term_factors(mul)
    fac_dum={};fac_repr={}
    for fac in args:
        fac_dum[fac]=_factor_dummies(fac);fac_repr[fac]=_key_text(fac)
    all_dums=sorted(set().union(*fac_dum.values()),key=S.default_sort_key)
    mask={d:_level_mask(d) for d in all_dums}
    dum_repr={d:_key_text(d) for d in all_dums}
    def _key(d):
        dumstruct=[fac for fac in fac_dum if d in fac_dum[fac]]
        other_dums=set().union(*[fac_dum[fac] for fac in dumstruct])
        other_dums.discard(d)
        masked_facs=[fac_repr[fac] for fac in dumstruct]
        for d2 in other_dums:masked_facs=[fac.replace(dum_repr[d2],mask[d2]) for fac in masked_facs]
        all_masked=[fac.replace(dum_repr[d],mask[d]) for fac in masked_facs]
        masked_facs=dict(zip(dumstruct,masked_facs))
        if _has_dups(all_masked):
            all_masked.sort()
            return mask[d],tuple(all_masked)
        keydict=dict(zip(dumstruct,all_masked))
        dumstruct.sort(key=lambda x:keydict[x])
        all_masked.sort()
        pos_val=[]
        for fac in dumstruct:
            if isinstance(fac,AntiSymmetricTensor):
                if d in fac.upper:pos_val.append('u')
                if d in fac.lower:pos_val.append('l')
            elif isinstance(fac,_LADDER):pos_val.append('u' if fac.is_creator else 'l')
            elif isinstance(fac,NO):
                for op in fac.ops:
                    if d in op.free_symbols:pos_val.append('u' if op.is_creator else 'l')
            else:
                facpos=-1
                while True:
                    facpos=masked_facs[fac].find(dum_repr[d],facpos+1)
                    if facpos==-1:break
                    pos_val.append(facpos)
        return mask[d],tuple(all_masked),pos_val[0],pos_val[-1]
    dumkey={d:_key(d) for d in all_dums}
    result=sorted(all_dums,key=lambda x:dumkey[x])
    if _has_dups(dumkey.values()):
        unordered={}
        for d,k in dumkey.items():unordered.setdefault(k,[]).append(d)
        unordered=[set(unordered[k]) for k in sorted(unordered) if len(unordered[k])>=2]
        result=_determine_ambiguous(mul,result,unordered)
    return result


def _subs_sequence(term,pairs):
    for old,new in pairs:term=term.subs(old,new)
    return term


_symbol_counter=[0]


def _determine_ambiguous(term,ordered,ambiguous_groups):
    all_ambiguous=set()
    for dummies in ambiguous_groups:all_ambiguous|=dummies
    all_ordered=set(ordered)-all_ambiguous
    if not all_ordered:
        group=[d for d in ordered if d in ambiguous_groups[0]]
        d=group[0];all_ordered.add(d);ambiguous_groups[0].remove(d)
    stored_counter=_symbol_counter[0];subslist=[]
    for d in [d for d in ordered if d in all_ordered]:
        subslist.append((d,S.Symbol('_]"]_%i'%_symbol_counter[0])));_symbol_counter[0]+=1
    newterm=_subs_sequence(term,subslist)
    neworder=_get_ordered_dummies(newterm)
    _symbol_counter[0]=stored_counter
    for group in ambiguous_groups:
        ordered_group=[d for d in neworder if d in group];ordered_group.reverse()
        result=[]
        for d in ordered:result.append(ordered_group.pop() if d in group else d)
        ordered=result
    return ordered


def _expanded_terms(expr):
    if isinstance(expr,OperatorExpr):
        return [t*_prod(word) for word,coefficient in expr._terms.items() for t in S.Add.make_args(S.expand(coefficient))]
    return list(S.Add.make_args(S.expand(S.sympify(expr))))


def substitute_dummies(expr,new_indices=False,pretty_indices={}):
    """Rename the dummies of every term canonically so equivalent terms combine, as in the reference."""
    if new_indices:
        letters_above=pretty_indices.get('above','');letters_below=pretty_indices.get('below','');letters_general=pretty_indices.get('general','')
        def _letter(letters,prefix,number):return letters[number] if number<len(letters) else '%s_%i'%(prefix,number-len(letters))
    aboves=[];belows=[];generals=[]
    dummies=[s for s in expr.free_symbols if _is_dummy(s)] if isinstance(expr,OperatorExpr) else list(S.sympify(expr).atoms(S.Dummy))
    if not new_indices:dummies=sorted(dummies,key=S.default_sort_key)
    counts={'a':0,'i':0,'p':0}
    for d in dummies:
        if _assumed(d,'above_fermi'):kind,group,letters='a',aboves,letters_above if new_indices else ''
        elif _assumed(d,'below_fermi'):kind,group,letters='i',belows,letters_below if new_indices else ''
        else:kind,group,letters='p',generals,letters_general if new_indices else ''
        if new_indices:
            group.append(S.Dummy(_letter(letters,kind,counts[kind]),**d.assumptions0));counts[kind]+=1
        else:group.append(d)
    new_terms=[]
    for term in _expanded_terms(expr):
        i=iter(belows);a=iter(aboves);p=iter(generals)
        subsdict={}
        for d in _get_ordered_dummies(term):
            if _assumed(d,'below_fermi'):subsdict[d]=next(i)
            elif _assumed(d,'above_fermi'):subsdict[d]=next(a)
            else:subsdict[d]=next(p)
        subslist=[];final_subs=[]
        for k,v in subsdict.items():
            if k==v:continue
            if v in subsdict:
                if subsdict[v] in subsdict:
                    x=S.Dummy('x');subslist.append((k,x));final_subs.append((x,v))
                else:final_subs.insert(0,(k,v))
            else:subslist.append((k,v))
        subslist.extend(final_subs)
        new_terms.append(_subs_sequence(term,subslist))
    return _sum(new_terms)


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
