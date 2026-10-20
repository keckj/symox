"""Einstein tensor expressions with native scalar coefficients."""
import itertools
import operator

from ._indexing import (_Object,TensorIndexType,TensorIndex,TensorSymmetry,
                        _IndexStructure,_index_key,_name_key,MAX_SLOTS)


def _scalar(value):
    from .. import sympify
    from .._core import Expr
    value=sympify(value)
    if not isinstance(value,Expr) or value.is_commutative is not True:
        raise TypeError('Tensor coefficients must be commutative native scalars')
    return value


def _split(expression):
    from .. import S
    if isinstance(expression,Tensor):return S.One,[expression]
    if isinstance(expression,TensMul):return expression.coeff,list(expression._factors)
    if isinstance(expression,TensExpr):raise NotImplementedError('Expand this tensor expression before taking its factors')
    return _scalar(expression),[]


def _normalize(factors):
    indices=[index for factor in factors for index in factor.indices]
    if len(indices)>MAX_SLOTS:raise NotImplementedError('Tensor products are limited to 32 slots')
    structure=_IndexStructure(indices)
    occupied={_name_key(index) for index,_ in structure.free};counters={};replacements={}
    for upper,lower in structure.dum:
        index=indices[upper];typ=index.tensor_index_type
        counter=counters.get(typ,0)
        while (typ,typ.dummy_name+'_'+str(counter)) in occupied:counter+=1
        name=typ.dummy_name+'_'+str(counter);counters[typ]=counter+1
        new=TensorIndex(name,typ)
        replacements[index]=new;replacements[-index]=-new
    return [Tensor._raw(factor.component,tuple(replacements.get(i,i) for i in factor.indices),type(factor)) for factor in factors]


def _product(coeff,factors,canonical=False):
    from .. import S
    coeff=_scalar(coeff)
    if coeff==0:return S.Zero
    if not factors:return coeff
    factors=_normalize(factors)
    if coeff==1 and len(factors)==1:
        result=factors[0];result._canonical=canonical;return result
    result=object.__new__(TensMul)
    result._args=tuple(([coeff] if coeff!=1 else [])+factors)
    result._coefficient=coeff;result._factors=tuple(factors);result._canonical=canonical
    return result


def _fresh_dummies(factors,occupied):
    structure=_IndexStructure([i for factor in factors for i in factor.indices])
    replacements={}
    for upper,lower in structure.dum:
        index=structure.indices[upper];counter=0
        while (index.tensor_index_type,'_tensor_internal_'+str(counter)) in occupied:counter+=1
        new=TensorIndex('_tensor_internal_'+str(counter),index.tensor_index_type)
        occupied.add(_name_key(new));replacements[index]=new;replacements[-index]=-new
    return [Tensor._raw(factor.component,tuple(replacements.get(i,i) for i in factor.indices),type(factor)) for factor in factors]


class TensorHead(_Object):
    __slots__=()
    def __new__(cls,name,index_types,symmetry=None,comm=0):
        index_types=tuple(index_types)
        if len(index_types)>MAX_SLOTS or any(not isinstance(t,TensorIndexType) for t in index_types):raise TypeError('Tensor slots require index types and rank at most 32')
        symmetry=TensorSymmetry.no_symmetry(len(index_types)) if symmetry is None else symmetry
        if not isinstance(symmetry,TensorSymmetry) or symmetry.rank!=len(index_types):raise ValueError('Symmetry rank does not match the tensor head')
        # A slot symmetry can exchange only slots belonging to the same type.
        for permutation in symmetry.generators:
            if any(index_types[i]!=index_types[j] for i,j in enumerate(permutation.array_form[:-2])):raise ValueError('Slot symmetries must preserve index types')
        comm=operator.index(comm)
        if comm not in (0,1,2):raise NotImplementedError('Only commuting, anticommuting and noncommuting tensor groups are supported')
        result=object.__new__(cls);result._args=(str(name),index_types,symmetry,comm);return result
    @property
    def name(self):return self.args[0]
    @property
    def index_types(self):return list(self.args[1])
    @property
    def symmetry(self):return self.args[2]
    @property
    def comm(self):return self.args[3]
    @property
    def rank(self):return len(self.args[1])
    @property
    def types(self):return sorted(set(self.index_types))
    def commutes_with(self,other):
        if not isinstance(other,TensorHead):return 0
        if self.comm==0 or other.comm==0:return 0
        return 1 if self.comm==other.comm==1 else None
    def __call__(self,*indices,**kwargs):return Tensor(self,indices,**kwargs)
    def __repr__(self):return '%s(%s)' % (self.name,','.join(t.name for t in self.index_types))
    def __lt__(self,other):return (self.name,repr(self.args))<(other.name,repr(other.args))
    def to_sympy(self):
        from sympy.tensor.tensor import TensorHead
        return TensorHead(self.name,[t.to_sympy() for t in self.index_types],self.symmetry.to_sympy(),self.comm)


def tensor_heads(s,index_types,symmetry=None,comm=0):
    from .. import symbols
    names=symbols(s,seq=True) if isinstance(s,str) else tuple(s)
    result=[TensorHead(name,index_types,symmetry,comm) for name in names]
    return result[0] if len(result)==1 else result


class TensExpr(_Object):
    __slots__=('_canonical',)
    is_commutative=False
    _symox_tensexpr=True
    @classmethod
    def fromiter(cls,arguments,**kwargs):return cls(*tuple(arguments),**kwargs)
    def __new__(cls,*args):
        result=object.__new__(cls);result._args=tuple(args);result._canonical=False;return result
    def __add__(self,other):return TensAdd(self,other)
    def __radd__(self,other):return TensAdd(other,self)
    def __sub__(self,other):return TensAdd(self,-other)
    def __rsub__(self,other):return TensAdd(other,-self)
    def __mul__(self,other):return TensMul(self,other)
    def __rmul__(self,other):return TensMul(other,self)
    def __neg__(self):return TensMul(-1,self)
    def __truediv__(self,other):
        if isinstance(other,TensExpr):raise ValueError('Division by an indexed tensor is undefined')
        return TensMul(1/_scalar(other),self)
    def __pow__(self,power):
        from .. import S
        power=operator.index(power)
        if power<0:raise ValueError('Negative tensor powers are undefined')
        if power==0:return S.One
        if power>32:raise NotImplementedError('Tensor powers exceed 32 factors')
        return TensMul(*([self]*power))
    @property
    def is_canon_bp(self):return self._canonical
    @property
    def rank(self):return len(self.get_free_indices())
    @property
    def ext_rank(self):return len(self.get_indices())
    @property
    def free_args(self):return sorted(self.get_free_indices(),key=_index_key)
    @property
    def free_indices(self):return set(self.get_free_indices())
    @property
    def index_types(self):return [i.tensor_index_type for i in self.get_indices()]
    @property
    def types(self):return sorted(set(self.index_types))
    @property
    def _index_structure(self):return _IndexStructure(self.get_indices())
    @property
    def free(self):return list(self._index_structure.free)
    @property
    def dum(self):return list(self._index_structure.dum)
    def get_indices(self):raise NotImplementedError('Use a concrete tensor expression')
    def get_free_indices(self):return self._index_structure.get_free_indices()
    @property
    def free_in_args(self):return self._argument_positions()[0]
    @property
    def dum_in_args(self):return self._argument_positions()[1]
    def _argument_positions(self):
        arguments=self.args if isinstance(self,TensMul) else (self,)
        positions=[]
        for component,arg in enumerate(arguments):
            if isinstance(arg,Tensor):positions.extend((slot,component) for slot in range(len(arg.indices)))
        structure=self._index_structure
        free=[(index,*positions[slot]) for index,slot in structure.free]
        dum=[(*positions[upper][:1],*positions[lower][:1],positions[upper][1],positions[lower][1]) for upper,lower in structure.dum]
        return free,dum
    def canon_bp(self):
        from ._canonical import canon_bp
        return canon_bp(self)
    def canonicalize(self):return self.canon_bp()
    def simplify(self,**kwargs):return self.canon_bp()
    def expand(self,**kwargs):return self
    def equals(self,other):return self==other or canon_equal(self,other)
    def contract_metric(self,g):
        from ._canonical import contract_metric
        return contract_metric(self,g)
    def contract_delta(self,delta):return self.contract_metric(delta)
    def substitute_indices(self,*pairs):
        from ._canonical import substitute_indices
        return substitute_indices(self,*pairs)
    def fun_eval(self,*pairs):
        from .tensor import deprecate_fun_eval
        deprecate_fun_eval();return self.substitute_indices(*pairs)
    def __call__(self,*indices):
        from .tensor import deprecate_call
        deprecate_call()
        if len(indices)!=self.rank:raise ValueError('Wrong number of free indices')
        return self.substitute_indices(*zip(self.free_args,indices))
    def _set_indices(self,*indices,is_canon_bp=False):
        coeff,factors=_split(self)
        if len(indices)!=sum(len(t.indices) for t in factors):raise ValueError('Index count mismatch')
        replaced=[];offset=0
        for tensor in factors:
            count=len(tensor.indices);new=indices[offset:offset+count]
            if any(index.tensor_index_type!=typ for index,typ in zip(new,tensor.component.index_types)):raise ValueError('Index types do not match tensor slots')
            replaced.append(Tensor._raw(tensor.component,new,type(tensor)));offset+=count
        return _product(coeff,replaced,is_canon_bp)
    def perm2tensor(self,g,is_canon_bp=False):
        from .tensor import perm2tensor
        return perm2tensor(self,g,is_canon_bp)
    def replace_with_arrays(self,replacement_dict,indices=None):
        from ._components import replace_with_arrays
        return replace_with_arrays(self,replacement_dict,indices)
    def matches(self,expr,repl_dict=None,old=False):
        from ._matching import matches
        return matches(self,expr,repl_dict,old)
    def match(self,pattern,old=False):return pattern.matches(self,old=old)


def canon_equal(left,right):
    from ._canonical import canon_bp
    try:return canon_bp(TensAdd(left,-right))==0
    except (ValueError,TypeError):return False


class Tensor(TensExpr):
    __slots__=()
    def __new__(cls,tensor_head,indices,*,is_canon_bp=False,**kwargs):
        if kwargs:raise TypeError('Unknown tensor options')
        if not isinstance(tensor_head,TensorHead):raise TypeError('A tensor requires a TensorHead')
        indices=tuple(indices)
        if len(indices)!=tensor_head.rank:raise ValueError('Wrong number of tensor indices')
        parsed=[]
        from .. import sympify
        for index,typ in zip(indices,tensor_head.index_types):
            if not isinstance(index,TensorIndex):
                index=sympify(index)
                if index.is_Symbol:index=TensorIndex(index,typ)
                elif index.is_Mul and index.as_coeff_Mul()[0]==-1 and index.as_coeff_Mul()[1].is_Symbol:index=-TensorIndex(index.as_coeff_Mul()[1],typ)
                else:raise TypeError('Tensor slots require typed indices or symbols')
            if index.tensor_index_type!=typ:raise ValueError('Index type does not match this tensor slot')
            parsed.append(index)
        result=Tensor._raw(tensor_head,parsed,cls)
        normalized=_normalize([result])[0];normalized._canonical=is_canon_bp;return normalized
    @staticmethod
    def _raw(head,indices,cls=None):
        result=object.__new__(Tensor if cls is None else cls)
        result._args=(head,tuple(indices));result._canonical=False;return result
    @property
    def component(self):return self.args[0]
    @property
    def components(self):return [self.component]
    @property
    def indices(self):return self.args[1]
    @property
    def coeff(self):
        from .. import S
        return S.One
    @property
    def nocoeff(self):return self
    def get_indices(self):return list(self.indices)
    def __repr__(self):return '%s(%s)' % (self.component.name,', '.join(map(str,self.indices)))
    __str__=__repr__
    def to_sympy(self):
        from sympy.tensor.tensor import Tensor
        return Tensor(self.component.to_sympy(),[i.to_sympy() for i in self.indices])


class TensMul(TensExpr):
    __slots__=('_coefficient','_factors')
    def __new__(cls,*args,**kwargs):
        from .. import S
        canonical=kwargs.pop('is_canon_bp',False);kwargs.pop('evaluate',None)
        if kwargs:raise TypeError('Unknown tensor product options')
        if any(isinstance(arg,TensAdd) for arg in args):
            choices=[arg.args if isinstance(arg,TensAdd) else (arg,) for arg in args]
            import math
            if math.prod(map(len,choices))>10000:raise NotImplementedError('Tensor expansion exceeds 10000 terms')
            return TensAdd(*(TensMul(*row) for row in itertools.product(*choices)))
        occupied={_name_key(i) for arg in args if isinstance(arg,TensExpr) for i in arg.get_indices()}
        coefficient=S.One;factors=[]
        for arg in args:
            coeff,part=_split(arg);coefficient*=coeff
            factors.extend(_fresh_dummies(part,occupied) if part else [])
        return _product(coefficient,factors,canonical)
    @property
    def coeff(self):return self._coefficient
    @property
    def nocoeff(self):
        from .. import S
        return _product(S.One,list(self._factors))
    @property
    def components(self):return [t.component for t in self._factors]
    def get_indices(self):return [i for factor in self._factors for i in factor.indices]
    def __repr__(self):return '*'.join(map(str,self.args))
    __str__=__repr__
    def to_sympy(self):
        from sympy.tensor.tensor import TensMul
        return TensMul(*(arg.to_sympy() for arg in self.args)).doit()


class TensAdd(TensExpr):
    __slots__=()
    def __new__(cls,*args,**kwargs):
        from .. import S,simplify
        kwargs.pop('evaluate',None)
        if kwargs:raise TypeError('Unknown tensor sum options')
        expanded=[term for arg in args for term in (arg.args if isinstance(arg,TensAdd) else (arg,))]
        if len(expanded)>10000:raise NotImplementedError('Tensor sum exceeds 10000 terms')
        terms={};scalars=S.Zero;signature=None
        for term in expanded:
            if not isinstance(term,TensExpr):
                scalar=_scalar(term)
                if scalar!=0:
                    if signature not in (None,frozenset()):raise ValueError('Cannot add scalars to a tensor with free indices')
                    signature=frozenset();scalars+=scalar
                continue
            free=frozenset(term.get_free_indices())
            if signature is not None and signature!=free:raise ValueError('All tensor sum terms must have the same free indices')
            signature=free;coeff=term.coeff;body=term.nocoeff
            terms[body]=terms.get(body,S.Zero)+coeff
        result=[]
        if scalars!=0:result.append(scalars)
        for body,coeff in terms.items():
            coeff=simplify(coeff)
            if coeff!=0:result.append(_product(coeff,_split(body)[1]))
        if not result:return S.Zero
        if len(result)==1:return result[0]
        value=object.__new__(cls);value._args=tuple(sorted(result,key=str));value._canonical=False;return value
    @property
    def coeff(self):raise ValueError('A tensor sum has no common scalar coefficient')
    @property
    def nocoeff(self):return self
    def get_free_indices(self):
        return next((arg.get_free_indices() for arg in self.args if isinstance(arg,TensExpr)),[])
    def get_indices(self):
        result=[]
        for arg in self.args:
            if isinstance(arg,TensExpr):
                result.extend(i for i in arg.get_indices() if i not in result)
        return result
    def __repr__(self):return ' + '.join(map(str,self.args))
    __str__=__repr__
    def to_sympy(self):
        from sympy.tensor.tensor import TensAdd
        return TensAdd(*(arg.to_sympy() for arg in self.args)).doit()


class TensorElement(TensExpr):
    __slots__=()
    def __new__(cls,expr,index_map):
        from .. import sympify
        if not isinstance(expr,TensExpr):raise TypeError('TensorElement requires a tensor expression')
        mapping=dict(index_map);free=expr.get_free_indices();names={i.name:i for i in free}
        parsed={names.get(str(index),index):sympify(value) for index,value in mapping.items()}
        parsed={i:v for i,v in parsed.items() if i in free}
        if not parsed:return expr
        for index,value in parsed.items():
            if value.is_integer is False or value.is_negative is True:raise ValueError('A component position must be a nonnegative integer')
            dim=index.tensor_index_type.dim
            if dim.is_Integer and value.is_Integer and value>=dim:raise ValueError('Component index exceeds its dimension')
        result=object.__new__(cls);result._args=(expr,tuple(sorted(parsed.items(),key=lambda pair:_index_key(pair[0]))));result._canonical=False;return result
    @property
    def expr(self):return self.args[0]
    @property
    def index_map(self):return dict(self.args[1])
    @property
    def coeff(self):
        from .. import S
        return S.One
    @property
    def nocoeff(self):return self
    def get_free_indices(self):return [i for i in self.expr.get_free_indices() if i not in self.index_map]
    def get_indices(self):return self.get_free_indices()
    def __repr__(self):return 'TensorElement(%s, %s)' % (self.expr,self.index_map)
    def to_sympy(self):
        from sympy.tensor.tensor import TensorElement
        return TensorElement(self.expr.to_sympy(),{i.to_sympy():v.to_sympy() for i,v in self.index_map.items()})
