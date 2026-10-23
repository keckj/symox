"""Angular-momentum eigenstates, ladder actions, rotations and two-spin coupling."""
import symox as S
from ._core import (Operator, OrthogonalKet, OrthogonalBra, Atom, OperatorExpr,
    TensorProductAtom, TensorProduct, hbar, represent)
from ..wigner import clebsch_gordan, wigner_d, wigner_d_small


def m_values(j):
    j=S.sympify(j)
    if not (2*j).is_Integer or j<0:raise ValueError('spin must be an explicit nonnegative integer or half-integer')
    dimension=int(2*j+1);return dimension,[j-i for i in range(dimension)]


class _SpinState:
    axis='z'
    def _spin_init(self,j,m):
        self.j,self.m=S.sympify(j),S.sympify(m)
        if self.j.is_number and self.m not in m_values(self.j)[1]:raise ValueError('magnetic quantum number lies outside the spin multiplet')
    def _basis_axis(self,basis):
        if basis is None:return self.axis
        if isinstance(basis,type) and issubclass(basis,_SpinOperator):basis=basis.axis
        elif isinstance(basis,(_SpinOperator,_SpinState)):basis=basis.axis
        if basis not in ('x','y','z'):raise NotImplementedError('cannot represent a spin state in the basis %s'%(basis,))
        return basis
    def _represent_base(self,**options):
        # Reference SpinState._represent_base: column of D^j_{m',m}(alpha, beta, gamma).
        alpha,beta,gamma=(S.sympify(options.get(k,0)) for k in ('alpha','beta','gamma'))
        _,values=m_values(self.j)
        return S.Matrix([WignerD(self.j,mval,self.m,alpha,beta,gamma).doit() for mval in values])
    def _represent_ket(self,**options):
        angles=_BASIS_ANGLES[self.axis,self._basis_axis(options.pop('basis',None))]
        return self._represent_base(**angles,**options)
    def _represent(self,**options):
        # Bras use the conjugate transpose of the dual ket representation.
        if self.is_bra:return self._dual()._represent(**options).conjugate().T
        return self._represent_ket(**options)
    def rewrite(self,*args,**options):
        # Reference SpinState._eval_rewrite_as_J{x,y,z}: expand the state in
        # the eigenbasis of another axis.
        if not args:return self
        rule=args[-1]
        name=rule if isinstance(rule,str) else getattr(rule,'__name__',None)
        if name not in ('Jx','Jy','Jz'):return self
        axis=name[1]
        coupled=isinstance(self,_CoupledSpinState)
        if coupled:evect=(_COUPLED_BRA if self.is_bra else _COUPLED_KET)[axis]
        else:evect=(_BRA if self.is_bra else _KET)[axis]
        return self._rewrite_basis({'x':Jx,'y':Jy,'z':Jz}[axis],evect,**options)
    def _rewrite_basis(self,basis,evect,**options):
        j=self.j
        coupled=isinstance(self,_CoupledSpinState)
        if not j.is_number:
            # The reference returns Sum(WignerD*state, (mi, -j, j)), which
            # the native Sum cannot hold with quantum states.
            raise NotImplementedError('rewriting a spin state of symbolic j (a Sum over m values) is not supported')
        if coupled:
            start=j**2 if j==int(j) else (2*j-1)*(2*j+1)/4
            extra=(self.jn,tuple(self.jcoupling))
        else:
            start=0;extra=()
        rep_options={k:v for k,v in options.items() if k!='coupled'}
        vect=represent(self,basis=basis,**rep_options)
        result=S.Integer(0)
        for i in range(int(2*j+1)):
            result+=vect[int(start)+i]*evect(j,j-i,*extra)
        if coupled and options.get('coupled') is False:
            return uncouple(result)
        return result
class JzKet(OrthogonalKet,_SpinState):
    def __init__(self,j,m):self._spin_init(j,m);super().__init__(self.j,self.m)
    def _dual(self):return _BRA[self.axis](self.j,self.m)
    def _represent(self,**options):return _SpinState._represent(self,**options)
class JzBra(OrthogonalBra,_SpinState):
    def __init__(self,j,m):self._spin_init(j,m);super().__init__(self.j,self.m)
    def _dual(self):return _KET[self.axis](self.j,self.m)
    def _represent(self,**options):return _SpinState._represent(self,**options)
    def _inner_with(self,ket):
        if isinstance(ket,_SpinState):
            if self.j!=ket.j:return S.Integer(0)
            return (self._represent(basis="z")*ket._represent(basis="z"))[0]
        return super()._inner_with(ket)
class JxKet(JzKet):axis='x'
class JxBra(JzBra):axis='x'
class JyKet(JzKet):axis='y'
class JyBra(JzBra):axis='y'
_KET={'x':JxKet,'y':JyKet,'z':JzKet};_BRA={'x':JxBra,'y':JyBra,'z':JzBra}
_P=S.pi;_BASIS_ANGLES={('z','z'):{},('z','x'):{'beta':3*_P/2},('z','y'):{'alpha':3*_P/2,'beta':_P/2,'gamma':_P/2},
    ('x','x'):{},('x','y'):{'alpha':3*_P/2},('x','z'):{'beta':_P/2},
    ('y','y'):{},('y','x'):{'gamma':_P/2},('y','z'):{'alpha':3*_P/2,'beta':-_P/2,'gamma':_P/2}}


class _SpinOperator(Operator):
    axis='z';is_hermitian=True
    def _apply_to(self,state,**options):
        if not isinstance(state,_SpinState) or not state.is_ket:return None
        if self.axis=='2':return hbar**2*state.j*(state.j+1)*state
        if self.axis==state.axis:return hbar*state.m*state
        if state.axis=='z' and self.axis in ('+','-'):
            step=1 if self.axis=='+' else -1
            if state.m==step*state.j:return S.Integer(0)
            return hbar*S.sqrt(state.j*(state.j+1)-state.m*(state.m+step))*JzKet(state.j,state.m+step)
        if state.axis=='z' and self.axis in ('x','y'):
            up=Jplus._apply_to(state);down=Jminus._apply_to(state)
            return (up+down)/2 if self.axis=='x' else (up-down)/(2*S.I)
        # Reference _apply_op rewrites the state in the operator basis, applies
        # it there and rewrites the result back in the state's own basis.
        _,values=m_values(state.j)
        vector=self._represent(j=state.j)*state._represent(basis="z")
        change=S.Matrix.hstack(*(_KET[state.axis](state.j,m)._represent(basis="z") for m in values))
        coefficients=(change.conjugate().T*vector).applyfunc(S.simplify)
        return sum((coefficients[i]*_KET[state.axis](state.j,m) for i,m in enumerate(values)),S.Integer(0))
    def _dual(self):return Jminus if self.axis=='+' else Jplus if self.axis=='-' else self
    def _commutator(self,other):
        if not isinstance(other,_SpinOperator):return None
        if self.axis=='2' or other.axis=='2' or self.axis==other.axis:return S.Integer(0)
        cyclic={('x','y'):Jz,('y','z'):Jx,('z','x'):Jy}
        if (self.axis,other.axis) in cyclic:return S.I*hbar*cyclic[self.axis,other.axis]
        if (other.axis,self.axis) in cyclic:return -S.I*hbar*cyclic[other.axis,self.axis]
        if self.axis=='+' and other.axis=='-':return 2*hbar*Jz
        if self.axis=='-' and other.axis=='+':return -2*hbar*Jz
        return None
    def _represent(self,**options):
        # Reference spin operators are only representable in the Jz basis.
        basis=options.get('basis')
        if basis is not None:
            axis=basis.axis if isinstance(basis,(_SpinOperator,_SpinState)) or (isinstance(basis,type) and issubclass(basis,_SpinOperator)) else basis
            if axis!='z':raise NotImplementedError('spin operators are represented in the Jz basis only')
        j=S.sympify(options.get('j',S.Rational(1,2)));dimension,values=m_values(j)
        if self.axis=='2':return hbar**2*j*(j+1)*S.eye(dimension)
        if self.axis=='z':return hbar*S.diag(*values)
        up=S.zeros(dimension)
        for column,m in enumerate(values):
            if column>0:up[column-1,column]=hbar*S.sqrt(j*(j+1)-m*(m+1))
        if self.axis=='+':return up
        if self.axis=='-':return up.T
        return (up+up.T)/2 if self.axis=='x' else (up-up.T)/(2*S.I)
class JzOp(_SpinOperator):axis='z'
class J2Op(_SpinOperator):axis='2'
class JxOp(_SpinOperator):axis='x'
class JyOp(_SpinOperator):axis='y'
class JplusOp(_SpinOperator):axis='+';is_hermitian=False
class JminusOp(_SpinOperator):axis='-';is_hermitian=False
Jz=JzOp('Jz');J2=J2Op('J2');Jx=JxOp('Jx');Jy=JyOp('Jy');Jplus=JplusOp('J+');Jminus=JminusOp('J-')


class WignerD(S.Function):
    """Held Wigner D-function D^j_{m,m'}(alpha, beta, gamma) (a scalar)."""
    _symox_user_hooks=True
    _symox_simplify_doit=True
    is_commutative=True
    def __new__(cls,*args,**hints):
        if len(args)!=6:raise ValueError('6 parameters expected, got %s'%(args,))
        obj=super().__new__(cls,*map(S.sympify,args),evaluate=False)
        return obj.doit() if hints.get('evaluate',False) else obj
    @property
    def j(self):return self.args[0]
    @property
    def m(self):return self.args[1]
    @property
    def mp(self):return self.args[2]
    @property
    def alpha(self):return self.args[3]
    @property
    def beta(self):return self.args[4]
    @property
    def gamma(self):return self.args[5]
    def _sympystr(self,printer,*args):return 'WignerD(%s)'%', '.join(printer._print(a) for a in self.args)
    def doit(self,**hints):
        # Port of the reference WignerD._eval_wignerd (Varshalovich conventions).
        j,m,mp,alpha,beta,gamma=self.args
        if alpha==0 and beta==0 and gamma==0:return S.KroneckerDelta(m,mp)
        if not j.is_number:raise ValueError('j parameter must be numerical to evaluate, got %s'%j)
        r=S.Integer(0)
        if beta==S.pi/2:
            for k in range(int(2*j+1)):
                if k>j+mp or k>j-m or k<mp-m:continue
                r+=S.Integer(-1)**k*S.binomial(j+mp,k)*S.binomial(j-mp,k+m-mp)
            r*=S.Integer(-1)**(m-mp)/2**j*S.sqrt(S.factorial(j+m)*S.factorial(j-m)/(S.factorial(j+mp)*S.factorial(j-mp)))
        else:
            _,mvals=m_values(j)
            for mpp in mvals:
                r+=WignerD(j,m,mpp,0,S.pi/2,0).doit()*(S.cos(-mpp*beta)+S.I*S.sin(-mpp*beta))*WignerD(j,mpp,-mp,0,S.pi/2,0).doit()
            r=r*S.I**(2*j-m-mp)*S.Integer(-1)**(2*m)
            r=S.simplify(r)
        r*=S.exp(-S.I*m*alpha)*S.exp(-S.I*mp*gamma)
        return r
class Rotation(Operator):
    is_unitary=True
    def _atom_repr(self):return 'R(%s,%s,%s)'%(self.alpha,self.beta,self.gamma)
    def _eval_inverse(self):return Rotation(-self.gamma,-self.beta,-self.alpha)
    def __init__(self,alpha,beta,gamma):self.alpha,self.beta,self.gamma=map(S.sympify,(alpha,beta,gamma));super().__init__(self.alpha,self.beta,self.gamma)
    @staticmethod
    def D(j,m,mp,alpha,beta,gamma):return WignerD(j,m,mp,alpha,beta,gamma)
    @staticmethod
    def d(j,m,mp,beta):return WignerD(j,m,mp,0,beta,0)
    def _represent(self,**options):
        # Reference Rotation._represent_base: held D^j_{m,m'}(alpha, beta, gamma)
        # entries, evaluated only with the doit option.
        j=S.sympify(options.get('j',S.Rational(1,2)));_,values=m_values(j);evaluate=bool(options.get('doit'))
        def entry(p,q):
            value=WignerD(j,values[p],values[q],self.alpha,self.beta,self.gamma)
            return value.doit() if evaluate else value
        return S.Matrix(len(values),len(values),entry)
    def _apply_to(self,state,**options):
        # Reference _apply_operator_uncoupled: sum_mp D^j_{m,mp} |j,mp> in the state's basis.
        if not isinstance(state,JzKet) or not state.is_ket or isinstance(state,JzKetCoupled) or not state.j.is_number:return None
        values=m_values(state.j)[1]
        return sum((WignerD(state.j,state.m,mp,self.alpha,self.beta,self.gamma).doit()*type(state)(state.j,mp) for mp in values),S.Integer(0))
    def _dual(self):return Rotation(-self.gamma,-self.beta,-self.alpha)


def _build_coupled(jcoupling,length):
    # Reference _build_coupled: the spaces and j values of each coupling step.
    n_list=[[n+1] for n in range(length)];coupled_jn=[];coupled_n=[]
    for n1,n2,j_new in jcoupling:
        coupled_jn.append(j_new);coupled_n.append((n_list[n1-1],n_list[n2-1]))
        n_sort=sorted(n_list[n1-1]+n_list[n2-1]);n_list[n_sort[0]-1]=n_sort
    return coupled_n,coupled_jn


def _coupled_arguments(j,m,jn,jcoupling):
    # Reference CoupledSpinState.__new__ validation and default coupling scheme.
    if len(jcoupling)==0:
        jcoupling=[(1,n,S.Add(*[jn[i] for i in range(n)])) for n in range(2,len(jn))]
        jcoupling.append((1,len(jn),j))
    elif len(jcoupling)==1:jcoupling=jcoupling[0]
    else:raise TypeError('CoupledSpinState only takes 3 or 4 arguments, got: %s'%(len(jcoupling)+3))
    if not isinstance(jn,(list,tuple,S.Tuple)):raise TypeError('jn must be Tuple, list or tuple, got %s'%jn.__class__.__name__)
    if not isinstance(jcoupling,(list,tuple,S.Tuple)):raise TypeError('jcoupling must be Tuple, list or tuple, got %s'%jcoupling.__class__.__name__)
    if not all(isinstance(term,(list,tuple,S.Tuple)) for term in jcoupling):raise TypeError('All elements of jcoupling must be list, tuple or Tuple')
    if not len(jn)-1==len(jcoupling):raise ValueError('jcoupling must have length of %d, got %d'%(len(jn)-1,len(jcoupling)))
    if not all(len(x)==3 for x in jcoupling):raise ValueError('All elements of jcoupling must have length 3')
    j=S.sympify(j);m=S.sympify(m);jn=tuple(S.sympify(ji) for ji in jn)
    jcoupling=tuple((S.sympify(n1),S.sympify(n2),S.sympify(ji)) for n1,n2,ji in jcoupling)
    if any(2*ji!=int(2*ji) for ji in jn if ji.is_number):raise ValueError('All elements of jn must be integer or half-integer, got: %s'%(S.Tuple(*jn),))
    if any(n1!=int(n1) or n2!=int(n2) for n1,n2,_ in jcoupling):raise ValueError('Indices in jcoupling must be integers')
    if any(n1<1 or n2<1 or n1>len(jn) or n2>len(jn) for n1,n2,_ in jcoupling):raise ValueError('Indices must be between 1 and the number of coupled spin spaces')
    if any(2*ji!=int(2*ji) for _,_,ji in jcoupling if ji.is_number):raise ValueError('All coupled j values in coupling scheme must be integer or half-integer')
    jcoupling=tuple((int(n1),int(n2),ji) for n1,n2,ji in jcoupling)
    coupled_n,coupled_jn=_build_coupled(jcoupling,len(jn));jvals=list(jn)
    for n,(n1,n2) in enumerate(coupled_n):
        j1=jvals[min(n1)-1];j2=jvals[min(n2)-1];j3=coupled_jn[n]
        if S.sympify(j1).is_number and S.sympify(j2).is_number and S.sympify(j3).is_number:
            if j1+j2<j3:raise ValueError('All couplings must have j1+j2 >= j3, in coupling number %d got j1,j2,j3: %d,%d,%d'%(n+1,j1,j2,j3))
            if abs(j1-j2)>j3:raise ValueError('All couplings must have |j1+j2| <= j3, in coupling number %d got j1,j2,j3: %d,%d,%d'%(n+1,j1,j2,j3))
        jvals[min(n1+n2)-1]=j3
    if len(jcoupling)>0 and jcoupling[-1][2]!=j:raise ValueError('Last j value coupled together must be the final j of the state')
    return j,m,jn,jcoupling


class _CoupledSpinState:
    def _coupled_init(self,j,m,jn,jcoupling):
        j,m,self.jn,self.jcoupling=_coupled_arguments(j,m,jn,jcoupling)
        self._spin_init(j,m)
        Atom.__init__(self,self.j,self.m,S.Tuple(*self.jn),S.Tuple(*(S.Tuple(*t) for t in self.jcoupling)))
    @property
    def coupling(self):return S.Tuple(*(S.Tuple(*t) for t in self.jcoupling))
    @property
    def coupled_jn(self):return _build_coupled(self.jcoupling,len(self.jn))[1]
    @property
    def coupled_n(self):return _build_coupled(self.jcoupling,len(self.jn))[0]
    def _coupled_label(self):
        label=['%s'%self.j,'%s'%self.m]+['j%d=%s'%(i,ji) for i,ji in enumerate(self.jn,start=1)]
        for ji,(n1,n2) in zip(self.coupled_jn[:-1],self.coupled_n[:-1]):
            label.append('j(%s)=%s'%(','.join(str(i) for i in sorted(n1+n2)),ji))
        return ','.join(label)
    def _subs_parameters(self,mapping,**kwargs):
        sub=lambda v:S.sympify(v).subs(mapping,**kwargs)
        return type(self)(sub(self.j),sub(self.m),[sub(j) for j in self.jn],[(a,b,sub(c)) for a,b,c in self.jcoupling])
    def _inner_with(self,ket):
        # Reference coupled states define no inner products; <..|..> stays held.
        from ._core import InnerProduct
        return InnerProduct(self,ket)


class JzKetCoupled(_CoupledSpinState,JzKet):
    def __init__(self,j,m,jn,*jcoupling):self._coupled_init(j,m,jn,jcoupling)
    def _dual(self):return _COUPLED_BRA[self.axis](self.j,self.m,self.jn,self.jcoupling)
    @classmethod
    def uncoupled_class(cls):return _KET[cls.axis]
    def _atom_repr(self):return '|%s>'%self._coupled_label()
    def _represent_ket(self,**options):
        # Reference CoupledSpinState._represent_coupled_base: the j block of the
        # direct sum of the multiplets of total spin sum(jn), sum(jn)-1, ...
        total=S.Add(*self.jn)
        if not self.j.is_number:raise ValueError('State must not have symbolic j value to represent')
        if not total.is_number:raise ValueError('State must not have symbolic j values to represent')
        dimension=sum(range(int(2*total+1),0,-2));result=S.zeros(dimension,1)
        start=int(self.j**2 if self.j.is_integer else (2*self.j-1)*(1+2*self.j)/4)
        block=_SpinState._represent_ket(self,**options)
        for i in range(block.rows):result[start+i,0]=block[i,0]
        return result
class JzBraCoupled(_CoupledSpinState,JzBra):
    def __init__(self,j,m,jn,*jcoupling):self._coupled_init(j,m,jn,jcoupling)
    def _dual(self):return _COUPLED_KET[self.axis](self.j,self.m,self.jn,self.jcoupling)
    @classmethod
    def uncoupled_class(cls):return _BRA[cls.axis]
    def _represent(self,**options):return self._dual()._represent(**options).conjugate().T
    def _atom_repr(self):return '<%s|'%self._coupled_label()
class JxKetCoupled(JzKetCoupled):axis='x'
class JxBraCoupled(JzBraCoupled):axis='x'
class JyKetCoupled(JzKetCoupled):axis='y'
class JyBraCoupled(JzBraCoupled):axis='y'
_COUPLED_KET={'x':JxKetCoupled,'y':JyKetCoupled,'z':JzKetCoupled};_COUPLED_BRA={'x':JxBraCoupled,'y':JyBraCoupled,'z':JzBraCoupled}
for _cls in (JzKet,JxKet,JyKet):_cls.coupled_class=classmethod(lambda cls:_COUPLED_KET[cls.axis])
for _cls in (JzBra,JxBra,JyBra):_cls.coupled_class=classmethod(lambda cls:_COUPLED_BRA[cls.axis])


def _confignum_to_difflist(config_num,diff,list_len):
    diff_list=[]
    for n in range(list_len):
        prev_diff=diff;rem_spots=list_len-n-1;rem_configs=S.binomial(diff+rem_spots-1,diff)
        while config_num>=rem_configs:
            config_num-=rem_configs;diff-=1;rem_configs=S.binomial(diff+rem_spots-1,diff)
        diff_list.append(prev_diff-diff)
    return diff_list


def _cg(j1,m1,j2,m2,j3,m3):return clebsch_gordan(j1,j2,j3,m1,m2,m3)


def _uncouple(state,jn,jcoupling_list):
    # Port of the reference _uncouple.
    if isinstance(state,_CoupledSpinState):
        jn=state.jn;coupled_n=state.coupled_n;coupled_jn=state.coupled_jn;evect=state.uncoupled_class()
    elif isinstance(state,_SpinState):
        if jn is None:raise ValueError('Must specify j-values for coupled state')
        if not isinstance(jn,(list,tuple)):raise TypeError('jn must be list or tuple')
        jn=tuple(map(S.sympify,jn))
        if jcoupling_list is None:jcoupling_list=[(1,1+i,S.Add(*[jn[k] for k in range(i+1)])) for i in range(1,len(jn))]
        if not isinstance(jcoupling_list,(list,tuple)):raise TypeError('jcoupling must be a list or tuple')
        if not len(jcoupling_list)==len(jn)-1:raise ValueError('Must specify 2 fewer coupling terms than the number of j values')
        coupled_n,coupled_jn=_build_coupled(jcoupling_list,len(jn));evect=type(state)
    else:raise TypeError('state must be a spin state')
    j,m=state.j,state.m;coupling_list=[];j_list=list(jn)
    for j3,(n1,n2) in zip(coupled_jn,coupled_n):
        j1=j_list[n1[0]-1];j2=j_list[n2[0]-1];coupling_list.append((n1,n2,j1,j2,S.sympify(j3)));j_list[min(n1+n2)-1]=j3
    if not (j.is_number and m.is_number and all(x.is_number for x in jn)):
        # The reference returns a held Sum of CG coefficients times states.
        raise NotImplementedError('uncoupling symbolic spins (a Sum over m values) is not supported')
    diff_max=[2*x for x in jn];diff=S.Add(*jn)-m;n=len(jn);result=S.Integer(0)
    for config_num in range(int(S.binomial(diff+n-1,diff))):
        diff_list=_confignum_to_difflist(config_num,diff,n)
        if any(d>p for d,p in zip(diff_list,diff_max)):continue
        coefficient=S.Integer(1)
        for j1_n,j2_n,j1,j2,j3 in coupling_list:
            m1=S.Add(*[jn[x-1]-diff_list[x-1] for x in j1_n]);m2=S.Add(*[jn[x-1]-diff_list[x-1] for x in j2_n])
            coefficient*=_cg(j1,m1,j2,m2,j3,m1+m2)
        if coefficient:result+=coefficient*TensorProduct(*[evect(ji,ji-d) for ji,d in zip(jn,diff_list)])
    return result


def _couple(tp,jcoupling_list):
    # Port of the reference _couple.
    states=tp.factors;coupled_evect=_COUPLED_KET[states[0].axis] if states[0].is_ket else _COUPLED_BRA[states[0].axis]
    if jcoupling_list is None:jcoupling_list=[(1,n+1) for n in range(1,len(states))]
    if not len(jcoupling_list)==len(states)-1:raise TypeError('jcoupling_list must be length %d, got %d'%(len(states)-1,len(jcoupling_list)))
    if not all(len(coupling)==2 for coupling in jcoupling_list):raise ValueError('Each coupling must define 2 spaces')
    if any(n1==n2 for n1,n2 in jcoupling_list):raise ValueError('Spin spaces cannot couple to themselves')
    if all(S.sympify(n1).is_number and S.sympify(n2).is_number for n1,n2 in jcoupling_list):
        j_test=[0]*len(states)
        for n1,n2 in jcoupling_list:
            if j_test[n1-1]==-1 or j_test[n2-1]==-1:raise ValueError("Spaces coupling j_n's are referenced by smallest n value")
            j_test[max(n1,n2)-1]=-1
    jn=[state.j for state in states];mn=[state.m for state in states]
    coupling_list=[];n_list=[[i+1] for i in range(len(states))]
    for n1,n2 in jcoupling_list:
        j1_n=list(n_list[n1-1]);j2_n=list(n_list[n2-1]);coupling_list.append((j1_n,j2_n));n_list[min(n1,n2)-1]=sorted(j1_n+j2_n)
    if not all(state.j.is_number and state.m.is_number for state in states):
        # The reference returns a held Sum over the coupled j values.
        raise NotImplementedError('coupling symbolic spins (a Sum over j values) is not supported')
    diff_max=[S.Add(*[jn[n-1]-mn[n-1] for n in coupling[0]+coupling[1]]) for coupling in coupling_list]
    result=S.Integer(0)
    for diff in range(int(diff_max[-1])+1):
        n=len(coupling_list)
        for config_num in range(int(S.binomial(diff+n-1,diff))):
            diff_list=_confignum_to_difflist(config_num,diff,n)
            if any(d>mx for d,mx in zip(diff_list,diff_max)):continue
            cg_terms=[];coupled_j=list(jn);jcoupling=[]
            for (j1_n,j2_n),coupling_diff in zip(coupling_list,diff_list):
                j1=coupled_j[min(j1_n)-1];j2=coupled_j[min(j2_n)-1];j3=j1+j2-coupling_diff
                coupled_j[min(j1_n+j2_n)-1]=j3
                m1=S.Add(*[mn[x-1] for x in j1_n]);m2=S.Add(*[mn[x-1] for x in j2_n]);m3=m1+m2
                cg_terms.append((j1,m1,j2,m2,j3,m3));jcoupling.append((min(j1_n),min(j2_n),j3))
            if any(abs(t[5])>t[4] for t in cg_terms):continue
            if any(t[0]+t[2]<t[4] for t in cg_terms):continue
            if any(abs(t[0]-t[2])>t[4] for t in cg_terms):continue
            coefficient=S.Mul(*[_cg(*t) for t in cg_terms])
            if coefficient:result+=coefficient*coupled_evect(j3,m3,jn,jcoupling)
    return result


def _map_atoms(expr,function):
    if isinstance(expr,OperatorExpr):
        result=S.Integer(0)
        for word,c in expr._terms.items():
            term=c
            for atom in word:term=term*function(atom)
            result+=term
        return result
    return expr


def uncouple(expr,jn=None,jcoupling_list=None):
    """Reference uncouple: rewrite every spin state in uncoupled tensor products."""
    def visit(atom):
        if isinstance(atom,_SpinState):return _uncouple(atom,jn,jcoupling_list)
        if isinstance(atom,TensorProductAtom):return TensorProduct(*(uncouple(f,jn,jcoupling_list) for f in atom.factors))
        return atom
    if isinstance(expr,Atom) and not isinstance(expr,TensorProductAtom):return visit(expr)
    return _map_atoms(expr,visit)


def couple(expr,jcoupling_list=None):
    """Reference couple: couple every tensor product of spin states."""
    def visit(atom):
        if not isinstance(atom,TensorProductAtom):return atom
        if not all(isinstance(state,_SpinState) and not isinstance(state,_CoupledSpinState) for state in atom.factors):return atom
        if not all(type(state) is type(atom.factors[0]) for state in atom.factors):raise TypeError('All states must be the same basis')
        return _couple(atom,jcoupling_list)
    if isinstance(expr,TensorProductAtom):return visit(expr)
    return _map_atoms(expr,visit)

__all__=['m_values','Jplus','Jminus','Jx','Jy','Jz','J2','Rotation','WignerD','JxKet','JxBra','JyKet','JyBra','JzKet','JzBra','JzOp','J2Op','JxKetCoupled','JxBraCoupled','JyKetCoupled','JyBraCoupled','JzKetCoupled','JzBraCoupled','couple','uncouple']
