"""Native special-function constructors, metadata and exact helper formulas.

Hyper-representation continuation formulas follow SymPy 1.14 (BSD-3-Clause;
see LICENSE-SYMPY). No SymPy code is imported during native computation.
"""
import functools
import math
import operator
from . import _core


def _s():
    import symox
    return symox


def intlike(n):
    try:
        result=int(n)
    except (TypeError,ValueError,OverflowError,AttributeError):
        return False
    return not bool(n-result)


def multigamma(x,p):
    s=_s();x,p=s.sympify(x),s.sympify(p)
    if p.is_positive is False or p.is_integer is False:
        raise ValueError('Order parameter p must be a positive integer')
    prefactor=s.pi**(p*(p-1)/4)
    if p.is_Integer:
        if int(p)>1024:raise NotImplementedError('explicit multigamma products are limited to dimension 1024')
        return prefactor*s.Mul(*(s.gamma(x+s.Rational(1-k,2)) for k in range(1,int(p)+1)))
    k=s.Dummy('k',integer=True)
    return prefactor*s.Product(s.gamma(x+(1-k)/2),(k,1,p))


class TupleArg(tuple):
    def __new__(cls,*args):
        return super().__new__(cls,(_s().sympify(a) for a in args))
    @property
    def args(self):return tuple(self)
    def _symox_(self):return _s()._tuple(self)
    def as_leading_term(self,*x,logx=None,cdir=0):
        return type(self)(*(v.as_leading_term(*x,logx=logx,cdir=cdir) for v in self))
    def limit(self,x,xlim,dir='+'):
        return type(self)(*(v.limit(x,xlim,dir=dir) for v in self))
    def subs(self,*args,**kwargs):return type(self)(*(v.subs(*args,**kwargs) for v in self))
    def to_sympy(self):
        from sympy.functions.special.hyper import TupleArg
        return TupleArg(*(v.to_sympy() for v in self))


_BESSEL={'besselj','bessely','besseli','besselk','hankel1','hankel2','jn','yn','hn1','hn2'}
_AIRY={'airyai','airybi','airyaiprime','airybiprime'}
_MATHIEU={'mathieuc','mathieus','mathieucprime','mathieusprime'}
_REPS={}


class _HeadMeta(type):
    def __instancecheck__(cls,value):
        return isinstance(value,_core.Expr) and value.class_name in cls._heads
    def __call__(cls,*args,**kwargs):
        if cls.__dict__.get('_abstract',False):
            raise TypeError('%s is an abstract function base'%cls.__name__)
        if kwargs.keys()-{'evaluate'}:raise TypeError('unsupported function option')
        if kwargs.get('evaluate') is False and cls.__name__ in _MATHIEU:
            raise NotImplementedError('held Mathieu constructors with evaluate=False are not implemented')
        if len(args)!=cls._arity:raise TypeError('%s takes %d arguments'%(cls.__name__,cls._arity))
        return _core._apply(cls.__name__,*args)


class BesselBase(metaclass=_HeadMeta):
    _abstract=True;_heads=_BESSEL
    @property
    def order(self):return self.args[0]
    @property
    def argument(self):return self.args[1]
    def fdiff(self,argindex=2):return _fdiff(self,argindex)


class SphericalBesselBase(BesselBase):
    _abstract=True;_heads={'jn','yn','hn1','hn2'}


class SphericalHankelBase(SphericalBesselBase):
    _abstract=True;_heads={'hn1','hn2'}


class AiryBase(metaclass=_HeadMeta):
    _abstract=True;_heads=_AIRY
    def fdiff(self,argindex=1):return _fdiff(self,argindex)


class MathieuBase(metaclass=_HeadMeta):
    _abstract=True;_heads=_MATHIEU;unbranched=True
    def fdiff(self,argindex=1):return _fdiff(self,argindex)


class TupleParametersBase(metaclass=_HeadMeta):
    _abstract=True;_heads={'hyper','meijerg'};is_commutative=True
    def fdiff(self,argindex=3):return _fdiff(self,argindex)


def _make_class(name,base,arity):
    return _HeadMeta(name,(base,),{'_heads':{name},'_arity':arity,'__module__':__name__})


hn1=_make_class('hn1',SphericalHankelBase,2)
hn2=_make_class('hn2',SphericalHankelBase,2)
for _name in sorted(_MATHIEU):globals()[_name]=_make_class(_name,MathieuBase,3)


def assume_integer_order(fn):
    @functools.wraps(fn)
    def guarded(self,nu,z):
        if _s().sympify(nu).is_integer is True:return fn(self,nu,z)
        return None
    return guarded


def jn_zeros(n,k,method='sympy',dps=15):
    """Positive spherical-Bessel zeros by native evaluation and bisection."""
    n,k,dps=operator.index(n),operator.index(k),operator.index(dps)
    if n<0 or k<0:raise ValueError('order and zero count must be nonnegative integers')
    if method not in ('sympy','scipy'):raise NotImplementedError('unknown zero-finding method')
    if not 1<=dps<=15:raise NotImplementedError('native Bessel zero finding currently supports 1 through 15 digits')
    if n>128 or k>10000:raise NotImplementedError('native Bessel zeros are limited to order 128 and 10000 zeros')
    roots=[]
    # j_n has no positive zero below sqrt(n*(n+1)). Starting near
    # the turning point avoids underflow and finds the first root in order.
    left=max(0.5,math.sqrt(n*(n+1)))
    def value(x):
        z=_core._numeric('jn',float(n),x)
        if z is None or not math.isfinite(z.real):raise ArithmeticError('native Bessel evaluation failed')
        return z.real
    fl=value(left)
    for _ in range(k):
        right=left+math.pi/4;fr=value(right)
        for _scan in range(1000):
            if fl*fr<=0:break
            left,fl=right,fr;right+=math.pi/4;fr=value(right)
        else:raise ArithmeticError('could not bracket the next spherical Bessel zero')
        for _bisect in range(80):
            middle=(left+right)/2
            if middle in (left,right):break
            fm=value(middle)
            if fl*fm<=0:right=middle
            else:left,fl=middle,fm
        root=(left+right)/2
        roots.append(_s().Float(format(root,'.15g'),dps))
        left=root+math.pi/8;fl=value(left)
    return roots


def _params(values):
    s=_s()
    if isinstance(values,_core.Expr) and values.class_name=='Tuple':values=values.args
    values=[s.sympify(v) for v in values]
    if any(v.is_finite is False for v in values):raise ValueError('special-function parameters must be finite')
    return tuple(sorted(values,key=s.default_sort_key))


class _MeijerMeta(_HeadMeta):
    def __call__(cls,*args,**kwargs):
        if kwargs.keys()-{'evaluate'}:raise TypeError('unsupported Meijer G option')
        if len(args)==5:args=((args[0],args[1]),(args[2],args[3]),args[4])
        if len(args)!=3:raise TypeError('meijerg expects three or five arguments')
        try:
            aa,bb,z=args
            if isinstance(aa,_core.Expr):aa=aa.args
            if isinstance(bb,_core.Expr):bb=bb.args
            if len(aa)!=2 or len(bb)!=2:raise TypeError('parameter groups require two tuples')
            an,ao=_params(aa[0]),_params(aa[1]);bm,bo=_params(bb[0]),_params(bb[1])
        except (TypeError,AttributeError):raise TypeError('invalid Meijer G parameter groups') from None
        for a in an:
            for b in bm:
                delta=a-b
                if delta.is_Integer and delta>0:raise ValueError('upper and lower poles must not overlap')
        s=_s()
        return _core._apply('meijerg',s._tuple([s._tuple(an),s._tuple(ao)]),s._tuple([s._tuple(bm),s._tuple(bo)]),z)


class meijerg(TupleParametersBase,metaclass=_MeijerMeta):
    _heads={'meijerg'}


def appellf1(a,b1,b2,c,x,y):
    s=_s();a,b1,b2,c,x,y=map(s.sympify,(a,b1,b2,c,x,y))
    if (s.default_sort_key(b1),s.default_sort_key(x))>(s.default_sort_key(b2),s.default_sort_key(y)):
        b1,b2,x,y=b2,b1,y,x
    return _core._apply('appellf1',a,b1,b2,c,x,y)


def marcumq(m,a,b):return _core._apply('marcumq',m,a,b)


class HyperRep(metaclass=_HeadMeta):
    _abstract=True;_heads=set()
    @classmethod
    def _expr_small(cls,*args):return _rep_formula(cls.__name__,*args)
    @classmethod
    def _expr_small_minus(cls,*args):return _rep_formula(cls.__name__,*args,minus=True)
    @classmethod
    def _expr_big(cls,*args):return _rep_formula(cls.__name__,*args,big=True)
    @classmethod
    def _expr_big_minus(cls,*args):return _rep_formula(cls.__name__,*args,big=True,minus=True)


for _suffix in ('power1','power2','log1','log2','atanh','asin1','asin2','sqrts1','sqrts2','cosasin','sinasin'):
    _name='HyperRep_'+_suffix
    _REPS[_name]=_make_class(_name,HyperRep,2 if _suffix in ('power1','power2','sqrts1','sqrts2','cosasin','sinasin') else 1)
    globals()[_name]=_REPS[_name]
HyperRep._heads=set(_REPS)


def _rep_formula(name,*args,big=False,minus=False):
    s=_s();args=tuple(map(s.sympify,args));n=args[-1] if big else s.Integer(0)
    values=args[:-1] if big else args
    z=values[-1];a=values[0] if len(values)==2 else None
    f=name.removeprefix('HyperRep_')
    if big and not n.is_Integer:
        raise NotImplementedError('continuation formulas require a concrete integer winding number')
    odd=bool(int(n)%2) if big else False
    half=s.Rational(1,2);ii=s.I;pi=s.pi
    root=s.sqrt(z)
    if f=='power1':
        if not big or a.is_integer is True:return (1+z if minus else 1-z)**a
        return (1+z)**a*s.exp(2*n*pi*ii*a) if minus else (z-1)**a*s.exp((2*n-1)*pi*ii*a)
    if f=='power2':
        if not big:return 2**(2*a-1)*(1+s.sqrt(1+z if minus else 1-z))**(1-2*a)
        if minus:
            sign=-1 if odd else 1
            return sign*2**(2*a-1)*(s.sqrt(1+z)+sign)**(1-2*a)*s.exp(-2*pi*ii*a*n)
        sign=1 if odd else -1;branch=n-1 if odd else n
        return 2**(2*a-1)*(1+sign*ii*s.sqrt(z-1))**(1-2*a)*s.exp(-2*branch*pi*ii*a)
    if f=='log1':
        if not big:return s.log(1+z if minus else 1-z)
        return s.log(1+z)+2*n*pi*ii if minus else s.log(z-1)+(2*n-1)*pi*ii
    if f in ('atanh','asin1','asin2'):
        if not big:
            if z==0:return s.Integer(1)
            value=(s.atan(root) if f=='atanh' else s.asinh(root))/root if minus else (s.atanh(root) if f=='atanh' else s.asin(root))/root
        elif f=='atanh':
            value=(s.atan(root)-(pi if odd else 0))/root if minus else (s.acoth(root)+(-1 if odd else 1)*ii*pi/2)/root
        else:
            value=(-1)**n*(s.asinh(root)+n*pi*ii)/root if minus else (-1)**n*((half-n)*pi+ii*s.acosh(root))/root
        if f=='asin2':
            value/= _rep_formula('HyperRep_power1',half,z,*((n,) if big else ()),big=big,minus=minus)
        return value
    if f in ('sqrts1','sqrts2'):
        if minus:
            angle=2*a*s.atan(root)-(2*pi*a if big and odd else 0)
            factor=(1+z)**a*(s.exp(2*pi*ii*n*a) if big else 1)
            return factor*s.cos(angle) if f=='sqrts1' else factor*root*s.sin(angle)
        if not big:
            left=(1-root)**(2*a);right=(1+root)**(2*a)
        elif odd:
            left=(root-1)**(2*a)*s.exp(2*pi*ii*a*n)
            right=(root+1)**(2*a)*s.exp(2*pi*ii*a*(n-1))
        else:
            left=(root-1)**(2*a)*s.exp(2*pi*ii*a*(n-1))
            right=(root+1)**(2*a)*s.exp(2*pi*ii*a*n)
        return (left+right)/2 if f=='sqrts1' else root*(left-right)/2
    if f=='log2':
        if not big:return s.log((1+s.sqrt(1+z if minus else 1-z))/2)
        if minus:return pi*ii*n+s.log((s.sqrt(1+z)+(-1 if odd else 1))/2)
        return (n-half)*pi*ii+s.log(root/2)+(-1 if odd else 1)*ii*s.asin(1/root)
    if f in ('cosasin','sinasin'):
        if not big:
            if minus:
                angle=2*a*s.asinh(root)
                return s.cosh(angle) if f=='cosasin' else -root/s.sqrt(1+z)*s.sinh(angle)
            angle=2*a*s.asin(root)
            return s.cos(angle) if f=='cosasin' else root/s.sqrt(1-z)*s.sin(angle)
        angle=2*a*s.asinh(root)+2*a*pi*ii*n if minus else 2*a*s.acosh(root)+a*pi*ii*(2*n-1)
        return s.cosh(angle) if f=='cosasin' else -s.sinh(angle)/s.sqrt(1+1/z if minus else 1-1/z)
    raise ValueError('unknown hypergeometric representative')


def _rep_rewrite(expr,rule):
    cls=_REPS[expr.class_name]
    x,n=expr.args[-1].extract_branch_factor(allow_half=True)
    minus=not n.is_Integer
    if minus:n-= _s().Rational(1,2)
    values=expr.args[:-1]+(x,)
    small=_rep_formula(expr.class_name,*values,minus=minus)
    if rule=='nonrepsmall':return small
    big=_rep_formula(expr.class_name,*values,n,big=True,minus=minus)
    return small if big==small else _s().Piecewise((big,_s().Abs(x)>1),(small,True))


def _fdiff(expr,argindex=None):
    name=expr.class_name
    if argindex is None:
        argindex=2 if name in _BESSEL or name=='marcumq' else 3 if name in ('hyper','meijerg') else 5 if name=='appellf1' else 1
    index=operator.index(argindex)
    if index<1:raise ValueError('argument index must be positive')
    value=_core._special_partial(expr,index-1)
    if value is None:raise NotImplementedError('partial derivative of %s with respect to argument %d is not implemented'%(name,index))
    return value


def _property(expr,name):
    head=expr.class_name
    if name=='order' and head in _BESSEL:return expr.args[0]
    if name=='argument':
        if head in _BESSEL:return expr.args[1]
        if head in ('hyper','meijerg') or head in _MATHIEU:return expr.args[2]
        if head in _AIRY:return expr.args[0]
    if head=='meijerg':
        an,ao=expr.args[0].args;bm,bo=expr.args[1].args
        groups={'an':an.args,'aother':ao.args,'bm':bm.args,'bother':bo.args,
                'ap':an.args+ao.args,'bq':bm.args+bo.args}
        if name in groups:return TupleArg(*groups[name])
        if name=='nu':return sum(groups['bq'],_s().Integer(0))-sum(groups['ap'],_s().Integer(0))
        if name=='delta':return len(groups['bm'])+len(groups['an'])-_s().Rational(len(groups['ap'])+len(groups['bq']),2)
    if head=='hyper':
        if name in ('ap','bq'):return TupleArg(*expr.args[0 if name=='ap' else 1].args)
        if name=='eta':return sum(expr.ap,_s().Integer(0))-sum(expr.bq,_s().Integer(0))
    raise AttributeError('%s has no %s property'%(head,name))


def _integrand(expr,t):
    if expr.class_name!='meijerg':raise AttributeError('integrand is available for Meijer G functions')
    s=_s();t=s.sympify(t)
    return expr.argument**t*s.Mul(*(s.gamma(b-t) for b in expr.bm))*s.Mul(*(s.gamma(1-a+t) for a in expr.an))/s.Mul(*(s.gamma(1-b+t) for b in expr.bother))/s.Mul(*(s.gamma(a-t) for a in expr.aother))


def _period(expr):
    if expr.class_name!='meijerg':raise AttributeError('get_period is available for Meijer G functions')
    def period(values):
        denominators=[]
        for i,value in enumerate(values):
            if not value.is_Rational or any((value-other).is_Integer for other in values[i+1:]):return None
            denominators.append(int(value.q))
        return math.lcm(*denominators) if denominators else 1
    a,b=period(expr.an),period(expr.bm)
    p,q=len(expr.ap),len(expr.bq)
    selected=[a,b] if p==q else [b] if p<q else [a]
    return _s().oo if None in selected else 2*_s().pi*math.lcm(*selected)


def install():
    s=_s();Expr=_core.Expr
    for name in ('order','argument','an','aother','bm','bother','ap','bq','nu','delta','eta'):
        previous=getattr(Expr,name,None)
        def getter(expr,_name=name,_previous=previous):
            try:return _property(expr,_name)
            except AttributeError:
                if _previous is not None:return _previous.__get__(expr,type(expr))
                raise
        setattr(Expr,name,property(getter))
    Expr.fdiff=_fdiff;Expr.integrand=_integrand;Expr.get_period=_period
    old_rewrite=Expr.rewrite
    def rewrite(expr,*args,deep=True,**kwargs):
        if not args:
            return old_rewrite(expr,*args,deep=deep,**kwargs)
        rule=getattr(args[-1],'__name__',str(args[-1]));patterns=args[:-1]
        if rule not in ('nonrep','nonrepsmall','Integral','Sum','besseli'):
            return old_rewrite(expr,*args,deep=deep,**kwargs)
        if patterns and isinstance(patterns[0],(list,tuple,set,frozenset)):patterns=patterns[0]
        names={getattr(p,'__name__',str(p)) for p in patterns}
        if rule not in ('nonrep','nonrepsmall'):
            expr=old_rewrite(expr,*args,deep=deep,**kwargs)
        def visit(e):
            values=tuple(visit(a) for a in e.args) if deep else e.args
            if values!=e.args:e=e.func(*values)
            if rule in ('nonrep','nonrepsmall') and e.class_name in _REPS and (not names or e.class_name in names or 'HyperRep' in names):
                return _rep_rewrite(e,rule)
            if e.class_name=='marcumq' and (not names or 'marcumq' in names):
                m,a,b=e.args
                if rule=='Integral':
                    t=kwargs.get('x',s.Dummy('x'))
                    return a**(1-m)*s.Integral(t**m*s.exp(-(t*t+a*a)/2)*s.besseli(m-1,a*t),(t,b,s.oo))
                if rule=='Sum':
                    k=kwargs.get('k',s.Dummy('k',integer=True))
                    return s.exp(-(a*a+b*b)/2)*s.Sum((a/b)**k*s.besseli(k,a*b),(k,1-m,s.oo))
                if rule=='besseli' and a==b and m.is_Integer and m>=1:
                    if int(m)>1024:raise NotImplementedError('explicit Marcum Q Bessel sums are limited to order 1024')
                    return s.Rational(1,2)+s.exp(-a*a)*(s.besseli(0,a*a)/2+sum((s.besseli(i,a*a) for i in range(1,int(m))),s.Integer(0)))
            return e
        return visit(expr)
    Expr.rewrite=rewrite
