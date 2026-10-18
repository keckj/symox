"""Precision-aware numeric protocol over the native evaluator and mpmath.

Internal results carry (real_mpf, imaginary_mpf, real_accuracy, imag_accuracy).
Precision arguments in this module are bits; public Expr.evalf uses decimals.
"""
import math
import mpmath
from mpmath.libmp import prec_to_dps,mpf_pos,mpf_sum,mpf_abs,fzero,fone,finf,fninf,fnan
import symox as S
from .. import PrecisionExhausted,N
INF=float('inf');MINUS_INF=-INF


def _context(prec):
    context=mpmath.mp.clone();context.prec=int(prec);return context

def bitcount(n):return abs(int(n)).bit_length()
def fastlog(x):return MINUS_INF if not x or iszero(x) else x[2]+x[3]
def iszero(value,scaled=False):
    if scaled:return bool(value and isinstance(value[0],list) and value[1]==value[-1]==1)
    return not value or value[1]==0 and value[-1]==0

def scaled_zero(mag,sign=1):
    if isinstance(mag,tuple) and iszero(mag,scaled=True):return (mag[0][0],*mag[1:])
    if not isinstance(mag,int) or sign not in (-1,1):raise ValueError('scaled zero requires an integer magnitude and sign +/-1')
    return ([int(sign<0)],1,mag,1),-1

def pure_complex(expr,or_real=False):
    expr=S.sympify(expr);a,b=expr.as_real_imag()
    if (a.is_Number or a==0) and (b.is_Number or b==0) and (b!=0 or or_real):return a,b
    return None

def quad_to_mpmath(q,ctx=None):
    if q==S.zoo:raise NotImplementedError('complex infinity has no mpmath scalar representation')
    context=ctx or mpmath.mp;real,imag,_,_=q
    return context.make_mpc((real or fzero,imag)) if imag else context.make_mpf(real or fzero)

def finalize_complex(real,imag,prec):
    real=None if real==fzero else real;imag=None if imag==fzero else imag
    if real is None and imag is None:raise ValueError('complex zero has unknown relative accuracy')
    if real is None:return None,imag,None,prec
    if imag is None:return real,None,prec,None
    delta=fastlog(real)-fastlog(imag)
    return real,imag,prec+min(delta,0),prec+min(-delta,0)

def complex_accuracy(result):
    if result==S.zoo:return INF
    real,imag,ra,ia=result
    if not real:return ia if imag else INF
    if not imag:return ra
    return max(fastlog(real),fastlog(imag))-max(fastlog(real)-ra,fastlog(imag)-ia)

def chop_parts(value,prec):
    if value==S.zoo:return value
    real,imag,ra,ia=value
    if real and real not in (finf,fninf,fnan) and fastlog(real)<-prec+4:real=ra=None
    if imag and imag not in (finf,fninf,fnan) and fastlog(imag)<-prec+4:imag=ia=None
    return real,imag,ra,ia

def check_target(expr,result,prec):
    if complex_accuracy(result)<prec:raise PrecisionExhausted('requested accuracy was not reached for '+str(expr))


def _numeric(expr,prec,options):
    context=_context(prec+20);digits=prec_to_dps(prec+20)+8
    expr=S.sympify(expr)
    if options.get('subs'):expr=expr.subs(options['subs'],simultaneous=True)
    if expr==S.zoo:raise NotImplementedError('complex infinity is not a finite numeric result')
    value=expr.evalf(digits)
    if getattr(value,'is_Integral',False):return quad_to_mpmath(do_integral(value,prec,options),context)
    if getattr(value,'is_Sum',False) or getattr(value,'is_Product',False):
        return _sum_or_product(value,prec,options,product_=getattr(value,'is_Product',False))
    real,imag=value.as_real_imag()
    if real.free_symbols or imag.free_symbols:raise NotImplementedError('numerical evaluation requires values for all free symbols')
    try:
        a=context.mpf(str(real.evalf(digits)));b=context.mpf(str(imag.evalf(digits)))
    except (TypeError,ValueError):raise NotImplementedError('native numerical evaluation did not resolve '+str(expr)) from None
    return context.mpc(a,b) if b else a


def _quad(value,prec):
    context=_context(prec);real,imag=context.re(value),context.im(value)
    a=None if not real else mpf_pos(real._mpf_,prec)
    b=None if not imag else mpf_pos(imag._mpf_,prec)
    return a,b,prec if a is not None else None,prec if b is not None else None

def evalf(expr,prec,options):
    value=S.sympify(expr)
    if value==S.zoo:return value
    result=_quad(_numeric(value,prec,dict(options)),prec)
    return chop_parts(result,prec) if options.get('chop') else result

def as_mpmath(expr,prec,options):return quad_to_mpmath(evalf(expr,prec,options),_context(prec))
def evalf_subs(prec,subs):
    return {symbol:value.evalf(prec_to_dps(prec)+5) if hasattr(value,'evalf') else value for symbol,value in subs.items()}
def get_abs(expr,prec,options):return evalf(S.Abs(expr),prec,options)
def get_complex_part(expr,no,prec,options):return evalf(S.re(expr) if no==0 else S.im(expr),prec,options)
def get_integer_part(expr,no,options,return_ints=False):
    # no=-1 is floor, no=1 is ceiling; both components are rounded separately.
    context=_context(options.get('maxprec',333));value=_numeric(expr,context.prec,options)
    function=context.floor if no<0 else context.ceil
    a,b=int(function(context.re(value))),int(function(context.im(value)))
    return (a,b) if return_ints else _quad(context.mpc(a,b),context.prec)

def add_terms(terms,prec,target_prec):
    terms=[(value,accuracy) for value,accuracy in terms if not iszero(value)]
    if not terms:return None,None
    if len(terms)==1:return terms[0]
    result=mpf_sum([scaled_zero(v) if iszero(v,scaled=True) else v for v,a in terms],prec)
    error=max(fastlog(v)-a for v,a in terms)+math.ceil(math.log2(len(terms)))
    if iszero(result):return scaled_zero(int(error))
    return mpf_pos(result,target_prec),min(target_prec,int(fastlog(result)-error))


def do_integral(expr,prec,options):
    context=_context(prec+20)
    if len(expr.limits)!=1 or len(expr.limits[0])!=3:raise NotImplementedError('numeric integral requires one definite interval')
    variable,a,b=expr.limits[0];lower,upper=as_mpmath(a,prec+20,options),as_mpmath(b,prec+20,options)
    from .numbers import sympify_mpmath
    def integrand(value):return _numeric(expr.function.subs(variable,sympify_mpmath(value)),prec+20,options)
    result=context.quad(integrand,[lower,upper])
    return _quad(result,prec)


def _sum_or_product(expr,prec,options,product_=False):
    context=_context(prec+20)
    if len(expr.limits)!=1 or len(expr.limits[0])!=3:raise NotImplementedError('numeric sums require one definite index range')
    variable,a,b=expr.limits[0];a=as_mpmath(a,prec+20,options);b=as_mpmath(b,prec+20,options)
    from .numbers import sympify_mpmath
    function=lambda index:_numeric(expr.function.subs(variable,sympify_mpmath(index)),prec+20,options)
    if context.isfinite(a) and context.isfinite(b):
        if b-a>1000000:raise ValueError('explicit numeric summation exceeds one million terms')
        values=(function(context.mpf(i)) for i in range(int(a),int(b)+1))
        return context.fprod(values) if product_ else context.fsum(values)
    return context.nprod(function,[a,b]) if product_ else context.nsum(function,[a,b])


def hypsum(expr,n,start,prec):return _sum_or_product(S.Sum(expr,(n,start,S.oo)),prec,{})
def check_convergence(numer,denom,n):
    numerator,denominator=S.Poly(numer,n),S.Poly(denom,n)
    rate=denominator.degree()-numerator.degree()
    if rate:return rate,None,None
    constant=denominator.LC()/numerator.LC()
    if S.Abs(constant)!=1:return 0,constant,None
    if numerator.degree()==0:return 0,constant,0
    return 0,constant,(denominator.all_coeffs()[1]-numerator.all_coeffs()[1])/denominator.LC()


def _evaluation_handler(name):
    def handler(expr,prec,options):return evalf(expr,prec,options)
    handler.__name__='evalf_'+name
    return handler
for _name in 'abs add alg_num atan ceiling exp float floor im integer integral log mul piecewise pow prod rational re sum symbol trig'.split():
    globals()['evalf_'+_name]=_evaluation_handler(_name)


class EvalfMixin:
    def evalf(self,n=15,subs=None,maxn=100,chop=False,strict=False,quad=None,verbose=False):
        evaluator=getattr(self,'_eval_evalf',None)
        if evaluator is not None:
            result=evaluator(max(1,math.ceil(n*math.log2(10))))
            if result is not None:return result
        from symox import sympify
        expression=sympify(self)
        if expression is self:raise NotImplementedError('a custom numerical object must implement _eval_evalf or a scalar conversion')
        return expression.evalf(n,subs=subs,maxn=maxn,chop=chop,strict=strict)
    n=evalf

__all__=['PrecisionExhausted','N','evalf','EvalfMixin','add_terms','as_mpmath','bitcount','check_convergence',
 'check_target','chop_parts','complex_accuracy','do_integral','evalf_subs','fastlog','finalize_complex','get_abs',
 'get_complex_part','get_integer_part','hypsum','iszero','pure_complex','quad_to_mpmath','scaled_zero']
__all__ += ['evalf_'+name for name in 'abs add alg_num atan ceiling exp float floor im integer integral log mul piecewise pow prod rational re sum symbol trig'.split()]
