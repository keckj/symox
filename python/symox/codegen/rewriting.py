"""Cost-directed expression optimizations for stable numerical code generation."""
import symox as S
from .cfunctions import exp2,expm1,log2,log1p
from .._kind import matrix_tag


class Optimization:
    def __init__(self,cost_function=None,priority=1):self.cost_function,self.priority=cost_function,priority
    def __call__(self,expr):raise NotImplementedError('an optimization must implement its transform')
    def cheapest(self,*args):return min(args,key=self.cost_function)
class ReplaceOptim(Optimization):
    def __init__(self,query,value,**kwargs):super().__init__(**kwargs);self.query,self.value=query,value
    def __call__(self,expr):
        return expr.replace(self.query,self.value)
class _Transform(Optimization):
    def __init__(self,transform,**kwargs):super().__init__(**kwargs);self.transform=transform
    def __call__(self,expr):return self.transform(expr)

def optimize(expr,optimizations):
    # sympy.codegen.rewriting.optimize: highest priority first, a cost
    # function keeps the cheaper of the old and the rewritten expression.
    for optimization in sorted(optimizations,key=lambda op:op.priority,reverse=True):
        result=optimization(expr)
        expr=result if optimization.cost_function is None else optimization.cheapest(expr,result)
    return expr

def _head(expr,name):return getattr(expr,'class_name','')==name

exp2_opt=ReplaceOptim(lambda p:p.is_Pow and p.base==2,lambda p:exp2(p.exp))

_u=S.Wild('u',properties=[lambda x:not x.is_number and not x.is_Add])
_v=S.Wild('v');_w=S.Wild('w');_n=S.Wild('n',properties=[lambda x:x.is_number])

sinc_opt1=ReplaceOptim(S.sin(_w)/_w,S.sinc(_w))
sinc_opt2=ReplaceOptim(S.sin(_n*_w)/_w,_n*S.sinc(_n*_w))
sinc_opts=(sinc_opt1,sinc_opt2)

log2_opt=ReplaceOptim(_v*S.log(_w)/S.log(2),_v*log2(_w),cost_function=lambda expr:expr.count(
    lambda e:(e.is_Pow and e.exp.is_negative) or ((_head(e,'log') or _head(e,'log2')) and not e.args[0].is_number)))
log2const_opt=ReplaceOptim(S.log(2)*log2(_w),S.log(_w))

logsumexp_2terms_opt=ReplaceOptim(
    lambda l:_head(l,'log') and l.args[0].is_Add and len(l.args[0].args)==2 and all(_head(t,'exp') for t in l.args[0].args),
    lambda l:S.Max(*[e.args[0] for e in l.args[0].args])+log1p(S.exp(S.Min(*[e.args[0] for e in l.args[0].args]))))


def _sift(seq,keyfunc):
    true,false=[],[]
    for item in seq:(true if keyfunc(item) else false).append(item)
    return true,false


class FuncMinusOneOptim(ReplaceOptim):
    """Rewrite ``c*f(x) + n`` with ``c*f_m_1(x)`` when it cancels the number term."""
    def __init__(self,func,func_m_1,opportunistic=True):
        weight=10
        super().__init__(lambda e:e.is_Add,self.replace_in_Add,
                         cost_function=lambda expr:expr.count_ops()-weight*expr.count(func_m_1))
        self.func,self.func_m_1,self.opportunistic=func,func_m_1,opportunistic
    def _group_Add_terms(self,add):
        numbers,non_num=_sift(add.args,lambda arg:arg.is_number)
        numsum=sum(numbers)
        terms_with_func,other=_sift(non_num,lambda arg:arg.has(self.func))
        return numsum,terms_with_func,other
    def replace_in_Add(self,e):
        numsum,terms_with_func,other_non_num_terms=self._group_Add_terms(e)
        if numsum==0:return e
        substituted,untouched=[],[]
        for with_func in terms_with_func:
            if with_func.is_Mul:
                func,coeff=_sift(with_func.args,lambda arg:arg.func==self.func)
                if len(func)==1 and len(coeff)==1:func,coeff=func[0],coeff[0]
                else:coeff=None
            elif with_func.func==self.func:func,coeff=with_func,S.Integer(1)
            else:coeff=None
            if coeff is not None and coeff.is_number and S.sign(coeff)==-S.sign(numsum):
                do_substitute=abs(coeff+numsum)<abs(numsum) if self.opportunistic else coeff+numsum==0
                if do_substitute:
                    numsum+=coeff;substituted.append(coeff*self.func_m_1(*func.args));continue
            untouched.append(with_func)
        return e.func(numsum,*substituted,*untouched,*other_non_num_terms)
    def __call__(self,expr):
        alt1=super().__call__(expr);alt2=super().__call__(expr.factor())
        return self.cheapest(alt1,alt2)


expm1_opt=FuncMinusOneOptim(S.exp,expm1)
cosm1_opt=FuncMinusOneOptim(S.cos,lambda x:-2*S.sin(x/2)**2)
log1p_opt=ReplaceOptim(lambda e:_head(e,'log'),
    lambda l:S.expand_log(l.replace(S.log,lambda arg:S.log(arg.factor()))).replace(S.log(_u+1),log1p(_u)))

def create_expand_pow_optimization(limit,*,base_req=lambda b:b.is_symbol):
    return ReplaceOptim(
        lambda e:e.is_Pow and base_req(e.base) and e.exp.is_Integer and abs(e.exp)<=limit,
        lambda p:(S.UnevaluatedExpr(S.Mul(*([p.base]*+p.exp),evaluate=False)) if p.exp>0 else
                  1/S.UnevaluatedExpr(S.Mul(*([p.base]*-p.exp),evaluate=False))))

def _logadd(expr,base2=False):
    def visit(e):
        if base2:
            numerator,denominator=S.fraction(e)
            if denominator!=S.log(2) or not _head(numerator,'log') or len(numerator.args)!=1:return e
            summand=numerator.args[0]
            if not summand.is_Add or len(summand.args)!=2:return e
            a,b=summand.args
            if not (a.is_Pow and a.base==2 and b.is_Pow and b.base==2):return e
            x,y=a.exp,b.exp
            if x.is_real is not True or y.is_real is not True:return e
            return S.Max(x,y)+log1p(2**(-S.Abs(x-y)))/S.log(2)
        if not _head(e,'log') or len(e.args)!=1 or not e.args[0].is_Add or len(e.args[0].args)!=2:return e
        a,b=e.args[0].args
        if not (_head(a,'exp') and _head(b,'exp')):return e
        x,y=a.args[0],b.args[0]
        if x.is_real is False or y.is_real is False:return e
        return S.Max(x,y)+log1p(S.exp(-S.Abs(x-y)))
    return S.bottom_up(expr,visit)
logaddexp_opt=_Transform(_logadd);logaddexp2_opt=_Transform(lambda e:_logadd(e,True))

def _powm1(e):
    def visit(node):
        if not node.is_Add:return node
        if S.Integer(-1) not in node.args:return node
        power=next((a for a in node.args if a.is_Pow and a.base.is_positive is True),None)
        return node if power is None else node-power+1+expm1(power.exp*S.log(power.base))
    return S.bottom_up(e,visit)
powm1_opt=_Transform(_powm1)

def _matinv(expr):
    from .matrix_nodes import MatrixSolve
    if matrix_tag(expr)=='MatMul':
        args=expr.args
        # The reference rewrites only a full-rank MatrixSymbol times a vector.
        if len(args)==2 and matrix_tag(args[0])=='Inverse' and args[1].shape[1]==1:
            inner=args[0].args[0]
            if isinstance(inner,S.MatrixSymbol) and bool(S.ask(S.Q.fullrank(inner))):return MatrixSolve(inner,args[1])
    return expr
matinv_opt=_Transform(_matinv)
optims_c99=(expm1_opt,log1p_opt,exp2_opt,log2_opt,log2const_opt)
optims_numpy=optims_c99+(logaddexp_opt,logaddexp2_opt)+sinc_opts
optims_scipy=(cosm1_opt,powm1_opt)
__all__=['Optimization','ReplaceOptim','FuncMinusOneOptim','optimize','create_expand_pow_optimization',
 'expm1_opt','cosm1_opt','exp2_opt','log1p_opt','log2_opt','log2const_opt','logaddexp_opt',
 'logaddexp2_opt','logsumexp_2terms_opt','powm1_opt','sinc_opt1','sinc_opt2','sinc_opts','matinv_opt',
 'optims_c99','optims_numpy','optims_scipy']
