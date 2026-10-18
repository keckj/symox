"""Cost-directed expression optimizations for stable numerical code generation."""
import symox as S
from .cfunctions import exp2,expm1,log2,log1p


class Optimization:
    def __init__(self,cost_function=None,priority=1):self.cost_function,self.priority=cost_function,priority
    def __call__(self,expr):raise NotImplementedError('an optimization must implement its transform')
class ReplaceOptim(Optimization):
    def __init__(self,query,value,**kwargs):super().__init__(**kwargs);self.query,self.value=query,value
    def __call__(self,expr):
        if callable(self.query) and not isinstance(self.query,type) and not hasattr(self.query,'_symox_head_name'):
            return S.bottom_up(expr,lambda node:self.value(*node.args) if self.query(node) else node)
        return expr.replace(self.query,self.value)
class _Transform(Optimization):
    def __init__(self,transform,**kwargs):super().__init__(**kwargs);self.transform=transform
    def __call__(self,expr):return self.transform(expr)

def optimize(expr,optimizations):
    for optimization in sorted(optimizations,key=lambda op:op.priority,reverse=True):
        result=optimization(expr);cost=optimization.cost_function
        if cost is None or cost(result)<cost(expr):expr=result
    return expr

def _head(expr,name):return getattr(expr,'class_name','')==name

def _minus_one(expr,source,target,opportunistic):
    def visit(node):
        if not node.is_Add:return node
        terms=list(node.args)
        for i,term in enumerate(terms):
            coefficient,base=term.as_coeff_Mul()
            if not _head(base,source):continue
            constant=sum((t for t in terms if t.is_number),S.Integer(0))
            if constant==-coefficient or opportunistic and constant.is_negative:
                rest=sum((t for j,t in enumerate(terms) if j!=i),S.Integer(0))
                return rest+coefficient+coefficient*target(*base.args)
        return node
    return S.bottom_up(expr,visit)
class FuncMinusOneOptim(Optimization):
    def __init__(self,func,func_m_1,opportunistic=True):
        super().__init__();self.func,self.func_m_1,self.opportunistic=func,func_m_1,opportunistic
    def __call__(self,expr):return _minus_one(expr,getattr(self.func,'__name__',str(self.func)),self.func_m_1,self.opportunistic)

def create_expand_pow_optimization(limit,*,base_req=lambda base:True):
    def transform(expr):
        def expand(node):
            if node.is_Pow and node.exp.is_Integer and 1<abs(int(node.exp))<=limit and base_req(node.base):
                power=int(node.exp);base=node.base if power>0 else S.Pow(node.base,-1,evaluate=False)
                return S.Mul(*([base]*abs(power)),evaluate=False)
            return node
        return S.bottom_up(expr,expand)
    return _Transform(transform)

expm1_opt=FuncMinusOneOptim(S.exp,expm1)
cosm1_opt=FuncMinusOneOptim(S.cos,lambda x:-2*S.sin(x/2)**2)
exp2_opt=ReplaceOptim(lambda e:e.is_Pow and e.base==2,lambda base,exp:exp2(exp))
log1p_opt=ReplaceOptim(lambda e:_head(e,'log') and len(e.args)==1 and e.args[0].is_Add and S.Integer(1) in e.args[0].args,lambda arg:log1p(arg-1))

def _log2(expr):
    def visit(e):
        numerator,denominator=S.fraction(e)
        if _head(numerator,'log') and denominator==S.log(2):return log2(numerator.args[0])
        return e
    return S.bottom_up(expr,visit)
log2_opt=_Transform(_log2);log2const_opt=_Transform(_log2)

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
logaddexp_opt=_Transform(_logadd);logaddexp2_opt=_Transform(lambda e:_logadd(e,True));logsumexp_2terms_opt=_Transform(_logadd)

def _powm1(e):
    def visit(node):
        if not node.is_Add:return node
        if S.Integer(-1) not in node.args:return node
        power=next((a for a in node.args if a.is_Pow and a.base.is_positive is True),None)
        return node if power is None else node-power+1+expm1(power.exp*S.log(power.base))
    return S.bottom_up(e,visit)
powm1_opt=_Transform(_powm1)

def _sinc(e):
    def visit(node):
        numerator,denominator=S.fraction(node)
        if _head(numerator,'sin'):
            arg=numerator.args[0];ratio=S.cancel(arg/denominator)
            if ratio.is_number:return ratio*S.sinc(arg)
        return node
    return S.bottom_up(e,visit)
sinc_opt1=_Transform(_sinc);sinc_opt2=_Transform(_sinc)

def _matinv(expr):
    from .matrix_nodes import MatrixSolve
    if getattr(expr,'kind',None)=='MatMul':
        args=expr.args
        if len(args)==2 and getattr(args[0],'kind',None)=='Inverse':return MatrixSolve(args[0].args[0],args[1])
    return expr
matinv_opt=_Transform(_matinv)
__all__=['Optimization','ReplaceOptim','FuncMinusOneOptim','optimize','create_expand_pow_optimization',
 'expm1_opt','cosm1_opt','exp2_opt','log1p_opt','log2_opt','log2const_opt','logaddexp_opt',
 'logaddexp2_opt','logsumexp_2terms_opt','powm1_opt','sinc_opt1','sinc_opt2','matinv_opt']
