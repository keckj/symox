"""Bounded polynomial and small-term approximations for numerical expressions."""
from itertools import product
import symox as S
from .rewriting import Optimization


def _interval(expr,bounds):
    if expr in bounds:return tuple(map(float,bounds[expr]))
    if expr.is_number:
        try:value=float(expr.evalf());return value,value
        except (ValueError,TypeError):return None
    intervals=[_interval(a,bounds) for a in expr.args]
    if any(v is None for v in intervals):return None
    if expr.is_Add:return sum(a for a,b in intervals),sum(b for a,b in intervals)
    if expr.is_Mul:
        low=high=1.
        for a,b in intervals:values=[low*a,low*b,high*a,high*b];low,high=min(values),max(values)
        return low,high
    if expr.is_Pow and expr.exp.is_Integer:
        a,b=intervals[0];power=int(expr.exp)
        if power<0 and a<=0<=b:return None
        values=[a**power,b**power]
        if power>0 and power%2==0 and a<=0<=b:values.append(0.)
        return min(values),max(values)
    import math
    if expr.class_name in ('sin','cos'):return -1.,1.
    if expr.class_name=='exp':return math.exp(intervals[0][0]),math.exp(intervals[0][1])
    return None


class SumApprox(Optimization):
    def __init__(self,bounds,reltol,**kwargs):super().__init__(**kwargs);self.bounds,self.reltol=dict(bounds),float(reltol)
    def __call__(self,expr):
        def simplify(node):
            if not node.is_Add:return node
            interval=_interval(node,self.bounds)
            if interval is None or interval[0]<=0<=interval[1]:return node
            minimum=min(abs(v) for v in interval);budget=self.reltol*minimum;terms=[]
            for term in node.args:
                bound=_interval(term,self.bounds)
                magnitude=max(map(abs,bound)) if bound is not None else float('inf')
                if magnitude<=budget:budget-=magnitude
                else:terms.append(term)
            return S.Add(*terms)
        return S.bottom_up(expr,simplify)


class SeriesApprox(Optimization):
    def __init__(self,bounds,reltol,max_order=4,n_point_checks=4,**kwargs):
        super().__init__(**kwargs);self.bounds=dict(bounds);self.reltol=float(reltol);self.max_order=int(max_order);self.n_point_checks=int(n_point_checks)
    def __call__(self,expr):
        variables=sorted(expr.free_symbols,key=str)
        if len(variables)!=1 or variables[0] not in self.bounds:return expr
        x=variables[0];a,b=map(float,self.bounds[x]);count=max(2,self.n_point_checks)
        points=[a+(b-a)*i/(count-1) for i in range(count)]
        for order in range(1,self.max_order+1):
            try:
                polynomial=S.series(expr,x,0,order).removeO()
                if not polynomial.is_polynomial(x):continue
                error=0.
                for point in points:
                    exact=complex(expr.subs(x,point).evalf());approx=complex(polynomial.subs(x,point).evalf())
                    error=max(error,abs(approx-exact)/max(abs(exact),1e-300))
                if error<=self.reltol:return polynomial
            except (ValueError,TypeError,NotImplementedError):continue
        return expr
__all__=['SumApprox','SeriesApprox']
