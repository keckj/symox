"""Bounded polynomial and small-term approximations for numerical expressions."""
import math
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
    """The reference SumApprox: drop Add terms bounded far below the others."""
    def __init__(self,bounds,reltol,**kwargs):super().__init__(**kwargs);self.bounds,self.reltol=bounds,reltol
    def __call__(self,expr):
        return expr.factor().replace(self.query,lambda arg:self.value(arg))
    def query(self,expr):return expr.is_Add
    def value(self,add):
        from symox import Interval
        from symox.calculus.singularities import is_increasing,is_decreasing
        for term in add.args:
            if term.is_number or term in self.bounds or len(term.free_symbols)!=1:continue
            fs,=term.free_symbols
            if fs not in self.bounds:continue
            intrvl=Interval(*self.bounds[fs])
            if is_increasing(term,intrvl,fs):
                self.bounds[term]=(term.subs({fs:self.bounds[fs][0]}),term.subs({fs:self.bounds[fs][1]}))
            elif is_decreasing(term,intrvl,fs):
                self.bounds[term]=(term.subs({fs:self.bounds[fs][1]}),term.subs({fs:self.bounds[fs][0]}))
            else:return add
        if all(term.is_number or term in self.bounds for term in add.args):
            bounds=[(term,term) if term.is_number else self.bounds[term] for term in add.args]
            largest_abs_guarantee=0
            for lo,hi in bounds:
                if lo<=0<=hi:continue
                largest_abs_guarantee=max(largest_abs_guarantee,min(abs(lo),abs(hi)))
            new_terms=[]
            for term,(lo,hi) in zip(add.args,bounds):
                if max(abs(lo),abs(hi))>=largest_abs_guarantee*self.reltol:new_terms.append(term)
            return add.func(*new_terms)
        return add


class SeriesApprox(Optimization):
    """Replace one-variable function calls by a truncated series (sympy.codegen.approximations)."""
    def __init__(self,bounds,reltol,max_order=4,n_point_checks=4,**kwargs):
        super().__init__(**kwargs)
        self.bounds=bounds;self.reltol=reltol;self.max_order=max_order
        if n_point_checks%2==1:
            raise ValueError("Checking the solution at expansion point is not helpful")
        self.n_point_checks=n_point_checks
        self._prec=math.ceil(-math.log10(self.reltol))
    def __call__(self,expr):
        return expr.factor().replace(self.query,lambda arg:self.value(arg))
    def query(self,expr):
        # The reference's UndefinedFunction instance test never holds for an
        # applied function, so every one-argument function call qualifies.
        return expr.is_Function and len(expr.args)==1
    def value(self,fexpr):
        free_symbols=fexpr.free_symbols
        if len(free_symbols)!=1:return fexpr
        symb,=free_symbols
        if symb not in self.bounds:return fexpr
        lo,hi=self.bounds[symb]
        x0=(lo+hi)/2
        cheapest=None
        for n in range(self.max_order+1,0,-1):
            fseri=fexpr.series(symb,x0=x0,n=n).removeO()
            n_ok=True
            for idx in range(self.n_point_checks):
                x=lo+idx*(hi-lo)/(self.n_point_checks-1)
                val=fseri.xreplace({symb:x});ref=fexpr.xreplace({symb:x})
                if abs((1-val/ref).evalf(self._prec))>self.reltol:
                    n_ok=False;break
            if n_ok:cheapest=fseri
            else:break
        return fexpr if cheapest is None else cheapest
__all__=['SumApprox','SeriesApprox']
