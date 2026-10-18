"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Complement DisjointUnion EmptySet FiniteSet Intersection Interval ProductSet Set SymmetricDifference Union UniversalSet imageset')

from .._kind import Kind
from symox import Set,Union,Intersection

class SetKind(Kind):
    _instances={}
    def __new__(cls,element_kind=None):
        if element_kind not in cls._instances:
            obj=super().__new__(cls);obj.element_kind=element_kind;cls._instances[element_kind]=obj
        return cls._instances[element_kind]
    def __repr__(self):return 'SetKind()' if self.element_kind is None else 'SetKind(%r)'%self.element_kind
    def __reduce__(self):return type(self),(self.element_kind,)


def _set_operation(x,y,operation):
    import symox as S
    from itertools import product
    if x==S.EmptySet or y==S.EmptySet:return S.EmptySet
    sx,sy=isinstance(x,S.Set),isinstance(y,S.Set)
    if not sx and not sy:return operation(S.sympify(x),S.sympify(y))
    if (not sx or x.is_finite_set is True) and (not sy or y.is_finite_set is True):
        return S.FiniteSet(*(operation(a,b) for a,b in product(x if sx else [S.sympify(x)],y if sy else [S.sympify(y)])))
    a,b=S.Dummy('a'),S.Dummy('b')
    if sx and sy:return S.imageset(S.Lambda((a,b),operation(a,b)),x,y)
    return S.imageset(S.Lambda(a,operation(a,y) if sx else operation(x,a)),x if sx else y)


def _interval_add(x,y,subtract=False):
    from symox import Interval,FiniteSet,Union
    if isinstance(x,Interval) and isinstance(y,Interval):
        return Interval(x.start-y.end if subtract else x.start+y.start,x.end-y.start if subtract else x.end+y.end,
                        left_open=x.left_open or (y.right_open if subtract else y.left_open),
                        right_open=x.right_open or (y.left_open if subtract else y.right_open))
    if isinstance(x,Union):return Union(*(_interval_add(a,y,subtract) for a in x.args))
    if isinstance(y,Union):return Union(*(_interval_add(x,b,subtract) for b in y.args))
    if isinstance(x,Interval) and not isinstance(y,Set):
        return Interval(x.start-y if subtract else x.start+y,x.end-y if subtract else x.end+y,x.left_open,x.right_open)
    return _set_operation(x,y,(lambda a,b:a-b) if subtract else (lambda a,b:a+b))


def set_add(x,y):return _interval_add(x,y)
def set_sub(x,y):return _interval_add(x,y,True)
def set_mul(x,y):
    from symox import Interval,Min,Max
    if isinstance(x,Interval) and not isinstance(y,Set) and getattr(y,'is_real',isinstance(y,(int,float))) is True:
        if y==0:
            from symox import FiniteSet
            return FiniteSet(0)
        if y>0:return Interval(x.start*y,x.end*y,x.left_open,x.right_open)
        if y<0:return Interval(x.end*y,x.start*y,x.right_open,x.left_open)
    if not isinstance(x,Set) and isinstance(y,Interval):return set_mul(y,x)
    if isinstance(x,Interval) and isinstance(y,Interval) and not (x.left_open or x.right_open or y.left_open or y.right_open):
        values=[a*b for a in (x.start,x.end) for b in (y.start,y.end)]
        if not any(getattr(v,'is_nan',False) or str(v)=='nan' for v in values):return Interval(Min(*values),Max(*values))
    return _set_operation(x,y,lambda a,b:a*b)
def set_div(x,y):return _set_operation(x,y,lambda a,b:a/b)
def set_pow(x,y):return _set_operation(x,y,lambda a,b:a**b)
def set_function(function,x):
    from symox import Dummy,Lambda,imageset
    variable=Dummy('x')
    return imageset(Lambda(variable,function(variable)),x) if isinstance(x,Set) else function(x)
def simplify_union(args):return Union(*args)
def simplify_intersection(args):return Intersection(*args)
def is_function_invertible_in_set(func,setv):
    from symox import Dummy,exp,log,is_increasing,is_decreasing
    if func in (exp,log):return True
    variable=Dummy('u',real=True);expression=func(variable)
    derivative=expression.diff(variable)
    if derivative.is_positive is True or derivative.is_negative is True:return True
    try:
        if is_increasing(expression,setv,variable) or is_decreasing(expression,setv,variable):
            # Non-strict monotonicity alone does not establish injectivity.
            if derivative!=0 and derivative.is_zero is False:return True
    except (ValueError,TypeError,NotImplementedError):pass
    return None

__all__ += ['SetKind','set_add','set_sub','set_mul','set_div','set_pow','set_function',
 'simplify_union','simplify_intersection','is_function_invertible_in_set']
