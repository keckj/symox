"""Additional special functions, with exact identities and mpmath evaluation."""
import operator
from functools import lru_cache
from . import _core
from ._compat import _ExprType


def _held(name,*args):
    return _core._apply(name,*args,evaluate=False)


def betainc(a,b,x1,x2):
    from . import sympify,beta,expand
    a,b,x1,x2 = map(sympify,(a,b,x1,x2))
    if x1 == x2:
        return sympify(0)
    if x1 == 0 and x2 == 1:
        return beta(a,b)
    if a.is_Integer and b.is_Integer and a > 0 and b > 0:
        from math import comb
        return expand(sum(((-1)**j*comb(int(b)-1,j)*(x2**(a+j)-x1**(a+j))/(a+j) for j in range(int(b))),sympify(0)))
    return _held('betainc',a,b,x1,x2)


def betainc_regularized(a,b,x1,x2):
    from . import sympify,beta
    a,b,x1,x2 = map(sympify,(a,b,x1,x2))
    if x1 == x2:return sympify(0)
    if x1 == 0 and x2 == 1:return sympify(1)
    value = betainc(a,b,x1,x2)
    return value/beta(a,b) if value.class_name != 'betainc' else _held('betainc_regularized',a,b,x1,x2)


def erf2inv(x,y):
    from . import erf,erfinv,sympify
    x,y = sympify(x),sympify(y)
    return x if y == 0 else erfinv(y+erf(x))


def Ynm_c(n,m,theta,phi):
    from . import Ynm,conjugate
    return conjugate(Ynm(n,m,theta,phi))


def Znm(n,m,theta,phi):
    from . import Ynm,sqrt,I,sympify
    n,m,theta,phi = map(sympify,(n,m,theta,phi))
    y = Ynm(n,m,theta,phi)
    if m.is_zero:return y
    if m.is_positive:return (y+Ynm_c(n,m,theta,phi))/sqrt(2)
    if m.is_negative:return (y-Ynm_c(n,m,theta,phi))/(sqrt(2)*I)
    return _held('Znm',n,m,theta,phi)


def chebyshevt_root(n,k):
    from . import cos,pi,sympify
    n,k = operator.index(n),operator.index(k)
    if n <= 0 or not 0 <= k < n:raise ValueError('require n > 0 and 0 <= k < n')
    return cos(pi*sympify(2*k+1)/(2*n))


def chebyshevu_root(n,k):
    from . import cos,pi,sympify
    n,k = operator.index(n),operator.index(k)
    if n <= 0 or not 0 <= k < n:raise ValueError('require n > 0 and 0 <= k < n')
    return cos(pi*sympify(k+1)/(n+1))


def hermite_prob(n,x):
    from . import sympify
    from .polys.orthopolys import hermite_prob_poly
    n,x = sympify(n),sympify(x)
    if n.is_integer is False or n.is_negative is True:raise ValueError('polynomial degree must be a nonnegative integer')
    if n.is_Integer:return hermite_prob_poly(int(n),x)
    return _held('hermite_prob',n,x)


def jacobi_normalized(n,a,b,x):
    from . import jacobi,gamma,sqrt,factorial,sympify
    n,a,b,x = map(sympify,(n,a,b,x))
    norm = 2**(a+b+1)*gamma(n+a+1)*gamma(n+b+1)/((2*n+a+b+1)*gamma(n+a+b+1)*factorial(n))
    return jacobi(n,a,b,x)/sqrt(norm)


def riemann_xi(s):
    from . import sympify,Rational
    s = sympify(s)
    if s in (0,1):return Rational(1,2)
    return _held('riemann_xi',s)


def stieltjes(n,a=1):
    from . import sympify,digamma
    n,a = sympify(n),sympify(a)
    if n.is_integer is False or n.is_negative is True:raise ValueError('Stieltjes index must be nonnegative and integral')
    if a.is_integer is True and a.is_nonpositive is True:raise ValueError('Stieltjes constant is singular at this argument')
    if n == 0:return -digamma(a)
    return _held('stieltjes',n,a)


class SingularityFunction(metaclass=_ExprType):
    _matches = staticmethod(lambda expr:expr.class_name == 'SingularityFunction')
    def __new__(cls,x,a,n):
        from . import sympify,oo
        x,a,n = map(sympify,(x,a,n))
        if n.is_integer is False:raise ValueError('singularity order must be integral')
        displacement = x-a
        if displacement.is_negative is True:return sympify(0)
        if displacement.is_positive is True:
            if n.is_nonnegative is True:return displacement**n
            if n.is_negative is True:return sympify(0)
        if displacement == 0:
            if n.is_positive is True:return sympify(0)
            if n == 0:return sympify(1)
            if n.is_negative is True:return oo
        return _held('SingularityFunction',x,a,n)


def singularityintegrate(f,x):
    from . import sympify,expand,Add,binomial
    f,x = sympify(f),sympify(x)
    if f.is_Add:
        terms = [singularityintegrate(t,x) for t in f.args]
        return None if any(t is None for t in terms) else Add(*terms)
    factors = f.args if f.is_Mul else (f,)
    sf = [v for v in factors if v.class_name == 'SingularityFunction']
    if len(sf) != 1:return None
    node = sf[0]
    variable,a,n = node.args
    if variable != x or n.is_Integer is not True:return None
    weight = f/node
    if weight.has(x) is False:
        return weight*SingularityFunction(x,a,n+1)/(n+1 if n >= 0 else 1)
    if n < 0:
        k = -int(n)-1
        return Add(*((-1)**j*binomial(k,j)*weight.diff(x,j).subs(x,a)*SingularityFunction(x,a,n+j+1) for j in range(k+1)))
    try:
        from . import Poly
        degree = int(Poly(weight,x).degree())
    except (ValueError,TypeError):return None
    from . import factorial
    return expand(Add(*(weight.diff(x,j).subs(x,a)/factorial(j)*SingularityFunction(x,a,n+j+1)/(n+j+1) for j in range(degree+1))))


_NUMERIC = {'betainc','betainc_regularized','riemann_xi','stieltjes'}
_SPECIAL = _NUMERIC | {'hermite_prob','SingularityFunction','Znm','ComplexRootOf','RootSum','Permanent'}


def _nodes(expr):
    seen,stack = set(),[expr]
    while stack:
        node = stack.pop()
        if not isinstance(node,_core.Expr) or node in seen:continue
        seen.add(node)
        yield node
        stack.extend(node.args)


@lru_cache(maxsize=8192)
def _special_nodes(expr):
    return tuple(node for node in _nodes(expr) if node.class_name in _SPECIAL)


def _numeric_value(node,digits):
    import mpmath
    mp = mpmath.mp.clone()
    from . import Float,I
    with mp.workdps(digits+12):
        def number(arg):
            re,im = arg.as_real_imag()
            return mp.mpc(str(re.evalf(digits+12)),str(im.evalf(digits+12))) if im != 0 else mp.mpf(str(re.evalf(digits+12)))
        args = [number(a) for a in node.args]
        name = node.class_name
        if name in ('betainc','betainc_regularized'):value = mp.betainc(*args,regularized=name.endswith('_regularized'))
        elif name == 'stieltjes':value = mp.stieltjes(int(args[0]),args[1])
        elif name == 'riemann_xi':
            z = args[0]
            value = mp.mpf('.5') if z in (0,1) else z*(z-1)*mp.gamma(z/2)*mp.zeta(z)/(2*mp.pi**(z/2))
        else:return node
        real,imag = mp.re(value),mp.im(value)
        return Float(mp.nstr(real,digits+5),digits)+I*Float(mp.nstr(imag,digits+5),digits)


def _partial(node,index):
    from . import diff,log,Integral,Dummy,gamma,zeta,pi,hermite_prob,SingularityFunction,DiracDelta
    name,args = node.class_name,node.args
    if name in ('betainc','betainc_regularized'):
        a,b,l,u = args
        if index >= 2:
            value = (-1 if index == 2 else 1)*args[index]**(a-1)*(1-args[index])**(b-1)
            if name.endswith('_regularized'):
                from . import beta
                value /= beta(a,b)
            return value
        t = Dummy('t')
        from . import beta,digamma
        value = Integral(t**(a-1)*(1-t)**(b-1)*log(t if index == 0 else 1-t),(t,l,u))
        if name.endswith('_regularized'):
            value = value/beta(a,b)-node*(digamma(args[index])-digamma(a+b))
        return value
    if name == 'hermite_prob' and index == 1:return args[0]*hermite_prob(args[0]-1,args[1])
    if name == 'SingularityFunction' and index in (0,1):
        x,a,n = args
        return (-1 if index else 1)*(n if n.is_positive is True else 1)*SingularityFunction(x,a,n-1)
    if name == 'riemann_xi':
        t = Dummy("xi_argument")
        return diff(t*(t-1)*gamma(t/2)*zeta(t)/(2*pi**(t/2)),t).subs(t,args[0])
    return None


def install(namespace):
    from functools import wraps
    Expr = _core.Expr
    native_heaviside = namespace['Heaviside']
    @wraps(native_heaviside)
    def Heaviside(arg,H0=None,*,evaluate=None):
        from . import sympify,Rational,Piecewise,Eq
        arg = sympify(arg)
        if H0 is None or sympify(H0) == Rational(1,2):return native_heaviside(arg,evaluate=evaluate)
        H0 = sympify(H0)
        if arg == 0:return H0
        if arg.is_nonzero is True:return native_heaviside(arg,evaluate=evaluate)
        return Piecewise((H0,Eq(arg,0)),(native_heaviside(arg,evaluate=evaluate),True))
    from . import _BuiltinFunctionHead
    namespace['Heaviside'] = _BuiltinFunctionHead('Heaviside', Heaviside)
    old_evalf,old_diff,old_doit,old_rewrite = Expr.evalf,Expr.diff,Expr.doit,Expr.rewrite
    @wraps(old_evalf)
    def evalf(self,n=15,*args,**kwargs):
        value = self
        substitutions = args[0] if args else kwargs.get('subs')
        # Ordinary scalar evaluation must reach the native substitution path
        # before arithmetic can round the replacement values. Keep the special
        # adapter when a special node is present or inserted by the mapping.
        if not _special_nodes(value) and not (
                isinstance(substitutions,dict) and
                any(_special_nodes(v) for v in substitutions.values())):
            return old_evalf(value,n,*args,**kwargs)
        if substitutions:value = value.subs(substitutions)
        replacements = {}
        for node in reversed(_special_nodes(value)):
            replaced = node.xreplace(replacements) if replacements else node
            if node.class_name in _NUMERIC and not replaced.free_symbols:
                replacements[node] = _numeric_value(replaced,int(n))
            elif node.class_name == 'ComplexRootOf':
                from ._public_roots import numeric_root
                replacements[node] = numeric_root(node,int(n))
        if replacements:value = value.xreplace(replacements)
        return old_evalf(value,n,*args,**kwargs)
    Expr.evalf = evalf
    @wraps(old_diff)
    def derivative(self,*symbols,**kwargs):
        if not any(node.class_name not in ('ComplexRootOf','Permanent','RootSum') for node in _special_nodes(self)):
            return old_diff(self,*symbols,**kwargs)
        if kwargs.get('evaluate') is False:return old_diff(self,*symbols,**kwargs)
        if len(symbols) == 1 and hasattr(symbols[0],'shape'):return old_diff(self,*symbols,**kwargs)
        from . import sympify,Dummy,Add,Derivative
        specs = []
        i = 0
        while i < len(symbols):
            spec = symbols[i]
            if isinstance(spec,(list,tuple)) or getattr(spec,'is_Tuple',False):
                variable,count = spec.args if getattr(spec,'is_Tuple',False) else spec
            else:
                variable,count = spec,1
                if i+1 < len(symbols) and isinstance(symbols[i+1],int):i += 1;count = symbols[i]
            specs.extend([sympify(variable)]*operator.index(count));i += 1
        if not specs:
            free = self.free_symbols
            if len(free) != 1:return old_diff(self,*symbols,**kwargs)
            specs = list(free)
        value = self
        for variable in specs:
            memo = {}
            def d(node):
                if node in memo:return memo[node]
                if not node.has(variable):return sympify(0)
                if node.class_name in _SPECIAL:
                    pieces = []
                    for index,arg in enumerate(node.args):
                        if not arg.has(variable):continue
                        partial = _partial(node,index)
                        if partial is None:
                            t = Dummy('partial_argument')
                            partial = Derivative(node.xreplace({arg:t}),t,evaluate=False).subs(t,arg)
                        pieces.append(partial*d(arg))
                    result = Add(*pieces)
                elif node.is_Add:result = Add(*(d(a) for a in node.args))
                elif node.is_Mul:result = Add(*(d(a)*node/a for a in node.args))
                elif node.is_Pow:
                    from . import log
                    a,b = node.args
                    result = b*a**(b-1)*d(a)+node*log(a)*d(b)
                else:
                    children = [a for a in node.args if isinstance(a,Expr) and any(v.class_name in _SPECIAL for v in _nodes(a))]
                    replacements = {a:Dummy('arg') for a in children}
                    transformed = node.xreplace(replacements)
                    inverse = {v:k for k,v in replacements.items()}
                    result = old_diff(transformed,variable).xreplace(inverse)
                    result += Add(*(old_diff(transformed,t).xreplace(inverse)*d(a) for a,t in replacements.items()))
                memo[node] = result
                return result
            value = d(value)
        return value
    Expr.diff = derivative
    original_fdiff = Expr.fdiff
    @wraps(original_fdiff)
    def fdiff(self,argindex=None):
        if self.class_name in _SPECIAL:
            index = 0 if argindex is None else operator.index(argindex)-1
            if not 0 <= index < len(self.args):raise ValueError('function argument index is out of range')
            partial = _partial(self,index)
            if partial is not None:return partial
        return original_fdiff(self,argindex)
    Expr.fdiff = fdiff
    old_root_diff = namespace['diff']
    from ._dispatch_arguments import FirstArgument
    first_diff=FirstArgument(old_root_diff)
    @wraps(old_root_diff)
    def diff(*args,**kwargs):
        expr=first_diff.get(args,kwargs)
        if isinstance(expr,Expr) and _special_nodes(expr):
            expr,rest,options=first_diff.split(args,kwargs)
            return expr.diff(*rest,**options)
        return old_root_diff(*args,**kwargs)
    namespace['diff'] = diff
    @wraps(old_doit)
    def doit(self,*args,**kwargs):
        if self.class_name == 'RootSum':
            from ._public_roots import root_sum_doit
            return root_sum_doit(self)
        if self.class_name == 'Permanent':
            matrix = self.args[0]
            if matrix.kind == 'ImmutableMatrix':return matrix.as_explicit().per()
        return old_doit(self,*args,**kwargs)
    Expr.doit = doit
    @wraps(old_rewrite)
    def rewrite(self,*args,**kwargs):
        if len(args) != 1:return old_rewrite(self,*args,**kwargs)
        target = args[0]
        from . import Piecewise,DiracDelta,Heaviside
        name = target if isinstance(target,str) else getattr(target,'__name__',str(target))
        if self.class_name == 'riemann_xi' and name in ('zeta','gamma'):
            from . import gamma,zeta,pi
            s = self.args[0]
            return s*(s-1)*gamma(s/2)*zeta(s)/(2*pi**(s/2))
        if self.class_name == 'SingularityFunction':
            x,a,n = self.args
            if name == 'Piecewise':
                from . import oo,Eq
                return Piecewise((oo,Eq(x,a)),(0,True)) if n.is_negative is True else Piecewise(((x-a)**n,x>=a),(0,True))
            if name in ('Heaviside','DiracDelta') and n.is_Integer:
                if n < 0:
                    return DiracDelta(x-a) if n == -1 else DiracDelta(x-a,-int(n)-1)
                return (x-a)**n*Heaviside(x-a,1)
        return old_rewrite(self,target,**kwargs)
    Expr.rewrite = rewrite
    original_subs,original_xreplace = Expr.subs,Expr.xreplace
    constructors = {name:globals()[name] for name in _NUMERIC|{'SingularityFunction','hermite_prob','Znm'}}
    def reevaluate(value):
        if not isinstance(value,Expr):return value
        replacements = {}
        for node in reversed(_special_nodes(value)):
            constructor = constructors.get(node.class_name)
            if constructor is None:continue
            updated = original_xreplace(node,replacements) if replacements else node
            evaluated = constructor(*updated.args)
            if evaluated != node:replacements[node] = evaluated
        return original_xreplace(value,replacements) if replacements else value
    @wraps(original_subs)
    def subs(self,*args,**kwargs):
        result = original_subs(self,*args,**kwargs)
        return reevaluate(result) if _special_nodes(self) else result
    @wraps(original_xreplace)
    def xreplace(self,rule):
        result = original_xreplace(self,rule)
        return reevaluate(result) if _special_nodes(self) else result
    Expr.subs,Expr.xreplace = subs,xreplace
    old_N = namespace['N']
    def N(x,n=15,**kwargs):
        from . import sympify
        if not isinstance(x,Expr) and hasattr(x,'applyfunc'):
            if getattr(x,'is_NDimArray',False):x=x.as_immutable()
            return x.applyfunc(lambda entry:N(entry,n,**kwargs))
        value = sympify(x)
        return value.evalf(n,**kwargs) if any(v.class_name in _NUMERIC|{'ComplexRootOf'} for v in _special_nodes(value)) else old_N(x,n,**kwargs)
    namespace['N'] = N
