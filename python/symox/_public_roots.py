"""Exact algebraic root objects and symmetric sums over polynomial roots.

Real numerical evaluation uses rational Sturm intervals. Complex evaluation
refines square-free factors using mpmath; complex boxes are not certified.
"""
import operator
from . import _core
from ._compat import _ExprType

_roots = {}
_sums = {}


def _polynomial(f,x=None):
    from . import Poly,Symbol,sympify
    if isinstance(f,Poly):
        if x is None:x = f.gens[0]
        f = f.as_expr()
    f = sympify(f)
    if x is None:
        symbols = sorted(f.free_symbols,key=str)
        if len(symbols) != 1:raise ValueError('a univariate polynomial is required')
        x = symbols[0]
    x = sympify(x)
    polynomial = Poly(f,x)
    if len(polynomial.gens) != 1 or polynomial.degree() < 1:raise ValueError('polynomial degree must be positive')
    if any(c.is_Rational is not True for c in polynomial.all_coeffs()):
        raise NotImplementedError('algebraic root isolation currently requires rational coefficients')
    return polynomial


def _mp_roots(polynomial,digits):
    import mpmath
    mp = mpmath.mp.clone()
    from mpmath.libmp.libhyper import NoConvergence
    values = []
    with mp.workdps(digits+25):
        for factor,multiplicity in polynomial.factor_list()[1]:
            coeffs = [mp.mpf(str(c.p))/int(c.q) for c in factor.all_coeffs()]
            try:
                roots = mp.polyroots(coeffs,maxsteps=1000,extraprec=max(80,digits),error=False)
            except NoConvergence:
                seeds = factor.nroots()
                def f(z):return mp.polyval(coeffs,z)
                roots = []
                for seed in seeds:
                    re,im = seed.as_real_imag()
                    z = mp.mpc(float(re),float(im))
                    roots.append(mp.findroot(f,(z,z+mp.mpc('0.001','0.001')),solver='secant',maxsteps=1000))
            values.extend((root,multiplicity) for root in roots)
        threshold = mp.power(10,-digits-5)
        reals,complexes = [],[]
        for value,multiplicity in values:
            bucket = reals if abs(mp.im(value)) < threshold else complexes
            bucket.extend([mp.re(value) if bucket is reals else value]*multiplicity)
        reals.sort()
        # Round only for ordering roots with theoretically equal real parts.
        complexes.sort(key=lambda z:(mp.nint(mp.re(z)*mp.power(10,digits//2)),mp.im(z)))
        return reals+complexes


class ComplexRootOf(metaclass=_ExprType):
    _matches = staticmethod(lambda expr:expr.class_name == 'ComplexRootOf')
    def __new__(cls,f,x,index=None,radicals=False,expand=True):
        if index is None:index,x = x,None
        index = operator.index(index)
        p = _polynomial(f,x)
        degree = int(p.degree())
        if index < 0:index += degree
        if not 0 <= index < degree:raise IndexError('root index out of range')
        from . import Tuple,Integer
        # Rational roots are always returned exactly, including repeated roots.
        candidates = None
        try:
            roots = p.all_roots(radicals=radicals)
            if len(roots) == degree:
                candidates = sorted(roots,key=lambda r:(0,float(r.evalf(18))) if r.is_real is True else (1,float(r.as_real_imag()[0].evalf(18)),float(r.as_real_imag()[1].evalf(18))))
        except (TypeError,ValueError,NotImplementedError):pass
        if candidates is not None and (radicals or candidates[index].is_Rational):return candidates[index]
        intervals = p.intervals()
        real_count,position = sum(m for bounds,m in intervals),0
        for (a,b),multiplicity in intervals:
            if position <= index < position+multiplicity and a == b:return a
            position += multiplicity
        node = _core._undefined_function('ComplexRootOf',(Tuple(*p.monic().all_coeffs()),Integer(index)),
                                        {'real':True,'algebraic':True} if index < real_count else {'algebraic':True},None)
        _roots[node] = (p,index)
        return node


CRootOf = ComplexRootOf


class RootOf(metaclass=_ExprType):
    _matches = ComplexRootOf._matches
    def __new__(cls,f,x,index=None,radicals=True,expand=True):
        return ComplexRootOf(f,x,index,radicals=radicals,expand=expand)


rootof = RootOf


def _metadata(node):
    if node in _roots:return _roots[node]
    from . import Poly,Dummy,Add
    coefficients,index = node.args
    x = Dummy('x')
    coefficients = coefficients.args
    p = Poly(Add(*(c*x**(len(coefficients)-i-1) for i,c in enumerate(coefficients))),x)
    return p,int(index)


def numeric_root(node,digits):
    from . import Float,Rational,I
    p,index = _metadata(node)
    intervals = p.intervals()
    position = 0
    for bounds,multiplicity in intervals:
        if position <= index < position+multiplicity:
            a,b = bounds
            if a == b:return a.evalf(digits)
            while a == 0 or b == 0 or max(abs(a),abs(b)) > 2*min(abs(a),abs(b)):
                a,b = p.refine_root(a,b,steps=32)
                if a == b:return a.evalf(digits)
            eps = min(abs(a),abs(b))*Rational(1,10**(digits+8))
            a,b = p.refine_root(a,b,eps=eps)
            return ((a+b)/2).evalf(digits)
        position += multiplicity
    import mpmath
    mp = mpmath.mp.clone()
    value = _mp_roots(p,digits+8)[index]
    with mp.workdps(digits+12):
        return Float(mp.nstr(mp.re(value),digits+4),digits)+I*Float(mp.nstr(mp.im(value),digits+4),digits)


class RootSum(metaclass=_ExprType):
    _matches = staticmethod(lambda expr:expr.class_name == 'RootSum')
    def __new__(cls,expr,func=None,x=None,auto=True,quadratic=False):
        from . import Lambda,Dummy,Tuple,Matrix,eye,zeros,sympify
        p = _polynomial(expr,x)
        variable = p.gens[0]
        if func is None:func = Lambda(variable,variable)
        if not isinstance(func,Lambda):raise ValueError('root sum requires a univariate Lambda')
        variables = func.variables
        variables = variables.args if getattr(variables,'is_Tuple',False) else tuple(variables)
        if len(variables) != 1:raise ValueError('root sum requires one bound variable')
        t = variables[0]
        body = func.expr
        if auto:
            # The trace of r(C), C the companion matrix, equals sum r(root).
            numerator,denominator = body.as_numer_denom()
            from . import Poly
            try:
                a,b = Poly(numerator,t),Poly(denominator,t)
            except (ValueError,TypeError):pass
            else:
                n = int(p.degree())
                coefficients = p.monic().all_coeffs()
                C = zeros(n)
                for i in range(1,n):C[i,i-1] = 1
                for i in range(n):C[i,n-1] = -coefficients[n-i]
                def matrix_polynomial(poly):
                    value = zeros(n)
                    for coefficient in poly.all_coeffs():value = value*C+coefficient*eye(n)
                    return value
                return (matrix_polynomial(a)*matrix_polynomial(b).inv()).trace().simplify()
        node = _core._apply('RootSum',Tuple(*p.monic().all_coeffs()),func,evaluate=False)
        _sums[node] = p,func
        return node


def root_sum_doit(node):
    from . import Dummy,Poly,Add
    if node in _sums:p,function = _sums[node]
    else:
        coefficients,function = node.args
        coefficients = coefficients.args
        x = Dummy('x')
        p = Poly(Add(*(c*x**(len(coefficients)-i-1) for i,c in enumerate(coefficients))),x)
    return Add(*(function(RootOf(p,i)) for i in range(int(p.degree()))))


def install():
    Expr = _core.Expr
    original_is_number = Expr.is_number
    Expr.is_number = property(lambda self:True if self.class_name == 'ComplexRootOf' else original_is_number.__get__(self,Expr))
    previous_poly = getattr(Expr,'poly',None)
    def poly(self):
        if self.class_name == 'ComplexRootOf':return _metadata(self)[0]
        if previous_poly is not None:return previous_poly.__get__(self,Expr)
        raise AttributeError('poly')
    Expr.poly = property(poly)
    previous_index = getattr(Expr,'index',None)
    def index(self):
        if self.class_name == 'ComplexRootOf':return _metadata(self)[1]
        if previous_index is not None:return previous_index.__get__(self,Expr)
        raise AttributeError('index')
    Expr.index = property(index)
