"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'factor_nc factor_terms gcd_terms')

import symox as _S


def decompose_power(expr):
    expr=_S.sympify(expr);base,exponent=expr.as_base_exp()
    if exponent.is_Rational:return base**_S.Rational(1,exponent.q),int(exponent.p)
    coefficient,tail=exponent.as_coeff_Mul()
    if coefficient.is_Rational:return base**(tail/coefficient.q),int(coefficient.p)
    return expr,1

def decompose_power_rat(expr):
    expr=_S.sympify(expr);base,exponent=expr.as_base_exp()
    return (base,exponent) if exponent.is_Rational else decompose_power(expr)


class Factors:
    def __init__(self,factors=None):
        if isinstance(factors,Factors):factors=factors.factors
        elif factors is None:factors={}
        elif not isinstance(factors,dict):
            expr=_S.sympify(factors)
            if expr==1:factors={}
            else:
                factors={}
                for base,power in expr.as_powers_dict().items():
                    if base.is_Rational and base.q!=1:
                        if base.p!=1:factors[_S.Integer(base.p)]=factors.get(_S.Integer(base.p),0)+power
                        factors[_S.Integer(base.q)]=factors.get(_S.Integer(base.q),0)-power
                    else:factors[base]=factors.get(base,0)+power
        self.factors={_S.sympify(base):_S.sympify(power) for base,power in factors.items() if power!=0 and base!=1}
        self.gens=frozenset(self.factors)
    def as_expr(self):return _S.Mul(*(base**power for base,power in self.factors.items()))
    def mul(self,other):
        factors=dict(self.factors)
        for base,power in other.factors.items():factors[base]=factors.get(base,0)+power
        return type(self)(factors)
    def quo(self,other):return self.mul(other.pow(-1))
    def pow(self,exponent):return type(self)({base:power*exponent for base,power in self.factors.items()})
    def gcd(self,other):return type(self)({base:min(self.factors.get(base,0),other.factors.get(base,0)) for base in self.gens|other.gens})
    def lcm(self,other):return type(self)({base:max(self.factors.get(base,0),other.factors.get(base,0)) for base in self.gens|other.gens})
    def normal(self,other):
        common=self.gcd(other)
        return self.quo(common),other.quo(common)
    def div(self,other):return self.normal(other)
    def rem(self,other):return self.div(other)[1]
    __mul__=mul
    __truediv__=quo
    __pow__=pow
    __divmod__=div
    def __eq__(self,other):return isinstance(other,Factors) and self.factors==other.factors
    def __hash__(self):return hash(frozenset(self.factors.items()))
    def __repr__(self):return 'Factors(%r)'%self.factors


class Term:
    def __init__(self,term,numer=None,denom=None):
        term=_S.sympify(term)
        if numer is None and denom is None:
            if term.is_commutative is not True:raise ValueError('a term requires commutative factors')
            self.coeff,rest=term.as_coeff_Mul();a,b=rest.as_numer_denom();self.numer,self.denom=Factors(a),Factors(b)
        else:self.coeff,self.numer,self.denom=term,Factors(numer),Factors(denom)
    def as_expr(self):return self.coeff*self.numer.as_expr()/self.denom.as_expr()
    def mul(self,other):
        a,b=self.numer.mul(other.numer).normal(self.denom.mul(other.denom))
        return type(self)(self.coeff*other.coeff,a,b)
    def inv(self):return type(self)(1/self.coeff,self.denom,self.numer)
    def quo(self,other):return self.mul(other.inv())
    def pow(self,n):return type(self)(self.coeff**n,self.numer.pow(n),self.denom.pow(n)) if n>=0 else self.inv().pow(-n)
    def gcd(self,other):return type(self)(_S.gcd(self.coeff,other.coeff),self.numer.gcd(other.numer),self.denom.gcd(other.denom))
    def lcm(self,other):return type(self)(_S.lcm(self.coeff,other.coeff),self.numer.lcm(other.numer),self.denom.lcm(other.denom))
    __mul__=mul
    __truediv__=quo
    __pow__=pow
    def __eq__(self,other):return isinstance(other,Term) and self.coeff==other.coeff and self.numer==other.numer and self.denom==other.denom
    def __hash__(self):return hash((self.coeff,self.numer,self.denom))
    def __repr__(self):return 'Term(%r, %r, %r)'%(self.coeff,self.numer,self.denom)

__all__ += ['decompose_power','decompose_power_rat','Factors','Term']
