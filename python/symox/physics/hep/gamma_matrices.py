"""Four-dimensional Clifford algebra over native indexed tensors."""
import symox as S
from ...tensor.tensor import TensorIndexType, TensorHead, Tensor, TensExpr, TensMul, TensAdd, contract_metric, canon_bp

LorentzIndex = TensorIndexType('LorentzIndex', dim=4, dummy_name='L')
GammaMatrix = TensorHead('GammaMatrix', [LorentzIndex], comm=2)
_metric = LorentzIndex.metric


def _parts(expression):
    if isinstance(expression, Tensor):return S.Integer(1), [expression]
    if isinstance(expression, TensMul):return expression.coeff, list(expression._factors)
    return S.sympify(expression), []
def _product(coefficient, factors):return TensMul(coefficient, *factors) if factors else coefficient

def extract_type_tens(expression, component):
    """Separate factors of one tensor head, preserving their relative order."""
    if isinstance(expression, TensAdd):raise ValueError('extract one tensor monomial at a time')
    coefficient, factors = _parts(expression)
    return _product(1,[f for f in factors if f.component==component]), _product(coefficient,[f for f in factors if f.component!=component])

def _trace(indices):
    if not indices:return S.Integer(4)
    if len(indices)%2:return S.Integer(0)
    result = S.Integer(0)
    for j in range(1,len(indices)):
        result += (-1)**(j+1)*_metric(indices[0],indices[j])*_trace(indices[1:j]+indices[j+1:])
    return result

def gamma_trace(t):
    if isinstance(t, S.Matrix):return t.trace()
    if isinstance(t, TensAdd):return sum((gamma_trace(term) for term in t.args),S.Integer(0))
    coefficient, factors = _parts(t)
    indices = [f.indices[0] for f in factors if f.component==GammaMatrix]
    other = [f for f in factors if f.component!=GammaMatrix]
    result = _product(coefficient, other)*_trace(indices)
    return canon_bp(contract_metric(result,_metric)) if isinstance(result,TensExpr) else S.simplify(result)

def _reduce(expression):
    if isinstance(expression,TensAdd):return sum((_reduce(term) for term in expression.args),S.Integer(0))
    coefficient, factors = _parts(expression)
    gammas = [i for i,f in enumerate(factors) if f.component==GammaMatrix]
    for k, i in enumerate(gammas):
        mu = factors[i].indices[0]
        for j in gammas[k+1:]:
            if factors[j].indices[0] != -mu:continue
            # Anticommute the first contracted gamma towards its partner.
            middle = [p for p in gammas if i<p<j]
            if not middle:
                remaining=factors[:i]+factors[i+1:j]+factors[j+1:]
                return _reduce(_product(4*coefficient,remaining))
            nextpos=middle[0];nu=factors[nextpos].indices[0]
            exchanged=list(factors);exchanged[i],exchanged[nextpos]=exchanged[nextpos],exchanged[i]
            reduced=[f for pos,f in enumerate(factors) if pos not in (i,nextpos)]+[_metric(mu,nu)]
            contraction=contract_metric(_product(2*coefficient,reduced),_metric)
            return _reduce(contraction)-_reduce(_product(coefficient,exchanged))
    return _product(coefficient,factors)

def kahane_simplify(expression):
    """Eliminate contracted gamma pairs using the four-dimensional Clifford relation."""
    if isinstance(expression,S.Matrix):return expression
    result = _reduce(expression)
    if isinstance(result,TensExpr):return result
    return result*S.eye(4)

def simplify_gamma_expression(expression):return kahane_simplify(expression)

def simplify_gpgp(ex, sort=True):
    """Reduce adjacent equal slashed-vector pairs; retain other ordered factors."""
    if isinstance(ex,TensAdd):return sum((simplify_gpgp(t,sort) for t in ex.args),S.Integer(0))
    coefficient,factors=_parts(ex)
    gamma=[f for f in factors if f.component==GammaMatrix]
    vectors=[f for f in factors if f.component!=GammaMatrix]
    if len(gamma)==2 and len(vectors)==2 and vectors[0].component==vectors[1].component and vectors[0].component.rank==1:
        gi=[g.indices[0] for g in gamma];vi=[v.indices[0] for v in vectors]
        if set(gi)==set(-i for i in vi):
            from ...tensor.tensor import TensorIndex
            i=TensorIndex('p',LorentzIndex);head=vectors[0].component
            return coefficient*head(i)*head(-i)
    return _reduce(ex)

__all__=['extract_type_tens','gamma_trace','kahane_simplify','simplify_gamma_expression','simplify_gpgp']
