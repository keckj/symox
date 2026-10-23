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
    if isinstance(t, S.DenseMatrix):return t.trace()
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
    if isinstance(expression,S.DenseMatrix):return expression
    result = _reduce(expression)
    if isinstance(result,TensExpr):return result
    return result*S.eye(4)

def simplify_gamma_expression(expression):return kahane_simplify(expression)

def _sorted_components(ex):
    # reference TensMul._sort_args_for_sorted_components
    coefficient, cv = _parts(ex)
    cv = list(cv)
    sign = 1
    n = len(cv) - 1
    for i in range(n):
        for j in range(n, i, -1):
            c = cv[j-1].component.commutes_with(cv[j].component)
            if c not in (0, 1):
                continue
            typ1 = sorted({t.name for t in cv[j-1].component.index_types})
            typ2 = sorted({t.name for t in cv[j].component.index_types})
            if (typ1, cv[j-1].component.name) > (typ2, cv[j].component.name):
                cv[j-1], cv[j] = cv[j], cv[j-1]
                if c:
                    sign = -sign
    return sign*coefficient, cv


def simplify_gpgp(ex, sort=True):
    """Simplify products ``G(i)*p(-i)*G(j)*p(-j) -> p(i)*p(-i)``.

    Reference algorithm (SymPy 1.14), except that a lone pair is not
    multiplied by 4: ``p_slash*p_slash`` is ``p(i)*p(-i)`` times the
    identity, not its trace (the reference returns ``4*p(i)*p(-i)``).
    """
    if isinstance(ex,TensAdd):return sum((simplify_gpgp(t,sort) for t in ex.args),S.Integer(0))
    if not isinstance(ex,(TensMul,Tensor)):return ex
    from ...tensor.tensor import TensorIndex
    coefficient,factors=_sorted_components(ex) if sort else _parts(ex)
    while True:
        # positions of every factor's slots, as in the reference ``ex.dum``
        comp_map=[]
        for i,f in enumerate(factors):
            comp_map.extend([i]*len(f.indices))
        flat=[ix for f in factors for ix in f.indices]
        dum=[]
        for a,ia in enumerate(flat):
            if ia.is_up:
                for b,ib in enumerate(flat):
                    if ib==-ia:
                        dum.append((a,b,comp_map[a],comp_map[b]))
        dum.sort()
        pairs=[]
        for i,f in enumerate(factors):
            if f.component!=GammaMatrix:continue
            for dx in dum:
                if dx[2]==i:pos=dx[3]
                elif dx[3]==i:pos=dx[2]
                else:continue
                comp=factors[pos].component
                if comp.comm==0 and comp.rank==1:
                    pairs.append((i,pos))
        hit=None
        for k,ai in enumerate(pairs[:-1]):
            if ai[0]!=pairs[k+1][0]-1:continue
            if factors[ai[1]].component!=factors[pairs[k+1][1]].component:continue
            hit=k;break
        if hit is None:
            return _product(coefficient,factors)
        ai,aj=pairs[hit],pairs[hit+1]
        elim={ai[0],ai[1],aj[0],aj[1]}
        head=factors[ai[1]].component
        mu=TensorIndex('mu',GammaMatrix.index_types[0])
        rest=[f for j,f in enumerate(factors) if j not in elim]
        new=_product(coefficient,rest+[head(mu),head(-mu)])
        coefficient,factors=_parts(new)


__all__=['extract_type_tens','gamma_trace','kahane_simplify','simplify_gamma_expression','simplify_gpgp']
