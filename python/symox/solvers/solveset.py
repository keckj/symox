"""Set-valued solving through the native scalar equation engine."""
from .. import solveset, solveset_real, solveset_complex, linsolve, nonlinsolve

__all__=['solveset','solveset_real','solveset_complex','linsolve','nonlinsolve']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'NonlinearError linear_eq_to_matrix substitution')

import symox as S


def domain_check(f,symbol,p):
    p=S.sympify(p)
    if p.is_finite is False:return False
    try:
        for node in S.preorder_traversal(S.sympify(f)):
            value=node.subs(symbol,p)
            if value in (S.oo,-S.oo,S.zoo,S.nan) or value.is_finite is False:return False
    except (ValueError,ZeroDivisionError):return False
    return True


def invert_complex(f_x,y,x,domain=S.Complexes):
    f_x,y=S.sympify(f_x),S.sympify(y)
    if not getattr(x,'is_Symbol',False) or not f_x.has(x):raise ValueError('inversion requires a dependent symbol')
    if y.has(x):raise ValueError('the inverse value must be independent of the variable')
    return x,S.solveset(f_x-y,x,domain=domain)
def invert_real(f_x,y,x):return invert_complex(f_x,y,x,S.Reals)


def linear_coeffs(eq,*syms,dict=False):
    matrix,rhs=S.linear_eq_to_matrix([eq],list(syms))
    values=[matrix[0,i] for i in range(len(syms))]+[-rhs[0]]
    if dict:return {symbol:value for symbol,value in zip([*syms,S.Integer(1)],values) if value!=0}
    return values


def solve_decomposition(f,symbol,domain):return S.solveset(f,symbol,domain=domain)


def solvify(f,symbol,domain):
    solution=S.solveset(f,symbol,domain=domain)
    if solution==S.EmptySet:return []
    if isinstance(solution,S.ConditionSet):raise NotImplementedError('the equation did not yield an explicit solution set')
    if solution.is_finite_set is True:return list(solution)
    period=S.periodicity(f,symbol)
    if period is not None:
        principal=solution.intersect(S.Interval.Ropen(0,period))
        if principal.is_finite_set is True:return list(principal)
        roots=S.solve(f,symbol)
        values=[]
        for root in roots:
            if root.free_symbols<=S.sympify(f).free_symbols-{symbol} and domain_check(f,symbol,root):
                if S.Interval.Ropen(0,period).contains(root)==S.true:values.append(root)
        return values or None
    if isinstance(solution,S.Union):
        values=[value for part in solution.args if part.is_finite_set is True for value in part]
        return values or None
    return None

__all__ += ['domain_check','invert_complex','invert_real','linear_coeffs','solve_decomposition','solvify']
