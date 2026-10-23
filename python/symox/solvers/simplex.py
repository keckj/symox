"""Exact two-phase primal simplex with Bland pivots and arbitrary variable bounds."""
import symox as S


class InfeasibleLPError(ValueError):pass
class UnboundedLPError(ValueError):pass


def _sign(x):
    x=S.simplify(x)
    if x==0:return 0
    if x.is_positive is True:return 1
    if x.is_negative is True:return -1
    if x.is_number:return 1 if float(x.evalf())>0 else -1
    raise ValueError('linear programming requires numeric coefficients')


def _optimize(rows,basis,cost):
    while True:
        reduced=[S.simplify(c-sum(cost[b]*row[j] for b,row in zip(basis,rows))) for j,c in enumerate(cost)]
        entering=next((j for j,c in enumerate(reduced) if j not in basis and _sign(c)<0),None)
        if entering is None:break
        eligible=[i for i,row in enumerate(rows) if _sign(row[entering])>0]
        if not eligible:raise UnboundedLPError('\nObjective function can assume arbitrarily large values!')
        leaving=min(eligible,key=lambda i:(rows[i][-1]/rows[i][entering],basis[i]))
        pivot=rows[leaving][entering];rows[leaving]=[S.cancel(v/pivot) for v in rows[leaving]]
        for i in range(len(rows)):
            if i==leaving:continue
            coefficient=rows[i][entering]
            if coefficient!=0:rows[i]=[S.cancel(a-coefficient*b) for a,b in zip(rows[i],rows[leaving])]
        basis[leaving]=entering
    solution=[S.Integer(0)]*len(cost)
    for i,b in enumerate(basis):solution[b]=S.simplify(rows[i][-1])
    return S.simplify(sum(c*v for c,v in zip(cost,solution))),solution


def _simplex(cost,constraints):
    n=len(cost);rows=[];basis=[];artificial=set();width=n
    # Build a feasible artificial basis, reversing constraints with negative RHS.
    records=[]
    for coefficients,relation,rhs in constraints:
        coefficients=list(coefficients);rhs=S.sympify(rhs)
        if _sign(rhs)<0:
            coefficients=[-a for a in coefficients];rhs=-rhs;relation={'<=':'>=','>=':'<=','=':'='}[relation]
        extras={}
        if relation=='<=':extras[width]=S.Integer(1);basic=width;width+=1
        else:
            if relation=='>=':extras[width]=S.Integer(-1);width+=1
            extras[width]=S.Integer(1);basic=width;artificial.add(width);width+=1
        records.append((coefficients,extras,rhs));basis.append(basic)
    rows=[[*(coefficients+[S.Integer(0)]*(width-n)),rhs] for coefficients,_,rhs in records]
    for row,(_,extras,_) in zip(rows,records):
        for j,value in extras.items():row[j]=value
    phase1=[S.Integer(int(j in artificial)) for j in range(width)]
    optimum,_=_optimize(rows,basis,phase1)
    if _sign(optimum)!=0:raise InfeasibleLPError('\nThe constraint set is empty!')
    redundant=[]
    for i,b in enumerate(basis):
        if b not in artificial:continue
        pivot=next((j for j in range(width) if j not in artificial and j not in basis and rows[i][j]!=0),None)
        if pivot is None:redundant.append(i);continue
        factor=rows[i][pivot];rows[i]=[S.cancel(v/factor) for v in rows[i]]
        for k in range(len(rows)):
            if k!=i:
                factor=rows[k][pivot];rows[k]=[S.cancel(a-factor*b) for a,b in zip(rows[k],rows[i])]
        basis[i]=pivot
    for i in reversed(redundant):rows.pop(i);basis.pop(i)
    columns=[j for j in range(width) if j not in artificial];renumber={j:i for i,j in enumerate(columns)}
    rows=[[row[j] for j in columns]+[row[-1]] for row in rows];basis=[renumber[b] for b in basis]
    fullcost=[S.sympify(cost[j]) if j<n else S.Integer(0) for j in columns]
    optimum,solution=_optimize(rows,basis,fullcost)
    return optimum,solution[:n]


def _matrix_rows(matrix,n):
    if matrix is None:return []
    if isinstance(matrix,S.DenseMatrix):
        if matrix.cols!=n:raise ValueError('constraint matrix has the wrong number of columns')
        return [[S.sympify(matrix[i,j]) for j in range(n)] for i in range(matrix.rows)]
    rows=[list(map(S.sympify,row)) for row in matrix]
    if any(len(row)!=n for row in rows):raise ValueError('constraint matrix has the wrong number of columns')
    return rows


def _bounds(bounds,n):
    if bounds is None:return [(0,None)]*n
    if isinstance(bounds,dict):return [bounds.get(i,(0,None)) for i in range(n)]
    if len(bounds)==2 and all(v is None or not isinstance(v,(tuple,list)) for v in bounds):return [tuple(bounds)]*n
    if len(bounds)!=n:raise ValueError('bounds must contain one pair per variable')
    return list(bounds)


def linprog(c,A=None,b=None,A_eq=None,b_eq=None,bounds=None):
    c=list(map(S.sympify,list(c)));n=len(c);constraints=[]
    for matrix,rhs,relation in ((A,b,'<='),(A_eq,b_eq,'=')):
        rows=_matrix_rows(matrix,n);values=[] if rhs is None else list(map(S.sympify,list(rhs)))
        if len(rows)!=len(values):raise ValueError('constraint matrix and right-hand side sizes differ')
        constraints.extend((row,relation,value) for row,value in zip(rows,values))
    offsets=[];mapping=[];upper=[];columns=0
    for i,(lower,higher) in enumerate(_bounds(bounds,n)):
        lower=None if lower is None or lower==-S.oo else S.sympify(lower)
        higher=None if higher is None or higher==S.oo else S.sympify(higher)
        if lower is not None and higher is not None and _sign(higher-lower)<0:raise InfeasibleLPError('lower bound exceeds upper bound')
        if lower is not None:
            offsets.append(lower);mapping.append({columns:S.Integer(1)});columns+=1
            if higher is not None:upper.append((i,higher-lower))
        elif higher is not None:
            offsets.append(higher);mapping.append({columns:S.Integer(-1)});columns+=1
        else:
            offsets.append(S.Integer(0));mapping.append({columns:S.Integer(1),columns+1:S.Integer(-1)});columns+=2
    transformed=[]
    for row,relation,rhs in constraints:
        coeff=[S.Integer(0)]*columns
        for i,a in enumerate(row):
            for j,factor in mapping[i].items():coeff[j]+=a*factor
        transformed.append((coeff,relation,rhs-sum(a*v for a,v in zip(row,offsets))))
    for i,bound in upper:
        coeff=[S.Integer(0)]*columns
        for j,value in mapping[i].items():coeff[j]=value
        transformed.append((coeff,'<=',bound))
    objective=[S.Integer(0)]*columns
    for i,a in enumerate(c):
        for j,factor in mapping[i].items():objective[j]+=a*factor
    optimum,solution=_simplex(objective,transformed)
    original=[S.simplify(offset+sum(v*solution[j] for j,v in row.items())) for offset,row in zip(offsets,mapping)]
    return S.simplify(sum(a*v for a,v in zip(c,original))),original


def _symbolic(f,constr,maximize=False):
    f=S.sympify(f);constr=list(constr);symbols=set(f.free_symbols)
    for relation in constr:symbols.update(getattr(relation,'free_symbols',set()))
    symbols=sorted(symbols,key=str);rows=[];rhs=[];equal=[];equal_rhs=[]
    for relation in constr:
        if relation in (True,S.true):continue
        if relation in (False,S.false):raise InfeasibleLPError('\nInconsistent/False constraint')
        name=relation.class_name
        if name not in ('Equality','LessThan','GreaterThan'):raise TypeError('constraints must be non-strict linear equalities or inequalities')
        expr=relation.lhs-relation.rhs
        if name=='GreaterThan':expr=-expr
        matrix,vector=S.linear_eq_to_matrix([expr],symbols)
        (equal if name=='Equality' else rows).append(list(matrix))
        (equal_rhs if name=='Equality' else rhs).append(vector[0])
    matrix,constant=S.linear_eq_to_matrix([f],symbols);cost=list(matrix)
    if maximize:cost=[-a for a in cost]
    try:
        value,solution=linprog(cost,rows,rhs,equal,equal_rhs,[(None,None)]*len(symbols))
    except InfeasibleLPError:
        if _univariate_conflict(constr):raise InfeasibleLPError('\nInconsistent/False constraint') from None
        raise
    value=(-value if maximize else value)-constant[0]
    return S.simplify(value),dict(zip(symbols,solution))

def _univariate_conflict(constr):
    # SymPy reports contradictory single-variable bounds (x <= -1, x >= 0)
    # as an inconsistent constraint before running the simplex.
    bounds={}
    for relation in constr:
        name=getattr(relation,'class_name',None)
        if name not in ('Equality','LessThan','GreaterThan') or len(relation.free_symbols)!=1:continue
        x,=relation.free_symbols
        expr=relation.lhs-relation.rhs;a=expr.coeff(x);rest=S.expand(expr-a*x)
        if rest.has(x) or _sign(a)==0:continue
        value=-rest/a;low,high=bounds.get(x,(None,None))
        upper=name=='Equality' or (name=='LessThan')==(_sign(a)>0)
        lower=name=='Equality' or not upper
        if upper:high=value if high is None or _sign(value-high)<0 else high
        if lower:low=value if low is None or _sign(value-low)>0 else low
        bounds[x]=(low,high)
        if low is not None and high is not None and _sign(high-low)<0:return True
    return False

def lpmin(f,constr):return _symbolic(f,constr)
def lpmax(f,constr):return _symbolic(f,constr,True)

def show_linprog(c,A=None,b=None,A_eq=None,b_eq=None,bounds=None):
    c=list(c);variables=S.symbols('x:'+str(len(c)));conditions=[]
    for matrix,rhs,equality in ((A,b,False),(A_eq,b_eq,True)):
        rows=_matrix_rows(matrix,len(c));rhs=[] if rhs is None else list(rhs)
        for row,value in zip(rows,rhs):
            lhs=sum(a*x for a,x in zip(row,variables));conditions.append(S.Eq(lhs,value) if equality else lhs<=value)
    for x,(lower,upper) in zip(variables,_bounds(bounds,len(c))):
        if lower is not None:conditions.append(x>=lower)
        if upper is not None:conditions.append(x<=upper)
    return sum(S.sympify(a)*x for a,x in zip(c,variables)),conditions

__all__=['InfeasibleLPError','UnboundedLPError','linprog','lpmin','lpmax','show_linprog']
