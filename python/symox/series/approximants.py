"""Consecutive Pade approximants from a native formal continued fraction."""
import itertools


def approximants(l,X=None,simplify=False):
    from .. import S,Symbol,sympify,cancel,simplify as simp
    X=Symbol('x') if X is None else sympify(X)
    values=list(itertools.islice(iter(l),257))
    if len(values)>256:
        raise ValueError('Pade approximation accepts at most 256 coefficients')
    values=list(map(sympify,values))
    p1,q1,p2,q2=[S.One],[S.Zero],[S.Zero],[S.One]
    while values:
        first=next((i for i,value in enumerate(values) if value!=0),None)
        if first is None:
            return
        leading=values[first]
        reciprocal=[cancel(1/leading)]
        for degree in range(first+1,len(values)):
            coefficient=-sum((values[j+1]*reciprocal[first-j-1]
                              for j in range(first,degree)),S.Zero)
            reciprocal.append(cancel(coefficient/leading))
        values=reciprocal
        factor,values[0]=values[0],S.Zero
        p=[S.Zero]*max(len(p2),first+len(p1))
        q=[S.Zero]*max(len(q2),first+len(q1))
        for k,value in enumerate(p2):p[k]=factor*value
        for k,value in enumerate(q2):q[k]=factor*value
        for k,value in enumerate(p1):p[k+first]+=value
        for k,value in enumerate(q1):q[k+first]+=value
        p,q=[cancel(value) for value in p],[cancel(value) for value in q]
        while len(p)>1 and p[-1]==0:p.pop()
        while len(q)>1 and q[-1]==0:q.pop()
        p1,p2,q1,q2=p2,p,q2,q
        numerator=sum((value*X**k for k,value in enumerate(p)),S.Zero)
        denominator=sum((value*X**k for k,value in enumerate(q)),S.Zero)
        if denominator==0:
            return
        result=cancel(numerator/denominator)
        yield simp(result) if simplify else result


__all__=['approximants']
