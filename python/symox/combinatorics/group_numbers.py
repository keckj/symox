"""Exact order criteria and known counts of finite groups."""
from itertools import combinations
from math import gcd,prod
from operator import index
from .. import factorint


def _factors(n):
    n=index(n)
    if n<1:raise ValueError('group order must be a positive integer')
    return {int(p):int(e) for p,e in factorint(n).items()}
def _nilpotent(factors):
    return all(pow(q,k,p)!=1 for p in factors for q,e in factors.items() for k in range(1,e+1))
def is_nilpotent_number(n):return _nilpotent(_factors(n))
def is_abelian_number(n):
    factors=_factors(n);return max(factors.values(),default=0)<3 and _nilpotent(factors)
def is_cyclic_number(n):
    factors=_factors(n);return all(e==1 for e in factors.values()) and _nilpotent(factors)


def groups_count(n):
    factors=_factors(n);n=index(n)
    if n==1:return 1
    if len(factors)==1:
        p,e=next(iter(factors.items()))
        tables={2:(1,1,2,5,14,51,267,2328,56092,10494213,49487367289),3:(1,1,2,5,15,67,504,9310,1396077,5937876645)}
        if p in tables and e<len(tables[p]):return tables[p][e]
        if e<=4:return {1:1,2:2,3:5,4:15}[e]
        if e==5:return 61+2*p+2*gcd(p-1,3)+gcd(p-1,4)
        if e==6:return 3*p*p+39*p+344+24*gcd(p-1,3)+11*gcd(p-1,4)+2*gcd(p-1,5)
        if e==7:
            if p==5:return 34297
            return 3*p**5+12*p**4+44*p**3+170*p*p+707*p+2455+(4*p*p+44*p+291)*gcd(p-1,3)+(p*p+19*p+135)*gcd(p-1,4)+(3*p+31)*gcd(p-1,5)+4*gcd(p-1,7)+5*gcd(p-1,8)+gcd(p-1,9)
    if all(e==1 for e in factors.values()):
        primes=list(factors);total=0
        if len(primes)>24:raise ValueError('squarefree group count exceeds the subset enumeration limit')
        for k in range(len(primes)+1):
            for selected in combinations(primes,k):
                selected=set(selected)
                total+=prod((p**sum(q%p==1 for q in selected)-1)//(p-1) for p in primes if p not in selected)
        return total
    known={12:5,18:5,20:5,24:15,28:4,36:14,40:14,44:4,45:2,48:52,50:5,52:5,54:15,56:13,60:13,63:4,68:5,72:50,75:3,76:4,80:52,84:15,88:12,90:10,92:4}
    if n in known:return known[n]
    raise ValueError('the number of groups of this order is not implemented')

__all__=['groups_count','is_abelian_number','is_cyclic_number','is_nilpotent_number']
