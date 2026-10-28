"""Independent integer/Fraction rounding; no string-to-number conversion."""
from fractions import Fraction as F
from pathlib import Path
import json
import symox as sx
import sympy as sp


def dyad(value):
    if value.is_Float:
        sign,man,exp,_=value._mpf_
        return F((-1 if sign else 1)*man)*(F(2)**exp)
    assert value.is_Rational,repr(value)
    return F(int(value.p),int(value.q))


def components(value,mod):
    return dyad(mod.re(value)),dyad(mod.im(value))


def raw_components(value,mod):
    if mod is sp:
        def convert(t):
            sign,man,exp,_=t
            return F((-1 if sign else 1)*man)*(F(2)**exp)
        return tuple(convert(t) for t in value._mpc_)
    return components(value.as_expr(),mod)


def encode(q):return [q.numerator,q.denominator]
def pair(q):return [encode(x) for x in q]

def quantum(q,p):
    assert q
    a,b=abs(q.numerator),q.denominator
    e=a.bit_length()-b.bit_length()
    if (a < b*(1<<e)) if e>=0 else (a*(1<<-e)<b):e-=1
    return F(2)**(e-p+1)


def rounded(q,p):
    if not q:return F(0)
    scale=quantum(q,p);x=abs(q)/scale
    n,r=divmod(x.numerator,x.denominator)
    if 2*r>x.denominator or 2*r==x.denominator and n%2:n+=1
    return (-n if q<0 else n)*scale


def exact_op(a,b,op):
    ar,ai=a;br,bi=b
    if op=='/':
        d=br*br+bi*bi
        return (ar*br+ai*bi)/d,(ai*br-ar*bi)/d
    if op=='*':return ar*br-ai*bi,ar*bi+ai*br
    if op=='-':return ar-br,ai-bi
    if op=='+':return ar+br,ai+bi
    raise ValueError(op)


def trace(rows,k,mod,p):
    events=[]
    def operation(a,b,op,label):
        aa=raw_components(a,mod);bb=raw_components(b,mod)
        assert aa==components(k.to_sympy(a),mod) and bb==components(k.to_sympy(b),mod)
        exact=exact_op(aa,bb,op);expected=tuple(rounded(x,p) for x in exact)
        result={'/':lambda:a/b,'*':lambda:a*b,'-':lambda:a-b}[op]()
        actual=raw_components(result,mod)
        assert actual==components(k.to_sympy(result),mod)
        events.append(dict(label=label,op=op,a=pair(aa),b=pair(bb),exact=pair(exact),
                           rounded=pair(expected),actual=pair(actual),correct=actual==expected,
                           error_ulps=[encode((a-e)/quantum(e,p)) if e else None for a,e in zip(actual,expected)]))
        return result
    r=0;pivots=[]
    for col in range(len(rows[0])):
        pivot=next((i for i in range(r,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        inv=operation(k.one,rows[r][col],'/',f'pivot{col}.reciprocal')
        rows[r]=[operation(v,inv,'*',f'pivot{col}.normalize[{j}]') for j,v in enumerate(rows[r])]
        for i in range(len(rows)):
            if i==r or not rows[i][col]:continue
            factor=rows[i][col]
            rows[i]=[operation(a,operation(factor,b,'*',f'pivot{col}.eliminate[{i},{j}].product'),'-',f'pivot{col}.eliminate[{i},{j}].difference') for j,(a,b) in enumerate(zip(rows[i],rows[r]))]
        pivots.append(col);r+=1
        if r==len(rows):break
    return dict(events=events,pivots=pivots)


def run():
    report=[]
    for p in (140,211):
        k=sx.ComplexField(p);rk=sp.ComplexField(p)
        raw=[[(2,1),(3,-2),(5,3)],[(7,-1),(11,2),(17,1)]]
        rows=[[k.convert(sx.Rational(a,7)+sx.I*sx.Rational(b,11)) for a,b in row] for row in raw]
        stored=[[components(k.to_sympy(v),sx) for v in row] for row in rows]
        # Explicitly construct reference input from native stored dyades.
        rr=[[rk.convert(sp.Rational(a.numerator,a.denominator)+sp.I*sp.Rational(b.numerator,b.denominator)) for a,b in row] for row in stored]
        assert stored==[[components(rk.to_sympy(v),sp) for v in row] for row in rr]
        a=trace(rows,k,sx,p);b=trace(rr,rk,sp,p)
        assert [e['label'] for e in a['events']]==[e['label'] for e in b['events']]
        divergences=[dict(native=x,reference=y) for x,y in zip(a['events'],b['events']) if x['actual']!=y['actual']]
        record=dict(precision=p,stored_input=[[pair(v) for v in row] for row in stored],native=a,reference=b,
                    first_divergence=divergences[0] if divergences else None,
                    native_incorrect=[e for e in a['events'] if not e['correct']],
                    reference_incorrect=[e for e in b['events'] if not e['correct']])
        report.append(record)
    return report

if __name__=='__main__':
    report=run()
    print(json.dumps(report,indent=2))
