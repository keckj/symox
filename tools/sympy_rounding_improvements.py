"""Record proved differences instead of copying reference rounding bugs."""
import argparse
import json
from pathlib import Path
import sympy as sp
import symox as sx


def audit():
    if sp.__version__!='1.14.0':raise ValueError('Reference version must be 1.14.0')
    n=10**100
    cases=[
        ('sqrt(10**200+1)',n,'N^2 < N^2+1 < (N+1)^2 for N=10^100, hence N < sqrt(N^2+1) < N+1.'),
        ('-sqrt(10**200+1)',-n,'Negating N < sqrt(N^2+1) < N+1 and truncating toward zero gives -N.'),
        ('exp(10**(-1000))',1,'For epsilon=10^-1000, 0 < epsilon < 1/2 < log(2); therefore 1 < exp(epsilon) < 2.'),
        ('exp(-10**(-1000))',0,'The same epsilon gives 1/2 < exp(-epsilon) < 1; truncation toward zero is zero.'),
    ]
    records=[]
    for source,expected,proof in cases:
        actual=int(sx.sympify(source))
        assert actual==expected
        try:reference={'value':str(int(sp.sympify(source)))}
        except Exception as error:reference={'error':type(error).__name__,'message':str(error)}
        records.append({'operation':'int','expression':source,'native':str(actual),
                        'expected':str(expected),'reference':reference,'proof':proof,
                        'scope':'This exact input family; no universal int/evalf certificate.'})
    literal='1_234.567_890_123_456_789'
    native=sx.Float(literal,40)
    correct=sp.Float(literal.replace('_',''),40)
    reference=sp.Float(literal,40)
    assert native.to_sympy()==correct
    records.append({'operation':'Float','expression':literal,'decimal_precision':40,
        'native':str(native),'reference':str(reference),'expected':str(correct),
        'proof':'Underscores separate digits and do not change the decimal exponent. The canonical literal is 1234.567890123456789.',
        'scope':'Decimal separators with explicit precision; no global Float certificate.'})
    return {'reference_version':sp.__version__,'records':records}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    result=audit()
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(f"{len(result['records'])} proved differences recorded")
