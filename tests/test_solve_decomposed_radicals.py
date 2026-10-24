"""Complete quadratic preimages without polar expansion of quartic radicals."""
import json
import subprocess
import sys

import sympy as sp


def test_composed_quartic_roots_and_checked_system_complete():
    # A process boundary makes the former native expansion hang a bounded
    # regression failure; Python signal handlers cannot interrupt Rust here.
    program = '''
import json
import symox as s
x,y=s.symbols('x y')
p=x**8-2*x**4+x**2-3
out={}
for name,f in [('octic',p),('repeated',(x**4-2*x**2-3)**2),('quadratic_shift',((x**2+1)**4-2*(x**2+1)**2+(x**2+1)-3))]:
    roots=s.roots(f,x)
    out[name]=[[str(r),int(m)] for r,m in roots.items()]
out['system']=[[str(a),str(b)] for a,b in s.solve([x*x+y*y-4,y-x**4+1],[x,y])]
x,y=s.symbols('x y',real=True)
out['real_system']=[[str(a),str(b)] for a,b in s.solve([x*x+y*y-4,y-x**4+1],[x,y])]
x=s.Symbol('x',positive=True)
out['positive_system']=[[str(a),str(b)] for a,b in s.solve([x*x+y*y-4,y-x**4+1],[x,y])]
print(json.dumps(out))
'''
    child = subprocess.run([sys.executable,'-c',program],capture_output=True,text=True,timeout=60)
    assert child.returncode == 0, child.stderr
    out = json.loads(child.stdout)
    x = sp.Symbol('x')
    p = x**8-2*x**4+x**2-3
    shifted = (x**2+1)**4-2*(x**2+1)**2+(x**2+1)-3
    for name, polynomial, multiplicity in [('octic',p,1),('repeated',x**4-2*x**2-3,2),('quadratic_shift',shifted,1)]:
        roots = [(sp.sympify(r),m) for r,m in out[name]]
        assert len(roots) == sp.degree(polynomial,x) and all(m == multiplicity for _,m in roots),name
        assert not any(r.has(sp.Float) for r,_ in roots)
        expected = sp.nroots(polynomial,n=65,maxsteps=300)
        numeric = [r.evalf(70) for r,_ in roots]
        # Compare every reference root with a distinct returned root, in
        # addition to residuals: residual tests alone miss missing solutions.
        unmatched = numeric.copy()
        for ref in expected:
            best = min(range(len(unmatched)),key=lambda i:abs(complex((unmatched[i]-ref).evalf(60))))
            assert abs(complex((unmatched.pop(best)-ref).evalf(60))) < 1e-48
        assert not unmatched
        assert all(abs(complex(polynomial.subs(x,r).evalf(65))) < 1e-48 for r,_ in roots)
    for name,count in [('system',8),('real_system',2),('positive_system',1)]:
        rows = [(sp.sympify(a),sp.sympify(b)) for a,b in out[name]]
        assert len(rows) == count
        values = [a.evalf(65) for a,_ in rows]
        assert all(abs(complex(a-b)) > 1e-20 for i,a in enumerate(values) for b in values[i+1:])
        for a,b in rows:
            assert abs(complex((a*a+b*b-4).evalf(65))) < 1e-48
            assert abs(complex((b-a**4+1).evalf(65))) < 1e-48
        if name in ('real_system','positive_system'):
            assert all(abs(complex(v.evalf(65)).imag) < 1e-48 for row in rows for v in row)
        if name=='positive_system':
            assert float(sp.re(rows[0][0].evalf(65))) > 0
