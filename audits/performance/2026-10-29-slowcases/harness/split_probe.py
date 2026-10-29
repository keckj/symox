import os, sys, time
lib=os.environ.get("LIB","sympy")
exec("from %s import *" % lib)
exec("from %s.integrals.rationaltools import ratint_logpart, ratint_ratpart, log_to_real" % lib)
x=symbols('x'); t=Dummy('t')
k=int(sys.argv[1]); f=1/(x**5+x+1)**k
T0=time.perf_counter()
p,q=fraction(together(f)); P,Q=Poly(p,x),Poly(q,x)
g,h=ratint_ratpart(P,Q,x); P2,Q2=h.as_numer_denom()
L=ratint_logpart(Poly(P2,x),Poly(Q2,x),x,t)
res=g
for hh,qq in L:
    print('deg q', qq.degree(), 'factors', [ff.degree() for ff,_ in qq.factor_list()[1]], flush=True)
    for ff,_ in qq.factor_list()[1]:
        hr=Poly({m: Poly(c, t).rem(ff).as_expr() for m,c in hh.as_poly(x).as_dict().items()}, x) if False else None
        # reduce t-coefficients of hh modulo ff
        hx=Poly(hh.as_expr(), x)
        coeffs=[Poly(c,t).rem(ff).as_expr() for c in hx.all_coeffs()]
        hred=Poly(sum(c*x**i for i,c in enumerate(reversed(coeffs))), x, t)
        T=time.perf_counter(); r=log_to_real(hred, ff, x, t); print(' factor deg',ff.degree(),'lrt %.2f'%(time.perf_counter()-T), None if r is None else len(str(r)), flush=True)
        res+= r if r is not None else RootSum(ff, Lambda(t, t*log(hred.as_expr())))
print('total %.2f len %d'%(time.perf_counter()-T0, len(str(res))))
d=diff(res,x)-f
for pt in [Rational(1,3), Rational(-5,7), 2]:
    print(pt, N(d.subs(x,pt),30))
