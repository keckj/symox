def more2(prde, P, DEx, Poly, x, t):
    D1 = lambda: DEx(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    yield 'bd_prim_eq', lambda: prde.param_rischDE(P('-1', t), P('x', t), [(P('1', t), P('1', t)), (P('t', t), P('1', t))], D1())
    yield 'bd_prim_eq2', lambda: prde.param_rischDE(P('-2', t), P('x', t), [(P('t', t), P('1', t)), (P('x*t**2', t), P('1', t))], D1())
    yield 'bd_prim_eq3', lambda: prde.param_rischDE(P('-x - 1', t), P('x', t), [(P('t', t), P('1', t)), (P('1', t), P('x', t))], D1())
    yield 'li_prim_eq', lambda: prde.limited_integrate(P('t', t), P('1', t), [(P('1', t), P('x', t))], D1())
