"""Printing parity regressions for audit batch print-a1.

Expected strings are hard-coded from SymPy 1.14.
"""

import pytest

import symox as L


CASES = {
    'multiline_latex': (
        "x, y, alpha = L.symbols('x y alpha'); e = L.sin(alpha*y) + L.exp(L.I*alpha) - L.cos(L.log(y)); out = '\\n'.join([L.multiline_latex(x, e), L.multiline_latex(x, e, 2), L.multiline_latex(x, e, terms_per_line=2, environment='eqnarray', use_dots=True), L.multiline_latex(x, e, environment='IEEEeqnarray')])",
        '\\begin{align*}\nx = & e^{i \\alpha} \\\\\n& + \\sin{\\left(\\alpha y \\right)} \\\\\n& - \\cos{\\left(\\log{\\left(y \\right)} \\right)} \n\\end{align*}\n\\begin{align*}\nx = & e^{i \\alpha} + \\sin{\\left(\\alpha y \\right)} \\\\\n& - \\cos{\\left(\\log{\\left(y \\right)} \\right)} \n\\end{align*}\n\\begin{eqnarray}\nx & = & e^{i \\alpha} + \\sin{\\left(\\alpha y \\right)} \\dots\\nonumber\\\\\n& & - \\cos{\\left(\\log{\\left(y \\right)} \\right)} \n\\end{eqnarray}\n\\begin{IEEEeqnarray}{rCl}\nx & = & e^{i \\alpha} \\nonumber\\\\\n& & + \\sin{\\left(\\alpha y \\right)} \\nonumber\\\\\n& & - \\cos{\\left(\\log{\\left(y \\right)} \\right)} \n\\end{IEEEeqnarray}',
    ),
    'print_mathml': (
        "import io, contextlib; x = L.Symbol('x'); buf = io.StringIO()\nwith contextlib.redirect_stdout(buf):\n    L.print_mathml(x + 1); L.print_mathml(x + 1, printer='presentation')\nout = buf.getvalue()",
        '<apply>\n\t<plus/>\n\t<ci>x</ci>\n\t<cn>1</cn>\n</apply>\n\n<mrow>\n\t<mi>x</mi>\n\t<mo>+</mo>\n\t<mn>1</mn>\n</mrow>\n\n',
    ),
    'mathematica_code': (
        "x = L.Symbol('x'); out = L.mathematica_code(L.sin(x).series(x).removeO()) + ' ' + L.mathematica_code(-x/3 + x**2*L.Rational(2, 5))",
        '(1/120)*x^5 - 1/6*x^3 + x (2/5)*x^2 - 1/3*x',
    ),
    'ordinal': (
        "from importlib import import_module; o = import_module(L.__name__ + '.sets.ordinals'); w = o.omega; out = '%s %s %s' % (o.Ordinal(o.OmegaPower(w + 1, 1), o.OmegaPower(3, 2)), (w + 1) * w, w**2 * 3 + w)",
        'w**(w + 1) + w**3*2 w**2 w**2*3 + w',
    ),
    'rawlines': (
        "from importlib import import_module; m = import_module(L.__name__ + '.utilities.misc'); out = m.rawlines('this\\nis a test  \\n') + '|' + m.rawlines('a\\\\b\\nc') + '|' + m.strlines('this' * 20)",
        "(\n    'this\\n'\n    'is a test  \\n'\n    '\\n'\n)|(\n    'a\\\\b\\n'\n    'c'\n)|(\n'thisthisthisthisthisthisthisthisthisthisthisthisthisthisthisth'\n'isthisthisthisthis'\n)",
    ),
    'codegen_ast': (
        "from importlib import import_module; a = import_module(L.__name__ + '.codegen.ast'); x = L.Symbol('x'); v = a.Variable(x, type=a.real); out = '%r|%s|%r|%r|%r|%s' % (v, a.Variable.deduced(L.Symbol('r', real=True)), a.Type('float64'), a.intc, a.value_const, a.For(x, L.Range(3), [a.aug_assign(x, '+', 1)]))",
        "Variable(Symbol('x'), type=FloatBaseType(String('real')))|Variable(r, type=real)|Type(String('float64'))|IntBaseType(String('intc'))|Attribute(String('value_const'))|For(x, iterable=Range(0, 3, 1), body=CodeBlock(\n    AddAugmentedAssignment(x, 1)\n))",
    ),
    'sumapprox': (
        "from importlib import import_module; ap = import_module(L.__name__ + '.codegen.approximations'); rw = import_module(L.__name__ + '.codegen.rewriting'); x, y, z = L.symbols('x y z'); bounds = {x: (-1, 1), y: (1000, 2000), z: (-10, 3)}; e = 3*(x + y + L.exp(z)); out = ', '.join(str(rw.optimize(e, [ap.SumApprox(bounds, reltol=r)])) for r in (1e-3, 1e-2, 1e-1))",
        '3*(x + y + exp(z)), 3*y + 3*exp(z), 3*y',
    ),
    'reduce_order': (
        "import sys, importlib; h = importlib.import_module(L.__name__ + '.simplify.hyperexpand'); h = sys.modules[L.__name__ + '.simplify.hyperexpand']; out = '%s|%s|%s' % (h.reduce_order(h.Hyper_Function((3, 4), (5, 6)))[0], h.reduce_order(h.Hyper_Function((1, 2), (3, 4)))[0], h.reduce_order_meijer(h.G_Function([3, 4], [5, 6], [3, 4], [1, 2]))[0])",
        'Hyper_Function((3, 4), (5, 6))|Hyper_Function((1, 2), (3, 4))|G_Function((4, 3), (5, 6), (3, 4), (2, 1))',
    ),
    'hermite_prob': (
        "n, x = L.symbols('n x'); k = L.Symbol('k', integer=False); out = '%s %s %s %s %s %s %s' % (L.hermite_prob(n, -x), L.hermite_prob(n, 1 - x), L.hermite_prob(n, 0), L.hermite_prob(n, L.oo), L.hermite_prob(-1, x), L.hermite_prob(k, x), L.hermite_prob(3, x))",
        '(-1)**n*hermite_prob(n, x) (-1)**n*hermite_prob(n, x - 1) sqrt(pi)/gamma(1/2 - n/2) oo hermite_prob(-1, x) hermite_prob(k, x) x**3 - 3*x',
    ),
    'options': (
        "from importlib import import_module; po = import_module(L.__name__ + '.polys.polyoptions'); x, y, z = L.symbols('x y z'); out = '%s %s %s %s' % (po.Options((x, y, z), {'domain': 'ZZ'}), po.build_options((x,), {'domain': 'QQ', 'expand': False}), po.Options((x,), {'extension': [L.sqrt(2)]}), po.Options((), {'gens': (x,), 'field': True}))",
        "{'domain': ZZ, 'gens': (x, y, z), 'auto': False} {'domain': QQ, 'expand': False, 'gens': (x,), 'auto': False} {'extension': {sqrt(2)}, 'gens': (x,), 'domain': QQ<sqrt(2)>, 'auto': False} {'gens': (x,), 'field': True, 'auto': False}",
    ),
    'decompose': (
        "from importlib import import_module; ip = import_module(L.__name__ + '.integrals.intpoly'); x, y = L.symbols('x y'); out = '%s %s %s %s' % (ip.decompose(x**2 + x*y + x + y + x**3*y**2 + y**5), ip.decompose(2*x + 3), ip.decompose(x*(x + 1)), ip.decompose((x + y)**2))",
        '{1: x + y, 2: x**2 + x*y, 5: x**3*y**2 + y**5} {0: 3, 1: 2*x} {x + 1: x*(x + 1)} {2: (x + y)**2}',
    ),
    'best_origin': (
        "from importlib import import_module; ip = import_module(L.__name__ + '.integrals.intpoly'); x, y = L.symbols('x y'); l = L.Segment2D(L.Point(0, 3), L.Point(1, 1)); out = '%s %s %s' % (ip.best_origin((2, 1), 3, l, x**3*y**7), ip.best_origin((2, 0), 3, l, x + y), ip.best_origin((0, 1), 3, l, 5))",
        '(0, 3.0) (0, 3) (0, 3.0)',
    ),
    'disambiguate': (
        "from importlib import import_module; sm = import_module(L.__name__ + '.core.symbol'); y = L.Symbol('y'); out = '%s %s' % (sm.disambiguate(L.Symbol('_x'), L.Dummy('x'), L.Dummy('x')), sm.disambiguate(L.Symbol('x')/y, L.Dummy('x')/y))",
        '(x_2, x, x_1) (x_1/y, x/y)',
    ),
    'propkb': (
        "from importlib import import_module; inf = import_module(L.__name__ + '.logic.inference'); x, y, z = L.symbols('x y z'); kb = inf.PropKB(); kb.tell(x | y); kb.tell(y); r1 = str(kb.clauses); kb.retract(y); kb.tell((x >> y) & (z | ~x)); out = '%s %s %s' % (r1, kb.clauses, kb.ask(y))",
        '[y, x | y] [x | y, y | ~x, z | ~x] True',
    ),
    'kernS': (
        "from importlib import import_module; sp = import_module(L.__name__ + '.core.sympify'); out = ' '.join(L.srepr(sp.kernS(s)) for s in ['-(x + 1)', '-(x + 1)*y', 'x - (x + 1)', '2*(x + y)', '(1 - x)**2'])",
        "Mul(Integer(-1), Add(Symbol('x'), Integer(1))) Mul(Integer(-1), Symbol('y'), Add(Symbol('x'), Integer(1))) Add(Symbol('x'), Mul(Integer(-1), Add(Symbol('x'), Integer(1)))) Mul(Integer(2), Add(Symbol('x'), Symbol('y'))) Pow(Add(Integer(1), Mul(Integer(-1), Symbol('x'))), Integer(2))",
    ),
    'collect_abs': (
        "from importlib import import_module; rs = import_module(L.__name__ + '.simplify.radsimp'); x, y, z = L.symbols('x y z'); out = ' '.join(L.srepr(rs.collect_abs(e)) for e in [abs(1/x), abs(x)*abs(y)*z, abs(x)/abs(y)**L.Rational(1, 2), abs(x + 1)/abs(x**2 - 1)])",
        "Abs(Pow(Symbol('x'), Integer(-1))) Mul(Symbol('z'), Abs(Mul(Symbol('x'), Symbol('y')))) Abs(Mul(Symbol('x'), Pow(Symbol('y'), Rational(-1, 2)))) Abs(Mul(Add(Symbol('x'), Integer(1)), Pow(Add(Pow(Symbol('x'), Integer(2)), Integer(-1)), Integer(-1))))",
    ),
    'collect_sqrt': (
        "from importlib import import_module; rs = import_module(L.__name__ + '.simplify.radsimp'); a, b = L.symbols('a b'); r2, r3, r5 = [L.sqrt(i) for i in (2, 3, 5)]; out = '%s %s %s' % (rs.collect_sqrt(a*r2 + b*r2 + a*r3 + b*r5, evaluate=False), rs.collect_sqrt(a + b, evaluate=False), rs.collect_sqrt(2*r2*a + 4*r2*b))",
        '((sqrt(3)*a, sqrt(5)*b, sqrt(2)*(a + b)), 3) ((a + b,), 0) 2*sqrt(2)*(a + 2*b)',
    ),
    'function_range': (
        "from importlib import import_module; cu = import_module(L.__name__ + '.calculus.util'); x = L.Symbol('x'); out = '%s %s %s %s' % (cu.function_range(L.log(x), x, L.S.Reals), cu.function_range(x**3, x, L.S.Reals), cu.function_range(x, x, L.S.Reals), cu.function_range(L.sin(x), x, L.S.Reals))",
        'Interval(-oo, oo) Interval(-oo, oo) Interval(-oo, oo) Interval(-1, 1)',
    ),
    'type_paths': (
        "out = '%s %s' % (type(L.Tuple(1, 2)), type(L.AlgebraicNumber(L.sqrt(2))))",
        "<class 'sympy.core.containers.Tuple'> <class 'sympy.core.numbers.AlgebraicNumber'>",
    ),
    'eqs_to_ring': (
        "from importlib import import_module; ps = import_module(L.__name__ + '.polys.solvers'); a, x, y = L.symbols('a x y'); eqs, ring = ps.sympy_eqs_to_ring([x + 2*y + a, 2*x - y], [x, y]); out = '%s %s %s' % (eqs, ring, L.sring(x + a, x, field=True)[0])",
        '[x + 2*y + a, 2*x - y] ZZ(a)[x,y] Polynomial ring in x over ZZ(a) with lex order',
    ),
    'cse_release': (
        "from importlib import import_module; cm = import_module(L.__name__ + '.simplify.cse_main'); x, y = L.symbols('x y'); eqs = [(x + y - 1)**2, x, x + y, (x + y)/(2*x + 1) + (x + y - 1)**2, (2*x + 1)**(x + y)]; out = str(cm.cse_release_variables(*L.cse(eqs)))",
        '([(x0, x + y), (x1, (x0 - 1)**2), (x2, 2*x + 1), (_3, x0/x2 + x1), (_4, x2**x0), (x2, None), (_0, x1), (x1, None), (_2, x0), (x0, None), (_1, x)], (_0, _1, _2, _3, _4))',
    ),
    'conjugate': (
        "a, b, x, y, z, n, m = L.symbols('a b x y z n m'); out = ' '.join(str(L.conjugate(e)) for e in [L.beta(x, y), L.betainc(a, b, x, y), L.betainc_regularized(a, b, x, y), L.erfi(z), L.erf2(x, y), L.li(z), L.fresnels(z), L.airyai(z), L.airybiprime(z), L.lowergamma(a, z), L.uppergamma(a, z), L.jacobi(n, a, b, x), L.gegenbauer(n, a, x), L.legendre(n, x), L.assoc_legendre(n, m, x), L.hermite(n, x), L.hermite_prob(n, x), L.chebyshevt(n, x), L.laguerre(n, x), L.assoc_laguerre(n, a, x), L.sech(z), L.csch(z), L.li(-L.Symbol('p', positive=True)), L.lowergamma(a, 0)])",
        'beta(conjugate(x), conjugate(y)) betainc(conjugate(a), conjugate(b), conjugate(x), conjugate(y)) betainc_regularized(conjugate(a), conjugate(b), conjugate(x), conjugate(y)) erfi(conjugate(z)) erf2(conjugate(x), conjugate(y)) li(conjugate(z)) fresnels(conjugate(z)) airyai(conjugate(z)) airybiprime(conjugate(z)) lowergamma(conjugate(a), conjugate(z)) uppergamma(conjugate(a), conjugate(z)) jacobi(n, conjugate(a), conjugate(b), conjugate(x)) gegenbauer(n, conjugate(a), conjugate(x)) legendre(n, conjugate(x)) assoc_legendre(n, conjugate(m), conjugate(x)) hermite(n, conjugate(x)) hermite_prob(n, conjugate(x)) chebyshevt(n, conjugate(x)) laguerre(n, conjugate(x)) assoc_laguerre(n, conjugate(a), conjugate(x)) sech(conjugate(z)) csch(conjugate(z)) conjugate(li(-p)) 0',
    ),
    'ground_roots': (
        "x = L.Symbol('x'); out = '%s %s' % (L.ground_roots(x**6 - 4*x**4 + 4*x**3 - x**2), L.Poly(x**6 - 4*x**4 + 4*x**3 - x**2).ground_roots())",
        '{1: 2, 0: 2} {1: 2, 0: 2}',
    ),
    'pow_of_pow': (
        "s, x = L.symbols('s x'); k = L.Symbol('k', integer=True); h = L.Rational(1, 2)**s; out = '%s %s %s %s %s %s' % (1/h, h**2, (x**s)**k, (L.Rational(2, 3)**s)**-1, L.sqrt(h), h**k)",
        '2**s 2**(-2*s) x**(k*s) (2/3)**(-s) sqrt((1/2)**s) 2**(-k*s)',
    ),
    'lerchphi': (
        "z, s = L.symbols('z s'); out = str(L.expand_func(L.lerchphi(z, s, L.Rational(3, 2))))",
        '-2**s/z + 2**(s - 1)*(polylog(s, sqrt(z))/sqrt(z) - polylog(s, sqrt(z)*exp_polar(I*pi))/sqrt(z))/z',
    ),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_print_a1_parity(name):
    code, expected = CASES[name]
    ns = {"L": L}
    exec(code, ns)
    assert ns["out"].replace("symox", "sympy") == expected
