"""Printing semantics, real tool execution, XML and native dependency checks."""
import builtins
import ctypes
import io
import json
import math
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import pytest
import symox as sx
import sympy as sp
import symox.printing as p

x,y=sx.symbols('x y')


@pytest.mark.parametrize('name',['ccode','cxxcode','fcode','pycode','rust_code','jscode','julia_code','octave_code','rcode','maple_code','glsl_code'])
def test_code_entrypoints(name):
    text=getattr(p,name)(x*x+sx.sin(y),assign_to='result')
    assert isinstance(text,str) and 'result' in text and 'x' in text and 'sin' in text.lower()


def test_mathematica_code_has_no_assign_to():
    # SymPy 1.14: mathematica_code(expr, **settings) has no assign_to.
    assert p.mathematica_code(x*x+sx.sin(y))=='x^2 + Sin[y]'
    with pytest.raises(TypeError,match="Unknown setting 'assign_to'"):
        p.mathematica_code(x*x+sx.sin(y),assign_to='result')


@pytest.mark.parametrize('expr',[x+1,x*y,sx.Rational(2,3),sx.sqrt(x),sx.sin(x),sx.Eq(x,1),sx.pi,sx.E,sx.I,
                                sx.Matrix([[x,1],[0,y]]),sx.Integral(x,(x,0,1)),sx.Derivative(sx.sin(x),x,evaluate=False)])
def test_content_mathml_against_reference(expr):
    def canonical(text):return ET.tostring(ET.fromstring(text),encoding='unicode')
    assert canonical(p.mathml(expr))==canonical(sp.mathml(expr.to_sympy()))


@pytest.mark.parametrize('mode',['content','presentation'])
def test_mathml_escape_and_piecewise(mode):
    expression=sx.Symbol('a<&"')+sx.Piecewise((x,x>0),(0,True))
    assert ET.fromstring(p.mathml(expression,printer=mode)) is not None
    assert '&lt;' in p.mathml(expression,printer=mode)
    with pytest.raises(ValueError):p.mathml(x,printer='bad')


@pytest.mark.skipif(shutil.which('node') is None,reason='Node.js unavailable')
@pytest.mark.parametrize('expression',[sx.sin(x)+sx.cos(x),x*x+sx.Rational(2,3),sx.sqrt(x+2),
    sx.Piecewise((x*x,x>0),(-x,True)),sx.exp(x)/(1+x*x),sx.log(x+3),sx.Abs(x),sx.Max(x,2)])
def test_javascript_executes(expression):
    values=[-1.5,0.0,1.25,3.0]
    source='const values='+json.dumps(values)+'; console.log(JSON.stringify(values.map(x => '+p.jscode(expression)+')));'
    actual=json.loads(subprocess.check_output(['node','-e',source],text=True))
    expected=[float(expression.subs(x,value)) for value in values]
    assert actual==pytest.approx(expected,rel=2e-13,abs=2e-13)


def test_unknown_functions_and_user_mapping():
    expression=sx.Function('custom')(x)
    # The reference Maple printer allows unknown functions and has no
    # user_functions setting.
    assert p.maple_code(expression)=='custom(x)'
    with pytest.raises(TypeError):p.maple_code(expression,user_functions={'custom':'mapped'})
    for fn in [p.jscode,p.julia_code,p.octave_code,p.rcode,p.glsl_code]:
        with pytest.raises(NotImplementedError):fn(expression)
        assert 'mapped' in fn(expression,user_functions={'custom':'mapped'})
    # SymPy's MCodePrinter prints an unknown function as ``name[args]``.
    assert p.mathematica_code(expression)=='custom[x]'
    assert p.mathematica_code(expression,user_functions={'custom':'mapped'})=='mapped[x]'
    # SymPy 1.14 jscode does not rename reserved words.
    assert p.jscode(sx.Symbol('for'))=='for'
    # The reference JavaScript printer has no reserved-word list.
    assert p.jscode(sx.Symbol('for'),error_on_reserved=True)=='for'
    assert p.jscode(x,human=False)==(set(),set(),'x')


def test_symbolic_python_program_and_representations():
    expression=sx.Function('f')(sx.Symbol('if'))+x*x
    namespace={}
    exec(p.python(expression),namespace)
    assert namespace['e']==expression
    assert p.srepr(x+1)==sx.srepr(x+1)
    assert p.sstrrepr('hello')=="'hello'"
    assert p.srepr([x,1])==sp.srepr([sp.Symbol('x'),1])


@pytest.mark.parametrize('headings',[None,'automatic',[['a','b'],['c','d']]])
def test_table_form(headings):
    data=[[1,2],[0,4]]
    options={} if headings is None else {'headings':headings}
    actual=p.TableForm(data,**options)
    expected=sp.TableForm([[1,2],[0,4]],**options)
    assert str(actual)==str(expected)
    assert p.latex(actual)==sp.latex(expected)
    assert actual.as_matrix()==sx.Matrix(data)
    assert data==[[1,2],[0,4]]


def test_pretty_unicode_and_console_wrappers(capsys,monkeypatch):
    previous=p.pprint_use_unicode(False)
    try:
        assert all(ord(char)<128 for char in p.pretty(sx.sqrt(x),use_unicode=False))
        assert isinstance(p.pretty(sx.sqrt(x)),str)
    finally:p.pprint_use_unicode(previous)
    for printer,formatter in [('print_ccode','ccode'),('print_fcode','fcode'),('print_glsl','glsl_code'),
        ('print_jscode','jscode'),('print_latex','latex'),('print_maple_code','maple_code'),
        ('print_mathml','mathml'),('print_python','python'),('print_rcode','rcode')]:
        getattr(p,printer)(x+1)
        if printer=='print_mathml':
            # SymPy 1.14 prints pretty XML (toprettyxml) here.
            assert capsys.readouterr().out=='<apply>\n\t<plus/>\n\t<ci>x</ci>\n\t<cn>1</cn>\n</apply>\n\n'
            continue
        assert capsys.readouterr().out==getattr(p,formatter)(x+1)+'\n'
    p.print_tree(x+1,assumptions=False)
    assert 'Add:' in capsys.readouterr().out
    messages=[];monkeypatch.setattr('pydoc.pager',messages.append)
    p.pager_print(x+1)
    assert messages==[p.pretty(x+1)]


@pytest.mark.skipif(shutil.which('dot') is None,reason='Graphviz unavailable')
def test_dot_graph_and_depth():
    expression=sx.Symbol('a"b')+sx.Symbol('α')+x*x
    graph=p.dotprint(expression)
    result=subprocess.run(['dot','-Tsvg'],input=graph,text=True,capture_output=True,check=True)
    assert '<svg' in result.stdout and 'α' in result.stdout
    assert ' -> ' not in p.dotprint(expression,maxdepth=0)
    assert p.dotprint(x+1)==sp.dotprint(sp.Symbol('x')+1)


@pytest.mark.skipif(shutil.which('latex') is None or shutil.which('dvipng') is None,reason='TeX converters unavailable')
def test_preview_real_png_and_tex(tmp_path):
    result=io.BytesIO()
    tex=tmp_path/'formula.tex'
    p.preview(sx.sqrt(x)+1,viewer='BytesIO',outputbuffer=result,euler=False,outputTexFile=tex)
    assert result.getvalue().startswith(b'\x89PNG\r\n\x1a\n')
    assert tex.exists() and '\\sqrt' in tex.read_text()


def test_gtk_headless_and_multiline_latex(monkeypatch):
    import importlib
    module=importlib.import_module('symox.printing.gtk')
    assert module.print_gtk(x+1,start_viewer=False) is None
    monkeypatch.setattr(module.shutil,'which',lambda exe:None)
    with pytest.raises(ImportError,match='mathmlviewer'):module.print_gtk(x+1)
    text=p.multiline_latex(y,x*x+x+1,terms_per_line=1)
    assert text.startswith('\\begin{align*}') and text.endswith('\\end{align*}')
    assert text.count('\\\\')==2
    with pytest.raises(ValueError):p.multiline_latex(y,x,environment='bad')
    # SymPy 1.14 accepts terms_per_line=0.
    assert p.multiline_latex(y,x,terms_per_line=0)=='\\begin{align*}\ny = & x \n\\end{align*}'


def test_smtlib_executes_in_z3():
    z3=pytest.importorskip('z3')
    i=sx.Symbol('i',integer=True)
    text=p.smtlib_code([i>2,i<2])
    solver=z3.Solver();solver.add(z3.parse_smt2_string(text));assert solver.check()==z3.unsat
    text=p.smtlib_code([x>sx.Rational(1,3),x<2])
    solver=z3.Solver();solver.add(z3.parse_smt2_string(text));assert solver.check()==z3.sat
    text=p.smtlib_code(sx.Eq(sx.Piecewise((x,x>0),(0,True)),3))
    solver=z3.Solver();solver.add(z3.parse_smt2_string(text));assert solver.check()==z3.sat
    assert 'declare' not in p.smtlib_code(x+1,auto_declare=False,auto_assert=False)


def test_native_printing_without_sympy(monkeypatch):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='sympy' or name.startswith('sympy.'):raise AssertionError('imported SymPy')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    assert p.jscode(x*x)=='Math.pow(x, 2)'
    assert p.mathml(x+1).startswith('<apply>')
    assert p.smtlib_code(x>0).endswith('(assert (> x 0))')
    assert '1' in str(p.TableForm([[1,2]]))


def test_target_specific_functions_and_ambiguous_symbols():
    assert p.glsl_code(sx.atan2(x,y))=='atan(x, y)'
    assert p.julia_code(sx.atan2(x,y))=='atan2(x, y)'
    assert p.octave_code(sx.loggamma(x))=='gammaln(x)'
    assert p.maple_code(sx.Eq(x,y))=='x = y'
    # SymPy 1.14 prints Max through its Piecewise rewrite in GLSL.
    assert p.glsl_code(sx.Max(x,y,2))=='(((x <= 2 && y <= 2) ? (\n   2\n)\n: ((x >= y) ? (\n   x\n)\n: (\n   y\n))))'
    a=sx.Symbol('a',positive=True);b=sx.Symbol('a',negative=True)
    # SymPy 1.14 jscode does not check symbol ambiguity.
    assert p.jscode(a+b)=='a + a'
    with pytest.raises(ValueError,match='Ambiguous'):p.smtlib_code(a+b)


@pytest.mark.parametrize('language,compiler',[('fortran','gfortran'),('rust','rustc')])
def test_compiled_native_targets(language,compiler,tmp_path):
    executable=shutil.which(compiler)
    if executable is None:pytest.skip(compiler+' unavailable')
    expression=x*x+sx.sin(x)
    if language=='fortran':
        source='program demo\nimplicit none\nreal(kind=8) :: x,result\nx=0.75d0\n'+p.fcode(expression,assign_to='result',standard=95,source_format='free')+'\nprint *,result\nend program demo\n'
        path=tmp_path/'demo.f90'
    else:
        source='fn main() {let x:f64=0.75; println!("{:.17}", '+p.rust_code(expression)+');}\n'
        path=tmp_path/'demo.rs'
    path.write_text(source)
    output=tmp_path/('demo.exe' if sys.platform=='win32' else 'demo')
    subprocess.run([executable,str(path),'-o',str(output)],check=True,capture_output=True)
    value=float(subprocess.check_output([str(output)],text=True).strip().replace('D','E'))
    assert value==pytest.approx(0.75**2+math.sin(0.75),rel=1e-14)


# Expected strings are hard-coded from SymPy 1.14 mathematica_code.
MATHEMATICA_CASES = [
    ('f(x)',
     'f[x]'),
    ('f(x, y)',
     'f[x, y]'),
    ('g(x)*f(x + 1)',
     'f[x + 1]*g[x]'),
    ('f(x).diff(x, 3)',
     'Hold[D[f[x], {x, 3}]]'),
    ('f(x) + sin(x)',
     'f[x] + Sin[x]'),
    ('Derivative(f(x, y), (x, 2), y)',
     'Hold[D[f[x, y], {x, 2}, y]]'),
    ('sin(x) + cos(x) + tan(x) + cot(x) + sec(x) + csc(x)',
     'Sin[x] + Cos[x] + Tan[x] + Cot[x] + Csc[x] + Sec[x]'),
    ('asin(x) + acos(x) + atan(x) + acot(x) + asec(x) + acsc(x)',
     'ArcCos[x] + ArcCot[x] + ArcCsc[x] + ArcSec[x] + ArcSin[x] + ArcTan[x]'),
    ('sinh(x) + cosh(x) + tanh(x) + coth(x) + sech(x) + csch(x)',
     'Sinh[x] + Cosh[x] + Tanh[x] + Coth[x] + Csch[x] + Sech[x]'),
    ('asinh(x) + acosh(x) + atanh(x) + acoth(x) + asech(x) + acsch(x)',
     'ArcCosh[x] + ArcCoth[x] + ArcCsch[x] + ArcSech[x] + ArcSinh[x] + ArcTanh[x]'),
    ('atan2(y, x)',
     'ArcTan[x, y]'),
    ('sign(x) + floor(x) + ceiling(y)',
     'ceiling[y] + floor[x] + sign[x]'),
    ('Max(x, y, 2) + Min(x, y)',
     'Max[2, x, y] + Min[x, y]'),
    ('conjugate(x) + Abs(x)',
     'Conjugate[x] + Abs[x]'),
    ('exp(x) + log(x) + log(x, 2)',
     'Exp[x] + Log[x] + Log[x]/Log[2]'),
    ('sqrt(x) + x**Rational(1, 3)',
     'x^(1/3) + x^(1/2)'),
    ('erf(x) + erfc(x) + erfi(x) + erfinv(x) + erfcinv(x)',
     'Erf[x] + Erfc[x] + InverseErfc[x] + Erfi[x] + InverseErf[x]'),
    ('erf2(x, y)',
     'Erf[x, y]'),
    ('expint(x, y) + Ei(x) + Si(x) + Ci(x) + Shi(x) + Chi(x) + li(x)',
     'CoshIntegral[x] + CosIntegral[x] + ExpIntegralEi[x] + SinhIntegral[x] + SinIntegral[x] + ExpIntegralE[x, y] + LogIntegral[x]'),
    ('gamma(x) + loggamma(x) + uppergamma(x, y) + lowergamma(x, y) + polygamma(x, y) + beta(x, y)',
     'Beta[x, y] + Gamma[x] + LogGamma[x] + lowergamma[x, y] + PolyGamma[x, y] + Gamma[x, y]'),
    ('factorial(x) + factorial2(x) + subfactorial(x)',
     'Factorial[x] + Factorial2[x] + Subfactorial[x]'),
    ('binomial(x, y)',
     'Gamma[x + 1]/(Gamma[y + 1]*Gamma[x - y + 1])'),
    ('RisingFactorial(x, y) + FallingFactorial(x, y)',
     'FactorialPower[x, y] + Pochhammer[x, y]'),
    ('laguerre(x, y) + assoc_laguerre(x, y, z) + hermite(x, y) + chebyshevt(x, y) + chebyshevu(x, y) + legendre(x, y)',
     'LaguerreL[x, y, z] + ChebyshevT[x, y] + ChebyshevU[x, y] + HermiteH[x, y] + LaguerreL[x, y] + LegendreP[x, y]'),
    ('jacobi(x, y, z, k) + gegenbauer(x, y, z) + assoc_legendre(x, y, z)',
     'LegendreP[x, y, z] + GegenbauerC[x, y, z] + JacobiP[x, y, z, k]'),
    ('zeta(x) + zeta(x, y) + dirichlet_eta(x) + polylog(x, y) + lerchphi(x, y, z)',
     'LerchPhi[x, y, z] + DirichletEta[x] + PolyLog[x, y] + Zeta[x] + Zeta[x, y]'),
    ('besselj(x, y) + bessely(x, y) + besseli(x, y) + besselk(x, y) + hankel1(x, y) + hankel2(x, y) + jn(x, y) + yn(x, y)',
     'BesselI[x, y] + BesselJ[x, y] + BesselK[x, y] + BesselY[x, y] + HankelH1[x, y] + HankelH2[x, y] + SphericalBesselJ[x, y] + SphericalBesselY[x, y]'),
    ('airyai(x) + airybi(x) + airyaiprime(x) + airybiprime(x)',
     'AiryAi[x] + AiryAiPrime[x] + AiryBi[x] + AiryBiPrime[x]'),
    ('elliptic_k(x) + elliptic_e(x) + elliptic_f(x, y) + elliptic_pi(x, y)',
     'EllipticE[x] + EllipticE[x, y] + EllipticK[x] + EllipticPi[x, y]'),
    ('hyper((1, 2), (3,), z)',
     'HypergeometricPFQ[{1, 2}, {3}, z]'),
    ('meijerg(((1, 2), ()), ((3,), ()), z)',
     'MeijerG[{{1, 2}, {}}, {{3}, {}}, z]'),
    ('KroneckerDelta(x, y)',
     'KroneckerDelta[x, y]'),
    ('LambertW(x) + LambertW(x, 2)',
     'ProductLog[x] + ProductLog[2, x]'),
    ('Sum(x**k, (k, 0, n))',
     'Hold[Sum[x^k, {k, 0, n}]]'),
    ('Integral(sin(x), x)',
     'Hold[Integrate[Sin[x], x]]'),
    ('Integral(sin(x), (x, 0, pi))',
     'Hold[Integrate[Sin[x], {x, 0, Pi}]]'),
    ('Integral(x*y, (x, 0, 1), (y, 0, 2))',
     'Hold[Integrate[x*y, {x, 0, 1}, {y, 0, 2}]]'),
    ('Derivative(sin(x)*y, x, y, evaluate=False)',
     'Hold[D[y*Sin[x], x, y]]'),
    ('Matrix([[x, 1], [0, y]])',
     '{{x, 1}, {0, y}}'),
    ('ImmutableSparseMatrix(2, 3, {(0, 1): x, (1, 2): 3})',
     'SparseArray[{{1, 2} -> x, {2, 3} -> 3}, {2, 3}]'),
    ('Array([[1, x], [y, 2]])',
     '{{1, x}, {y, 2}}'),
    ('Eq(x, y)',
     'x == y'),
    ('Ne(x, y)',
     'x != y'),
    ('x < y',
     'x < y'),
    ('Ge(x, y)',
     'x >= y'),
    ('And(x > 0, Or(y < 1, Not(z > 2)))',
     '(z <= 2 || y < 1) && x > 0'),
    ('(x, y, 1)',
     '{x, y, 1}'),
    ('[x, [y, 2]]',
     '{x, {y, 2}}'),
    ('Tuple(x, y)',
     '{x, y}'),
    ('pi + E + I',
     'E + Pi + I'),
    ('oo',
     'Infinity'),
    ('-oo',
     '-Infinity'),
    ('zoo',
     'ComplexInfinity'),
    ('nan',
     'Indeterminate'),
    ('GoldenRatio + EulerGamma + Catalan',
     'EulerGamma + Catalan + GoldenRatio'),
    ('TribonacciConstant',
     '(1/3 + (1/3)*(19 - 3*33^(1/2))^(1/3) + (1/3)*(3*33^(1/2) + 19)^(1/3))'),
    ('-x/3 + x**2*Rational(2, 5)',
     '(2/5)*x^2 - 1/3*x'),
    ('x**y**z',
     'x^(y^z)'),
    ('1/(x*y)',
     '1/(x*y)'),
    ('(x + y)**-1',
     '(x + y)^(-1)'),
    ("Float('1.0e-30')",
     '1.00000000000000e-30'),
    ('2*I*x',
     '2*I*x'),
    ('A*B*x',
     '(x)*A**B'),
]


def _mathematica_namespace():
    namespace = dict(vars(sx))
    namespace.update(zip('xyznk', sx.symbols('x y z n k')))
    namespace.update(f=sx.Function('f'), g=sx.Function('g'))
    namespace.update(zip('AB', sx.symbols('A B', commutative=False)))
    return namespace


@pytest.mark.parametrize('source,expected', MATHEMATICA_CASES)
def test_mathematica_code_matches_sympy(source, expected):
    assert p.mathematica_code(eval(source, _mathematica_namespace())) == expected


def test_mathematica_code_settings_match_sympy():
    from symox.printing.codeprinter import PrintMethodNotImplementedError
    f = sx.Function('f')
    assert p.mathematica_code(f(x), user_functions={'f': 'F'}) == 'F[x]'
    assert p.mathematica_code(sx.sin(x), user_functions={'sin': 'MySin'}) == 'MySin[x]'
    rules = {'exp': [(lambda e: e.is_Symbol, 'EE'), (lambda e: True, 'Exp2')]}
    assert p.mathematica_code(sx.exp(x) + sx.exp(x + 1), user_functions=rules) == 'EE[x] + Exp2[x + 1]'
    assert p.mathematica_code(sx.sin(x) + sx.zoo, human=False) == (set(), set(), 'Sin[x] + ComplexInfinity')
    # MCodePrinter has no Piecewise printer: every Piecewise reaches the
    # unsupported ExprCondPair, with or without a final True branch.
    message = ("Unsupported by <class 'symox.printing.mathematica.MCodePrinter'>: "
               "<class 'symox.functions.elementary.piecewise.ExprCondPair'>\n"
               "Set the printer option 'strict' to False in order to generate partially printed code.")
    for pw in (sx.Piecewise((x, x > 0)), sx.Piecewise((x, x > 0), (-x, True))):
        with pytest.raises(PrintMethodNotImplementedError) as info:
            p.mathematica_code(pw)
        assert str(info.value) == message
    assert p.mathematica_code(sx.Piecewise((x, x > 0), (0, True)), strict=False) == (
        '(* Not supported in Wolfram Language: *)\n(* ExprCondPair *)\n(* ExprCondPair *)\n'
        'Piecewise[(x, x > 0), (0, True)]')


@pytest.mark.parametrize('name,comment', [('mathematica_code', '(* Not supported in Wolfram Language: *)\n(* Indexed *)\n'),
                                          ('maple_code', '# Not supported in maple:\n# Indexed\n')])
def test_indexed_unsupported_in_mathematica_and_maple(name, comment):
    # SymPy's Maple and Mathematica printers have no Indexed printer.
    from symox.printing.codeprinter import PrintMethodNotImplementedError
    expression = sx.IndexedBase('A')[sx.Idx('i')] + 1
    with pytest.raises(PrintMethodNotImplementedError, match="<class 'symox.tensor.indexed.Indexed'>"):
        getattr(p, name)(expression)
    assert getattr(p, name)(expression, strict=False) == comment + 'A[i] + 1'


PIECEWISE_DEFAULT_MESSAGE = ('All Piecewise expressions must contain an (expr, True) statement to be used '
                             'as a default condition. Without one, the generated expression may not evaluate '
                             'to anything under some condition.')


@pytest.mark.parametrize('name,settings', [
    ('ccode', {}), ('ccode', {'standard': 'C89'}), ('cxxcode', {}), ('fcode', {}), ('fcode', {'source_format': 'free'}),
    ('rust_code', {}), ('jscode', {}), ('julia_code', {}), ('julia_code', {'inline': False}), ('octave_code', {}),
    ('octave_code', {'inline': False}), ('glsl_code', {}), ('maple_code', {}), ('maple_code', {'inline': False})])
@pytest.mark.parametrize('assign_to', [None, 'z'])
def test_piecewise_without_default_raises_value_error(name, settings, assign_to):
    # SymPy 1.14 raises ValueError in these printers, also with strict=False.
    for pw in (sx.Piecewise((x, x > 0)), x + sx.Piecewise((x, x > 0), (y, x < -1))):
        for extra in ({}, {'strict': False}):
            with pytest.raises(ValueError) as info:
                getattr(p, name)(pw, assign_to=assign_to, **settings, **extra)
            assert str(info.value) == PIECEWISE_DEFAULT_MESSAGE


def test_piecewise_without_default_in_rcode_and_pycode():
    # SymPy 1.14 rcode and pycode print a missing default as NA / None.
    pw = sx.Piecewise((x, x > 0), (y, x < -1))
    assert p.rcode(pw) == 'ifelse(x > 0,x,ifelse(x < -1,y,NA))'
    assert p.rcode(pw, assign_to='z') == 'z = ifelse(x > 0,x,ifelse(x < -1,y,NA));'
    assert p.pycode(pw) == '((x) if (x > 0) else (y) if (x < -1) else None)'


@pytest.mark.parametrize('name,expected', [
    ('julia_code', 'if (x > 0)\n    z = x + 1\nelseif (x < -1)\n    z = y\nelse\n    z = x .* y\nend'),
    ('octave_code', 'if (x > 0)\n  z = x + 1;\nelseif (x < -1)\n  z = y;\nelse\n  z = x.*y;\nend')])
def test_piecewise_assignment_without_inline(name, expected):
    pw = sx.Piecewise((x + 1, x > 0), (y, x < -1), (x * y, True))
    assert getattr(p, name)(pw, assign_to='z', inline=False) == expected
