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


@pytest.mark.parametrize('name',['ccode','cxxcode','fcode','pycode','rust_code','jscode','julia_code','octave_code','rcode','maple_code','mathematica_code','glsl_code'])
def test_code_entrypoints(name):
    text=getattr(p,name)(x*x+sx.sin(y),assign_to='result')
    assert isinstance(text,str) and 'result' in text and 'x' in text and 'sin' in text.lower()


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
    for fn in [p.jscode,p.julia_code,p.octave_code,p.rcode,p.mathematica_code,p.glsl_code]:
        with pytest.raises(NotImplementedError):fn(expression)
        assert 'mapped' in fn(expression,user_functions={'custom':'mapped'})
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
    for options in [{'terms_per_line':0},{'environment':'bad'}]:
        with pytest.raises(ValueError):p.multiline_latex(y,x,**options)


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
