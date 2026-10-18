"""Printer metadata, composable layouts and optional backend contracts."""
import importlib
import inspect
import pytest
import symox as S
import sympy as R


@pytest.mark.parametrize('text',['a_x^1','var_sub1__sup_sub2','alpha12','_x','x__1__2_3','xhat',''])
def test_symbol_conventions_match_reference(text):
    from symox.printing.conventions import split_super_sub
    from sympy.printing.conventions import split_super_sub as reference
    assert split_super_sub(text)==reference(text)


def test_latex_translation_and_escaping():
    module=importlib.import_module('symox.printing.latex');reference=importlib.import_module('sympy.printing.latex')
    for text in ['alpha','alphahatdotprime','abc','xhat','beta','Omega']:
        assert module.translate(text)==reference.translate(text)
    text='a_b&c%$#{}~^\\'
    assert module.latex_escape(text)==reference.latex_escape(text)


def test_pretty_symbols_shapes_and_centering():
    native=importlib.import_module('symox.printing.pretty.pretty_symbology')
    reference=importlib.import_module('sympy.printing.pretty.pretty_symbology')
    old_native,old_reference=native.pretty_use_unicode(),reference.pretty_use_unicode()
    try:
        for enabled in [False,True]:
            native.pretty_use_unicode(enabled);reference.pretty_use_unicode(enabled)
            for name in ['alpha_1','x_12','beta','theta','x__2']:
                assert native.pretty_symbol(name)==reference.pretty_symbol(name)
            for shape in ['(',')','[',']','{','}','-']:
                for size in range(1,7):assert native.xobj(shape,size)==reference.xobj(shape,size)
            assert native.center('abc',6)==reference.center('abc',6)
            assert native.center_pad(3,6)==(' ','  ')
            if enabled:assert native.pretty_atom('Pi')=='π'
    finally:native.pretty_use_unicode(old_native);reference.pretty_use_unicode(old_reference)


def test_layout_composition_aligns_baselines_and_preserves_entries():
    from symox.printing.pretty.stringpict import stringPict,prettyForm
    from sympy.printing.pretty.stringpict import stringPict as Picture
    left=stringPict('a\nb',1);right=stringPict('x\ny\nz',1)
    expected=Picture('a\nb',1).right(' + ',Picture('x\ny\nz',1))
    assert left.right(' + ',right)==expected
    assert stringPict('x+3').below(stringPict.LINE,'3')[0].replace('─','-')=='x+3\n---\n 3 '
    quotient=prettyForm('x')/prettyForm('y')
    assert quotient.height()==3 and quotient.baseline==1
    power=quotient**prettyForm('2')
    assert power.baseline>quotient.baseline and '2' in str(power)
    assert prettyForm.apply('f',prettyForm('x')).s=='f(x)'


def test_printer_context_and_signatures():
    from symox.printing.printer import Printer,printer_context,print_function
    printer=Printer()
    with pytest.raises(RuntimeError):
        with printer_context(printer,a=1):
            assert printer._context['a']==1
            raise RuntimeError()
    assert printer._context=={}
    class CustomPrinter(Printer):_default_settings={'precision':17,'order':None}
    @print_function(CustomPrinter)
    def output(expr,**settings):return CustomPrinter(settings).doprint(expr)
    assert 'precision' in inspect.signature(output).parameters
    from symox.printing.str import StrReprPrinter
    assert StrReprPrinter().doprint(['x',S.Symbol('x')])=="['x', x]"
    from symox.printing.mathml import MathMLPrinter
    assert '<plus/>' in MathMLPrinter().doprint(S.Symbol('x')+1)


def test_function_metadata_is_connected_to_the_language_printers():
    from symox.printing.julia import julia_code,known_fcns_src1
    from symox.printing.octave import octave_code,known_fcns_src2
    x=S.Symbol('x')
    assert 'besselj' in known_fcns_src1 and julia_code(S.besselj(2,x))=='besselj(2, x)'
    assert known_fcns_src2['LambertW']=='lambertw' and octave_code(S.LambertW(x))=='lambertw(x)'
    with pytest.raises(NotImplementedError):julia_code(S.Function('besselj')(2,x))


def test_optional_aesara_broadcasting_contract():
    from symox.printing.aesaracode import dim_handling,aesara_code
    x,y=S.symbols('x y')
    assert dim_handling([x,y],dim=2)=={x:(False,False),y:(False,False)}
    assert dim_handling([x,y],dims={x:0,y:2})=={x:(True,True),y:(False,False)}
    assert dim_handling([x,y])=={}
    if importlib.util.find_spec('aesara') is None:
        with pytest.raises(ImportError):aesara_code(S.sin(x))


def test_testing_adapters_execute_real_assertions_and_doctests(tmp_path):
    from symox.testing.pytest import raises,warns_deprecated_sympy
    from symox.testing.runtests import sympytestfile,split_list,doctest,SymPyDocTestRunner,SymPyDocTestFinder
    from symox.utilities.exceptions import sympy_deprecation_warning
    with raises(ZeroDivisionError) as caught:1/0
    assert isinstance(caught.value,ZeroDivisionError)
    with warns_deprecated_sympy():sympy_deprecation_warning('old API',deprecated_since_version='0.1',active_deprecations_target='test')
    assert split_list(list(range(10)),'2/3')==[3,4,5]
    assert split_list(list(range(10)),'3/3')==[6,7,8,9]
    path=tmp_path/'example.txt';path.write_text('>>> from symox import Rational\n>>> Rational(1, 2) + Rational(1, 3)\n5/6\n')
    result=sympytestfile(str(path),module_relative=False)
    assert result.failed==0 and result.attempted==2
    assert doctest(str(path))


def test_interactive_transformations_preserve_strings_and_rationals():
    from symox.interactive.session import int_to_Integer,enable_automatic_symbols,enable_automatic_int_sympification
    import ast
    from types import SimpleNamespace
    namespace={'Integer':S.Integer}
    result=eval(int_to_Integer('1/2 + 0x10'),namespace)
    assert result==S.Rational(33,2)
    assert eval(int_to_Integer("'123'"),namespace)=='123'
    shell=SimpleNamespace(user_ns={},ast_transformers=[],input_transformers_cleanup=[])
    transformer=enable_automatic_symbols(shell)
    transformer.visit(ast.parse('unknown_x + 1'))
    assert shell.user_ns['unknown_x']==S.Symbol('unknown_x')
    enable_automatic_int_sympification(shell)
    assert shell.input_transformers_cleanup[0](['1/2\n'])


def test_python_session_source_transformations_are_connected(monkeypatch):
    import code
    import symox.interactive as interactive
    import sys
    original=sys.displayhook
    captured={}
    def interact(console,banner=''):
        console.runsource('answer = 1/2 + missing_name')
        captured.update(console.locals)
    monkeypatch.setattr(code.InteractiveConsole,'interact',interact)
    try:interactive.init_session(quiet=True,auto_symbols=True,auto_int_to_Integer=True)
    finally:sys.displayhook=original
    assert captured['answer']==S.Rational(1,2)+S.Symbol('missing_name')
