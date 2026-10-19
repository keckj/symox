"""Symbol printer hooks retain native output and compose through super()."""
import pytest
import symox as sx
import sympy as sp
from symox.printing.latex import LatexPrinter
from sympy.printing.latex import LatexPrinter as ReferencePrinter


@pytest.mark.parametrize('name', ['x','alpha','beta_1','x12','x_1_2','x__2','theta^2',
                                 'x_{raw}','alpha_beta','x_omega2','xhat','_x'])
@pytest.mark.parametrize('style', ['plain','bold'])
def test_symbol_hook_signature_and_name_components(name,style):
    actual=LatexPrinter()._print_Symbol(sx.Symbol(name),style=style)
    expected=ReferencePrinter()._print_Symbol(sp.Symbol(name),style=style)
    assert actual==expected


@pytest.mark.parametrize('name', ['alpha_beta','x^2','x__2'])
def test_disabled_splitting_and_explicit_symbol_names(name):
    for native,reference in [(LatexPrinter({'disable_split_super_sub':True}),
                              ReferencePrinter({'disable_split_super_sub':True}))]:
        assert native.doprint(sx.Symbol(name))==reference.doprint(sp.Symbol(name))
    symbol=sx.Symbol(name)
    assert LatexPrinter({'symbol_names':{symbol:r'CUSTOM'}})._print_Symbol(symbol,style='bold')=='CUSTOM'


def test_downstream_override_can_call_parent_symbol_hook():
    class Custom(LatexPrinter):
        def _print_Symbol(self,expr):
            return 'custom<' + super()._print_Symbol(expr) + '>'
    assert Custom().doprint(sx.Symbol('alpha'))==r'custom<\alpha>'


@pytest.mark.parametrize('mode', ['plain','inline','equation','equation*'])
@pytest.mark.parametrize('itex', [False,True])
def test_symbol_modes_are_applied_once(mode,itex):
    settings=dict(mode=mode,itex=itex)
    assert sx.latex(sx.Symbol('alpha'),**settings)==sp.latex(sp.Symbol('alpha'),**settings)
    assert sx.latex(sx.Symbol('alpha')+1,**settings)==sp.latex(sp.Symbol('alpha')+1,**settings)


@pytest.mark.parametrize('token',['omega2',r'x\_2','alphahatdotprime','xvec','alpha_beta'])
def test_translate_does_not_reparse_symbol_components(token):
    from symox.printing.latex import translate
    from sympy.printing.latex import translate as reference
    assert translate(token)==reference(token)
