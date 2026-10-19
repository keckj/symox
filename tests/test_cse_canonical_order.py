"""Canonical CSE uses public arguments and node counts without width cutoffs."""
import sympy as sp
import symox as sx


def compare(expressions,order='canonical'):
    replacements,reduced=sx.cse(expressions,order=order)
    expected=sp.cse([expression.to_sympy() for expression in expressions],order=order)
    assert [(a.to_sympy(),b.to_sympy()) for a,b in replacements]==expected[0]
    assert [value.to_sympy() for value in reduced]==expected[1]


def test_canonical_order_prefers_smaller_trees_before_sort_keys():
    a,b,c=sx.symbols('a b c');f=sx.Function('f');g=sx.Function('g')
    expressions=[f(a+b+c)+g(a+b+c)+sx.sin(a)*sx.cos(b),
                 sx.sin(a)+sx.cos(b)+f(a+b+c)]
    compare(expressions)


def test_canonical_order_remains_active_on_wide_sums():
    x=sx.Symbol('x');functions=[sx.Function('f'+str(i))(x) for i in range(72)]
    expression=sum(functions)
    compare([expression+sx.sin(expression),*(function+sx.cos(function) for function in functions)])


def test_storage_order_option_reconstructs_exact_input():
    x,y,z=sx.symbols('x y z')
    expressions=[(x+y)*z+(x+y)**2,sx.sin(x+y)]
    for order in ['none','canonical',None,'invalid']:
        replacements,reduced=sx.cse(expressions,order=order)
        for symbol,value in reversed(replacements):reduced=[expression.xreplace({symbol:value}) for expression in reduced]
        assert reduced==expressions
