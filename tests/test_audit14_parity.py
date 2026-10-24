"""audit14 (wave 14): sorted differentiation variables, the
Derivative(f(g(x)), g(x)) chain form, Derivative._eval_subs cases and
Subs derivatives. Expected values are hard-coded from SymPy 1.14."""
from symox import *
del test

x, y, z = symbols('x y z')
f, g = symbols('f g', cls=Function)


def test_sorted_variables():
    # Derivative.__new__ sorts with _sort_variable_count.
    assert str(diff(f(x, y), y, x)) == 'Derivative(f(x, y), x, y)'
    assert str(Derivative(f(x, y), y, x, evaluate=True)) == 'Derivative(f(x, y), x, y)'
    assert str(Derivative(Derivative(f(x, y), y), x).doit()) == 'Derivative(f(x, y), x, y)'
    assert Derivative(f(x, y, z), z, x, y, evaluate=True).variables == (x, y, z)
    assert str(diff(f(x, y)*g(z), z, y, x)) == 'Derivative(g(z), z)*Derivative(f(x, y), x, y)'


def test_function_argument_chain_rule():
    assert str(diff(f(g(x)), x)) == 'Derivative(f(g(x)), g(x))*Derivative(g(x), x)'
    assert str(Derivative(f(g(x)), g(x)).subs(g(x), x**2)) == 'Subs(Derivative(f(_xi), _xi), _xi, x**2)'
    assert str(diff(f(x, y, g(x)), y, g(x), x)) == (
        'Derivative(g(x), x)*Derivative(f(x, y, g(x)), y, (g(x), 2))'
        ' + Subs(Derivative(f(x, y, _xi), _xi, x, y), _xi, g(x))')
    assert str(diff(f(x, g(x)), x, 2)) == (
        '(Derivative(f(x, g(x)), (g(x), 2))*Derivative(g(x), x)'
        ' + Subs(Derivative(f(x, _xi), _xi, x), _xi, g(x)))*Derivative(g(x), x)'
        ' + Derivative(f(x, g(x)), g(x))*Derivative(g(x), (x, 2))'
        ' + Derivative(g(x), x)*Subs(Derivative(f(_xi_1, g(x)), _xi_1, g(x)), _xi_1, x)'
        ' + Subs(Derivative(f(_xi_1, g(x)), (_xi_1, 2)), _xi_1, x)')


def test_subs_derivatives():
    assert str(Subs(Derivative(f(x, y), y), y, z).doit()) == 'Derivative(f(x, z), z)'
    assert str(diff(Subs(f(x, y), y, z), z)) == 'Derivative(f(x, z), z)'
    assert str(diff(Subs(Derivative(f(x, y), y), y, g(x)), g(x))) == 'Derivative(f(x, g(x)), (g(x), 2))'


def test_function_argument_inside_other_argument():
    # Intentional correction: SymPy 1.14 gives Derivative(f(g(x), g(x)**2), g(x))
    # for the first slot, which also differentiates the second argument, and
    # the substitution f -> x*y, g -> sin then gives 5*sin(x)**2*cos(x)
    # instead of diff(sin(x)**3, x) = 3*sin(x)**2*cos(x).
    e = diff(f(g(x), g(x)**2), x)
    v = e.subs(f, Lambda((x, y), x*y)).doit().subs(g, Lambda(x, sin(x))).doit()
    assert simplify(v - diff(sin(x)**3, x)) == 0
    for point in (Rational(1, 3), Rational(7, 5)):
        assert abs(N(v.subs(x, point) - 3*sin(point)**2*cos(point), 30)) < 1e-25


def test_hyper_function_suitable_origin_and_has():
    from symox.simplify.hyperexpand import Hyper_Function
    assert Hyper_Function((S.Half,), (Rational(3, 2),))._is_suitable_origin() is True
    assert Hyper_Function((S.Half,), (S.Half,))._is_suitable_origin() is False
    assert Hyper_Function((S.Half,), (Rational(-1, 2),))._is_suitable_origin() is False
    assert Hyper_Function((S.Half,), (0,))._is_suitable_origin() is False
    assert Hyper_Function((S.Half,), (-1, 1,))._is_suitable_origin() is False
    assert Hyper_Function((S.Half, 0), (1,))._is_suitable_origin() is False
    a, b, c = symbols('a, b, c', cls=Dummy)
    h = Hyper_Function([2, -a], [b])
    assert h.has(a) and h.has(Tuple(b)) and not h.has(c)
