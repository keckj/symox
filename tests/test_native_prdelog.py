"""Native base-case is_log_deriv_k_t_radical_in_field,
parametric_log_deriv_heu and limited_integrate (prdelog), and the native
splitfactor precheck.

The cases are base-case calls recorded from risch_integrate on the audit
integrand corpus plus seeded random inputs
(audits/performance/2026-10-30-prdelog/scripts/make_tests.py); the expected
outputs (srepr, Poly domains and container types) are hard-coded from
SymPy 1.14 (sympy.integrals.prde). The native kernel must answer each case
(_prde_native.STATS) unless listed as a decline.
"""
import pytest

from symox import Symbol, Poly, Rational, S, srepr, sympify, symbols, log
from symox.integrals.risch import DifferentialExtension, splitfactor
from symox.integrals.prde import (is_log_deriv_k_t_radical_in_field,
                                  parametric_log_deriv_heu,
                                  parametric_log_deriv, limited_integrate)
from symox.integrals import _prde_native
from symox.integrals import _risch_splitfactor_native as _sf
from symox import _core

x, t = symbols('x t')
null = None  # the cases are JSON

CASES = [
    ["limited", [["x", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Rational(1, 2), Pow(Symbol('x'), Integer(2)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["1", "ZZ"], ["x", "ZZ"]], ["tuple", ["One", "Integer(1)"], ["Symbol", "Symbol('x')"]]],
    ["limited", [["0", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(0)", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["-1", "ZZ"], ["1", "ZZ"]], null],
    ["is_log_deriv", [["1", "ZZ"], ["1", "ZZ"]], null],
    ["is_log_deriv", [["1", "ZZ"], ["x + 1", "ZZ"]], ["tuple", ["One", "Integer(1)"], ["Add", "Add(Symbol('x'), Integer(1))"]]],
    ["limited", [["-x", "ZZ"], ["x + 1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["is_log_deriv", [["-2*x", "ZZ"], ["1", "ZZ"]], null],
    ["is_log_deriv", [["2*x", "ZZ"], ["1", "ZZ"]], null],
    ["limited", [["1", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Symbol('x')", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["-2", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Integer(2), Symbol('x'))", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["-3", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Integer(3), Symbol('x'))", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["6", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(6), Symbol('x'))", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["x**2", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Rational(1, 3), Pow(Symbol('x'), Integer(3)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["-2*x**2", "ZZ"], ["3", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Rational(2, 9), Pow(Symbol('x'), Integer(3)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(0)", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(1)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x**2", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(-1)", "ZZ"], ["Poly", "Symbol('x')", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x + 1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["limited", [["1", "ZZ"], ["x", "ZZ"], [[["1", "ZZ"], ["x + 1", "ZZ"]]]], null],
    ["limited", [["1", "ZZ"], ["1", "ZZ"], [[["2*x", "ZZ"], ["x**2 + 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Symbol('x')", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["x", "ZZ"], ["1", "ZZ"], [[["2*x", "ZZ"], ["x**2 + 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Rational(1, 2), Pow(Symbol('x'), Integer(2)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x**2 + 1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["limited", [["x", "ZZ"], ["x**2 + 1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["is_log_deriv", [["-2", "ZZ"], ["1", "ZZ"]], null],
    ["is_log_deriv", [["2", "ZZ"], ["1", "ZZ"]], null],
    ["limited", [["-1", "ZZ"], ["2*x + 2", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["limited", [["1", "ZZ"], ["1", "ZZ"], [[["-1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Symbol('x')", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["x", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x + 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Rational(1, 2), Pow(Symbol('x'), Integer(2)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["x**3", "ZZ"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Rational(1, 4), Pow(Symbol('x'), Integer(4)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["-x**3", "ZZ"], ["2", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Rational(1, 8), Pow(Symbol('x'), Integer(4)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x - 1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["limited", [["-2 + log(2)", "ZZ[log(2)]"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Symbol('x'), Add(Integer(-2), log(Integer(2))))", "ZZ[log(2)]"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x**2 + 2*x + 1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(-1)", "ZZ"], ["Poly", "Add(Symbol('x'), Integer(1))", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["x + 1", "ZZ"], ["x**2", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(-1)", "ZZ"], ["Poly", "Symbol('x')", "QQ"]], ["list", ["Poly", "Integer(1)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["x**2 + x", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], null],
    ["limited", [["-1", "ZZ"], ["2*x", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(0)", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Rational(-1, 2)", "QQ"]]]],
    ["is_log_deriv", [["1", "ZZ"], ["x**2", "ZZ"]], null],
    ["is_log_deriv", [["-1", "ZZ"], ["x**2", "ZZ"]], null],
    ["is_log_deriv", [["-96*x**4 - 172*x**3 + 4*x**2 + 109*x + 30", "QQ"], ["-16*x**5 - 40*x**4 - 12*x**3 + 26*x**2 + 12*x", "ZZ"]], ["tuple", ["Integer", "Integer(2)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(1, 2)), Integer(2)), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Rational(3, 2), Symbol('x'))), Integer(5)))"]]],
    ["is_log_deriv", [["-2*x**3/3 - x**2 + 2*x/3 - 2", "QQ"], ["3*x - 2", "ZZ"]], null],
    ["heu", [["-96*x**4 - 172*x**3 + 4*x**2 + 109*x + 30", "QQ"], ["-16*x**5 - 40*x**4 - 12*x**3 + 26*x**2 + 12*x", "ZZ"], ["-x - 3", "ZZ"], ["3/2 - x", "QQ"]], null],
    ["heu", [["32*x**6 + 272*x**5 + 292*x**4 - 242*x**3 - 283*x**2 + 123*x/2 + 45", "QQ"], ["16*x**6 + 16*x**5 - 48*x**4 - 44*x**3 + 27*x**2 + 18*x", "QQ"], ["-x - 3", "ZZ"], ["3/2 - x", "QQ"]], null],
    ["limited", [["-96*x**4 - 172*x**3 + 4*x**2 + 109*x + 30", "QQ"], ["-16*x**5 - 40*x**4 - 12*x**3 + 26*x**2 + 12*x", "ZZ"], [[["-2", "ZZ"], ["1 - 3*x/2", "QQ"]], [["-3*x**2 + 2*x/3", "QQ"], ["-3*x**2", "ZZ"]]]], null],
    ["limited", [["-2", "ZZ"], ["1/2", "QQ"], [[["-2", "ZZ"], ["1 - 3*x/2", "QQ"]], [["-3*x**2 + 2*x/3", "QQ"], ["-3*x**2", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Integer(4), Symbol('x'))", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["-2*x**3 + 11*x**2/3 + 4*x/3 - 6", "QQ"], ["2*x**3 - 3*x**2 - 4*x + 4", "ZZ"]], null],
    ["is_log_deriv", [["-x**2/2", "QQ"], ["-x - 4/3", "QQ"]], null],
    ["heu", [["-2*x**3 + 11*x**2/3 + 4*x/3 - 6", "QQ"], ["2*x**3 - 3*x**2 - 4*x + 4", "ZZ"], ["x", "ZZ"], ["-2", "ZZ"]], null],
    ["heu", [["2*x**4 + x**3 - 34*x**2/3 + 4*x/3 + 12", "QQ"], ["-4*x**3 + 6*x**2 + 8*x - 8", "ZZ"], ["x", "ZZ"], ["-2", "ZZ"]], null],
    ["limited", [["-2*x**3 + 11*x**2/3 + 4*x/3 - 6", "QQ"], ["2*x**3 - 3*x**2 - 4*x + 4", "ZZ"], []], null],
    ["limited", [["-3*x**2 - 2*x + 1", "ZZ"], ["x", "ZZ"], []], null],
    ["is_log_deriv", [["-112*x**3 - 112*x**2 + 88*x + 56", "ZZ"], ["-16*x**4 + 24*x**2 - 8", "ZZ"]], ["tuple", ["One", "Integer(1)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Integer(-1)), Integer(5)), Pow(Add(Symbol('x'), Integer(1)), Integer(-2)), Pow(Add(Pow(Symbol('x'), Integer(2)), Rational(-1, 2)), Integer(2)))"]]],
    ["is_log_deriv", [["x**2/3 - x + 2/3", "QQ"], ["x**3 + x**2 + x - 3", "ZZ"]], null],
    ["heu", [["-112*x**3 - 112*x**2 + 88*x + 56", "ZZ"], ["-16*x**4 + 24*x**2 - 8", "ZZ"], ["-1", "ZZ"], ["-4*x - 4", "ZZ"]], null],
    ["heu", [["464*x**4 + 896*x**3 + 72*x**2 - 576*x - 216", "ZZ"], ["64*x**5 + 64*x**4 - 96*x**3 - 96*x**2 + 32*x + 32", "ZZ"], ["-1", "ZZ"], ["-4*x - 4", "ZZ"]], null],
    ["limited", [["-112*x**3 - 112*x**2 + 88*x + 56", "ZZ"], ["-16*x**4 + 24*x**2 - 8", "ZZ"], [[["x**2/3 + 3*x/2 - 1", "QQ"], ["-3*x**2 - 4*x", "ZZ"]], [["-2", "ZZ"], ["2*x + 1", "ZZ"]]]], null],
    ["limited", [["-3*x**2 + 2*x - 1", "ZZ"], ["2", "ZZ"], [[["x**2/3 + 3*x/2 - 1", "QQ"], ["-3*x**2 - 4*x", "ZZ"]], [["-2", "ZZ"], ["2*x + 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Add(Mul(Integer(-1), Rational(1, 2), Pow(Symbol('x'), Integer(3))), Mul(Rational(1, 2), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Rational(1, 2), Symbol('x')))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["-2/3", "QQ"], ["2*x + 3", "ZZ"]], ["tuple", ["Integer", "Integer(3)"], ["Pow", "Pow(Add(Symbol('x'), Rational(3, 2)), Integer(-1))"]]],
    ["is_log_deriv", [["-1", "ZZ"], ["4*x", "ZZ"]], ["tuple", ["Integer", "Integer(4)"], ["Pow", "Pow(Symbol('x'), Integer(-1))"]]],
    ["heu", [["-2/3", "QQ"], ["2*x + 3", "ZZ"], ["-2*x**2 + 2*x + 3", "ZZ"], ["3*x/2 + 2", "QQ"]], ["tuple", ["Integer", "Integer(3)"], ["Zero", "Integer(0)"], ["Pow", "Pow(Add(Symbol('x'), Rational(3, 2)), Integer(-1))"]]],
    ["heu", [["-4*x**3 - 2*x**2 + 11*x + 23/3", "QQ"], ["3*x**2 + 17*x/2 + 6", "QQ"], ["-2*x**2 + 2*x + 3", "ZZ"], ["3*x/2 + 2", "QQ"]], ["tuple", ["Integer", "Integer(3)"], ["Integer", "Integer(3)"], ["Pow", "Pow(Add(Symbol('x'), Rational(3, 2)), Integer(-1))"]]],
    ["limited", [["-2/3", "QQ"], ["2*x + 3", "ZZ"], []], null],
    ["limited", [["-x - 3", "ZZ"], ["x**2 + 3", "ZZ"], []], null],
    ["is_log_deriv", [["2*x**6/9 - 5*x**5/9 + 118*x**4/81 - 1901*x**3/81 + 9*x**2 - 1603*x/27 + 247/9", "QQ"], ["-2*x**5/3 + 17*x**4/3 - 442*x**3/27 + 1097*x**2/27 - 106*x/9 - 16", "QQ"]], null],
    ["is_log_deriv", [["2*x + 4/3", "QQ"], ["-2*x**3 + 3*x**2 - x + 1", "ZZ"]], null],
    ["heu", [["2*x**6/9 - 5*x**5/9 + 118*x**4/81 - 1901*x**3/81 + 9*x**2 - 1603*x/27 + 247/9", "QQ"], ["-2*x**5/3 + 17*x**4/3 - 442*x**3/27 + 1097*x**2/27 - 106*x/9 - 16", "QQ"], ["-x**2/3 + 4*x/3 + 2/3", "QQ"], ["3", "ZZ"]], null],
    ["heu", [["2*x**7/3 - 23*x**6/3 + 973*x**5/27 - 2441*x**4/27 + 1921*x**3/27 + 2083*x**2/27 - 797*x/3 + 151/3", "QQ"], ["-2*x**5 + 17*x**4 - 442*x**3/9 + 1097*x**2/9 - 106*x/3 - 48", "QQ"], ["-x**2/3 + 4*x/3 + 2/3", "QQ"], ["3", "ZZ"]], null],
    ["limited", [["2*x**6/9 - 5*x**5/9 + 118*x**4/81 - 1901*x**3/81 + 9*x**2 - 1603*x/27 + 247/9", "QQ"], ["-2*x**5/3 + 17*x**4/3 - 442*x**3/27 + 1097*x**2/27 - 106*x/9 - 16", "QQ"], [[["-2", "ZZ"], ["-x**2 + x + 1", "ZZ"]]]], null],
    ["limited", [["2", "ZZ"], ["2*x + 4/3", "QQ"], [[["-2", "ZZ"], ["-x**2 + x + 1", "ZZ"]]]], null],
    ["is_log_deriv", [["-4*x**3 - 17*x**2/3 - 5*x - 1/2", "QQ"], ["4*x**3 + 17*x**2/3 + 5*x + 2", "QQ"]], null],
    ["is_log_deriv", [["3*x**2 - 3*x - 3", "ZZ"], ["-x**3 + 2*x**2 - x - 2", "ZZ"]], null],
    ["heu", [["-4*x**3 - 17*x**2/3 - 5*x - 1/2", "QQ"], ["4*x**3 + 17*x**2/3 + 5*x + 2", "QQ"], ["-4*x - 1", "ZZ"], ["3 - 2*x", "ZZ"]], null],
    ["heu", [["-8*x**4 - 82*x**3/3 - 98*x**2/3 - 27*x - 7/2", "QQ"], ["-8*x**4 + 2*x**3/3 + 7*x**2 + 11*x + 6", "QQ"], ["-4*x - 1", "ZZ"], ["3 - 2*x", "ZZ"]], null],
    ["limited", [["-4*x**3 - 17*x**2/3 - 5*x - 1/2", "QQ"], ["4*x**3 + 17*x**2/3 + 5*x + 2", "QQ"], []], null],
    ["limited", [["x**2/3 - 2*x - 4", "QQ"], ["-2*x**2/3 + x + 4", "QQ"], []], null],
    ["is_log_deriv", [["-4*x**4/3 - 106*x**3/27 + 209*x**2/18 + 10*x/9 - 173/27", "QQ"], ["2*x**4 - 65*x**3/18 + x**2/3 - 23*x/9 + 2/3", "QQ"]], null],
    ["is_log_deriv", [["4", "ZZ"], ["x + 1", "ZZ"]], ["tuple", ["One", "Integer(1)"], ["Pow", "Pow(Add(Symbol('x'), Integer(1)), Integer(4))"]]],
    ["heu", [["-4*x**4/3 - 106*x**3/27 + 209*x**2/18 + 10*x/9 - 173/27", "QQ"], ["2*x**4 - 65*x**3/18 + x**2/3 - 23*x/9 + 2/3", "QQ"], ["-2*x**2/3 + 4*x + 4", "QQ"], ["-2/3", "QQ"]], null],
    ["heu", [["-8*x**6/3 + 562*x**5/27 - 112*x**4/9 - 1636*x**3/81 - 713*x**2/27 - 428*x/27 + 778/81", "QQ"], ["-4*x**4/3 + 65*x**3/27 - 2*x**2/9 + 46*x/27 - 4/9", "QQ"], ["-2*x**2/3 + 4*x + 4", "QQ"], ["-2/3", "QQ"]], null],
    ["limited", [["-4*x**4/3 - 106*x**3/27 + 209*x**2/18 + 10*x/9 - 173/27", "QQ"], ["2*x**4 - 65*x**3/18 + x**2/3 - 23*x/9 + 2/3", "QQ"], [[["4*x**2 + 2*x + 4", "ZZ"], ["-3*x**2 + 2*x/3 - 1/3", "QQ"]], [["-x**2 + x + 1/3", "QQ"], ["4*x/3 + 1", "QQ"]]]], null],
    ["limited", [["3*x + 3/2", "QQ"], ["-1", "ZZ"], [[["4*x**2 + 2*x + 4", "ZZ"], ["-3*x**2 + 2*x/3 - 1/3", "QQ"]], [["-x**2 + x + 1/3", "QQ"], ["4*x/3 + 1", "QQ"]]]], ["tuple", ["tuple", ["Poly", "Add(Mul(Integer(-1), Rational(3, 2), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Rational(3, 2), Symbol('x')))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["3/2", "QQ"], ["3*x/2 - 1", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Add", "Add(Symbol('x'), Rational(-2, 3))"]]],
    ["is_log_deriv", [["x + 4/3", "QQ"], ["1/2 - x", "QQ"]], null],
    ["limited", [["3/2", "QQ"], ["3*x/2 - 1", "QQ"], [[["4/3 - 4*x/3", "QQ"], ["x**2 + 2*x/3 + 1", "QQ"]]]], null],
    ["limited", [["-x**3/2 + x**2 + 2*x", "QQ"], ["1", "ZZ"], [[["4/3 - 4*x/3", "QQ"], ["x**2 + 2*x/3 + 1", "QQ"]]]], ["tuple", ["tuple", ["Poly", "Add(Mul(Integer(-1), Rational(1, 8), Pow(Symbol('x'), Integer(4))), Mul(Rational(1, 3), Pow(Symbol('x'), Integer(3))), Pow(Symbol('x'), Integer(2)))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["2", "QQ"], ["4*x/3 - 2", "QQ"]], ["tuple", ["Integer", "Integer(2)"], ["Pow", "Pow(Add(Symbol('x'), Rational(-3, 2)), Integer(3))"]]],
    ["is_log_deriv", [["-x**3 - x**2/3 + 2*x + 3", "QQ"], ["-2*x**2 + 4*x/3 - 3", "QQ"]], null],
    ["is_log_deriv", [["-3*x/2 - 16", "QQ"], ["3*x**2/2 + 11*x + 12", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(4, 3)), Integer(-2)), Add(Symbol('x'), Integer(6)))"]]],
    ["is_log_deriv", [["-4*x", "ZZ"], ["2*x/3 - 2/3", "QQ"]], null],
    ["limited", [["1 - 4*x", "ZZ"], ["-1/3", "QQ"], [[["2 - 4*x/3", "QQ"], ["1 - x", "ZZ"]], [["1 - x/2", "QQ"], ["-4*x - 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Add(Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(3), Symbol('x')))", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["6*x - 10/3", "QQ"], ["6*x**2 - 17*x/3 + 4/3", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(-1, 2)), Integer(-1)), Pow(Add(Symbol('x'), Rational(-4, 9)), Integer(2)))"]]],
    ["is_log_deriv", [["7/9 - x/3", "QQ"], ["2*x**2/3 - 10*x/9 + 4/9", "QQ"]], ["tuple", ["Integer", "Integer(2)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Integer(-1)), Integer(4)), Pow(Add(Symbol('x'), Rational(-2, 3)), Integer(-5)))"]]],
    ["is_log_deriv", [["x/3 - 2/3", "QQ"], ["2*x**2 + x - 1", "ZZ"]], ["tuple", ["Integer", "Integer(6)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(-1, 2)), Integer(-1)), Pow(Add(Symbol('x'), Integer(1)), Integer(2)))"]]],
    ["limited", [["1 - 3*x**2", "ZZ"], ["-3", "ZZ"], [[["2*x + 1", "ZZ"], ["4*x**2 - 3*x + 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Add(Mul(Rational(1, 3), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Rational(1, 3), Symbol('x')))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["is_log_deriv", [["27*x**2 - 18*x - 8", "QQ"], ["9*x**3 - 27*x**2/2 + 14*x - 2", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(-1, 6)), Integer(-1)), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(-1), Rational(4, 3), Symbol('x')), Rational(4, 3)), Integer(2)))"]]],
    ["heu", [["27*x**2 - 18*x - 8", "QQ"], ["9*x**3 - 27*x**2/2 + 14*x - 2", "QQ"], ["4*x**2 - 3*x", "ZZ"], ["2", "ZZ"]], ["tuple", ["One", "Integer(1)"], ["Zero", "Integer(0)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(-1, 6)), Integer(-1)), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(-1), Rational(4, 3), Symbol('x')), Rational(4, 3)), Integer(2)))"]]],
    ["heu", [["108*x**5 - 243*x**4 + 579*x**3/2 - 96*x**2 - 18*x - 16", "QQ"], ["18*x**3 - 27*x**2 + 28*x - 4", "QQ"], ["4*x**2 - 3*x", "ZZ"], ["2", "ZZ"]], ["tuple", ["One", "Integer(1)"], ["Integer", "Integer(3)"], ["Mul", "Mul(Pow(Add(Symbol('x'), Rational(-1, 6)), Integer(-1)), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(-1), Rational(4, 3), Symbol('x')), Rational(4, 3)), Integer(2)))"]]],
    ["limited", [["2 - x", "ZZ"], ["4/3", "QQ"], [[["2*x - 3", "ZZ"], ["x**2 + 4*x - 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Add(Mul(Integer(-1), Rational(3, 8), Pow(Symbol('x'), Integer(2))), Mul(Rational(3, 2), Symbol('x')))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
    ["heu", [["82*x**3/3 + 149*x**2/9 + 19*x/9", "QQ"], ["-4*x**4 - 10*x**3/3 - 2*x**2/3", "QQ"], ["-3*x**2 - 2", "ZZ"], ["1/2", "QQ"]], ["tuple", ["Integer", "Integer(6)"], ["Zero", "Integer(0)"], ["Mul", "Mul(Pow(Symbol('x'), Integer(-19)), Pow(Add(Symbol('x'), Rational(1, 3)), Integer(-10)), Pow(Add(Symbol('x'), Rational(1, 2)), Integer(-12)))"]]],
    ["heu", [["36*x**6 + 30*x**5 + 30*x**4 + 101*x**3/3 + 221*x**2/18 + 19*x/18", "QQ"], ["-2*x**4 - 5*x**3/3 - x**2/3", "QQ"], ["-3*x**2 - 2", "ZZ"], ["1/2", "QQ"]], ["tuple", ["Integer", "Integer(6)"], ["Integer", "Integer(18)"], ["Mul", "Mul(Pow(Symbol('x'), Integer(-19)), Pow(Add(Symbol('x'), Rational(1, 3)), Integer(-10)), Pow(Add(Symbol('x'), Rational(1, 2)), Integer(-12)))"]]],
    ["limited", [["0", "QQ"], ["3*x/2 - 1", "QQ"], [[["-2*x**2 - x/2 - 2", "QQ"], ["2*x**2 + 3*x + 2/3", "QQ"]], [["3*x**2 + 1", "ZZ"], ["1 - 2*x", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Integer(0)", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["heu", [["6", "ZZ"], ["1 - 3*x", "ZZ"], ["2*x**2 - x - 1/2", "QQ"], ["1/2", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Zero", "Integer(0)"], ["Pow", "Pow(Add(Symbol('x'), Rational(-1, 3)), Integer(-2))"]]],
    ["heu", [["-12*x**3 + 10*x**2 + x + 2", "QQ"], ["1/2 - 3*x/2", "QQ"], ["2*x**2 - x - 1/2", "QQ"], ["1/2", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Integer", "Integer(2)"], ["Pow", "Pow(Add(Symbol('x'), Rational(-1, 3)), Integer(-2))"]]],
    ["heu", [["-16", "QQ"], ["4*x + 1/3", "QQ"], ["3*x**2 - 4*x/3 - 3", "QQ"], ["-x/2 - 2", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Zero", "Integer(0)"], ["Pow", "Pow(Add(Symbol('x'), Rational(1, 12)), Integer(-4))"]]],
    ["heu", [["36*x**3 - 13*x**2 - 88*x/3 + 29", "QQ"], ["-2*x**2 - 49*x/6 - 2/3", "QQ"], ["3*x**2 - 4*x/3 - 3", "QQ"], ["-x/2 - 2", "QQ"]], ["tuple", ["One", "Integer(1)"], ["Integer", "Integer(3)"], ["Pow", "Pow(Add(Symbol('x'), Rational(1, 12)), Integer(-4))"]]],
    ["limited", [["2", "ZZ"], ["-4/3", "QQ"], [[["-x**2 - x + 3", "ZZ"], ["-4*x**2 + 3*x - 1", "ZZ"]], [["x**2 + 2*x/3 + 3", "QQ"], ["x**2 + x + 1", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Rational(3, 2), Symbol('x'))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["heu", [["1", "ZZ"], ["-4*x - 1", "ZZ"], ["-x**2/2 + x + 3", "QQ"], ["-3*x/2 - 1", "QQ"]], ["tuple", ["Integer", "Integer(4)"], ["Zero", "Integer(0)"], ["Pow", "Pow(Add(Symbol('x'), Rational(1, 4)), Integer(-1))"]]],
    ["heu", [["2*x**3 - 7*x**2/2 - 29*x/2 - 4", "QQ"], ["6*x**2 + 11*x/2 + 1", "QQ"], ["-x**2/2 + x + 3", "QQ"], ["-3*x/2 - 1", "QQ"]], ["tuple", ["Integer", "Integer(4)"], ["Integer", "Integer(4)"], ["Pow", "Pow(Add(Symbol('x'), Rational(1, 4)), Integer(-1))"]]],
    ["limited", [["0", "ZZ"], ["-x", "ZZ"], [[["-1", "ZZ"], ["3 - 3*x", "ZZ"]], [["-2", "ZZ"], ["-x**2/2 + x + 4", "QQ"]]]], ["tuple", ["tuple", ["Poly", "Integer(0)", "ZZ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["1", "ZZ"], ["2", "ZZ"], [[["1", "ZZ"], ["-2*x**2 - 4*x + 3", "ZZ"]], [["2", "ZZ"], ["x - 3", "ZZ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Rational(1, 2), Symbol('x'))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"], ["Poly", "Integer(0)", "QQ"]]]],
    ["limited", [["-1/3", "QQ"], ["3", "ZZ"], [[["3", "ZZ"], ["4*x**2/3 + 4*x - 1", "QQ"]]]], ["tuple", ["tuple", ["Poly", "Mul(Integer(-1), Rational(1, 9), Symbol('x'))", "QQ"], ["Poly", "Integer(1)", "QQ"]], ["list", ["Poly", "Integer(0)", "QQ"]]]],
]

# Calls answered by the stock Python body (outside the native domain).
DECLINED = {
    # ZZ[log(2)] coefficients: not a rational base field.
    ("limited", '[["-2 + log(2)", "ZZ[log(2)]"], ["1", "ZZ"], [[["1", "ZZ"], ["x", "ZZ"]]]]'),
}


def _de():
    return DifferentialExtension(extension={'D': [Poly(1, x)]})


def _poly(p):
    return Poly(sympify(p[0], locals={'x': x}), x, domain=p[1])


def _enc(v):
    if v is None:
        return None
    if isinstance(v, Poly):
        return ['Poly', srepr(v.as_expr()), str(v.get_domain())]
    if isinstance(v, (tuple, list)):
        return [type(v).__name__] + [_enc(i) for i in v]
    return [type(v).__name__, srepr(v)]


_KEY = {'is_log_deriv': 'is_log_deriv_k_t_radical_in_field',
        'heu': 'parametric_log_deriv_heu', 'limited': 'limited_integrate'}


def _run(kind, args, DE):
    if kind == 'is_log_deriv':
        return is_log_deriv_k_t_radical_in_field(_poly(args[0]), _poly(args[1]), DE)
    if kind == 'heu':
        return parametric_log_deriv_heu(*[_poly(a) for a in args], DE)
    G = [(_poly(a), _poly(d)) for a, d in args[2]]
    return limited_integrate(_poly(args[0]), _poly(args[1]), G, DE)


@pytest.mark.parametrize('case', CASES, ids=[str(i) for i in range(len(CASES))])
def test_matches_sympy_reference(case):
    kind, args, expected = case
    import json
    declined = (kind, json.dumps(args)) in DECLINED
    key = _KEY[kind]
    before = _prde_native.STATS[key]
    assert _enc(_run(kind, args, _de())) == expected
    assert _prde_native.STATS[key] == before + (0 if declined else 1)


def test_reference_examples():
    # sympy/integrals/tests/test_prde.py::test_limited_integrate
    DE = _de()
    before = _prde_native.STATS['limited_integrate']
    G = [(Poly(x, x), Poly(x + 1, x))]
    assert limited_integrate(Poly(-(1 + x + 5*x**2 - 3*x**3), x),
                             Poly(1 - x - x**2 + x**3, x), G, DE) == \
        ((Poly(x**2 - x + 2, x), Poly(x - 1, x, domain='QQ')), [2])
    G = [(Poly(1, x), Poly(x, x))]
    assert limited_integrate(Poly(5*x**2, x), Poly(3, x), G, DE) == \
        ((Poly(5*x**3/9, x), Poly(1, x, domain='QQ')), [0])
    assert _prde_native.STATS['limited_integrate'] == before + 2
    # Base-case is_log_deriv_k_t_radical_in_field.
    before = _prde_native.STATS['is_log_deriv_k_t_radical_in_field']
    assert is_log_deriv_k_t_radical_in_field(Poly(5, x), Poly(2*x, x), DE) == (2, x**5)
    assert is_log_deriv_k_t_radical_in_field(Poly(3*x + 1, x), Poly(x**2 + x, x), DE) == \
        (1, x*(x + 1)**2)
    assert is_log_deriv_k_t_radical_in_field(Poly(1, x), Poly(x**2 - 2, x), DE) is None
    assert is_log_deriv_k_t_radical_in_field(Poly(1, x), Poly(x**2, x), DE) is None
    assert _prde_native.STATS['is_log_deriv_k_t_radical_in_field'] == before + 4


def test_parametric_log_deriv_delegates_to_native_heu():
    DE = _de()
    before = _prde_native.STATS['parametric_log_deriv_heu']
    assert parametric_log_deriv(Poly(5, x), Poly(2*x, x), Poly(1, x), Poly(1, x), DE) is None
    assert _prde_native.STATS['parametric_log_deriv_heu'] == before + 1


def test_explicit_dummies_and_cases_use_stock_body():
    DE = _de()
    keys = ('is_log_deriv_k_t_radical_in_field', 'parametric_log_deriv_heu')
    before = [_prde_native.STATS[k] for k in keys]
    z = Symbol('z')
    assert is_log_deriv_k_t_radical_in_field(Poly(5, x), Poly(2*x, x), DE, z=z) == (2, x**5)
    c1 = Symbol('c1')
    assert parametric_log_deriv_heu(Poly(5, x), Poly(2*x, x), Poly(1, x), Poly(1, x),
                                    DE, c1=c1) is None
    assert [_prde_native.STATS[k] for k in keys] == before


def test_non_base_extension_uses_stock_body():
    # sympy/integrals/tests/test_prde.py::test_parametric_log_deriv
    # (primitive case).
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    before = _prde_native.STATS['parametric_log_deriv_heu']
    assert parametric_log_deriv_heu(Poly(5*t**2 + t - 6, t), Poly(2*x*t**2, t),
                                    Poly(-1, t), Poly(x*t**2, t), DE) == (2, 6, t*x**5)
    assert _prde_native.STATS['parametric_log_deriv_heu'] == before


def test_limited_integrate_without_unknown_polynomials():
    # param_rischDE where param_poly_rischDE returns no polynomial (h == []):
    # C is eye(m) and limited_integrate returns None (SymPy 1.14).
    DE = _de()
    before = _prde_native.STATS['limited_integrate']
    G = [(Poly(Rational(1, 2), x), Poly(4*x**2 - 1, x)), (Poly(4, x), Poly(x + 1, x))]
    assert limited_integrate(Poly(-x**2/3 + 2*x/3 - Rational(2, 3), x), Poly(2*x**2, x),
                             G, DE) is None
    assert _prde_native.STATS['limited_integrate'] == before + 1


def test_splitfactor_precheck_filters_non_rational_domains():
    precheck = _core._risch_splitfactor_precheck
    one = Poly(1, x, domain='QQ')
    assert precheck(Poly(x**2 + x, x, domain='QQ'), one) is True
    assert precheck(Poly(x**2 + x, x, domain='ZZ'), one) is False
    assert precheck(Poly(x**2 + x, x, domain='QQ'), Poly(1, x, domain='ZZ')) is False
    assert precheck(Poly(x**2 + x, x, domain='QQ'), Poly(1, t, domain='QQ')) is False
    assert precheck(Poly(x**2 + log(2)*x, x), one) is False
    if _sf._ROOT is not None:
        assert _sf._ROOT['precheck'] is precheck
    # Results are unchanged on either side of the filter.
    DE = _de()
    assert splitfactor(Poly(x**2 + x, x, domain='ZZ'), DE) == \
        (Poly(x**2 + x, x, domain='ZZ'), Poly(1, x, domain='ZZ'))
    assert splitfactor(Poly(x**2 + x, x, domain='QQ'), DE) == \
        (Poly(x**2 + x, x, domain='QQ'), Poly(1, x, domain='QQ'))
