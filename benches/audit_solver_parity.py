#!/usr/bin/env python3
"""Focused solver parity audit against SymPy (not the repository test suite).

Checks whole finite solution lists, exact polynomial defining data, real-domain
filters, pole exclusions, and multiprecision LambertW values. Both libraries
must be installed. Run with --out to retain the machine-readable report.
"""
import argparse
import json
import pickle
import time
from pathlib import Path

import symox as sx
import sympy as sp


def normalized(value):
    value = value.to_sympy() if hasattr(value, 'to_sympy') else sp.sympify(value)
    x = sp.Symbol('x')
    return value.xreplace({root: sp.CRootOf(sp.Poly.from_list(root.poly.all_coeffs(), x), root.index)
                          for root in value.atoms(sp.CRootOf)})


def equivalent(a, b, approximate=False):
    if a == b:
        return True
    difference = sp.simplify(a-b)
    if difference == 0:
        return True
    return approximate and difference.is_number and abs(complex(difference.evalf(30))) < 1e-12


def same_solutions(a, b, approximate=False):
    if len(a) != len(b):
        raise AssertionError(f'root count: symox {len(a)}, SymPy {len(b)}; {a} versus {b}')
    remaining = list(b)
    for value in a:
        for i, candidate in enumerate(remaining):
            if equivalent(value, candidate, approximate):
                remaining.pop(i)
                break
        else:
            raise AssertionError(f'unmatched solution {value}; remaining reference {remaining}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--only', action='append', help='Run checks whose name contains this text; repeatable.')
    args = parser.parse_args()
    report = {'reference': sp.__version__, 'symox': sx.__version__, 'checks': []}

    def check(name, fn):
        if args.only and not any(pattern in name for pattern in args.only):
            return
        start = time.monotonic()
        try:
            details = fn()
            row = {'name': name, 'status': 'passed'}
            if details:
                row.update(details)
        except Exception as error:
            row = {'name': name, 'status': 'failed', 'error': f'{type(error).__name__}: {error}'}
        row['seconds'] = time.monotonic()-start
        report['checks'].append(row)
        print(name, row['status'], row.get('error', ''), flush=True)

    fixtures = [
        ('quintic', 'x**5-x-1', {}),
        ('septic', 'x**7-x-1', {}),
        ('quintic real', 'x**5-x-1', {'real': True}),
        ('quintic positive', 'x**5-x-1', {'positive': True}),
        ('quintic negative', 'x**5-x-1', {'negative': True}),
        ('quintic integer', 'x**5-x-1', {'integer': True}),
        ('multiple real roots', 'x**5-5*x+1', {'real': True}),
        ('repeated factors', '(x**5-x-1)**2*(x-2)', {}),
        ('mixed radicals and root objects', '(x**5-x-1)*(x**2-2)', {}),
        ('rational polynomial', 'x**5/3-x/3-sp.Rational(1,3)', {}),
        ('exact quadratic cancellation', 'x**2-10**20*x+1', {}),
        ('radical exclusion', 'sqrt(x)+1', {}),
        ('principal logarithm exclusion', 'log(x)-4*sp.I', {}),
        ('principal arcsine exclusion', 'sp.asin(x)-2', {}),
        ('principal logarithm valid value', 'log(x)-3*sp.I', {}),
        ('Lambert positive', 'x*exp(x)-1', {}),
        ('Lambert two real branches', 'x*exp(x)+sp.Rational(1,10)', {}),
        ('Lambert complex', 'x*exp(x)+1', {}),
        ('Lambert no real roots', 'x*exp(x)+1', {'real': True}),
        ('Lambert affine', '(2*x+3)*exp(4*x+5)-7', {}),
        ('Lambert power two', 'x**2*exp(x)-1', {}),
        ('Lambert power extra real branch', '2*x**2*exp(x)-1', {}),
        ('Lambert power three', 'x**3*exp(x)-1', {}),
        ('Lambert exponential plus linear', 'x+exp(x)', {}),
        ('Lambert logarithm', 'x+log(x)', {}),
        ('Lambert affine logarithm', '2*x+3*log(x)-1', {}),
        ('Lambert shifted logarithm', 'log(2*x+3)-x', {}),
        ('Lambert complex logarithm', 'log(x)-x', {}),
        ('Lambert float input', 'x*exp(x)-0.1', {}),
        ('finite solutions only', 'x**2*exp(x)', {}),
        ('exponential has no zero', 'exp(x)', {}),
    ]
    for name, source, assumptions in fixtures:
        def run(source=source, assumptions=assumptions):
            values = []
            for library in (sx, sp):
                x = library.Symbol('x', **assumptions)
                namespace = {'sp': library, 'x': x, 'sqrt': library.sqrt, 'exp': library.exp, 'log': library.log}
                expression = eval(source, {'__builtins__': {}}, namespace)
                values.append([normalized(v) for v in library.solve(expression, x)])
            if name == 'quintic integer':
                # x^5-x-1 is 1 mod 2 for both possible residue classes.
                # This certifies that SymPy's returned nonintegral CRootOf
                # cannot satisfy the integer assumption.
                assert values[0] == [], values[0]
                assert all((r**5-r-1) % 2 == 1 for r in range(2))
                if not values[1]:
                    return None
                return {'comparison': 'proved_stronger_domain_filter',
                        'proof': 'No solution modulo 2, hence no integer solution.',
                        'reference_count': len(values[1])}
            if name == 'Lambert power three':
                # W(z)*exp(W(z))=z gives (3*W(z))^3*exp(3*W(z))=27*z^3.
                # All three n-th roots are valid; SymPy lists only one here.
                assert len(values[0]) == 3
                for candidate in values[0]:
                    ws = candidate.atoms(sp.LambertW)
                    assert len(ws) == 1
                    w, = ws
                    assert candidate == 3*w
                    assert sp.expand(27*w.args[0]**3-1) == 0
                assert all(any(equivalent(a,b) for a in values[0]) for b in values[1])
                if len(values[0]) == len(values[1]):
                    return None
                return {'comparison': 'proved_additional_solutions',
                        'proof': 'W identity and exact cube roots of 1/27.',
                        'additional_count': len(values[0])-len(values[1])}
            same_solutions(*values, approximate='0.1' in source)
        check(name, run)

    def root_protocols():
        x = sx.Symbol('x')
        f = x**5-x-1
        roots = sx.solve(f, x)
        assert len(roots) == 5
        for i, root in enumerate(roots):
            assert root == sx.CRootOf(f, i), (root, sx.CRootOf(f, i))
            assert not root.free_symbols
            assert root.is_algebraic is True
            assert root.is_rational is False
            assert root.is_real is (i == 0)
            assert pickle.loads(pickle.dumps(root)) == root
            assert sx.sympify(root.to_sympy()) == root
    check('root identity, facts, pickle and SymPy interchange', root_protocols)

    def root_numerics():
        x = sx.Symbol('x')
        roots = sx.solve(x**5-x-1,x)
        reference = sp.Symbol('x')
        for i, root in enumerate(roots):
            actual = normalized(root.evalf(60))
            expected = sp.CRootOf(reference**5-reference-1,i).evalf(80)
            assert abs((actual-expected).evalf(70)) < sp.Pow(10,-56), (root,actual,expected)
    check('all quintic roots at 60 digits', root_numerics)

    def root_system_and_flags():
        x,y = sx.symbols('x y')
        equation = x**5-x-1
        roots = sx.solve(equation,x)
        assert sx.solve(equation,x,dict=True) == [{x:root} for root in roots]
        variables, rows = sx.solve(equation,x,set=True)
        assert variables == [x] and rows == {(root,) for root in roots}
        result = sx.solve([equation,y-x],[x,y])
        assert len(result) == 5
        assert set(result) == {(root,root) for root in roots}
    check('quintic in triangular system and output flags', root_system_and_flags)

    def denominator_exclusion():
        for library in (sx, sp):
            x = library.Symbol('x')
            p = x**5-x-1
            numerator = library.expand(p*(x-2))
            denominator = library.expand(p*(x-3))
            result = library.solve(numerator/denominator, x)
            assert result == [2], result
    check('algebraic poles excluded exactly', denominator_exclusion)

    def symbolic_parameter():
        for name in ('a', '_u'):
            results = []
            for library in (sx, sp):
                x, a = library.symbols('x '+name)
                results.append([normalized(v) for v in library.solve(x*library.exp(x)-a, x)])
            same_solutions(*results)
    check('symbolic parameter and temporary-symbol collision', symbolic_parameter)

    def complete_sets():
        x,a = sx.symbols('x a')
        for value in (sx.S.One,sx.E):
            result = sx.solveset(x*sx.exp(x)-value,x,domain=sx.S.Complexes)
            assert result.class_name == 'ConditionSet', result
        result = sx.solveset(x*sx.exp(x)+a,x,domain=sx.S.Reals)
        assert result.class_name == 'ConditionSet', result
        result = sx.solveset(x*sx.exp(x)+sx.Rational(1,10),x,domain=sx.S.Reals)
        assert result == sx.FiniteSet(sx.LambertW(-sx.Rational(1,10)),sx.LambertW(-sx.Rational(1,10),-1)),result
    check('solveset never reports partial Lambert families as complete', complete_sets)

    for digits in (50, 100, 300):
        for argument, branch in [('1', 0), ('-1/10', 0), ('-1/10', -1), ('-1', 0), ('-1', -1), ('1', 2)]:
            def numerical(digits=digits, argument=argument, branch=branch):
                actual = sx.LambertW(sx.sympify(argument), branch).evalf(digits)
                converted = normalized(actual)
                assert not converted.has(sp.LambertW), actual
                reference = sp.LambertW(sp.sympify(argument), branch).evalf(digits+20)
                error = abs((converted-reference).evalf(digits+10))
                scale = max(sp.S.One, abs(reference))
                assert error < scale*sp.Pow(10, -digits+3), (error, digits, actual)
            check(f'LambertW({argument},{branch}) at {digits} digits', numerical)

    stress = [
        ('1/2+I/10**20', 0), ('-3/10+I/10**20', 0),
        ('-3/10-I/10**20', 0), ('-1/5+I/5', 0),
        ('-3/10+I/10**20', -1), ('-3/10-I/10**20', -1),
        ('-exp(-1)+1/10**80', 0), ('-exp(-1)-1/10**80', 0),
        ('-exp(-1)+1/10**80', -1), ('-exp(-1)-1/10**80', -1),
        ('10**1000', 0), ('-1/10**1000', -1),
    ]
    for argument, branch in stress:
        def numerical_stress(argument=argument, branch=branch):
            actual = normalized(sx.LambertW(sx.sympify(argument), branch).evalf(100))
            assert not actual.has(sp.LambertW), actual
            reference = sp.LambertW(sp.sympify(argument), branch).evalf(240)
            # Compare components separately: a small imaginary component must
            # retain its requested significant digits too.
            for a,b in zip(actual.as_real_imag(),reference.as_real_imag()):
                error = abs((a-b).evalf(120))
                scale = abs(b) if b else sp.S.One
                assert error < scale*sp.Pow(10,-96), (argument,branch,error,a,b)
        check(f'branch/scale stress LambertW({argument},{branch})', numerical_stress)

    failed = [row for row in report['checks'] if row['status'] != 'passed']
    report['passed'] = len(report['checks'])-len(failed)
    report['failed'] = len(failed)
    if args.out:
        args.out.write_text(json.dumps(report, indent=2)+'\n')
    print(f"{report['passed']} passed, {report['failed']} failed", flush=True)
    raise SystemExit(bool(failed))


if __name__ == '__main__':
    main()
