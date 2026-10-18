"""Scope normalization, local complex regularity and asymptotic expansions.

All algebra, substitution, limits and finite series run through the native
expression engine. Undecidable local regularity retains a three-valued result.
"""
import operator
from functools import partial
from math import gcd

from . import _core


def canonical_variables(self):
    from . import Symbol
    bound = self.variables if self.is_Subs else getattr(self, 'bound_symbols', ())
    reserved = {v.name for v in self.free_symbols - set(bound)}
    replacements, index = {}, 0
    for symbol in bound:
        while '_' + str(index) in reserved:
            index += 1
        replacements[symbol] = Symbol('_' + str(index))
        index += 1
    return replacements


def as_dummy(self):
    from . import Dummy
    if self.is_Symbol:
        return Dummy(self.name)
    memo = {}

    def visit(node):
        if node in memo:
            return memo[node]
        children = {arg: visit(arg) for arg in node.args if arg.args}
        children = {old: new for old, new in children.items() if old != new}
        result = node.xreplace(children) if children else node
        bound = set(getattr(result, 'bound_symbols', ()))
        if bound:
            # A binding variable can also occur freely in its own endpoint.
            # Native subs masks only the free occurrences before xreplace
            # renames the bound body and the variable declarations together.
            masks = {v: Dummy() for v in bound & result.free_symbols}
            if masks:
                result = result.subs(masks, simultaneous=True)
            result = result.xreplace(canonical_variables(result))
            if masks:
                result = result.xreplace({v: k for k, v in masks.items()})
        memo[node] = result
        return result

    return visit(self)


class FromIter:
    """Use the concrete expression head when accessed through an instance."""
    def __get__(self, instance, owner):
        if instance is not None:
            return partial(_fromiter, instance.func)
        return partial(_generic_fromiter, owner)


def _fromiter(constructor, args, **assumptions):
    return constructor(*tuple(args), **assumptions)


def _generic_fromiter(owner, args, **assumptions):
    if assumptions:
        raise TypeError('Expr.fromiter does not accept constructor assumptions')
    return _core._apply('Expr', *tuple(args))


def extract_branch_factor(self, allow_half=False):
    from . import Add, Mul, I, pi, S, ceiling, exp_polar
    ordinary, exponents = [], []
    for factor in Mul.make_args(self):
        if factor.is_Pow:
            base, power = factor.args
            if base.is_Function and base.name == 'exp_polar' and not base.is_UndefinedFunction:
                exponents.append(base.args[0]*power)
                continue
        if factor.is_Function and factor.name == 'exp_polar' and not factor.is_UndefinedFunction:
            exponents.append(factor.args[0])
        else:
            ordinary.append(factor)
    if not exponents:
        return self, S.Zero
    exponent = Add(*exponents).expand()
    multiple, extra = S.Zero, []
    for term in Add.make_args(exponent):
        coefficient = term.as_coefficient(I*pi)
        if coefficient is None:
            extra.append(term)
        else:
            multiple += coefficient
    if multiple.is_number:
        constant, tail = multiple, ()
    else:
        constant, tail = multiple.as_coeff_add(*multiple.free_symbols)
    turns = ceiling(constant/2 - S.Half)
    residual = constant - 2*turns
    if allow_half:
        difference = residual.extract_additively(1)
        if difference is not None:
            turns += S.Half
            residual = difference
    remaining = I*pi*Add(residual, *tail) + Add(*extra)
    if remaining != 0:
        ordinary.append(exp_polar(remaining))
    return Mul(*ordinary), turns


def _group(values):
    values = list(values)
    if None in values or values.count(False) > 1:
        return None
    return not any(v is False for v in values)


def is_meromorphic(self, x, a):
    from . import S, I
    x, a = _core.sympify(x), _core.sympify(a)
    if not x.is_Symbol:
        raise TypeError('%s should be of symbol type' % x)
    memo = {}

    def finite(value):
        answer = value.is_finite
        if answer is None and (value.is_Add or value.is_Mul):
            if all(finite(arg) is True for arg in value.args):
                return True
        return answer

    def local(node):
        if node in memo:
            return memo[node]
        answer = evaluate(node)
        memo[node] = answer
        return answer

    def evaluate(node):
        if node == S.NaN:
            return None
        if node in (S.Infinity, S.NegativeInfinity):
            return False
        if node.is_Add or node.is_Mul:
            return _group(local(arg) for arg in node.args)
        if node.is_Symbol or node == x:
            return True
        if node.is_Pow:
            base, power = node.args
            base_regular = local(base)
            if power.is_Integer:
                return base_regular
            power_regular = local(power)
            if base_regular is False:
                return False if power_regular else None
            if base_regular is None:
                return None
            value = base.subs(x, a)
            if value.is_zero is True or value.is_finite is False:
                return False
            if value.is_zero is None or value.is_finite is None:
                return None
            if power_regular is not True:
                return power_regular
            return power.subs(x, a).is_finite
        if not node.is_Function or node.is_UndefinedFunction:
            return True if x not in node.free_symbols and not node.is_UndefinedFunction else None
        name, args = node.name, node.args
        if len(args) != 1 or local(args[0]) is not True:
            return None
        value = args[0].subs(x, a)
        if value == S.NaN:
            return None
        if name in ('exp', 'exp_polar', 'sin', 'cos', 'tan', 'cot', 'sec', 'csc',
                    'gamma', 'erf', 'erfc'):
            return finite(value)
        if name == 'log':
            if value.is_zero is True or value.is_finite is False:
                return False
            return True if value.is_zero is False and value.is_finite is True else None
        singular = {'asin': (S.One, -S.One, S.Zero), 'acos': (S.One, -S.One, S.Zero),
                    'atan': (I, -I), 'acot': (I, -I)}
        if name in singular:
            if value in singular[name]:
                return False
            if value.is_finite is False:
                return name in ('atan', 'acot')
            if value.is_number and value.is_finite is True:
                return True
        return None

    return local(self)


def _variable(expr, x):
    if x is not None:
        x = _core.sympify(x)
        if not x.is_Symbol:
            raise ValueError('Series variable must be a symbol')
        return x
    symbols = expr.free_symbols
    if len(symbols) > 1:
        raise ValueError('A series variable is required for multivariate expressions')
    from . import Symbol
    return next(iter(symbols)) if symbols else Symbol('x')


def aseries(self, x=None, n=6, bound=0, hir=False):
    """Expand at positive infinity using comparable real exponential scales."""
    from . import Add, Dummy, O, S, exp, log, simplify
    x = _variable(self, x)
    n, bound = operator.index(n), operator.index(bound)
    if n < 1 or n > 256 or bound < 0 or bound > 32:
        raise ValueError('aseries requires 1 <= n <= 256 and 0 <= bound <= 32')
    if bound and self.is_Function and not self.is_UndefinedFunction and self.name == 'exp' and self.args[0].limit(x, S.Infinity) == S.Infinity:
        exponent = self.args[0]
        factors = exponent.args if exponent.is_Mul else (exponent,)
        dominant = next((factor for factor in factors if factor.is_Function and not factor.is_UndefinedFunction and factor.name == 'exp'
                         and factor.args[0].limit(x, S.Infinity) == S.Infinity), S.One)
        coefficient = simplify(exponent/dominant)
        normalized = coefficient.series(x, S.Infinity, bound).removeO()
        pieces = [dominant*term for term in Add.make_args(normalized) if term != 0]
        residual = exponent-Add(*pieces)
        from . import Mul
        return Mul(*(exp(piece) for piece in pieces), exp(residual), evaluate=False)
    pending, seen, scales = [self], set(), []
    while pending:
        node = pending.pop()
        if node in seen:
            continue
        seen.add(node)
        if node.is_Function and not node.is_UndefinedFunction and node.name == 'exp' and x in node.free_symbols:
            if node.args[0].limit(x, S.Infinity) == S.NegativeInfinity:
                scales.append(node)
                continue
        pending.extend(node.args)
    if not scales:
        result = self.series(x, S.Infinity, n)
        def positive_direction(node):
            children = {arg: positive_direction(arg) for arg in node.args if arg.args}
            children = {a: b for a, b in children.items() if a != b}
            node = node.xreplace(children) if children else node
            if node == log(1/x):
                return -log(x)
            if node.is_Pow and node.args[0] == 1/x:
                return x**(-node.args[1])
            return node
        return positive_direction(result)
    def ratio_limit(exponent, base):
        ratio=simplify(exponent/base)
        return ratio if ratio.is_Rational else ratio.limit(x,S.Infinity)

    # Select the most rapidly vanishing comparable scale. Slower exponentials
    # remain coefficient functions and can be expanded recursively unless hir
    # requests their hierarchical form. A nonconstant exponent ratio may have
    # a rational limit; retain its exact residual in the coefficient.
    base = scales[0].args[0]
    for node in scales[1:]:
        ratio=ratio_limit(node.args[0],base)
        if ratio==S.Infinity:
            base=node.args[0]
        elif ratio!=0:
            if not (ratio.is_Rational and ratio.is_positive is True):
                raise NotImplementedError('Asymptotic exponential scales cannot be compared in a rational scale')
            normalized=simplify(node.args[0]/ratio)
            if simplify(normalized-base).limit(x,S.Infinity)==S.Infinity:
                base=normalized
    selected=[]
    for node in scales:
        ratio=ratio_limit(node.args[0],base)
        if ratio==0:continue
        if not ratio.is_Rational or ratio.is_positive is not True:
            raise NotImplementedError('Asymptotic exponential scales require rational leading ratios')
        selected.append((node,ratio))
    denominator = 1
    for _,ratio in selected:
        d = int(ratio.q)
        denominator = denominator*d // gcd(denominator, d)
    base /= denominator
    slot = Dummy('w', positive=True)
    rewritten = self.xreplace({node: slot**int(ratio*denominator)*exp(simplify(node.args[0]-ratio*denominator*base))
                              for node, ratio in selected})
    expansion = rewritten.series(slot, 0, n)
    remainder = expansion.getO()
    polynomial = expansion.removeO().expand()
    terms = polynomial.args if polynomial.is_Add else (polynomial,)
    powers = {}
    for term in terms:
        coefficient, power = term.as_coeff_exponent(slot)
        if slot in coefficient.free_symbols or not power.is_Integer:
            raise NotImplementedError('Asymptotic scale is not a Laurent monomial')
        powers[int(power)] = powers.get(int(power), S.Zero) + coefficient
    result, coefficient_order = S.Zero, False
    for power, coefficient in sorted(powers.items()):
        if not hir and x in coefficient.free_symbols:
            coefficient = aseries(coefficient, x, n, max(0, bound-1), False)
        result += coefficient * exp(base*power) if power else coefficient
        if coefficient.getO() is not None:
            coefficient_order = True
            break
    if remainder is not None and not coefficient_order:
        result += O(exp(base*n), (x, S.Infinity))
    return result
