"""Core internals imported by Devito and TorchInductor, compared with SymPy 1.14.

Each probe runs the same code against ``symox`` and the reference ``sympy``
package and compares the printed results.
"""
import importlib
import warnings

import pytest
import symox as S
import sympy as R


def _module(backend, name):
    return importlib.import_module(backend.__name__ + '.' + name)


def _same(probe):
    native, reference = probe(S), probe(R)
    assert native == reference
    return native


def test_assumption_rule_tables_match_reference():
    def probe(m):
        a = _module(m, 'core.assumptions')
        rules = a._assume_rules
        generated = a._generate_assumption_rules()
        # The tables are defaultdicts: lookups may add empty entries.
        def table(t):
            return sorted((k, sorted(v)) for k, v in t.items() if v)
        return (sorted(rules.defined_facts), sorted(a._assume_defined),
                table(rules.full_implications), table(rules.prereq), table(rules.beta_triggers),
                [(sorted(c), v) for c, v in rules.beta_rules],
                table(generated.full_implications) == table(rules.full_implications))
    assert _same(probe)[-1]


@pytest.mark.parametrize('facts', [
    {'integer': True, 'positive': True}, {'zero': True}, {'prime': True},
    {'real': True, 'positive': False, 'zero': False}, {'composite': False, 'even': True, 'positive': True},
    {'extended_real': True, 'finite': False}, {'imaginary': True}])
def test_std_fact_kb_deductions_match_reference(facts):
    def probe(m):
        kb = _module(m, 'core.assumptions').StdFactKB(facts)
        return sorted(kb.items()), sorted(kb.generator.items()), sorted(kb.copy().items())
    _same(probe)


@pytest.mark.parametrize('facts', [{'integer': True, 'rational': False},
                                   {'positive': True, 'negative': True}, {'zero': True, 'odd': True}])
def test_std_fact_kb_rejects_inconsistent_facts_like_reference(facts):
    def probe(m):
        a, f = _module(m, 'core.assumptions'), _module(m, 'core.facts')
        with pytest.raises(f.InconsistentAssumptions) as info:
            a.StdFactKB(facts)
        return isinstance(info.value, ValueError), str(info.value)
    _same(probe)


def test_class_assumption_preparation_and_ask_match_reference():
    def probe(m):
        a = _module(m, 'core.assumptions')
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            class Example(metaclass=a.ManagedProperties):
                is_positive = True
                def __init__(self):
                    self._assumptions = type(self).default_assumptions
                def _eval_is_integer(self):
                    return True
        value = Example()
        answers = [value.is_rational, value.is_even, value.is_negative, value.is_prime]
        return ([type(w.message).__name__ for w in caught],
                sorted(Example._explicit_class_assumptions.items()),
                sorted(Example.default_assumptions.items()), sorted(Example._prop_handler),
                answers, sorted(value._assumptions.items()),
                value._assumptions is Example.default_assumptions,
                a.as_property('real'))
    _same(probe)


@pytest.mark.parametrize('text', ['a & b', 'a | !b', '!a', 'a & b & !c', 'a | b | c'])
def test_logic_fromstring_and_simplification_match_reference(text):
    def probe(m):
        L = _module(m, 'core.logic')
        expr = L.Logic.fromstring(text)
        return (type(expr).__name__, str(expr), str(L.Not(expr)),
                str(L.And(expr, L.Not('a'))), str(L.Or(expr, 'a')),
                L.fuzzy_and([True, None]), L.fuzzy_or([False, None]), L.fuzzy_not(None))
    _same(probe)


def test_sympifyit_and_priority_decorators_match_reference():
    def probe(m):
        d = _module(m, 'core.decorators')

        class Probe:
            _op_priority = 1.0

            @d._sympifyit('other', NotImplemented)
            def __add__(self, other):
                return ('add', str(other))

            @d._sympifyit('other')
            def strict(self, other):
                return str(other)

            @d.call_highest_priority('__rmul__')
            def __mul__(self, other):
                return 'mine'

        class High:
            _op_priority = 100.0

            def __rmul__(self, other):
                return 'theirs'

        p = Probe()
        out = [p.__add__(2), p.__add__(object()) is NotImplemented, p * 2, p * High(), p.strict(2), p.__add__('x')]
        with pytest.raises(m.SympifyError):
            p.strict(object())
        return out
    _same(probe)


def test_class_ordering_and_basic_compare_match_reference():
    def probe(m):
        b = _module(m, 'core.basic')
        x, y = m.symbols('x y')
        pairs = [(m.Symbol, m.Integer), (m.Add, m.Mul), (m.Pow, m.Mul), (m.Integer, m.Rational)]
        values = [y, x, m.Integer(2), m.Rational(1, 3), x + y, x * y, x ** 2, m.sin(x), m.pi, m.Float(2.5), m.Tuple(x, 1)]
        return (list(b.ordering_of_classes),
                [b._cmp_name(p, q) for p, q in pairs],
                [str(v) for v in sorted(values, key=b._args_sortkey)],
                [x.compare(y), y.compare(x), (x + 1).compare(x * 2), m.Integer(2).compare(m.Rational(1, 2))])
    _same(probe)


def test_operation_categories_match_reference():
    def probe(m):
        o = _module(m, 'core.operations')
        x, y = m.symbols('x y')
        a, b = m.symbols('a b')
        values = [x + y, x * y, m.Max(x, y), m.Min(x, y), m.And(a, b), m.Or(a, b), x ** y, m.sin(x),
                  m.Union(m.Interval(0, 1), m.Interval(2, 3))]
        classes = [m.Add, m.Mul, m.Pow, m.Symbol, m.Union, m.Intersection]
        return ([(isinstance(v, o.AssocOp), isinstance(v, o.LatticeOp)) for v in values],
                [(issubclass(c, o.AssocOp), issubclass(c, o.LatticeOp)) for c in classes],
                sorted(map(str, m.Add.make_args(x + y + 1))), [str(v) for v in m.Mul.make_args(x + y)],
                sorted(map(str, o.LatticeOp.make_args(m.Max(x, y)))), sorted(map(str, o.LatticeOp.make_args(x))),
                issubclass(o.ShortCircuit, Exception), o.LatticeOp.is_commutative)
    _same(probe)


def test_global_parameters_follow_evaluate_like_reference():
    def probe(m):
        p = _module(m, 'core.parameters')
        out = [p.global_parameters.evaluate, p.global_parameters.distribute, p.global_parameters.exp_is_pow]
        with m.evaluate(False):
            out.append(p.global_parameters.evaluate)
            out.append(str(m.Add(m.Integer(1), m.Integer(2))))
        out.append(p.global_parameters.evaluate)
        with p.distribute(True):
            out.append(p.global_parameters.distribute)
        return out
    assert _same(probe) == [True, True, False, False, '1 + 2', True, True]


def test_node_counts_and_ordered_match_reference():
    def probe(m):
        s = _module(m, 'core.sorting')
        x, y = m.symbols('x y')
        items = [m.Derivative(m.sin(x), (x, 2)), m.nan, [x, [y, 1]], {x: y + 1}, x + y, m.Tuple(x, 2.5),
                 m.ImmutableMatrix([[x, 1], [2, y]]), m.ImmutableSparseMatrix(2, 2, {(0, 0): x}), m.Float(1.5)]
        exprs = [x ** 2, y, x, m.Integer(1), x * y, m.sin(x) + 1, m.Rational(1, 2)]
        return [s._nodes(i) for i in items], [str(e) for e in s.ordered(exprs)]
    _same(probe)


def test_boolean_atom_classes_match_reference():
    def probe(m):
        b = _module(m, 'logic.boolalg')
        x = m.Symbol('x')
        return ([isinstance(v, b.BooleanAtom) for v in (m.true, m.false, x, m.Integer(1), m.And(x, m.Symbol('y')))],
                isinstance(m.true, b.BooleanTrue), isinstance(m.false, b.BooleanTrue),
                isinstance(m.false, b.BooleanFalse), b.BooleanTrue() is m.true, b.BooleanFalse() is m.false,
                b.BooleanAtom.is_Boolean, b.BooleanAtom.is_Atom, b.BooleanAtom._op_priority)
    _same(probe)


def test_add_mul_construction_helpers_match_reference():
    def probe(m):
        add, mul = _module(m, 'core.add'), _module(m, 'core.mul')
        I = m.Integer
        x, y, z = m.symbols('x y z')
        out = []
        terms = [z, x * y, m.sin(x), I(2), x ** 2, y, m.Rational(1, 2)]
        add._addsort(terms)
        out.append(terms)
        factors = [z, x * y, m.sin(x), I(2), x ** 2, y]
        mul._mulsort(factors)
        out.append(factors)
        out.append(m.srepr(add._unevaluated_Add(m.Float(1.0), x, I(2)).args[1:]))
        out.append(str(add._unevaluated_Add(m.Float(1.0), x, I(2)).args[0]))
        out.append(m.srepr(add._unevaluated_Add(x + 1, x + 2)))
        out.append(m.srepr(mul._unevaluated_Mul(m.Float(3.0), x, I(2)).args[1:]))
        out.append(str(mul._unevaluated_Mul(m.Float(3.0), x, I(2)).args[0]))
        out.append(m.srepr(mul._unevaluated_Mul(m.sqrt(2), m.sqrt(3))))
        out.append(m.srepr(mul._unevaluated_Mul(I(2), x, m.Mul(3, y, evaluate=False))))
        for coeff, factors, clear, sign in [(m.S.Half, x + 2, True, False), (m.S.Half, x + 2, False, False),
                                            (m.S.Half, (x + 2) * y, False, False), (I(-1), x + y, True, False),
                                            (I(-1), x + y, True, True), (I(2), 3 * x, True, False),
                                            (x, I(3), True, False), (I(3), I(4), True, False),
                                            (I(2), m.S.One, True, False), (I(2), m.sin(x), True, False)]:
            out.append(m.srepr(mul._keep_coeff(coeff, factors, clear=clear, sign=sign)))
        out.append(m.srepr(m.Add._from_args([y, x, I(2)])))
        out.append(m.srepr(m.Mul._from_args((I(1), x))))
        out.append((m.Add._from_args([]), m.Mul._from_args([]), m.Add._from_args([x]),
                    m.Add.identity, m.Mul.identity, m.Add._op_priority, m.Mul._op_priority))
        return [str(o) for o in out]
    _same(probe)


def test_evalf_table_keys_match_reference():
    def probe(m):
        e = _module(m, 'core.evalf')
        e.evalf_table.clear()
        e._create_evalf_table()
        names = sorted(getattr(k, '__name__', str(k)) for k in e.evalf_table)
        x = m.Symbol('x')
        n = _module(m, 'core.numbers')
        constants = [str(e.evalf_table[k](None, 15, {})) for k in (n.Zero, n.Half, n.Pi, n.Exp1, n.NegativeOne, n.One)
                     if k in e.evalf_table]
        aliases = [e.evalf_table[k].__name__ for k in (m.Add, m.Mul, m.Pow)]
        return names, constants, aliases, str(m.sin(x).subs(x, 1).evalf(20))
    names, constants, _, _ = _same(probe)
    assert len(names) == 33 and len(constants) == 6


def test_downstream_module_paths_match_reference():
    def probe(m):
        poly = _module(m, 'polys.polytools')
        missing = [n for n in R.polys.polytools.__all__ if not hasattr(poly, n)]
        core = _module(m, 'core.core')
        registry = type('Probe', (core.Registry,), {})()
        registry.value = 3
        out = [missing, core.Registry.__slots__, registry.value]
        del registry.value
        out.append(hasattr(registry, 'value'))
        cache = _module(m, 'polys.domains.modularinteger')._modular_integer_cache
        cache.clear()
        out.append(type(cache).__name__)
        matrices = _module(m, 'matrices.matrices')
        x = m.Symbol('x')
        out.append(isinstance(m.Matrix([[1, x]]), matrices.MatrixBase))
        out.append(issubclass(matrices.NonSquareMatrixError, Exception))
        out.append(_module(m, 'polys.rootoftools').ComplexRootOf.clear_cache())
        for name in ('integers', 'miscellaneous'):
            mod = _module(m, 'functions.elementary.' + name)
            out.append([str(getattr(mod, f)(m.Rational(7, 2))) for f in
                        (('floor', 'ceiling', 'frac') if name == 'integers' else ('sqrt', 'cbrt'))])
        kind = _module(m, 'core.kind')
        out.append([str(k) for k in (kind.UndefinedKind, kind.NumberKind, kind.BooleanKind)])
        return out
    assert _same(probe)[0] == []
