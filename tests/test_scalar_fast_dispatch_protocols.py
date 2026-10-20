"""Common scalar dispatch retains foreign and context-sensitive protocols."""
import pytest
import symox as sx


def test_foreign_type_equality_is_not_used_to_classify_numeric_operands():
    class Meta(type):
        def __eq__(self, other):
            raise AssertionError('type equality is not scalar classification')
    class Foreign(metaclass=Meta):
        def _symox_(self):
            return sx.Integer(2)
    x = sx.Symbol('x')
    assert x+Foreign() == x+2
    assert x*Foreign() == 2*x


def test_numeric_subclass_reflected_priority_is_preserved():
    class PriorityInteger(int):
        _op_priority = 100
        def __rmul__(self, other):
            return ('priority', other)
        def __rpow__(self, other):
            return ('power', other)
    x = sx.Symbol('x')
    assert x*PriorityInteger(2) == ('priority', x)
    assert x**PriorityInteger(2) == ('power', x)


def test_held_arithmetic_and_indeterminate_division_remain_context_local():
    x = sx.Symbol('x')
    with sx.evaluate(False):
        assert (x+x).args == (x,x)
        assert (x**1).args == (x,sx.Integer(1))
    assert x+x == 2*x
    sx.seterr(divide=True)
    try:
        with pytest.raises(ValueError, match='0/0'):
            sx.Integer(0)/0
    finally:
        sx.seterr(divide=False)
    assert sx.Integer(0)/0 is sx.S.NaN
