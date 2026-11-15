# Subclassing

[README](../README.md)

Every symox expression is one native Rust object type, `_core.Expr`. SymPy's
class hierarchy is reproduced on top of it in two ways: SymPy-named classes
classify native nodes for `isinstance`, `type` and dispatch, and user
subclasses of `Symbol`, `Function`, `Expr` and related classes become real
subclasses of the native type whose instances travel through native algorithms
as opaque host handles. This page describes what each mechanism preserves and
where it differs from SymPy 1.14. Behaviour gaps are tracked in the
[parity contract](algorithmic-parity.md).

## SymPy-named classes

`symox._node_types` generates one nominal class per SymPy core class listed in
`symox._sympy_classes`, with SymPy's method resolution order, and
`type(expr)` returns the nominal class of a native node:

```python
from symox import symbols, sin, Add, Float, Integer, Rational, Function

x, y = symbols('x y')
assert type(sin(x)) == sin and isinstance(x + y, Add)
assert issubclass(Integer, Rational) and issubclass(sin, Function)
assert [k.__name__ for k in type(Float(1.5)).__mro__][:6] == [
    'Float', 'Number', 'AtomicExpr', 'Atom', 'Expr', 'Basic']
e = x*y + 1
assert type(e)(*e.args) == e and e.func(*e.args) == e
```

The public names imported from `symox` (`Add`, `Mul`, `Symbol`, `Function`,
`Eq`, ...) are constructor classes that route to native construction. A
nominal class compares equal to its public class and hashes like it, so
`type(sin(x)) == sin` and dictionaries keyed by public classes work, but
identity does not hold: `type(x + y) is Add` and `type(sin(x)) is sin` are
false, and `Add.__base__` is `object` where SymPy gives `AssocOp`. Code that
walks `__base__` or `__mro__` of the public classes (Devito does) sees a
different chain from SymPy's.

For an undefined function `f = Function('f')`, `f(x).func == f` holds, but
`f(x).func is f` and `type(f(x)) == f` are false: the factory object and the
class of its applications are distinct. Function classes defined with
`class f(Function)` do not have this gap.

Attributes assigned to a public class after import, such as
`Symbol.custom = 42`, are forwarded to the nominal class and seen by native
instances and by user subclasses, except for names the native type already
defines.

## User subclasses

A class deriving from `Symbol`, `Dummy`, `Function`, `Expr`, `Atom`,
`AtomicExpr`, `Number`, `NumberSymbol` or `Basic` is registered when it is
defined. Its instances are native expressions that carry a reference to the
Python object, so arithmetic, `subs`, `xreplace`, `diff`, `expand`, `solve`,
matrices and printing return instances of the user class, with their instance
attributes:

```python
import copy, pickle
from symox import Symbol, Expr, DiracDelta, S, latex, srepr, solve

class Const(Symbol):
    def __new__(cls, name, meaning=''):
        obj = super().__new__(cls, name)
        obj.meaning = meaning
        return obj

    def __getnewargs_ex__(self):
        return ((self.name, self.meaning), {})

c = Const('c', 'speed')
assert (c*x + y).subs(y, 1).atoms(Const).pop().meaning == 'speed'
assert (c*x**2 + y).diff(c) == x**2 and solve(c*x - 1, x) == [1/c]
assert pickle.loads(pickle.dumps(c*x)).atoms(Const).pop().meaning == 'speed'
assert copy.deepcopy(c).meaning == 'speed' and srepr(c) == "Const('c')"

class step(Function):
    is_real = True

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            return S.One if arg >= 0 else S.Zero

    def fdiff(self, argindex=1):
        return DiracDelta(self.args[0])

    def _eval_is_nonnegative(self):
        return True

    def _sympystr(self, printer):
        return 'H[%s]' % printer._print(self.args[0])

    def _latex(self, printer):
        return r'\theta\left(%s\right)' % printer._print(self.args[0])

s = step(x**2)
assert (step(3), step(-1)) == (1, 0) and str(s) == 'H[x**2]'
assert latex(s) == r'\theta\left(x^{2}\right)'
assert s.diff(x) == 2*x*DiracDelta(x**2) and s.is_nonnegative
assert step.nargs == {1} and type(s.subs(x, y)) is step

class Pair(Expr):
    def __new__(cls, a, b):
        return Expr.__new__(cls, a, b)

    def _eval_derivative(self, v):
        return Pair(self.args[0].diff(v), self.args[1].diff(v))

p = Pair(x**2, x*y)
assert p.diff(x) == Pair(2*x, y) and str(p + p) == '2*Pair(x**2, x*y)'
assert srepr(p) == "Pair(Pow(Symbol('x'), Integer(2)), Mul(Symbol('x'), Symbol('y')))"
```

All of these results agree with SymPy 1.14. As in SymPy, pickling rebuilds a
`Symbol` subclass from its name only, so state passed to `__new__` survives
`pickle` and `deepcopy` only with `__getnewargs_ex__` (or `__getnewargs__`).
Calling `step(x, y)` raises `TypeError: step takes exactly 1 argument (2
given)`, as in SymPy.

## Honored hooks

The hooks below are consulted by native algorithms through `symox._host`:

- construction: `eval` (and an overridden `__new__`), `nargs` or the arity of
  `eval`'s signature, `func`, `_hashable_content` for equality and hashing;
- calculus and substitution: `fdiff` (default: `Derivative` and `Subs` as in
  SymPy), `_eval_derivative`, `_eval_subs`, `free_symbols`;
- assumptions: `_eval_is_<fact>` methods and `is_<fact> = True/False` class
  attributes, `is_commutative` and `_eval_is_commutative`;
- evaluation and rewriting: `_eval_evalf`, `_as_mpf_val` for `Number`
  subclasses (giving `float()` and `evalf`), `_eval_expand_<hint>` (including
  user-defined hints), `_eval_rewrite_as_<name>` and `_eval_rewrite`,
  `_eval_simplify`, `_eval_power` on expression classes, and `__add__` and
  `__mul__` of `Number` subclasses when `Add` and `Mul` fold coefficients;
- printing and ordering: `_sympystr`, `_latex`, `_pretty` (single-line
  results), `_sympyrepr`, `class_key`, `_op_priority` and a `name` property.

`_print_<ClassName>` methods of a Python printer subclass are used for user
classes: a `StrPrinter` subclass with `_print_step` prints `step(x) + 1`
through that method. `Singleton` metaclass classes register lazily in `S`.

## Known limits

These differences were measured against SymPy 1.14:

- `__slots__` on a subclass of a symox class does not remove the instance
  `__dict__`; SymPy's slotted hierarchy does.
- A subclass of a built-in function head such as `sin` inherits native
  evaluation at floats: `class mysin(sin)` gives `mysin(Float(0.5))` equal to
  `0.479425538604203`, while SymPy keeps `mysin(0.5)`.
- An exception raised inside a hook is reported only when the native algorithm
  returns. Polynomial simplifications treat noncommutative host factors as
  commuting, and native `subs` does not match noncommutative subsequences.
- The expand hints `frac`, `numer`, `denom` and `modulus` use the native
  implementation instead of user `_eval_expand_*` methods.
- The native tree hash ignores `_hashable_content`; equality uses it, so the
  only effect is extra hash collisions.
- `Basic.__new__` sympifies its arguments (SymPy stores them unchanged), and a
  `Wild` subclass overriding `matches` is not consulted.
- Nested Python `_print_<ClassName>` printer methods are not consulted by the
  native printers (see [printing](printing.md) for printer extension).
- There is no expression cache: equal applications are equal but not
  identical (`step(x) is step(x)` and `sin(x) is sin(x)` are false; SymPy's
  cache makes both true). SymPy's shared fact base between symbols with equal
  assumptions is not reproduced either.
- `codegen.ast` nodes and `IndexedBase`, `Indexed` and `Idx` are not `Basic`.

`tests/test_subclass_fidelity.py` reproduces subclassing patterns from
EinsteinPy, optlang, Lcapy, Devito and TorchInductor and runs each against
both symox and SymPy. It does not replace those projects' own suites: the
[featured-project audit](downstream-integration.md) records all five as native
incompatibilities. Its recorded EinsteinPy failure, an attribute assigned in a
`Symbol` subclass `__new__`, is the `Const` pattern above, which works with
the current build; the audit has to be rerun before any of these projects
can be claimed to work.
