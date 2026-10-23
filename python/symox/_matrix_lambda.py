"""``Lambda`` with an explicit matrix body.

Native Lambda bodies are scalar expressions. The reference also accepts an
explicit matrix body, stored as an ImmutableDenseMatrix: calling it
substitutes the arguments simultaneously and returns the matrix.
"""


class MatrixLambda:
    __slots__ = ('signature', 'expr')
    is_Lambda = True
    is_identity = False
    is_Function = False

    def __init__(self, signature, expr):
        from . import ImmutableMatrix, Tuple, Symbol
        from ._core import BadSignatureError
        if isinstance(signature, (Tuple, tuple, list)):
            signature = tuple(signature)
        else:
            signature = (signature,)
        if not all(isinstance(s, Symbol) for s in signature):
            # Nested signatures are not supported for matrix bodies.
            for s in signature:
                if not isinstance(s, Symbol):
                    raise BadSignatureError('Lambda signature should be only tuples and symbols, not %s' % (s,))
        if len(set(signature)) != len(signature):
            raise BadSignatureError('Duplicate symbol in Lambda signature')
        object.__setattr__(self, 'signature', signature)
        object.__setattr__(self, 'expr', ImmutableMatrix(expr))

    def __setattr__(self, name, value):
        raise AttributeError('Lambda is immutable')

    @property
    def variables(self):
        return self.signature

    @property
    def args(self):
        from . import Tuple
        return (Tuple(*self.signature), self.expr)

    @property
    def nargs(self):
        from . import FiniteSet
        return FiniteSet(len(self.signature))

    @property
    def bound_symbols(self):
        return set(self.signature)

    @property
    def free_symbols(self):
        return set(self.expr.free_symbols) - set(self.signature)

    def __call__(self, *args):
        from ._core import BadArgumentsError
        if len(args) != len(self.signature):
            raise BadArgumentsError('%s takes exactly %d argument%s (%d given)' % (
                self, len(self.signature), '' if len(self.signature) == 1 else 's', len(args)))
        from . import sympify
        return self.expr.subs(dict(zip(self.signature, map(sympify, args))), simultaneous=True)

    def subs(self, *args, **kwargs):
        from . import sympify
        mapping = dict(args[0]) if len(args) == 1 else {args[0]: args[1]}
        mapping = {sympify(k): sympify(v) for k, v in mapping.items() if sympify(k) not in self.signature}
        if not mapping:
            return self
        return MatrixLambda(self.signature, self.expr.subs(mapping, **kwargs))

    def _canonical(self):
        from . import Dummy
        dummies = [Dummy('x%d' % i) for i in range(len(self.signature))]
        return self.expr.subs(dict(zip(self.signature, dummies)), simultaneous=True), dummies

    def __eq__(self, other):
        # Equal up to renaming of the bound variables, as the reference.
        if not isinstance(other, MatrixLambda) or len(other.signature) != len(self.signature):
            return False
        renamed = other.expr.subs(dict(zip(other.signature, self.signature)), simultaneous=True)
        return renamed == self.expr

    def __hash__(self):
        return hash((MatrixLambda, len(self.signature), self.expr.shape))

    def __str__(self):
        from .printing.str import sstr
        from . import Tuple
        sig = self.signature[0] if len(self.signature) == 1 else Tuple(*self.signature)
        return 'Lambda(%s, %s)' % (sstr(sig), sstr(self.expr))

    __repr__ = __str__

    def _latex(self, printer):
        from . import Tuple
        sig = self.signature[0] if len(self.signature) == 1 else Tuple(*self.signature)
        return r'\left( %s \mapsto %s \right)' % (printer._print(sig), printer._print(self.expr))

    def _sympy_(self):
        import sympy
        return sympy.Lambda(tuple(sympy.sympify(s) for s in self.signature), sympy.sympify(self.expr))
