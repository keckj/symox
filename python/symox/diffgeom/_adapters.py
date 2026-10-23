"""Native expression dispatch for differential tensors and point evaluation."""
from functools import wraps
from .._dispatch_arguments import FirstArgument


def _field_expression(expr):
    # Unpickles a callable scalar-field expression.
    import symox
    install()
    return symox._core.Expr.__new__(symox._diffgeom_callable_type,expr)


def install():
    import symox
    from ._tensors import _Tensor
    from ._coordinates import Point,_scalar_at,_atoms,BaseScalarField,CoordinateSymbol
    Expr = symox._core.Expr
    if getattr(Expr,'_diffgeom_dispatch_installed',False):
        return
    from ._coordinates import _field_atoms
    class FieldExpression(Expr):
        __slots__ = ()
        def __reduce_ex__(self,protocol):
            return (_field_expression,(Expr.__new__(Expr,self),))
        def __call__(self,*args,**kwargs):
            if len(args)==1 and isinstance(args[0],Point) and not kwargs:
                return _scalar_at(self,args[0])
            raise TypeError('scalar fields are evaluated at differential-geometric points')
    symox._diffgeom_callable_type=FieldExpression
    symox._diffgeom_field_atoms=_field_atoms
    for name in ('__add__','__sub__','__mul__','__truediv__'):
        original = getattr(Expr,name)
        def arithmetic(original,name):
            @wraps(original)
            def operation(self,other):
                if isinstance(other,_Tensor):
                    if name == '__add__':
                        return other.__radd__(self)
                    if name == '__sub__':
                        return other.__rsub__(self)
                    if name == '__mul__':
                        return other.__rmul__(self)
                    raise ValueError('division by a differential tensor is undefined')
                return original(self,other)
            return operation
        setattr(Expr,name,arithmetic(original,name))
    for name in ('rcall',):
        original = getattr(Expr,name)
        def evaluation(original):
            @wraps(original)
            def operation(self,*args,**kwargs):
                if len(args) == 1 and isinstance(args[0],Point) and not kwargs:
                    return _scalar_at(self,args[0])
                if self == 0 and any(isinstance(a,_Tensor) for a in args) and not kwargs:
                    return self
                return original(self,*args,**kwargs)
            return operation
        setattr(Expr,name,evaluation(original))
    original_atoms = Expr.atoms
    @wraps(original_atoms)
    def atoms(self,*types):
        special = tuple(t for t in types if t in (BaseScalarField,CoordinateSymbol))
        if special:
            result = {field for atom in self.free_symbols if isinstance((field := _atoms.get(atom)),special)}
            ordinary = tuple(t for t in types if t not in special)
            return result | (original_atoms(self,*ordinary) if ordinary else set())
        return original_atoms(self,*types)
    Expr.atoms = atoms
    for name in ('diff','simplify','expand'):
        original = getattr(symox,name)
        def transform(original,name):
            first=FirstArgument(original)
            @wraps(original)
            def operation(*args,**kwargs):
                expression=first.get(args,kwargs)
                if isinstance(expression,_Tensor):
                    expression,rest,options=first.split(args,kwargs)
                    return getattr(expression,name)(*rest,**options)
                return original(*args,**kwargs)
            return operation
        setattr(symox,name,transform(original,name))
    Expr._diffgeom_dispatch_installed = True
