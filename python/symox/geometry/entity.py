"""Exact geometry containers backed by native scalar expressions."""


class GeometryError(ValueError):
    pass


def scalar(value):
    from ..calculus._analysis import scalar as convert
    return convert(value)


def clean(value):
    from .. import cancel, simplify
    return simplify(cancel(value))


def zero(value):
    value = clean(value)
    return value == 0 or value.is_zero is True


def sign(value):
    from .. import oo
    if value == oo:
        return 1
    if value == -oo:
        return -1
    value = clean(value)
    if value == 0 or value.is_zero is True:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    raise NotImplementedError('Geometry requires a decidable real sign')


def compare(left, right):
    from .. import oo
    if left == right:
        return 0
    if left == oo or right == -oo:
        return 1
    if left == -oo or right == oo:
        return -1
    return sign(left - right)


class GeometryEntity:
    __slots__ = ('_args',)

    def __init__(self, *args, **kwargs):
        pass

    def __setattr__(self, name, value):
        raise AttributeError('Geometry objects are immutable')

    @property
    def args(self):
        return self._args

    @property
    def func(self):
        return type(self)

    @property
    def free_symbols(self):
        return set().union(*(arg.free_symbols for arg in self.args))

    @property
    def bounds(self):
        from .point import Point
        from .. import Min, Max
        if isinstance(self, Point):
            points = [self]
        elif hasattr(self, 'vertices'):
            points = self.vertices
        elif hasattr(self, 'points'):
            points = self.points
        else:
            raise NotImplementedError('A finite planar bounding box is not available for this entity')
        if any(len(point) != 2 for point in points):
            raise NotImplementedError('Geometry bounding boxes currently require planar entities')
        return (Min(*(point.x for point in points)), Min(*(point.y for point in points)),
                Max(*(point.x for point in points)), Max(*(point.y for point in points)))

    def __hash__(self):
        return hash((type(self), self.args))

    def __eq__(self, other):
        return type(self) is type(other) and self.args == other.args

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, ', '.join(map(str, self.args)))

    __str__ = __repr__

    def __reduce__(self):
        return type(self), self.args

    def subs(self, *args, **kwargs):
        return self.func(*(arg.subs(*args, **kwargs) for arg in self.args))

    def xreplace(self, rule):
        if self in rule:
            return rule[self]
        return self.func(*(arg.xreplace(rule) for arg in self.args))

    def evalf(self, n=15, **options):
        return self.func(*(arg.evalf(n, **options) for arg in self.args), evaluate=False)

    n = evalf

    def simplify(self, **kwargs):
        from .. import simplify
        return self.func(*(arg.simplify(**kwargs) if isinstance(arg, GeometryEntity)
                           else simplify(arg, **kwargs) for arg in self.args))

    def to_sympy(self):
        import sympy
        return getattr(sympy, type(self).__name__)(*(arg.to_sympy() for arg in self.args))

    def _sympy_(self):
        return self.to_sympy()

    def __contains__(self, other):
        result = self.contains(other)
        if result is None:
            raise NotImplementedError('Geometric containment is undecidable')
        return bool(result)

    def translate(self, x=0, y=0, z=0):
        return self.func(*(arg.translate(x, y, z) if isinstance(arg, GeometryEntity)
                           else arg for arg in self.args))

    def rotate(self, angle, pt=None):
        return self.func(*(arg.rotate(angle, pt) if isinstance(arg, GeometryEntity)
                           else arg for arg in self.args))

    def scale(self, x=1, y=1, pt=None):
        return self.func(*(arg.scale(x, y, pt=pt) if isinstance(arg, GeometryEntity)
                           else arg for arg in self.args))


def construct(cls, args):
    value = object.__new__(cls)
    object.__setattr__(value, '_args', tuple(args))
    return value


def unique(values):
    result = []
    for value in values:
        if value not in result:
            result.append(value)
    return result
