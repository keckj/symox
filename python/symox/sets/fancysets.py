"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'ComplexRegion Complexes ImageSet Integers Naturals Naturals0 Range Rationals Reals')
from symox import ComplexRegion


def normalize_theta_set(theta):
    from symox import Interval,FiniteSet,Union,EmptySet,pi,Mod,oo
    if theta==EmptySet:return theta
    if isinstance(theta,FiniteSet):return FiniteSet(*(Mod(angle,2*pi) for angle in theta))
    if isinstance(theta,Union):return Union(*(normalize_theta_set(piece) for piece in theta.args))
    if not isinstance(theta,Interval):raise NotImplementedError('angular normalization requires intervals or finite sets')
    width=theta.end-theta.start
    if width.is_infinite is True or width>2*pi:return Interval.Ropen(0,2*pi)
    if width==2*pi:
        if theta.left_open and theta.right_open:
            a=Mod(theta.start,2*pi)
            from symox import Complement
            return Complement(Interval.Ropen(0,2*pi),FiniteSet(a))
        return Interval.Ropen(0,2*pi)
    a=Mod(theta.start,2*pi);b=a+width
    if b<2*pi:return Interval(a,b,theta.left_open,theta.right_open)
    if b==2*pi:return Union(Interval.Ropen(a,2*pi) if not bool(theta.left_open) else Interval.open(a,2*pi),FiniteSet(0) if not theta.right_open else EmptySet)
    return Union(Interval(a,2*pi,theta.left_open,True),Interval(0,b-2*pi,False,theta.right_open))


class CartesianComplexRegion(ComplexRegion):
    def __init__(self,sets):super().__init__(sets,polar=False)
    args=property(lambda self:(self.sets,))


class PolarComplexRegion(ComplexRegion):
    def __init__(self,sets):
        from symox import Union,ProductSet
        pieces=sets.args if isinstance(sets,Union) else (sets,)
        normalized=Union(*(ProductSet(piece.args[0],normalize_theta_set(piece.args[1])) for piece in pieces))
        super().__init__(normalized,polar=True)
    args=property(lambda self:(self.sets,))

__all__ += ['CartesianComplexRegion','PolarComplexRegion','normalize_theta_set']
