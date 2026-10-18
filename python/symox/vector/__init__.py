"""Three-dimensional coordinate vector calculus using native Symox scalars."""
from ._basis import (Vector,VectorAdd,VectorMul,BaseVector,VectorZero,Cross,Dot,
                     cross,dot,VectorKind,Dyadic,DyadicAdd,DyadicMul,BaseDyadic,DyadicZero)
from .coordsysrect import CoordSys3D,BaseScalar
from .point import Point
from .orienters import AxisOrienter,BodyOrienter,SpaceOrienter,QuaternionOrienter
from .operators import Gradient,Divergence,Curl,Laplacian,gradient,curl,divergence
from .functions import (express,matrix_to_vector,laplacian,is_conservative,
                        is_solenoidal,scalar_potential,directional_derivative,
                        scalar_potential_difference)
from .deloperator import Del
from .implicitregion import ImplicitRegion
from .parametricregion import ParametricRegion,parametric_region_list
from .integrals import ParametricIntegral,vector_integrate
from ._adapters import install as _install
_install()

__all__ = ['Vector','VectorAdd','VectorMul','BaseVector','VectorZero','Cross','Dot',
           'cross','dot','VectorKind','Dyadic','DyadicAdd','DyadicMul','BaseDyadic',
           'DyadicZero','BaseScalar','Del','CoordSys3D','express','matrix_to_vector',
           'laplacian','is_conservative','is_solenoidal','scalar_potential',
           'directional_derivative','scalar_potential_difference','Point','AxisOrienter',
           'BodyOrienter','SpaceOrienter','QuaternionOrienter','Gradient','Divergence',
           'Curl','Laplacian','gradient','curl','divergence','ParametricRegion',
           'parametric_region_list','ImplicitRegion','ParametricIntegral','vector_integrate']
