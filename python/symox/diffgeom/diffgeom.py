"""Local differential geometry with native scalar and matrix operations."""
from ._coordinates import Manifold,Patch,CoordSystem,CoordinateSymbol,Point,BaseScalarField
from ._tensors import (BaseVectorField,Commutator,Differential,TensorProduct,WedgeProduct,
                       LieDerivative,contravariant_order,covariant_order,vectors_in_basis)
from ._calculus import (BaseCovarDerivativeOp,CovarDerivativeOp,intcurve_diffequ,intcurve_series,
                        twoform_to_matrix,metric_to_Christoffel_1st,metric_to_Christoffel_2nd,
                        metric_to_Riemann_components,metric_to_Ricci_components)

__all__ = ['BaseCovarDerivativeOp','BaseScalarField','BaseVectorField','Commutator',
    'contravariant_order','CoordSystem','CoordinateSymbol','CovarDerivativeOp',
    'covariant_order','Differential','intcurve_diffequ','intcurve_series','LieDerivative',
    'Manifold','metric_to_Christoffel_1st','metric_to_Christoffel_2nd',
    'metric_to_Ricci_components','metric_to_Riemann_components','Patch','Point',
    'TensorProduct','twoform_to_matrix','vectors_in_basis','WedgeProduct']
