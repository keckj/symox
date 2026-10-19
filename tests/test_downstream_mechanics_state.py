"""Kane's state vectors remain available to downstream equation consumers."""
import pytest
import symox as sx
import sympy as sp
import symox.physics.mechanics as sm
import sympy.physics.mechanics as pm


@pytest.mark.parametrize('count',[1,2])
def test_kane_state_vectors_match_reference(count):
    def system(library,mechanics):
        q=mechanics.dynamicsymbols('q:'+str(count))
        u=mechanics.dynamicsymbols('u:'+str(count))
        time=mechanics.dynamicsymbols._t
        kinematic=[qi.diff(time)-2*ui for qi,ui in zip(q,u)]
        return mechanics.KanesMethod(mechanics.ReferenceFrame('N'),q,u,kinematic)
    actual=system(sx,sm);expected=system(sp,pm)
    for name in ['_q','_u','_qdot','_udot','_qdep','_udep','_uaux']:
        value,reference=getattr(actual,name),getattr(expected,name)
        assert value.shape==reference.shape
        assert value.to_sympy()==reference
    assert actual._q is actual.q and actual._u is actual.u
    assert actual.kindiffdict()=={qd:2*u for qd,u in zip(actual._qdot,actual.u)}


def test_numeric_rotations_accept_rounding_residuals_like_reference():
    actual=sm.ReferenceFrame('N');rotated=sm.ReferenceFrame('A')
    reference=pm.ReferenceFrame('N');expected=pm.ReferenceFrame('A')
    angles=[float(sp.pi/2),float(sp.pi),0.0]
    rotated.orient_body_fixed(actual,angles,'XYZ')
    expected.orient_body_fixed(reference,angles,'XYZ')
    import numpy as np
    np.testing.assert_allclose(np.array(rotated.dcm(actual).tolist(),dtype=float),
                               np.array(expected.dcm(reference).tolist(),dtype=float),atol=1e-15)


@pytest.mark.parametrize('method',['LU','ADJ','BLOCK','CRAMER'])
def test_matrix_solver_dispatch_and_invalid_method_match_reference(method):
    matrix=sx.Matrix([[3,1],[1,2]]);rhs=sx.Matrix([2,4])
    actual=matrix.solve(rhs,method=method)
    assert actual.to_sympy()==matrix.to_sympy().solve(rhs.to_sympy(),method=method)
    assert matrix*actual==rhs
    with pytest.raises(ValueError):matrix.solve(rhs,method='UNKNOWN')
    with pytest.raises(ValueError):matrix.to_sympy().solve(rhs.to_sympy(),method='UNKNOWN')
