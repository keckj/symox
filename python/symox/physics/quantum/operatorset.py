"""Mappings between supported eigenstates and their commuting observables."""
import symox as S
from ._core import Operator,StateBase,Ket,Bra
from .cartesian import XOp,YOp,ZOp,PxOp,XKet,PxKet,XBra,PxBra,PositionKet3D,PositionBra3D
from .spin import JxOp,JyOp,JzOp,J2Op,JxKet,JyKet,JzKet,JxBra,JyBra,JzBra
from .sho1d import NumberOp,SHOKet,SHOBra


def operators_to_state(operators,**options):
    ops=operators if isinstance(operators,(set,frozenset)) else {operators}
    classes={o if isinstance(o,type) else type(o) for o in ops}
    mapping={XOp:(XKet,('x',)),PxOp:(PxKet,('px',)),NumberOp:(SHOKet,(S.Symbol('n',integer=True,nonnegative=True),)),Operator:(Ket,('psi',))}
    if len(classes)==1:
        cls=next(iter(classes))
        if cls in mapping:
            state,args=mapping[cls];return state(*args)
    if classes=={XOp,YOp,ZOp}:return PositionKet3D('x','y','z')
    if J2Op in classes and len(classes)==2:
        axis=next(c for c in classes if c is not J2Op);state={JxOp:JxKet,JyOp:JyKet,JzOp:JzKet}.get(axis)
        if state:return state(S.Symbol('j'),S.Symbol('m'))
    return None

def state_to_operators(state,**options):
    cls=state if isinstance(state,type) else type(state)
    mapping={XKet:XOp,XBra:XOp,PxKet:PxOp,PxBra:PxOp,SHOKet:NumberOp,SHOBra:NumberOp,Ket:Operator,Bra:Operator}
    if cls in mapping:return mapping[cls]('O') if mapping[cls] is Operator else mapping[cls]()
    if cls in (PositionKet3D,PositionBra3D):return {XOp('X'),YOp('Y'),ZOp('Z')}
    spin={JxKet:JxOp,JxBra:JxOp,JyKet:JyOp,JyBra:JyOp,JzKet:JzOp,JzBra:JzOp}
    if cls in spin:return {spin[cls]('J'),J2Op('J2')}
    return None
__all__=['operators_to_state','state_to_operators']
