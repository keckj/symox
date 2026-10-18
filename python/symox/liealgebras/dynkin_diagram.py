"""Text Dynkin diagrams."""
from .cartan_type import CartanType


def DynkinDiagram(cartan_type):
    return CartanType(cartan_type).dynkin_diagram()
