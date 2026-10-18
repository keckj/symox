"""Cartan matrix constructor."""
from .cartan_type import CartanType


def CartanMatrix(cartan_type):
    return CartanType(cartan_type).cartan_matrix()
