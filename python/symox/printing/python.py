"""Reconstruct native symbolic expressions as executable Python programs."""
from .repr import srepr

def python(expr,**settings):
    return 'from symox import *\ne = '+srepr(expr,**settings)

def print_python(expr,**settings):
    print(python(expr,**settings))
