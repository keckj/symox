"""Coefficient typing for the native domain matrix interfaces."""
from typing import Union
from fractions import Fraction
from ... import DomainElement
RingElement=Union[DomainElement,int,Fraction]
__all__=['RingElement']
