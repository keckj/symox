"""The nabla operator and its dot/cross/gradient spelling."""
from .operators import gradient,divergence,curl


class Del:
    def __repr__(self):
        return 'delop'

    def __call__(self,field,doit=False):
        return gradient(field,doit=doit)

    gradient = __call__

    def dot(self,field,doit=False):
        return divergence(field,doit=doit)

    __and__ = dot

    def cross(self,field,doit=False):
        return curl(field,doit=doit)

    __xor__ = cross
