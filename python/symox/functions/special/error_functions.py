"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Chi Ci E1 Ei Li Shi Si erf erf2 erf2inv erfc erfcinv erfi erfinv expint fresnelc fresnels li')

from ...core.function import DefinedFunction

class TrigonometricIntegral(DefinedFunction):
    _matches=staticmethod(lambda expr:expr.class_name in ('Si','Ci','Shi','Chi'))
class FresnelIntegral(DefinedFunction):
    _matches=staticmethod(lambda expr:expr.class_name in ('fresnelc','fresnels'))


def real_to_real_as_real_imag(self,deep=True,**hints):
    from symox import Integer,I
    argument=self.args[0]
    if argument.is_extended_real is True:return (self.expand(deep=deep,**hints) if deep else self),Integer(0)
    real,imag=(argument.expand(deep=deep,**hints) if deep else argument).as_real_imag()
    positive,negative=self.func(real+I*imag),self.func(real-I*imag)
    return (positive+negative)/2,(positive-negative)/(2*I)

__all__ += ['TrigonometricIntegral','FresnelIntegral','real_to_real_as_real_imag']
