"""Expression strings for generated Python numerical callables."""
from .pycode import PythonCodePrinter,MpmathPrinter,ArrayPrinter
class LambdaPrinter(PythonCodePrinter):pass
class NumPyPrinter(ArrayPrinter):pass
class IntervalPrinter(MpmathPrinter):_module='mpmath.iv'
def lambdarepr(expr,**settings):return LambdaPrinter(settings).doprint(expr)
__all__=['LambdaPrinter','NumPyPrinter','MpmathPrinter','PythonCodePrinter','IntervalPrinter','lambdarepr']
