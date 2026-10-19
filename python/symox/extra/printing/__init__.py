"""Optional Symox printing extensions."""
__all__ = ['AdaptiveCppCodePrinter', 'JAXPrinter', 'JaxPrinter', 'PyTorchPrinter', 'SYCLCodePrinter', 'SlangCodePrinter', 'SlangPrinter', 'SyclCodePrinter', 'TorchPrinter', 'adaptivecppcode', 'jaxcode', 'slangcode', 'syclcode', 'torchcode']

def __getattr__(name):
    if name not in __all__:
        raise AttributeError(name)
    from .. import __getattr__ as resolve
    return resolve(name)
