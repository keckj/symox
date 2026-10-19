"""Preserve a wrapped operation's first-argument keyword during dispatch."""
from inspect import Parameter,signature


class FirstArgument:
    def __init__(self,function):
        try:self.signature=signature(function)
        except (TypeError,ValueError):self.signature=None
        parameter=next(iter(self.signature.parameters.values()),None) if self.signature else None
        self.name=(parameter.name if parameter is not None and parameter.kind in
                   (Parameter.POSITIONAL_ONLY,Parameter.POSITIONAL_OR_KEYWORD) else None)

    def get(self,args,kwargs):
        return args[0] if args else kwargs.get(self.name)

    def split(self,args,kwargs):
        # Validate before invoking callbacks or mutating specialized objects.
        if self.signature is not None:self.signature.bind(*args,**kwargs)
        if args:return args[0],args[1:],kwargs
        options=dict(kwargs)
        return options.pop(self.name),(),options
