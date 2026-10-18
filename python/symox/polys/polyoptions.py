"""Validated public polynomial option construction."""
from .polytools import Options as _Options,_build_options,_DEFAULTS,_FLAGS
from .polyerrors import OptionError,FlagError


class Options(_Options):
    def __init__(self,gens,args,flags=None,strict=False):
        args = dict(args)
        unknown = set(args)-set(_DEFAULTS)-{'order'}
        if unknown:raise OptionError('unknown polynomial option: '+', '.join(sorted(unknown)))
        if flags is not None:
            forbidden = set(args)&set(_FLAGS)-set(flags)
            if forbidden:raise FlagError('flags are not allowed: '+', '.join(sorted(forbidden)))
        if args.get('field') and args.get('modulus') is not None:raise OptionError('field and modulus cannot be combined')
        super().__init__(_build_options(tuple(gens),args))
        if len(set(self.gens)) != len(self.gens):raise OptionError('generators must be distinct')
        modulus = self.get('modulus')
        if modulus is not None:
            from .. import GF
            if int(modulus) < 2:raise OptionError('modulus must be at least two')
            self['domain'] = GF(int(modulus))
        elif isinstance(self.get('domain'),str):
            from .. import ZZ,QQ,RR,CC
            domains = {'ZZ':ZZ,'QQ':QQ,'RR':RR,'CC':CC}
            if self['domain'] not in domains:raise OptionError('unrecognized domain specification')
            self['domain'] = domains[self['domain']]
    def clone(self,updates=None):
        data = dict(self);data.update(updates or {})
        return Options((),data)


__all__ = ['Options']
