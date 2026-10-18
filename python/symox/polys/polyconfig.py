"""Context-local polynomial configuration with literal-only environment parsing."""
from ast import literal_eval
from contextlib import contextmanager
from contextvars import ContextVar
import os

_default_config={'USE_COLLINS_RESULTANT':False,'USE_SIMPLIFY_GCD':True,'USE_HEU_GCD':True,
 'USE_IRREDUCIBLE_IN_FACTOR':False,'USE_CYCLOTOMIC_FACTOR':True,'EEZ_RESTART_IF_NEEDED':True,
 'EEZ_NUMBER_OF_CONFIGS':3,'EEZ_NUMBER_OF_TRIES':5,'EEZ_MODULUS_STEP':2,
 'GF_IRRED_METHOD':'rabin','GF_FACTOR_METHOD':'zassenhaus','GROEBNER':'buchberger'}
_config=ContextVar('polynomial_configuration',default=dict(_default_config))

def configure():
    values=dict(_default_config)
    for key in values:
        value=os.getenv('SYMOX_'+key,os.getenv('SYMPY_'+key))
        if value is not None:
            try:value=literal_eval(value)
            except (ValueError,SyntaxError):pass
            values[key]=value
    _config.set(values)

def setup(key,value=None):
    key=key.upper();values=dict(_config.get());values[key]=_default_config[key] if value is None else value;_config.set(values)
def query(key):return _config.get().get(key.upper())
@contextmanager
def using(**kwargs):
    values=dict(_config.get());values.update({k.upper():v for k,v in kwargs.items()});token=_config.set(values)
    try:yield
    finally:_config.reset(token)
configure()
__all__=['configure','query','setup','using']
