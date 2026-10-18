from .. import capture,filldedent
__all__ = ['capture','filldedent']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.ntheory._intfunc', 'as_int')

import os as _os
import re as _re
import sys as _sys
from functools import wraps as _wraps
from shutil import which as _which

class Undecidable(ValueError):pass

def debug(*args):
    if _os.getenv('SYMOX_DEBUG',_os.getenv('SYMPY_DEBUG','False')).lower() in ('true','1'):
        print(*args,file=_sys.stderr)
def debugf(string,args):debug(string%args)
def debug_decorator(func):
    @_wraps(func)
    def wrapped(*args,**kwargs):
        result=func(*args,**kwargs);debug(func.__name__,args,kwargs,'->',result);return result
    return wrapped

def find_executable(executable,path=None):return _which(executable,path=path)
def func_name(x,short=False):
    name=getattr(x,'class_name',None) or getattr(getattr(x,'func',x),'__name__',type(x).__name__)
    aliases={'GreaterThan':'Ge','StrictGreaterThan':'Gt','LessThan':'Le','StrictLessThan':'Lt','Equality':'Eq','Unequality':'Ne'}
    return aliases.get(name,name) if short else name

def ordinal(num):
    suffix='th' if 10<=abs(int(num))%100<=20 else {1:'st',2:'nd',3:'rd'}.get(abs(int(num))%10,'th')
    return str(num)+suffix

def replace(string,*reps):
    replacements=dict(reps[0]) if len(reps)==1 and isinstance(reps[0],dict) else dict(reps)
    if not replacements:return string
    pattern=_re.compile('|'.join(map(_re.escape,sorted(replacements,key=lambda k:-len(k)))))
    return pattern.sub(lambda match:replacements[match.group()],string)

def translate(s,a,b=None,c=None):
    if a is None:
        if c is not None:raise ValueError('third translation argument is invalid when the mapping is None')
        mapping={};delete=b or ''
    elif isinstance(a,dict):mapping=dict(a);delete=b or ''
    else:
        if b is None or len(a)!=len(b):raise ValueError('translation strings must have the same length')
        mapping=dict(zip(a,b));delete=c or ''
    if delete:s=s.translate(str.maketrans('','',delete))
    return replace(s,mapping)

def rawlines(s):
    if '\n' not in s:return repr(s)
    # Adjacent quoted literals preserve whitespace, escapes and line endings.
    return '(\n'+''.join('    '+repr(line)+'\n' for line in s.splitlines(keepends=True))+')'

def strlines(s,c=64,short=False):
    if not isinstance(s,str):raise ValueError('expected a string')
    if '\n' in s:return rawlines(s)
    if c<8:raise ValueError('line width must be at least eight characters')
    chunks=[];current=''
    for char in s:
        if len(repr(current+char))>c-4 and current:chunks.append(repr(current));current=''
        current+=char
    if current or not chunks:chunks.append(repr(current))
    return chunks[0] if short and len(chunks)==1 else '(\n'+'\n'.join(chunks)+'\n)'
__all__ += ['Undecidable','debug','debugf','debug_decorator','find_executable','func_name','ordinal','replace','translate','rawlines','strlines']
