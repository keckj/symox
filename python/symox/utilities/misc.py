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
def _dense_matrix_name(x):
    # The mutable and immutable dense matrices share the native class, which
    # SymPy names MutableDenseMatrix and ImmutableDenseMatrix.
    from symox import _core,ImmutableMatrix
    if x is _core.Matrix:return 'MutableDenseMatrix'
    if x is ImmutableMatrix:return 'ImmutableDenseMatrix'
    if isinstance(x, _core.Matrix):return 'ImmutableDenseMatrix' if x._immutable else 'MutableDenseMatrix'
    return None

def func_name(x,short=False):
    name=_dense_matrix_name(x) or getattr(x,'class_name',None) or getattr(getattr(x,'func',x),'__name__',type(x).__name__)
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
    # The reference rendering: dedent() for plain text, adjacent literals otherwise.
    lines=s.split('\n')
    if len(lines)==1:return repr(lines[0])
    triple=["'''" in s,'"""' in s]
    if any(li.endswith(' ') for li in lines) or '\\' in s or all(triple):
        rv=[];trailing=s.endswith('\n');last=len(lines)-1
        for i,li in enumerate(lines):
            rv.append(repr(li+'\n') if i!=last or trailing else repr(li))
        return '(\n    %s\n)'%'\n    '.join(rv)
    rv='\n    '.join(lines)
    if triple[0]:return 'dedent("""\\\n    %s""")'%rv
    return "dedent('''\\\n    %s''')"%rv

def strlines(s,c=64,short=False):
    if not isinstance(s,str):raise ValueError('expecting string input')
    if '\n' in s:return rawlines(s)
    q='"' if repr(s).startswith('"') else "'"
    q=(q,)*2
    if '\\' in s:
        m='(\nr%s%%s%s\n)'%q;j='%s\nr%s'%q;c-=3
    else:
        m='(\n%s%%s%s\n)'%q;j='%s\n%s'%q;c-=2
    out=[]
    while s:
        out.append(s[:c]);s=s[c:]
    if short and len(out)==1:return (m%out[0]).splitlines()[1]
    return m%j.join(out)
__all__ += ['Undecidable','debug','debugf','debug_decorator','find_executable','func_name','ordinal','replace','translate','rawlines','strlines']
