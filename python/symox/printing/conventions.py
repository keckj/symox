"""Shared rules for symbol names and derivative notation."""
import re

def split_super_sub(text):
    text=str(text)
    if not text:return '',[],[]
    markers=list(re.finditer(r'__|_|\^',text[1:]))
    positions=[(match.start()+1,match.end()+1,match.group()) for match in markers]
    name=text[:positions[0][0]] if positions else text
    supers=[];subs=[]
    for j,(start,end,marker) in enumerate(positions):
        part=text[end:positions[j+1][0] if j+1<len(positions) else len(text)]
        (subs if marker=='_' else supers).append(part)
    match=re.fullmatch(r'([^\W\d_]+)(\d+)',name)
    if match:name=match.group(1);subs.insert(0,match.group(2))
    return name,supers,subs

def requires_partial(expr):
    if getattr(expr,'is_Derivative',False):return requires_partial(expr.expr)
    return sum(symbol.is_integer is not True for symbol in expr.free_symbols)>1
__all__=['split_super_sub','requires_partial']
