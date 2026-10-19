"""Unicode and ASCII mathematical glyphs with native symbol formatting."""
import unicodedata
from ._symbology_data import *
from ._symbology_data import _xobj_unicode,_xobj_ascii,_xsym
from ..conventions import split_super_sub
from .pretty import pprint_use_unicode as pretty_use_unicode,pprint_try_use_unicode as pretty_try_use_unicode


def U(name):
    try:return unicodedata.lookup(name)
    except KeyError:return None

def g(l):return U('GREEK SMALL LETTER '+l.upper())
def G(l):return U('GREEK CAPITAL LETTER '+l.upper())
def VF(txt):return U('VULGAR FRACTION '+txt)
def xstr(*args):return str(*args)
def line_width(text):return sum(0 if unicodedata.combining(c) else 2 if unicodedata.east_asian_width(c) in ('W','F') else 1 for c in str(text))
def center_pad(wstring,wtarget,fillchar=' '):
    if wstring>wtarget:raise ValueError('target width is smaller than the string')
    delta=wtarget-wstring;left=delta//2
    return fillchar*left,fillchar*(delta-left)
def center(string,width,fillchar=' '):
    left,right=center_pad(line_width(string),width,fillchar)
    return left+string+right

def xobj(symb,length):
    if length<=0:raise ValueError('a spatial glyph must have positive length')
    table=_xobj_unicode if pretty_use_unicode() else _xobj_ascii
    shape=table[symb];single=None
    if isinstance(shape,tuple) and isinstance(shape[0],tuple):shape,single=shape
    if isinstance(shape,tuple):
        extension=shape[0];top=shape[1] if len(shape)>1 else extension;bottom=shape[2] if len(shape)>2 else extension;middle=shape[3] if len(shape)>3 else None
    else:extension=top=bottom=shape;middle=None
    if length==1:return single or extension
    if middle is not None and length%2==0:length+=1
    result=[extension]*length;result[0]=top;result[-1]=bottom
    if middle is not None:result[length//2]=middle
    return result

def vobj(symb,height):return '\n'.join(xobj(symb,height))
def hobj(symb,width):return ''.join(xobj(symb,width))
def xsym(sym):return _xsym[sym][1 if pretty_use_unicode() else 0]
def pretty_atom(atom_name,default=None,printer=None):
    if pretty_use_unicode():
        if atom_name=='ImaginaryUnit' and printer is not None and printer._settings.get('imaginary_unit')=='j':return U('DOUBLE-STRUCK ITALIC SMALL J')
        return atoms_table[atom_name]
    if default is not None:return default
    raise KeyError('atom has only a Unicode representation')
def is_subscriptable_in_unicode(subscript):return all(char in sub for char in str(subscript))
def pretty_symbol(symb_name,bold_name=False):
    if not pretty_use_unicode():return str(symb_name)
    from symox import Symbol,_core
    text=_core.pretty(Symbol(str(symb_name)),use_unicode=True,wrap_line=False)
    if not bold_name:return text
    def bold(char):
        name=unicodedata.name(char,'').replace('LATIN ','').replace('GREEK ','').replace('LETTER ','')
        return U('MATHEMATICAL BOLD '+name) or char
    return ''.join(bold(c) for c in text)
def annotated(letter):
    if letter not in ('F','G'):raise KeyError(letter)
    if pretty_use_unicode():
        return (2,0,2,0,'┌─\n├─\n╵') if letter=='F' else (3,0,3,1,'╭─╮\n│ ╮\n╰─╯')
    return (3,0,3,0,' _\n|_\n|') if letter=='F' else (3,0,3,1,' __\n/ _\n\\_|')

__all__=['U','G','g','VF','greek_letters','greek_unicode','digit_2txt','symb_2txt','sub','sup','root','frac',
 'atoms_table','pretty_use_unicode','pretty_try_use_unicode','xstr','xobj','vobj','hobj','xsym','pretty_atom',
 'pretty_symbol','annotated','center_pad','center','is_subscriptable_in_unicode','line_width']
