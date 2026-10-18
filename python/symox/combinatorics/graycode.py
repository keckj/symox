"""Binary reflected Gray codes with exact rank and unrank."""
import operator


def _rank(bits):
    value=int(bits,2);rank=0
    while value:
        rank^=value;value>>=1
    return rank


class GrayCode:
    def __init__(self,n,*args,**kwargs):
        self.n=operator.index(n)
        if self.n<1:
            raise ValueError('Gray code width must be positive')
        if args or kwargs.keys()-{'start','rank'} or 'start' in kwargs and 'rank' in kwargs:
            raise ValueError('supply either start or rank')
        self.selections=1<<self.n
        if 'start' in kwargs:
            start=str(kwargs['start'])
            if len(start)>self.n or not start or set(start)-{'0','1'}:
                raise ValueError('invalid Gray bit string')
            self._rank=_rank(start)
        else:
            self._rank=operator.index(kwargs.get('rank',0))%self.selections
        self._skip=False

    @property
    def rank(self):
        return self._rank

    @property
    def current(self):
        return self.unrank(self.n,self.rank)

    @classmethod
    def unrank(cls,n,rank):
        n=operator.index(n);rank=operator.index(rank)
        if n<1 or not 0<=rank<1<<n:
            raise ValueError('Gray rank out of range')
        return format(rank^(rank>>1),'0%db'%n)

    def next(self,delta=1):
        return type(self)(self.n,rank=(self.rank+operator.index(delta))%self.selections)

    def generate_gray(self,**hints):
        if hints.keys()-{'start','rank'}:
            raise TypeError('unsupported Gray code hint')
        start=type(self)(self.n,**hints).rank if hints else self.rank
        try:
            for i in range(start,self.selections):
                self._rank=i
                if self._skip:
                    self._skip=False
                else:
                    yield self.current
        finally:
            self._rank=0

    def skip(self):
        self._skip=True


def _bits(value):
    value=''.join(map(str,value))
    if not value or set(value)-{'0','1'}:raise ValueError('expected a nonempty binary string')
    return value

def bin_to_gray(bin_list):
    bits=_bits(bin_list);return bits[0]+''.join(str(int(a)^int(b)) for a,b in zip(bits,bits[1:]))
def gray_to_bin(gray_list):
    bits=_bits(gray_list);result=[];bit=0
    for value in bits:bit^=int(value);result.append(str(bit))
    return ''.join(result)
def random_bitstring(n):
    from ..core.random import rng
    n=operator.index(n)
    if n<0:raise ValueError('bit string length must be nonnegative')
    return ''.join(str(rng.randrange(2)) for _ in range(n))
def get_subset_from_bitstring(super_set,bitstring):
    values=list(super_set);bits=_bits(bitstring)
    if len(values)!=len(bits):raise ValueError('bit string and superset lengths differ')
    return [value for value,bit in zip(values,bits) if bit=='1']
def graycode_subsets(gray_code_set):
    values=list(gray_code_set)
    if not values:yield [];return
    for bits in GrayCode(len(values)).generate_gray():yield get_subset_from_bitstring(values,bits)
