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
