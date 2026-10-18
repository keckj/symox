"""Baseline-aware mathematical text pictures and composable arithmetic layouts."""
import shutil
from .pretty_symbology import center,line_width,vobj,hobj,xsym,pretty_use_unicode


def _picture(value):return value if isinstance(value,stringPict) else stringPict(str(value))


class stringPict:
    LINE='line'
    def __init__(self,s,baseline=0):
        self.s=str(s);self.picture=self.equalLengths(self.s.splitlines());self.baseline=int(baseline);self.binding=None
    @staticmethod
    def equalLengths(lines):
        lines=list(lines) or [''];width=max(map(line_width,lines))
        return [center(line,width) for line in lines]
    def height(self):return len(self.picture)
    def width(self):return line_width(self.picture[0])
    @staticmethod
    def next(*args):
        values=[_picture(value) for value in args]
        if not values:return '',0
        baseline=max(v.baseline for v in values);below=max(v.height()-v.baseline for v in values);rows=[]
        for value in values:
            rows.append([' '*value.width()]*(baseline-value.baseline)+value.picture+[' '*value.width()]*(below-(value.height()-value.baseline)))
        return '\n'.join(''.join(parts) for parts in zip(*rows)),baseline
    def right(self,*args):return self.next(self,*args)
    def left(self,*args):return self.next(*args,self)
    @staticmethod
    def stack(*args):
        values=[value if value==stringPict.LINE else _picture(value) for value in args]
        width=max((v.width() for v in values if v!=stringPict.LINE),default=1)
        values=[stringPict(hobj('-',max(width,1))) if v==stringPict.LINE else v for v in values]
        rows=[center(row,width) for value in values for row in value.picture]
        baseline=values[0].height()+values[1].baseline if len(values)>1 else values[0].baseline if values else 0
        return '\n'.join(rows),baseline
    def below(self,*args):return self.stack(self,*args)[0],self.baseline
    def above(self,*args):
        text,base=self.stack(*args,self)
        return text,len(text.splitlines())-self.height()+self.baseline
    def parens(self,left='(',right=')',ifascii_nougly=False):
        if ifascii_nougly and not pretty_use_unicode():return self.s,self.baseline
        a=stringPict(vobj(left,self.height()),self.baseline) if left else ''
        b=stringPict(vobj(right,self.height()),self.baseline) if right else ''
        return self.next(a,self,b)
    def leftslash(self):
        height=self.height();slash='\n'.join(' '*(height-i-1)+'/'+' '*i for i in range(height))
        return self.next(stringPict(slash,self.baseline),self)
    def root(self,n=None):
        radical='√' if pretty_use_unicode() else '\\/'
        top=hobj('_',max(1,self.width()))
        body=stringPict(*self.above(top))
        text,baseline=self.next(radical,body)
        if n is not None and str(n)!='2':text,baseline=self.next(stringPict(str(n),baseline),stringPict(text,baseline))
        return text,baseline
    def terminal_width(self):return shutil.get_terminal_size((80,24)).columns
    def render(self,*args,**kwargs):
        if not kwargs.get('wrap_line',True):return str(self)
        width=kwargs.get('num_columns') or self.terminal_width()
        if self.width()<=width:return str(self)
        if width<4:raise ValueError('line wrapping needs at least four columns')
        chunks=[]
        for offset in range(0,self.width(),width):chunks.append('\n'.join(row[offset:offset+width] for row in self.picture))
        return '\n\n'.join(chunks)
    def __str__(self):return '\n'.join(self.picture)
    def __repr__(self):return 'stringPict(%r, %d)'%(str(self),self.baseline)
    def __eq__(self,other):return self.picture==other.picture if isinstance(other,stringPict) else str(self)==other if isinstance(other,str) else False
    def __getitem__(self,index):return self.picture[index]
    def __len__(self):return len(self.s)


class prettyForm(stringPict):
    ATOM,FUNC,DIV,POW,MUL,ADD,NEG,OPEN=range(8)
    def __init__(self,s,baseline=0,binding=0,unicode=None):super().__init__(s,baseline);self.binding=binding;self._unicode=unicode or s
    unicode=property(lambda self:self._unicode)
    def _bound(self,limit):return stringPict(*self.parens()) if self.binding>limit else self
    def __add__(self,*others):
        values=[self._bound(self.NEG)]
        for other in others:
            other=other if isinstance(other,prettyForm) else prettyForm(str(other))
            values.extend([' ' if other.binding==self.NEG else ' + ',other._bound(self.NEG)])
        return prettyForm(*self.next(*values),binding=self.ADD)
    def __truediv__(self,den,slashed=False):
        den=den if isinstance(den,prettyForm) else prettyForm(str(den))
        if slashed:return prettyForm(*self.next(self,'/',den),binding=self.DIV)
        return prettyForm(*self.stack(self._bound(self.DIV),self.LINE,den._bound(self.DIV)),binding=self.DIV)
    def __mul__(self,*others):
        values=[self._bound(self.MUL)]
        for other in others:
            other=other if isinstance(other,prettyForm) else prettyForm(str(other))
            values.extend([xsym('*'),other._bound(self.MUL)])
        return prettyForm(*self.next(*values),binding=self.MUL)
    def __pow__(self,exponent):
        exponent=_picture(exponent);base=self._bound(self.FUNC)
        left=[' '*base.width()]*exponent.height()+base.picture
        right=exponent.picture+[' '*exponent.width()]*base.height()
        text='\n'.join(a+b for a,b in zip(left,right));baseline=exponent.height()+base.baseline
        return prettyForm(text,baseline,binding=self.POW)
    @staticmethod
    def apply(function,*args):
        arguments=stringPict(*stringPict.next(*[part for i,arg in enumerate(args) for part in ([', ',arg] if i else [arg])]))
        return prettyForm(*stringPict.next(function,stringPict(*arguments.parens())),binding=prettyForm.FUNC)

__all__=['stringPict','prettyForm']
