"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'AtomicExpr Expr UnevaluatedExpr factor_terms')


def unchanged(func,*args):
    result=func(*args)
    same_head=result.func==func or (not getattr(result,'is_UndefinedFunction',False) and getattr(func,'_symox_head_name',None)==getattr(result,'class_name',object()))
    return same_head and result.args==args

class ExprBuilder:
    def __init__(self,op,args=None,validator=None,check=True):
        if not callable(op):raise TypeError('expression builder operator must be callable')
        self.op,self.args,self.validator=op,list(args or ()),validator
        if check:self.validate()
    def _arguments(self):return [arg.build() if isinstance(arg,ExprBuilder) else arg for arg in self.args]
    def validate(self):
        if self.validator:self.validator(*self._arguments())
    def build(self,check=True):
        args=self._arguments()
        if check and self.validator:self.validator(*args)
        return self.op(*args)
    def append_argument(self,arg,check=True):
        self.args.append(arg)
        if check:self.validate()
    def __getitem__(self,index):return self.op if index==0 else self.args[index-1]
    def search_element(self,elem):
        for i,arg in enumerate(self.args):
            if arg is elem:return (i,)
            if isinstance(arg,ExprBuilder):
                index=arg.search_element(elem)
                if index is not None:return (i,)+index
        return None
    def __repr__(self):return str(self.build())
__all__ += ['ExprBuilder','unchanged']
