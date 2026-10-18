"""Native rust code generation entry point."""
from symox._core import rust_code
__all__ = ['rust_code']

from .codeprinter import CodePrinter
class RustCodePrinter(CodePrinter):
    language='Rust';printmethod='_rust_code';_style='rust'
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self._integer_mode=False
        self._types={}
    def _typed_print(self,expr,integer):
        old=self._integer_mode
        self._integer_mode=integer
        try:return self._print(expr)
        finally:self._integer_mode=old
    def _is_integer_type(self,type):
        from ..codegen.ast import IntBaseType
        return isinstance(type,IntBaseType)
    def _print_Integer(self,expr):return str(int(expr))+('' if self._integer_mode else '.0')
    def _print_Symbol(self,expr):
        name=expr.name
        if not self._integer_mode and self._is_integer_type(self._types.get(expr)):
            return '('+name+' as f64)'
        return name
    def _declaration(self,variable,parameter=False):
        from ..codegen.ast import Variable
        if not isinstance(variable,Variable):variable=Variable(variable)
        self._types[variable.symbol]=variable.type
        old=self._integer_mode
        self._integer_mode=self._is_integer_type(variable.type)
        try:return super()._declaration(variable,parameter)
        finally:self._integer_mode=old
    def _print_FunctionDefinition(self,node):
        old_types,old_mode=self._types,self._integer_mode
        self._types={v.symbol:v.type for v in node.parameters}
        self._integer_mode=self._is_integer_type(node.return_type)
        self._return_integer=self._integer_mode
        try:return super()._print_FunctionDefinition(node)
        finally:self._types,self._integer_mode=old_types,old_mode
    def _print_Return(self,node):
        return 'return '+self._typed_print(node.return_,getattr(self,'_return_integer',False))+';'
    def _print_Assignment(self,node):
        kind=self._is_integer_type(self._types.get(node.lhs))
        return self._typed_print(node.lhs,True)+' = '+self._typed_print(node.rhs,kind)+';'
    def _print_AugmentedAssignment(self,node):
        kind=self._is_integer_type(self._types.get(node.lhs))
        return self._typed_print(node.lhs,True)+' '+node.binop+'= '+self._typed_print(node.rhs,kind)+';'
    _print_AddAugmentedAssignment=_print_SubAugmentedAssignment=_print_MulAugmentedAssignment=_print_DivAugmentedAssignment=_print_ModAugmentedAssignment=_print_AugmentedAssignment
    def _print_For(self,node):
        from ..codegen.ast import integer
        from ._nodes import indent
        iterator=node.iterable
        if not hasattr(iterator,'start'):raise TypeError('Rust loops require a range')
        self._types[node.target]=integer
        target=self._typed_print(node.target,True)
        a,b,c=[self._typed_print(v,True) for v in (iterator.start,iterator.stop,iterator.step)]
        body=indent(self._print(node.body))
        if iterator.step==1:return 'for '+target+' in '+a+'..'+b+' {\n'+body+'\n}'
        return '{\n    let mut '+target+' = '+a+';\n    while '+target+(' < ' if iterator.step>0 else ' > ')+b+' {\n'+indent(body)+'\n        '+target+' += '+c+';\n    }\n}'
    def _print_CodeBlock(self,node):return '\n'.join(self._print(arg) for arg in node.args)
    def _print_Pow(self,expr):
        b,e=expr.args
        if self._integer_mode:
            if not e.is_Integer or e < 0:raise TypeError('An integer Rust power needs a nonnegative integer exponent')
            return '('+self._print(b)+').pow('+str(int(e))+')'
        return '('+self._print(b)+').'+('powi' if e.is_Integer else 'powf')+'('+(str(int(e)) if e.is_Integer else self._print(e))+')'
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Function',False):
            name={'Abs':'abs','ceiling':'ceil','log':'ln'}.get(expr.class_name,expr.class_name)
            if len(expr.args)==1 and name in known_functions:
                return '('+self._print(expr.args[0])+').'+name+'()'
        from symox import _core
        return _core.rust_code(expr)
class TypeCast:
    def __init__(self,value,type):self.value,self.type=value,type
    @property
    def args(self):return self.value,self.type
    def _rust_code(self,printer):return '('+printer._print(self.value)+' as '+str(self.type)+')'
def float_floor(value):
    from symox import floor
    return floor(value)
def float_ceiling(value):
    from symox import ceiling
    return ceiling(value)
known_functions={name:name for name in 'sin cos tan asin acos atan sinh cosh tanh exp ln log10 log2 sqrt cbrt abs floor ceil'.split()}
__all__ += ['RustCodePrinter','TypeCast','float_floor','float_ceiling','known_functions']
