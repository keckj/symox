"""Executable statement rendering shared by compatible language printers."""
import json


def indent(text):return '\n'.join('    '+line if line else '' for line in str(text).splitlines())
def _none(value):return value is None or type(value).__name__=='NoneToken'


class NodePrinter:
    _style='c'
    def _print_String(self,node):return str(node.text)
    def _print_QuotedString(self,node):return json.dumps(node.text)
    def _print_NoneToken(self,node):return {'python':'None','rust':'()','fortran':''}.get(self._style,'NULL')
    def _print_Variable(self,node):return self._print(node.symbol)
    _print_Pointer=_print_Variable
    def _print_Type(self,node):
        name=str(node.name)
        if self._style=='python':return {'integer':'int','real':'float','complex':'complex','bool':'bool'}.get(name,'float')
        if self._style=='fortran':return {'integer':'integer','intc':'integer(c_int)','real':'real(kind=8)','float32':'real(kind=4)','float64':'real(kind=8)','float80':'real(kind=16)','complex':'complex(kind=8)','complex64':'complex(kind=4)','complex128':'complex(kind=8)','bool':'logical','untyped':'real(kind=8)'}.get(name,'integer(kind=8)' if name.endswith('64') else 'integer')
        if self._style=='rust':return {'integer':'i32','intc':'i32','real':'f64','float32':'f32','float64':'f64','complex':'num_complex::Complex64','bool':'bool','untyped':'f64','void':'()','int8':'i8','int16':'i16','int32':'i32','int64':'i64','uint8':'u8','uint16':'u16','uint32':'u32','uint64':'u64'}.get(name,name)
        return {'integer':'int','intc':'int','real':'double','float32':'float','float64':'double','float80':'long double','complex':'double complex','complex64':'float complex','complex128':'double complex','bool':'bool','untyped':'double','void':'void','int8':'int8_t','int16':'int16_t','int32':'int32_t','int64':'int64_t','uint8':'uint8_t','uint16':'uint16_t','uint32':'uint32_t','uint64':'uint64_t'}.get(name,name)
    _print_IntBaseType=_print_SignedIntType=_print_UnsignedIntType=_print_FloatBaseType=_print_FloatType=_print_ComplexBaseType=_print_ComplexType=_print_Type
    def _declaration(self,variable,parameter=False):
        from ..codegen.ast import Variable,untyped
        if not isinstance(variable,Variable):variable=Variable(variable)
        symbol=self._print(variable.symbol);typecode=self._print(variable.type);value=variable.value
        attrs={a.name for a in variable.attrs};pointer=type(variable).__name__=='Pointer'
        if self._style=='python':return symbol if parameter else symbol+' = '+('None' if _none(value) else self._print(value))
        if self._style=='rust':
            if parameter:return symbol+': '+('*mut ' if pointer else '')+typecode
            return 'let '+('' if 'value_const' in attrs else 'mut ')+symbol+': '+typecode+('' if _none(value) else ' = '+self._print(value))+';'
        if self._style=='fortran':
            extra=[]
            for attr in variable.attrs:
                if attr.name=='dimension':extra.append('dimension('+', '.join(map(str,attr.parameters))+')')
                elif attr.name.startswith('intent_'):extra.append('intent('+attr.name[7:]+')')
                elif attr.name=='intent':extra.append('intent('+str(attr.parameters[0])+')')
                elif attr.name in ('allocatable','pointer','target','optional','value','parameter'):extra.append(attr.name)
            text=typecode+(', '+', '.join(extra) if extra else '')+' :: '+symbol
            return text+(' = '+self._print(value) if not _none(value) else '')
        prefix=('const ' if 'value_const' in attrs else '')+typecode+' '+('*' if pointer else '')+symbol
        return prefix+(' = '+self._print(value) if not _none(value) else '')+('' if parameter else ';')
    def _print_Declaration(self,node):return self._declaration(node.variable)
    def _print_Element(self,node):
        symbol=self._print(node.symbol);indices=node.indices
        if self._style=='fortran':return symbol+'('+', '.join(self._print(i) for i in indices)+')'
        if self._style=='python':return symbol+'['+', '.join(self._print(i) for i in indices)+']'
        if not _none(node.strides):index=sum(i*s for i,s in zip(indices,node.strides))
        elif len(indices)==1:index=indices[0]
        else:raise ValueError('flattened array code requires strides for multidimensional elements')
        if not _none(node.offset):index+=node.offset
        return symbol+'['+self._print(index)+']'
    def _print_FunctionCall(self,node):return str(node.name)+'('+', '.join(self._print(v) for v in node.function_args)+')'
    def _print_AssignmentBase(self,node):return self._print_Assignment(node)
    def _print_AugmentedAssignment(self,node):
        lhs,rhs=self._print(node.lhs),self._print(node.rhs)
        if self._style=='fortran':return lhs+' = '+lhs+' '+node.binop+' '+rhs
        return lhs+' '+node.binop+'= '+rhs+('' if self._style=='python' else ';')
    _print_AddAugmentedAssignment=_print_SubAugmentedAssignment=_print_MulAugmentedAssignment=_print_DivAugmentedAssignment=_print_ModAugmentedAssignment=_print_AugmentedAssignment
    def _print_Return(self,node):return 'return '+self._print(node.return_)+('' if self._style in ('python','fortran') else ';')
    def _print_BreakToken(self,node):return {'python':'break','fortran':'exit'}.get(self._style,'break;')
    def _print_ContinueToken(self,node):return {'python':'continue','fortran':'cycle'}.get(self._style,'continue;')
    def _print_While(self,node):
        condition=self._print(node.condition);body=indent(self._print(node.body))
        if self._style=='python':return 'while '+condition+':\n'+(body or '    pass')
        if self._style=='fortran':return 'do while ('+condition+')\n'+body+'\nend do'
        return 'while ('+condition+') {\n'+body+'\n}'
    def _print_For(self,node):
        target=self._print(node.target);iterator=node.iterable;body=indent(self._print(node.body))
        if hasattr(iterator,'start') and hasattr(iterator,'stop'):
            start,stop,step=iterator.start,iterator.stop,iterator.step
        else:
            if self._style!='python':raise TypeError('generated loops require a range')
            return 'for '+target+' in '+self._print(iterator)+':\n'+body
        a,b,c=map(lambda v:self._print(v) if hasattr(v,'args') else str(v),(start,stop,step))
        if self._style=='python':return 'for '+target+' in range('+a+', '+b+', '+c+'):\n'+(body or '    pass')
        if self._style=='fortran':return 'do '+target+' = '+a+', '+self._print(stop-step)+', '+c+'\n'+body+'\nend do'
        if self._style=='rust':
            if step==1:return 'for '+target+' in '+a+'..'+b+' {\n'+body+'\n}'
            return '{\n    let mut '+target+' = '+a+';\n    while '+target+(' < ' if step>0 else ' > ')+b+' {\n'+indent(body)+'\n        '+target+' += '+c+';\n    }\n}'
        return 'for ('+target+' = '+a+'; '+target+(' < ' if step>0 else ' > ')+b+'; '+target+' += '+c+') {\n'+body+'\n}'
    def _print_Scope(self,node):
        body=self._print(node.body)
        return body if self._style in ('python','fortran') else '{\n'+indent(body)+'\n}'
    def _prototype(self,node):
        parameters=', '.join(self._declaration(v,True) for v in node.parameters)
        if self._style=='python':return 'def '+str(node.name)+'('+parameters+')'
        if self._style=='rust':return 'fn '+str(node.name)+'('+parameters+') -> '+self._print(node.return_type)
        return self._print(node.return_type)+' '+str(node.name)+'('+parameters+')'
    def _print_FunctionPrototype(self,node):return self._prototype(node)+(':' if self._style=='python' else ';')
    def _print_FunctionDefinition(self,node):
        if self._style=='fortran':
            name=str(node.name);params=', '.join(self._print(getattr(v,'symbol',v)) for v in node.parameters)
            # Fortran assigns a function result instead of returning a value.
            from ..codegen.ast import Return,Assignment
            statements=[Assignment(name,a.return_) if isinstance(a,Return) else a for a in node.body.args]
            declarations='\n'.join(self._declaration(v,True) for v in node.parameters)
            return self._print(node.return_type)+' function '+name+'('+params+')\n'+indent(declarations)+'\n'+indent('\n'.join(self._print(a) for a in statements))+'\nend function '+name
        body=indent(self._print(node.body));head=self._prototype(node)
        return head+(':\n'+(body or '    pass') if self._style=='python' else ' {\n'+body+'\n}')
    def _print_Print(self,node):
        args=', '.join(self._print(v) for v in node.print_args)
        if self._style=='python':return 'print('+args+('' if _none(node.file) else ', file='+str(node.file.name))+')'
        if self._style=='fortran':return 'print *, '+args
        fmt=json.dumps(node.format_string.text if hasattr(node.format_string,'text') else str(node.format_string)) if not _none(node.format_string) else json.dumps(' '.join(['%g']*len(node.print_args))+'\n')
        if self._style=='rust':return 'println!('+json.dumps(' '.join(['{}']*len(node.print_args)))+(', '+args if args else '')+');'
        return ('printf('+fmt if _none(node.file) else 'fprintf('+str(node.file.name)+', '+fmt)+(', '+args if args else '')+');'
    def _print_RuntimeError_(self,node):return json.dumps(str(node.message))
    def _print_Raise(self,node):
        message=self._print(node.exception)
        if self._style=='python':return 'raise RuntimeError('+message+')'
        if self._style=='rust':return 'panic!('+message+');'
        if self._style=='fortran':return 'error stop '+message
        return 'fprintf(stderr, "%s\\n", '+message+'); abort();'

    def _print_PreIncrement(self,node):return '++'+self._print(node.arg)
    def _print_PreDecrement(self,node):return '--'+self._print(node.arg)
    def _print_PostIncrement(self,node):return self._print(node.arg)+'++'
    def _print_PostDecrement(self,node):return self._print(node.arg)+'--'
    def _print_CommaOperator(self,node):return '('+', '.join(self._print(a) for a in node.args)+')'
    def _print_Label(self,node):return str(node.name)+':'+('' if _none(node.body) else '\n'+self._print(node.body))
    def _print_goto(self,node):return 'goto '+str(node.label)+';'
    def _print_sizeof(self,node):return 'sizeof('+self._print(node.arg)+')'
    def _print_alignof(self,node):return '_Alignof('+self._print(node.arg)+')'
    def _print_struct(self,node):return 'struct '+('' if _none(node.name) else str(node.name))+' {\n'+indent('\n'.join(self._declaration(getattr(v,'variable',v)) for v in node.declarations))+'};'
    def _print_union(self,node):return self._print_struct(node).replace('struct ','union ',1)
    def _print_using(self,node):return 'using '+str(node.alias)+' = '+self._print(node.type)+';'
    def _print_ArrayConstructor(self,node):return '['+', '.join(self._print(v) for v in node.elements)+']'
    def _print_Do(self,node):
        return ('do concurrent (' if node.concurrent else 'do ')+self._print(node.counter)+(' = ')+self._print(node.first)+(':' if node.concurrent else ', ')+self._print(node.last)+(':' if node.concurrent else ', ')+self._print(node.step)+(')' if node.concurrent else '')+'\n'+indent(self._print(node.body))+'\nend do'
    def _print_ImpliedDoLoop(self,node):return '('+self._print(node.expr)+', '+self._print(node.counter)+' = '+self._print(node.first)+', '+self._print(node.last)+', '+self._print(node.step)+')'
    def _print_FortranReturn(self,node):return 'return'+('' if _none(node.return_value) else ' '+self._print(node.return_value))
    def _print_GoTo(self,node):return 'go to '+('('+', '.join(map(str,node.labels))+'), '+self._print(node.expr) if not _none(node.expr) else str(node.labels[0]))
    def _print_Module(self,node):
        return 'module '+str(node.name)+'\n'+indent('\n'.join(self._print(d) for d in node.declarations))+'\ncontains\n'+indent('\n'.join(self._print(d) for d in node.definitions))+'\nend module '+str(node.name)
    def _print_Program(self,node):return 'program '+str(node.name)+'\n'+indent(self._print(node.body))+'\nend program '+str(node.name)
    def _print_Subroutine(self,node):
        head='subroutine '+str(node.name)+'('+', '.join(self._print(getattr(v,'symbol',v)) for v in node.parameters)+')'
        declarations='\n'.join(self._declaration(v,True) for v in node.parameters)
        return head+'\n'+indent(declarations)+'\n'+indent(self._print(node.body))+'\nend subroutine '+str(node.name)
    def _print_SubroutineCall(self,node):return 'call '+str(node.name)+'('+', '.join(self._print(v) for v in node.subroutine_args)+')'
    def _print_use(self,node):
        text='use '+str(node.namespace)
        if not _none(node.only):text+=', only: '+', '.join(self._print(v) if hasattr(v,'args') else str(v) for v in node.only)
        if not _none(node.rename):text+=', '+', '.join(self._print(v) for v in node.rename)
        return text
    def _print_use_rename(self,node):return str(node.local)+' => '+str(node.original)
    def _print_sum_(self,node):return 'sum('+self._print(node.array)+('' if _none(node.dim) else ', dim='+self._print(node.dim))+('' if _none(node.mask) else ', mask='+self._print(node.mask))+')'
    def _print_product_(self,node):return self._print_sum_(node).replace('sum(','product(',1)
