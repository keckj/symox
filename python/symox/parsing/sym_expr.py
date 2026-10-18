"""Numeric C/Fortran source adapters to the local program AST.

The parser accepts declarations, assignments, numeric functions, counted loops
and while loops. Unsupported constructs raise SyntaxError rather than being
silently discarded. This is a deliberately smaller grammar than Clang/LFortran.
"""
from pathlib import Path
import re
import symox as S
from ..codegen import ast as A

_TYPES={'int':A.intc,'integer':A.integer,'long':A.int64,'float':A.float32,'double':A.float64,'real':A.real,'logical':A.bool_,'bool':A.bool_,'void':A.void}
_IDENTIFIER=r'[A-Za-z_]\w*'


def _split(text,delimiter=','):
    result=[];level=0;start=0
    for i,char in enumerate(text):
        if char in '([':level+=1
        elif char in ')]':level-=1
        elif char==delimiter and not level:result.append(text[start:i].strip());start=i+1
    result.append(text[start:].strip());return result


def _balanced(text,start,left='(',right=')'):
    if text[start]!=left:raise SyntaxError('expected '+left)
    depth=0
    for i in range(start,len(text)):
        if text[i]==left:depth+=1
        elif text[i]==right:
            depth-=1
            if not depth:return text[start+1:i],i+1
    raise SyntaxError('unclosed '+left)


def _expression(text,symbols):
    text=text.strip()
    for operator,constructor in [('||',S.Or),('&&',S.And)]:
        if operator in text:return constructor(*(_expression(part,symbols) for part in text.split(operator)))
    if text.startswith('!'):return S.Not(_expression(text[1:],symbols))
    if any(token in text for token in ['?','[',']','{','}','^']):raise SyntaxError('unsupported numeric expression '+text)
    text=re.sub(r'(?<=\d)[fF]\b','',text)
    return S.parse_expr(text,local_dict={**symbols,'true':S.true,'false':S.false})


def _declare(typecode,text,symbols,initialize=False):
    declarations=[]
    for part in _split(text):
        match=re.fullmatch('('+_IDENTIFIER+r')(?:\s*=\s*(.+))?',part,re.S)
        if not match:raise SyntaxError('unsupported declaration '+part)
        name,value=match.groups();symbol=symbols.setdefault(name,S.Symbol(name))
        initial=_expression(value,symbols) if value is not None else S.Integer(0) if initialize else A.none
        if initial is not A.none and typecode in (A.float32,A.float64,A.real) and initial.is_number:initial=S.Float(str(initial))
        declarations.append(A.Declaration(A.Variable(symbol,typecode,initial)))
    return declarations


def _statement(text,symbols):
    text=text.strip()
    declaration=re.fullmatch(r'(?:const\s+)?(int|integer|long|float|double|real|bool|logical)\s+(.+)',text,re.S)
    if declaration:return _declare(_TYPES[declaration[1]],declaration[2],symbols)
    if text.startswith('return '):return [A.Return(_expression(text[7:],symbols))]
    match=re.fullmatch('('+_IDENTIFIER+r')\s*(\+=|-=|\*=|/=|=)\s*(.+)',text,re.S)
    if match:
        name,operator,value=match.groups();symbol=symbols.setdefault(name,S.Symbol(name));rhs=_expression(value,symbols)
        return [A.Assignment(symbol,rhs) if operator=='=' else A.aug_assign(symbol,operator[0],rhs)]
    match=re.fullmatch(r'(?:\+\+|--)?('+_IDENTIFIER+r')(?:\+\+|--)?',text)
    if match and ('++' in text or '--' in text):
        symbol=symbols.setdefault(match[1],S.Symbol(match[1]));return [A.Assignment(symbol,symbol+(1 if '++' in text else -1))]
    call=re.fullmatch('('+_IDENTIFIER+r')\((.*)\)',text,re.S)
    if call:return [A.FunctionCall(call[1],tuple(_expression(part,symbols) for part in _split(call[2]) if part))]
    raise SyntaxError('unsupported numeric statement '+text)


def _c_block(source,symbols):
    statements=[];position=0
    while position<len(source):
        while position<len(source) and source[position].isspace():position+=1
        if position==len(source):break
        tail=source[position:]
        function=re.match(r'(int|long|float|double|bool|void)\s+('+_IDENTIFIER+r')\s*\(',tail)
        if function:
            typecode,name=function[1],function[2];params,end=_balanced(source,position+function.end()-1)
            while end<len(source) and source[end].isspace():end+=1
            body,next_position=_balanced(source,end,'{','}');local=dict(symbols);parameters=[]
            for param in _split(params):
                if not param or param=='void':continue
                declaration=re.fullmatch(r'(int|long|float|double|bool)\s+('+_IDENTIFIER+')',param)
                if not declaration:raise SyntaxError('unsupported function parameter '+param)
                symbol=local.setdefault(declaration[2],S.Symbol(declaration[2]));parameters.append(A.Variable(symbol,_TYPES[declaration[1]]))
            statements.append(A.FunctionDefinition(_TYPES[typecode],name,parameters,A.CodeBlock(*_c_block(body,local))));position=next_position;continue
        loop=re.match(r'(for|while)\s*\(',tail)
        if loop:
            header,end=_balanced(source,position+loop.end()-1)
            while end<len(source) and source[end].isspace():end+=1
            body,next_position=_balanced(source,end,'{','}')
            if loop[1]=='while':node=A.While(_expression(header,symbols),A.CodeBlock(*_c_block(body,symbols)))
            else:
                parts=_split(header,';')
                if len(parts)!=3:raise SyntaxError('for requires initializer, comparison and increment')
                initial=_statement(parts[0],symbols)
                if isinstance(initial[0],A.Declaration):
                    variable=initial[0].variable;start=variable.value;symbol=variable.symbol
                    statements.append(A.Declaration(A.Variable(symbol,variable.type)))
                elif isinstance(initial[0],A.Assignment):symbol,start=initial[0].lhs,initial[0].rhs
                else:raise SyntaxError('unsupported for initializer')
                relation=_expression(parts[1],symbols)
                if not getattr(relation,'is_Relational',False) or relation.lhs!=symbol:raise SyntaxError('counted for requires a direct counter comparison')
                increment=parts[2].replace(' ','')
                if increment in (str(symbol)+'++','++'+str(symbol)):step=S.Integer(1)
                elif increment in (str(symbol)+'--','--'+str(symbol)):step=S.Integer(-1)
                else:
                    inc=re.fullmatch(re.escape(str(symbol))+r'([+-])=(.+)',increment)
                    if not inc:raise SyntaxError('unsupported for increment')
                    step=_expression(inc[2],symbols)*(1 if inc[1]=='+' else -1)
                if step.is_positive is True and relation.rel_op not in ('<','<=') or step.is_negative is True and relation.rel_op not in ('>','>='):raise SyntaxError('counter step and comparison directions differ')
                stop=relation.rhs+(step if relation.rel_op in ('<=','>=') else 0)
                node=A.For(symbol,S.Range(start,stop,step),A.CodeBlock(*_c_block(body,symbols)))
            statements.append(node);position=next_position;continue
        end=source.find(';',position)
        if end<0:raise SyntaxError('C statements require a terminating semicolon')
        if source[position:end].strip():statements.extend(_statement(source[position:end],symbols))
        position=end+1
    return statements


def _fortran_block(lines,symbols,position=0,terminators=()):
    statements=[]
    while position<len(lines):
        line=lines[position].strip();lower=line.lower()
        if any(lower.startswith(end) for end in terminators):return statements,position+1
        function=re.fullmatch(r'(integer|real|double precision|logical)(?:\([^)]*\))?\s+function\s+('+_IDENTIFIER+r')\s*\(([^)]*)\)(?:\s+result\(('+_IDENTIFIER+r')\))?',line,re.I)
        if function:
            typecode,name,params,result=function.groups();typecode='double' if typecode.lower()=='double precision' else typecode.lower()
            local=dict(symbols);names=[part.strip() for part in params.split(',') if part.strip()]
            for arg in names:local[arg]=S.Symbol(arg)
            body,position=_fortran_block(lines,local,position+1,('end function',))
            parameters=[];kept=[];declarations={str(node.variable.symbol):node.variable for node in body if isinstance(node,A.Declaration)}
            for arg in names:parameters.append(A.Variable(local[arg],declarations[arg].type if arg in declarations else A.real))
            kept=[node for node in body if not isinstance(node,A.Declaration) or str(node.variable.symbol) not in names]
            output=result or name;value=local.setdefault(output,S.Symbol(output))
            if output not in declarations:kept.insert(0,A.Declaration(A.Variable(value,_TYPES[typecode],S.Integer(0))))
            kept.append(A.Return(value));statements.append(A.FunctionDefinition(_TYPES[typecode],name,parameters,A.CodeBlock(*kept)));continue
        loop=re.fullmatch(r'do\s+('+_IDENTIFIER+r')\s*=\s*(.+)',line,re.I)
        if loop:
            parts=_split(loop[2]);symbol=symbols.setdefault(loop[1],S.Symbol(loop[1]))
            if len(parts) not in (2,3):raise SyntaxError('Fortran do requires first, last and optional step')
            first,last=map(lambda text:_expression(text,symbols),parts[:2]);step=_expression(parts[2],symbols) if len(parts)==3 else S.Integer(1)
            body,position=_fortran_block(lines,symbols,position+1,('end do','enddo'))
            statements.append(A.For(symbol,S.Range(first,last+step,step),A.CodeBlock(*body)));continue
        declaration=re.fullmatch(r'(integer|real|double precision|logical)(?:\([^)]*\))?(?:\s*,[^:]+)?\s*::\s*(.+)',line,re.I)
        if declaration:
            typecode=declaration[1].lower();statements.extend(_declare(_TYPES['double' if typecode=='double precision' else typecode],declaration[2],symbols,initialize=True))
        elif lower not in ('implicit none',''):statements.extend(_statement(line,symbols))
        position+=1
    if terminators:raise SyntaxError('unclosed Fortran block')
    return statements,position


class SymPyExpression:
    def __init__(self,source_code=None,mode=None):
        self._expr=[]
        if source_code is not None:
            if mode is None:
                if not isinstance(source_code,(list,tuple)):raise ValueError('source text requires a language mode')
                self._expr=list(source_code)
            else:self.convert_to_expr(source_code,mode)
    def convert_to_expr(self,src_code,mode):
        source=str(src_code)
        if '\n' not in source and len(source)<240:
            path=Path(source)
            if path.is_file():source=path.read_text()
        mode=mode.lower()
        if mode in ('c','c99'):
            source=re.sub(r'/\*.*?\*/|//[^\n]*','',source,flags=re.S);result=_c_block(source,{})
        elif mode in ('f','fortran','f95'):
            lines=[line.split('!',1)[0].strip() for line in source.splitlines()];result,_=_fortran_block(lines,{})
        else:raise NotImplementedError('unsupported source language '+mode)
        self._expr=result
    def return_expr(self):return list(self._expr)
    def convert_to_c(self):
        from ..printing.c import C99CodePrinter
        printer=C99CodePrinter();return [printer.doprint(expr) for expr in self._expr]
    def convert_to_fortran(self):
        from ..printing.fortran import FCodePrinter
        printer=FCodePrinter();return [printer.doprint(expr) for expr in self._expr]
    def convert_to_python(self):
        from ..printing.pycode import PythonCodePrinter
        printer=PythonCodePrinter();return [printer.doprint(expr) for expr in self._expr]

__all__=['SymPyExpression']
