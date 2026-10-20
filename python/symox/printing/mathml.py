"""Content and presentation MathML emitted directly from native expression trees."""
from xml.etree.ElementTree import Element,SubElement,tostring
from ._legacy import LegacyPrinter as Printer


def _tag(name,*children,text=None,**attributes):
    node=Element(name,attributes)
    if text is not None:node.text=str(text)
    node.extend(children)
    return node


def _pairs(expr):
    args=expr.args
    if args and isinstance(args[0],(tuple,list)):return args
    return list(zip(args[::2],args[1::2]))


class MathMLPrinterBase(Printer):
    presentation=False
    def doprint(self,expr):
        return tostring(self._node(expr),encoding='unicode',short_empty_elements=True).replace(' />','/>')
    def _node(self,expr):
        from symox import sympify,Rational
        if isinstance(expr,(list,tuple)):
            return _tag('list' if not self.presentation else 'mrow',*(self._node(v) for v in expr))
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):
            if not self.presentation:return _tag('matrix',*(_tag('matrixrow',*(self._node(v) for v in row)) for row in expr.tolist()))
            rows=[_tag('mtr',*(_tag('mtd',self._node(v)) for v in row)) for row in expr.tolist()]
            return _tag('mrow',_tag('mo',text='['),_tag('mtable',*rows),_tag('mo',text=']'))
        if not hasattr(expr,'class_name'):expr=sympify(expr)
        if self.presentation:return self._presentation(expr)
        name=expr.class_name
        if expr.is_Symbol:return _tag('ci',text=self._settings.get('symbol_names',{}).get(expr,expr.name))
        if expr.is_Integer or expr.is_Float:return _tag('cn',text=str(expr))
        if expr.is_Rational:return _tag('apply',_tag('divide'),_tag('cn',text=expr.p),_tag('cn',text=expr.q))
        constants={'Pi':'pi','Exp1':'exponentiale','ImaginaryUnit':'imaginaryi','Infinity':'infinity','BooleanTrue':'true','BooleanFalse':'false','NaN':'notanumber'}
        if name in constants:return _tag(constants[name])
        if name=='NegativeInfinity':return _tag('apply',_tag('minus'),_tag('infinity'))
        if not expr.args:return _tag('ci',text=str(expr))
        if expr.is_Add:return _tag('apply',_tag('plus'),*(self._node(v) for v in expr.as_ordered_terms()))
        if expr.is_Mul:return _tag('apply',_tag('times'),*(self._node(v) for v in expr.as_ordered_factors()))
        if expr.is_Pow:
            b,e=expr.args
            if e==Rational(1,2):return _tag('apply',_tag('root'),self._node(b))
            return _tag('apply',_tag('power'),self._node(b),self._node(e))
        if name=='Piecewise':
            return _tag('piecewise',*(_tag('otherwise',self._node(v)) if str(c) in ('True','true') else _tag('piece',self._node(v),self._node(c)) for v,c in _pairs(expr)))
        if name in ('Integral','Sum','Product'):
            body=expr.args[0]
            result=self._node(body)
            for bound in expr.args[1:]:
                var,*limits=bound.args
                parts=[_tag('bvar',self._node(var))]
                if len(limits)==2:parts += [_tag('lowlimit',self._node(limits[0])),_tag('uplimit',self._node(limits[1]))]
                elif len(limits)==1:parts += [_tag('uplimit',self._node(limits[0]))]
                result=_tag('apply',_tag({'Integral':'int','Sum':'sum','Product':'product'}[name]),*parts,result)
            return result
        if name=='Derivative':
            bounds=[]
            for var,count in expr.variable_count:
                nodes=[self._node(var)]
                if count!=1:nodes.append(_tag('degree',_tag('cn',text=count)))
                bounds.append(_tag('bvar',*nodes))
            return _tag('apply',_tag('diff'),*bounds,self._node(expr.expr))
        mapping={'Equality':'eq','Unequality':'neq','StrictLessThan':'lt','LessThan':'leq','StrictGreaterThan':'gt','GreaterThan':'geq',
                 'And':'and','Or':'or','Not':'not','Xor':'xor','Implies':'implies','Abs':'abs','log':'ln','ceiling':'ceiling','Min':'min','Max':'max'}
        head=mapping.get(name,name if name in ('sin','cos','tan','sinh','cosh','tanh','exp','floor','factorial','conjugate') else None)
        return _tag('apply',_tag(head) if head else _tag('ci',text=name),*(self._node(a) for a in expr.args))

    def _presentation(self,expr):
        from symox import Rational
        name=expr.class_name
        if expr.is_Symbol:return _tag('mi',text=self._settings.get('symbol_names',{}).get(expr,expr.name))
        if expr.is_Integer or expr.is_Float:return _tag('mn',text=str(expr))
        if expr.is_Rational:return _tag('mfrac',_tag('mn',text=expr.p),_tag('mn',text=expr.q))
        constants={'Pi':'π','Exp1':'e','ImaginaryUnit':'i','Infinity':'∞','BooleanTrue':'true','BooleanFalse':'false','NaN':'NaN'}
        if name in constants:return _tag('mi',text=constants[name])
        if expr.is_Pow:
            b,e=expr.args
            if e==Rational(1,2):return _tag('msqrt',self._node(b))
            return _tag('msup',self._node(b),self._node(e))
        if expr.is_Add or expr.is_Mul:
            if expr.is_Mul:
                numerator,denominator=expr.as_numer_denom()
                if denominator!=1:return _tag('mfrac',self._node(numerator),self._node(denominator))
            args=expr.as_ordered_terms() if expr.is_Add else expr.as_ordered_factors()
            children=[]
            for arg in args:
                if children:children.append(_tag('mo',text='+' if expr.is_Add else '·'))
                node=self._node(arg)
                if expr.is_Mul and arg.is_Add:node=_tag('mrow',_tag('mo',text='('),node,_tag('mo',text=')'))
                children.append(node)
            return _tag('mrow',*children)
        relations={'Equality':'=','Unequality':'≠','StrictLessThan':'<','LessThan':'≤','StrictGreaterThan':'>','GreaterThan':'≥'}
        if name in relations:return _tag('mrow',self._node(expr.lhs),_tag('mo',text=relations[name]),self._node(expr.rhs))
        children=[]
        for arg in expr.args:
            if children:children.append(_tag('mo',text=','))
            children.append(self._node(arg))
        return _tag('mrow',_tag('mi',text=name),_tag('mo',text='('),*children,_tag('mo',text=')'))


class MathMLContentPrinter(MathMLPrinterBase):pass
class MathMLPresentationPrinter(MathMLPrinterBase):presentation=True


def mathml(expr,printer='content',**settings):
    if printer not in ('content','presentation'):raise ValueError('printer must be content or presentation')
    return (MathMLContentPrinter if printer=='content' else MathMLPresentationPrinter)(settings).doprint(expr)


def print_mathml(expr,printer='content',**settings):
    print(mathml(expr,printer=printer,**settings))


class MathMLPrinter(MathMLPrinterBase):
    """Extensible content MathML renderer over native scalar expressions."""
