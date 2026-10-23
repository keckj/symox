"""Scalar expression emitters for JavaScript, Julia, Octave, R and CAS targets."""
from ._legacy import LegacyCXX17CodePrinter as CXX17CodePrinter
from .precedence import PRECEDENCE

_CONSTANTS={
 'javascript':{'Pi':'Math.PI','Exp1':'Math.E','Infinity':'Infinity','NegativeInfinity':'-Infinity','NaN':'NaN'},
 'julia':{'Pi':'pi','Exp1':'exp(1)','ImaginaryUnit':'im','Infinity':'Inf','NegativeInfinity':'-Inf','NaN':'NaN'},
 'octave':{'Pi':'pi','Exp1':'exp(1)','ImaginaryUnit':'1i','Infinity':'inf','NegativeInfinity':'-inf','NaN':'NaN'},
 'r':{'Pi':'pi','Exp1':'exp(1)','ImaginaryUnit':'1i','Infinity':'Inf','NegativeInfinity':'-Inf','NaN':'NaN'},
 'maple':{'Pi':'Pi','Exp1':'exp(1)','ImaginaryUnit':'I','Infinity':'infinity','NegativeInfinity':'-infinity','NaN':'undefined'},
 'mathematica':{'Pi':'Pi','Exp1':'E','ImaginaryUnit':'I','Infinity':'Infinity','NegativeInfinity':'-Infinity','NaN':'Indeterminate'},
 'glsl':{'Pi':'3.141592653589793','Exp1':'2.718281828459045'},
}
_COMMON=set('sin cos tan asin acos atan atan2 sinh cosh tanh asinh acosh atanh exp log floor ceiling Abs sign Min Max'.split())
_EXTRA={'julia':{'gamma','loggamma','factorial','binomial','erf','erfc'},
        'octave':{'gamma','loggamma','factorial','erf','erfc','cot','sec','csc'},
        'r':{'gamma','loggamma','factorial','choose'},
        'maple':{'gamma','loggamma','factorial','binomial','erf','erfc','cot','sec','csc'},
        'mathematica':{'gamma','loggamma','factorial','binomial','erf','erfc','cot','sec','csc'},
        'javascript':set(),'glsl':set()}
_MAP={'Abs':'abs','ceiling':'ceil','Min':'min','Max':'max'}
_MMA={'Abs':'Abs','ceiling':'Ceiling','Min':'Min','Max':'Max','loggamma':'LogGamma',
      'asin':'ArcSin','acos':'ArcCos','atan':'ArcTan','atan2':'ArcTan',
      'asinh':'ArcSinh','acosh':'ArcCosh','atanh':'ArcTanh'}
_RESERVED=set('if else for while do break continue return function class var let const new in import export default switch case true false null'.split())


class ScalarCodePrinter(CXX17CodePrinter):
    target='javascript'
    def __init__(self,settings=None,**kwargs):
        super().__init__(settings,**kwargs)
        self._user_functions=dict(self._settings.get('user_functions',{}))
        self.known_functions={}
        module_name={'javascript':'jscode','r':'rcode'}.get(self.target,self.target)
        if module_name in ('jscode','rcode','julia','octave','mathematica'):
            from importlib import import_module
            module=import_module('symox.printing.'+module_name)
            self.known_functions.update(getattr(module,'known_functions',{}))
            self.known_functions.update({name:name for name in getattr(module,'known_fcns_src1',())})
            self.known_functions.update(getattr(module,'known_fcns_src2',{}))
        self.known_functions.update(self._user_functions)
        if self.target=='julia' and 'atan2' not in self._user_functions:self.known_functions['atan2']='atan'

    def doprint(self,expr,assign_to=None):
        from symox import sympify
        from symox.codegen.ast import Assignment,CodeBlock,Comment
        if not isinstance(expr,(Assignment,CodeBlock,Comment)) and not hasattr(expr,'args'):
            expr=sympify(expr)
        emitted={}
        for symbol in getattr(expr,'free_symbols',()):
            name=self._print(symbol)
            if name in emitted and emitted[name]!=symbol:raise ValueError('Ambiguous emitted symbol name: '+name)
            emitted[name]=symbol
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):
            rows=[[self._print(value) for value in row] for row in expr.tolist()]
            if self.target=='mathematica':text='{'+', '.join('{'+', '.join(row)+'}' for row in rows)+'}'
            elif self.target=='javascript':text='['+', '.join('['+', '.join(row)+']' for row in rows)+']'
            elif self.target in ('julia','octave'):text='['+'; '.join(' '.join(row) for row in rows)+']'
            elif self.target=='r':text='matrix(c('+', '.join(value for row in rows for value in row)+'), nrow='+str(expr.rows)+', byrow=TRUE)'
            else:raise NotImplementedError('Matrix output is unsupported for '+self.target)
        else:text=self._print(expr)
        if assign_to is not None:
            lhs=assign_to if isinstance(assign_to,str) else self._print(assign_to)
            text=str(lhs)+(' <- ' if self.target=='r' else ' = ')+text
            if self.target in ('javascript','octave','glsl','maple'):text+=';'
        if not self._settings.get('human',True):return set(),set(),text
        return text

    def _print_Symbol(self,expr):
        name=expr.name
        if not name.isidentifier():raise ValueError('Invalid identifier: '+repr(name))
        if name in _RESERVED:
            if self._settings.get('error_on_reserved',False):raise ValueError('Reserved identifier: '+name)
            name+=self._settings.get('reserved_word_suffix','_')
        return name

    def _print_Dummy(self,expr):
        return self._print_Symbol(expr)+'_'+str(expr.dummy_index)

    def _print_Rational(self,expr):
        if self.target in ('r','glsl'):return f'{expr.p}.0/{expr.q}.0'
        if self.target=='julia':return f'{expr.p}//{expr.q}'
        return f'{expr.p}/{expr.q}'

    def _print_Float(self,expr):return str(expr)

    def _print_Pow(self,expr):
        from symox import Rational
        base,exponent=expr.args
        b=self._print(base);e=self._print(exponent)
        if self.target=='octave' and exponent==-1:
            return '1./'+self.parenthesize(base,PRECEDENCE['Pow'])
        if exponent==Rational(1,2) and self.target!='mathematica':
            return ('Math.' if self.target=='javascript' else '')+'sqrt('+b+')'
        if self.target=='javascript':return 'Math.pow('+b+', '+e+')'
        if self.target=='glsl':
            if exponent.is_Integer:e=str(int(exponent))+'.0'
            return 'pow('+b+', '+e+')'
        power={'julia':' .^ ','octave':'.^'}.get(self.target,'^')
        return self.parenthesize(base,PRECEDENCE['Pow'])+power+self.parenthesize(exponent,PRECEDENCE['Pow'])

    def _print_Mul(self,expr):
        # Retain the native ordered numerator/denominator decomposition and
        # choose scalar or elementwise operators for each target.
        args=expr.as_ordered_factors() if self.order not in ('none','old') else expr.args
        numerator=[];denominator=[];sign=''
        for arg in args:
            if self.target=='octave' and arg.is_Number and arg.is_negative:
                sign='' if sign else '-';arg=-arg
                if arg==1:continue
            if getattr(arg,'is_Pow',False) and arg.exp.is_Rational and arg.exp.is_negative:
                denominator.append(arg.base**(-arg.exp))
            else:numerator.append(arg)
        multiply={'octave':'.*','julia':' .* '}.get(self.target,'*')
        divide={'octave':'./','julia':' ./ '}.get(self.target,'/')
        top=multiply.join(self.parenthesize(v,PRECEDENCE['Mul'],strict=True) for v in numerator) or '1'
        if not denominator:return sign+top
        bottom=multiply.join(self.parenthesize(v,PRECEDENCE['Mul'],strict=True) for v in denominator)
        if len(denominator)>1:bottom='('+bottom+')'
        return sign+top+divide+bottom

    def _print_Assignment(self,expr):
        return self.doprint(expr.rhs,expr.lhs)

    def _print_CodeBlock(self,expr):return '\n'.join(self._print(a) for a in expr.args)
    def _print_Comment(self,expr):
        if self.target in ('julia','r','maple'):return '# '+expr.text
        if self.target=='octave':return '% '+expr.text
        if self.target=='mathematica':return '(* '+expr.text+' *)'
        return '// '+expr.text

    def _print_Piecewise(self,expr):
        from .mathml import _pairs
        pairs=_pairs(expr)
        if not pairs or str(pairs[-1][1]) not in ('True','true'):
            raise ValueError('Piecewise code needs a final True branch')
        default=self._print(pairs[-1][0])
        if self.target=='mathematica':
            return 'Piecewise[{'+', '.join('{'+self._print(v)+', '+self._print(c)+'}' for v,c in pairs[:-1])+'}, '+default+']'
        for value,condition in reversed(pairs[:-1]):
            v,c=self._print(value),self._print(condition)
            if self.target in ('javascript','glsl','julia'):default='(('+c+') ? ('+v+') : ('+default+'))'
            else:default='piecewise('+c+', '+v+', '+default+')' if self.target=='maple' else 'ifelse('+c+', '+v+', '+default+')'
        return default

    def _relation(self,expr,op):
        if self.target=='maple':op={'==':'=','!=':'<>'}.get(op,op)
        return self.parenthesize(expr.lhs,PRECEDENCE['Relational'])+' '+op+' '+self.parenthesize(expr.rhs,PRECEDENCE['Relational'])
    def _print_Equality(self,expr):return self._relation(expr,'==')
    def _print_Unequality(self,expr):return self._relation(expr,'!=')
    def _print_StrictLessThan(self,expr):return self._relation(expr,'<')
    def _print_LessThan(self,expr):return self._relation(expr,'<=')
    def _print_StrictGreaterThan(self,expr):return self._relation(expr,'>')
    def _print_GreaterThan(self,expr):return self._relation(expr,'>=')
    def _print_BooleanTrue(self,expr):return {'r':'TRUE','mathematica':'True'}.get(self.target,'true')
    def _print_BooleanFalse(self,expr):return {'r':'FALSE','mathematica':'False'}.get(self.target,'false')
    def _print_And(self,expr):return '('+' && '.join(self._print(a) for a in expr.args)+')'
    def _print_Or(self,expr):return '('+' || '.join(self._print(a) for a in expr.args)+')'
    def _print_Not(self,expr):return '!('+self._print(expr.args[0])+')'

    def emptyPrinter(self,expr):
        name=getattr(expr,'class_name',type(expr).__name__)
        if name in _CONSTANTS[self.target]:return _CONSTANTS[self.target][name]
        if name in ('isnan','isinf'):
            # The C99 classification BooleanFunctions print as calls in
            # every reference language printer (round parentheses included).
            return name+'('+', '.join(self._print(a) for a in expr.args)+')'
        if not getattr(expr,'is_Function',False):
            raise NotImplementedError(f'{name} is unsupported by the {self.target} printer')
        if self.target=='mathematica' and name in ('Sqrt','Cbrt','exp2','log2','log10','expm1','log1p','fma','hypot'):
            # MCodePrinter._print_Function prints a function it does not know
            # as ``Name[args]``.
            return name+'['+', '.join(self._print(a) for a in expr.args)+']'
        target=(self._user_functions if getattr(expr,'is_UndefinedFunction',False) else self.known_functions).get(name)
        if target is None and (getattr(expr,'is_UndefinedFunction',False) or name not in _COMMON|_EXTRA[self.target]):
            if name in ('cot','sec','csc') and self.target in ('javascript','glsl'):
                fn={'cot':'tan','sec':'cos','csc':'sin'}[name]
                return '1.0/'+('Math.' if self.target=='javascript' else '')+fn+'('+self._print(expr.args[0])+')'
            if self._settings.get('allow_unknown_functions',False) or not self._settings.get('strict',True):target=name
            else:raise NotImplementedError(f'{name} requires a user_functions mapping for {self.target}')
        if target is not None and not isinstance(target,str):
            target=next((fn for predicate,fn in target if predicate(*expr.args)),None)
            if target is None:raise ValueError('No matching user function for '+name)
        if target is None:
            if self.target=='mathematica':target=_MMA.get(name,name[0].upper()+name[1:])
            elif self.target=='maple':target={'asin':'arcsin','acos':'arccos','atan':'arctan','ceiling':'ceil','Abs':'abs','gamma':'GAMMA','loggamma':'lnGAMMA'}.get(name,_MAP.get(name,name))
            elif self.target=='r':target={'loggamma':'lgamma','ceiling':'ceiling','Min':'min','Max':'max'}.get(name,_MAP.get(name,name))
            elif self.target=='octave':target={'loggamma':'gammaln'}.get(name,_MAP.get(name,name))
            elif self.target in ('glsl','julia'):target={'atan2':'atan'}.get(name,_MAP.get(name,name))
            else:target=_MAP.get(name,name)
            if self.target=='javascript':target='Math.'+target
        args=list(expr.args)
        if self.target=='mathematica' and name=='atan2':args.reverse()
        left,right=('[',']') if self.target=='mathematica' else ('(',')')
        rendered=[self._print(a) for a in args]
        if self.target in ('glsl','octave') and name in ('Min','Max') and len(rendered)>2:
            result=rendered[0]
            for value in rendered[1:]:result=target+'('+result+', '+value+')'
            return result
        return target+left+', '.join(rendered)+right


class JavascriptCodePrinter(ScalarCodePrinter):target='javascript';printmethod='_javascript'
class JuliaCodePrinter(ScalarCodePrinter):target='julia';printmethod='_julia'
class OctaveCodePrinter(ScalarCodePrinter):target='octave';printmethod='_octave'
class RCodePrinter(ScalarCodePrinter):target='r';printmethod='_rcode'
class MapleCodePrinter(ScalarCodePrinter):target='maple';printmethod='_maple'
class MCodePrinter(ScalarCodePrinter):
    target='mathematica';printmethod='_mcode'
    def _print_Mul(self,expr):
        # The reference CodePrinter._print_Mul, under MCodePrinter's
        # commutative/noncommutative split.
        from symox import Mul,Pow,S
        from symox.core.mul import _keep_coeff
        from .precedence import precedence
        PREC=precedence(expr)
        c,nc=expr.args_cnc()
        expr=Mul(*c);prec=precedence(expr)
        coeff,rest=expr.as_coeff_Mul()
        if coeff<0:expr=_keep_coeff(-coeff,rest);sign='-'
        else:sign=''
        a=[];b=[];pow_paren=[]
        args=expr.as_ordered_factors() if self.order not in ('old','none') else Mul.make_args(expr)
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp!=-1:b.append(Pow(item.base,-item.exp,evaluate=False))
                else:
                    if len(item.args[0].args)!=1 and isinstance(item.base,Mul):pow_paren.append(item)
                    b.append(Pow(item.base,-item.exp))
            else:a.append(item)
        a=a or [S.One]
        if len(a)==1 and sign=='-':a_str=[self.parenthesize(a[0],0.5*(PRECEDENCE['Pow']+PRECEDENCE['Mul']))]
        else:a_str=[self.parenthesize(x,prec) for x in a]
        b_str=[self.parenthesize(x,prec) for x in b]
        for item in pow_paren:
            if item.base in b:b_str[b.index(item.base)]='(%s)'%b_str[b.index(item.base)]
        if not b:res=sign+'*'.join(a_str)
        elif len(b)==1:res=sign+'*'.join(a_str)+'/'+b_str[0]
        else:res=sign+'*'.join(a_str)+'/(%s)'%'*'.join(b_str)
        if nc:res+='*'+'**'.join(self.parenthesize(x,PREC) for x in nc)
        return res
class GLSLPrinter(ScalarCodePrinter):target='glsl';printmethod='_glsl'
