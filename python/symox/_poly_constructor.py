"""Public Poly input normalization, followed by native polynomial arithmetic."""
import operator
import re

from ._core import Poly
from .polys.polyerrors import (OptionError,GeneratorsError,GeneratorsNeeded,
                              PolynomialError,MultivariatePolynomialError)


_BOOLEAN={'all','auto','composite','expand','field','formal','frac','gaussian',
          'greedy','include','polys','series','split','strict','symmetric'}
_KNOWN=_BOOLEAN|{'domain','extension','gen','gens','method','modulus','order',
                 'sort','symbols','wrt'}
_EXCLUDES={
    'field':{'domain','split','gaussian'},
    'greedy':{'domain','split','gaussian','extension','modulus','symmetric'},
    'composite':{'domain','split','gaussian','extension','modulus','symmetric'},
    'domain':{'field','greedy','split','gaussian','extension'},
    'split':{'field','greedy','domain','gaussian','extension','modulus','symmetric'},
    'gaussian':{'field','greedy','domain','split','extension','modulus','symmetric'},
    'extension':{'greedy','domain','split','gaussian','modulus','symmetric'},
    'modulus':{'greedy','split','domain','gaussian','extension'},
    'symmetric':{'greedy','domain','split','gaussian','extension'},
}


def _options(gens,domain,modulus,extension,kwargs):
    args=dict(kwargs or {})
    if 'opt' in args and len(args)==1 and not gens:
        args=dict(args['opt'])
    defaults=dict(args.pop('defaults',{}))
    for key,value in [('domain',domain),('modulus',modulus),('extension',extension)]:
        if value is not None:args[key]=value
    for key,value in defaults.items():
        if key not in args and not any(key in _EXCLUDES.get(option,())
                                       for option,v in args.items() if v is not None):
            args[key]=value
    for key in args:
        if key not in _KNOWN:raise OptionError("'%s' is not a valid option"%key)
    args={key:value for key,value in args.items() if value is not None}
    for key in args.keys()&_BOOLEAN:
        if args[key] not in (True,False):
            raise OptionError("'%s' must have a boolean value assigned, got %s"%(key,args[key]))
        args[key]=bool(args[key])
    if gens and args.get('gens',()):
        raise OptionError("both '*gens' and keyword argument 'gens' supplied")
    for key in args:
        for excluded in _EXCLUDES.get(key,()):
            if excluded in args:
                raise OptionError("'%s' option is not allowed together with '%s'"%(key,excluded))
    if 'symmetric' in args and 'modulus' not in args:
        raise OptionError("'symmetric' option is only allowed together with 'modulus'")
    if 'split' in args:raise NotImplementedError("'split' option is not implemented yet")
    if 'order' in args:raise NotImplementedError("'order' keyword is not implemented yet")
    if 'method' in args:
        if not isinstance(args['method'],str):raise OptionError('expected a string')
        args['method']=args['method'].lower()
    if 'symbols' in args:
        try:args['symbols']=iter(args['symbols'])
        except TypeError:raise OptionError('expected an iterator or iterable container') from None
    if 'extension' in args and args['extension'] is False:
        raise OptionError("'False' is an invalid argument for 'extension'")
    if args.get('extension')==1:args['extension']=True
    from . import sympify
    from ._core import Expr
    wanted=gens or args.get('gens',())
    if isinstance(wanted,Expr):wanted=(wanted,)
    elif len(wanted)==1 and isinstance(wanted[0],(list,tuple)):wanted=wanted[0]
    if wanted==(None,):wanted=()
    wanted=tuple(sympify(g) for g in wanted)
    if len(set(wanted))!=len(wanted):raise GeneratorsError('duplicated generators: %s'%(wanted,))
    if any(getattr(g,'is_commutative',True) is False for g in wanted):
        raise GeneratorsError('non-commutative generators: %s'%(wanted,))
    if 'sort' in args or 'wrt' in args:_sort((),args)
    return wanted,args


def _sort(gens,args):
    sort=args.get('sort',())
    if isinstance(sort,str):sort=[name.strip() for name in sort.split('>')]
    elif hasattr(sort,'__getitem__'):sort=list(map(str,sort))
    else:raise OptionError("invalid argument for 'sort' option")
    wrt=args.get('wrt',())
    from ._core import Expr
    if isinstance(wrt,Expr):wrt=[str(wrt)]
    elif isinstance(wrt,str):
        if wrt.strip().endswith(','):raise OptionError('Bad input: missing parameter.')
        wrt=re.split(r'\s*,\s*|\s+',wrt.strip()) if wrt.strip() else []
    elif hasattr(wrt,'__getitem__'):wrt=list(map(str,wrt))
    else:raise OptionError("invalid argument for 'wrt' option")
    positions={name:i+1 for i,name in enumerate(sort)}
    conventional={name:124+i for i,name in enumerate('xyz')}
    conventional.update({name:216+i for i,name in enumerate('pqrstuvw')})
    conventional.update({name:301+i for i,name in enumerate('abcdefghijklmno')})
    def key(gen):
        text=str(gen)
        if text in wrt:return (-len(wrt)+wrt.index(text),text,0)
        name,index=re.match(r'^(.*?)(\d*)$',text).groups()
        return (positions.get(name,conventional.get(name,1000)),name,int(index or 0))
    return tuple(sorted(gens,key=key))


def _from_terms(terms,gens,domain,modulus,extension):
    from . import Add,Mul,Dummy
    slots=gens if all(getattr(g,'is_Symbol',False) for g in gens) else tuple(Dummy('g') for g in gens)
    expression=Add(*(coefficient*Mul(*(g**e for g,e in zip(slots,monomial)))
                     for monomial,coefficient in terms.items()))
    result=Poly._from_normalized(expression,slots,domain,modulus,extension)
    if slots==gens:return result
    try:return result._with_gens(*gens)
    except ValueError as error:
        if 'non-symbol generators over this coefficient domain' not in str(error):raise
        # Preserve the established direct-expression path for structured
        # coefficient domains whose rename operation is still restricted.
        expression=Add(*(coefficient*Mul(*(g**e for g,e in zip(gens,monomial)))
                         for monomial,coefficient in terms.items()))
        return Poly._from_normalized(expression,gens,domain,modulus,extension)


def _data(rep,gens,domain,modulus,extension):
    from . import sympify
    from ._core import DomainElement,Domain
    def coefficient(value):
        # Domain elements deliberately reject implicit scalar/float coercion.
        # Their parent provides the exact expression embedding.
        return value.parent().to_sympy(value) if isinstance(value,DomainElement) else sympify(value)
    if not gens:
        raise GeneratorsNeeded("Cannot initialize from '%s' without generators"%
                               ('dict' if isinstance(rep,dict) else 'list'))
    if isinstance(rep,dict):
        terms={}
        for monomial,value in rep.items():
            if not isinstance(monomial,tuple):monomial=(monomial,)
            if len(monomial)!=len(gens):raise PolynomialError('invalid monomial for generators')
            try:monomial=tuple(operator.index(e) for e in monomial)
            except TypeError:raise PolynomialError('monomial exponents must be nonnegative integers') from None
            if any(e<0 for e in monomial):raise PolynomialError('monomial exponents must be nonnegative integers')
            terms[monomial]=coefficient(value)
    else:
        if len(gens)!=1:raise MultivariatePolynomialError("'list' representation not supported")
        values=[coefficient(value) for value in rep]
        if (modulus is None and extension is None and
                (domain is None or isinstance(domain,(str,Domain)) and str(domain) in ('ZZ','QQ')) and
                all(value.is_Rational for value in values)):
            return Poly._from_dense_rational(values,gens[0],None if domain is None else str(domain))
        terms={(len(values)-1-i,):value for i,value in enumerate(values)}
    return _from_terms(terms,gens,domain,modulus,extension)


def _unexpanded(rep,gens,domain,modulus,extension,args):
    """Read coefficient/generator factors before the backend distributes them."""
    from . import Add,Mul,S,exp
    def is_coefficient(factor):
        if factor.is_Number:return True
        if domain is not None and not isinstance(domain,str):return factor in domain
        if extension is True:return factor.is_algebraic is True
        if args.get('greedy') is False:return factor.is_number
        return factor==S.ImaginaryUnit
    terms=[];found=set()
    for term in Add.make_args(rep):
        coefficient=[];powers={}
        for factor in Mul.make_args(term):
            if is_coefficient(factor):
                coefficient.append(factor);continue
            base,exponent=factor.as_base_exp()
            if base==S.Exp1 and not exponent.is_Number:
                multiplier,argument=exponent.as_coeff_Mul(rational=True)
                base,exponent=exp(argument),multiplier
            if exponent.is_Rational:
                base,exponent=base**(S.One/exponent.q),exponent.p
            elif not exponent.is_Integer:
                base,exponent=factor,S.One
            if exponent<0:base,exponent=1/base,-exponent
            powers[base]=powers.get(base,0)+int(exponent);found.add(base)
        terms.append((Mul(*coefficient),powers))
    if not gens:gens=_sort(found,args)
    if not gens:raise GeneratorsNeeded('Cannot initialize a polynomial without generators')
    indices={g:i for i,g in enumerate(gens)};data={}
    for coefficient,powers in terms:
        monomial=[0]*len(gens)
        for base,exponent in powers.items():
            if base in indices:monomial[indices[base]]+=exponent
            elif any(base.has(g) for g in gens):
                raise PolynomialError('%s contains an element of the set of generators'%base)
            else:coefficient*=base**exponent
        monomial=tuple(monomial);data[monomial]=data.get(monomial,0)+coefficient
    return _from_terms(data,gens,domain,modulus,extension)


def construct(rep,gens,domain=None,modulus=None,extension=None,kwargs=None):
    gens,args=_options(gens,domain,modulus,extension,kwargs)
    domain=args.get('domain');modulus=args.get('modulus');extension=args.get('extension')
    if args.get('gaussian'):
        from . import QQ_I
        domain=QQ_I
    if args.get('symmetric') is False:
        raise NotImplementedError('nonsymmetric finite-field polynomial coefficients are not implemented')
    if getattr(rep,'is_Matrix',False):
        raise TypeError('a noncommutative matrix cannot be a scalar polynomial')
    from .utilities.iterables import iterable
    if isinstance(rep,dict) or (not isinstance(rep,Poly) and iterable(rep,exclude=str)):
        result=_data(rep,gens,domain,modulus,extension)
    else:
        from . import sympify
        if not isinstance(rep,Poly):rep=sympify(rep,evaluate=type(rep) is not str)
        if getattr(rep,'is_Equality',False):rep=rep.lhs-rep.rhs
        special_inference=not gens and (args.get('greedy') is False or extension is True)
        if (args.get('expand') is False or special_inference) and not isinstance(rep,Poly):
            if args.get('expand',True):
                from . import expand
                rep=expand(rep)
            result=_unexpanded(rep,gens,domain,modulus,extension,args)
        else:result=Poly._from_normalized(rep,gens,domain,modulus,extension)
        if not gens and ('sort' in args or 'wrt' in args):
            result=result.reorder(*_sort(result.gens,args))
    if args.get('field'):result=result.to_field()
    return result
