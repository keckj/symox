"""Result transport and independent checks, always outside timed regions."""
import math


def encode(value, limit=2_000_000):
    # Importing SymPy is intentionally delayed until all native timings finish.
    import sympy as sp
    if value is None or isinstance(value,(bool,str)):
        return {"kind":"literal","value":value}
    if isinstance(value,int):
        return {"kind":"expr","value":sp.srepr(sp.Integer(value))}
    if isinstance(value,float):
        return {"kind":"expr","value":sp.srepr(sp.Float(value,17))}
    if isinstance(value,complex):
        return {"kind":"expr","value":sp.srepr(sp.Float(value.real,17)+sp.I*sp.Float(value.imag,17))}
    if type(value).__module__.startswith("mpmath"):
        if hasattr(value,"_mpf_"):
            return {"kind":"expr","value":sp.srepr(sp.Float(value._mpf_,precision=max(53,value._mpf_[3])))}
        if hasattr(value,"_mpc_"):
            real,imag=value._mpc_
            expression=sp.Float(real,precision=max(53,real[3]))+sp.I*sp.Float(imag,precision=max(53,imag[3]))
            return {"kind":"expr","value":sp.srepr(expression)}
    if hasattr(value,"shape") and hasattr(value,"dtype"):
        return {"kind":"array","shape":list(value.shape),"values":value.tolist()}
    if hasattr(value,"rows") and hasattr(value,"cols") and hasattr(value,"tolist"):
        return {"kind":"matrix","shape":[int(value.rows),int(value.cols)],
                "values":[encode(v,limit) for row in value.tolist() for v in row]}
    if isinstance(value,dict):
        return {"kind":"dict","items":[[encode(k,limit),encode(v,limit)] for k,v in value.items()]}
    if isinstance(value,(list,tuple,set)):
        return {"kind":"sequence","values":[encode(v,limit) for v in value]}
    if hasattr(value,"array_form"):
        return {"kind":"permutation","values":list(value.array_form)}
    if hasattr(value,"polys") and hasattr(value,"order"):
        return {"kind":"basis","values":[encode(v,limit) for v in value.polys],
                "order":str(value.order),"domain":str(value.domain)}
    if hasattr(value,"gens") and hasattr(value,"as_expr"):
        return {"kind":"poly","expr":encode(value.as_expr(),limit),
                "gens":[encode(g,limit) for g in value.gens],"domain":str(value.domain)}
    if type(value).__module__.startswith("symox"):
        converter=getattr(value,"to_sympy",None)
        if converter is None:
            raise TypeError(f"No explicit result exporter for {type(value).__name__}")
        value=converter()
    if not isinstance(value,sp.Basic):
        raise TypeError(f"Unsupported result transport: {type(value).__name__}")
    text=sp.srepr(value)
    if len(text)>limit:
        raise ValueError("Result exceeds serialization size budget")
    return {"kind":"expr","value":text}


def decode(data):
    import sympy as sp
    kind=data["kind"]
    if kind in ("literal","permutation"):return data.get("value",data.get("values"))
    if kind=="expr":return sp.sympify(data["value"])
    if kind=="matrix":return sp.Matrix(*data["shape"],[decode(v) for v in data["values"]])
    if kind=="sequence":return [decode(v) for v in data["values"]]
    if kind=="dict":return {decode(k):decode(v) for k,v in data["items"]}
    if kind=="array":return data
    if kind=="poly":return decode(data["expr"])
    if kind=="basis":return [decode(v) for v in data["values"]]
    raise ValueError(kind)


def held_counts(value):
    import sympy as sp
    from collections import Counter
    counts=Counter()
    if isinstance(value,dict):parts=[*value.keys(),*value.values()]
    elif isinstance(value,(list,tuple)):parts=value
    elif isinstance(value,sp.MatrixBase):parts=list(value)
    else:parts=[value]
    names={"Integral","Derivative","Sum","Product","Limit","ConditionSet",
           "LaplaceTransform","InverseLaplaceTransform","FourierTransform",
           "InverseFourierTransform","Determinant"}
    for part in parts:
        if isinstance(part,(dict,list,tuple,sp.MatrixBase)):
            counts.update(held_counts(part));continue
        if isinstance(part,sp.Basic):
            for node in sp.preorder_traversal(part):
                if type(node).__name__ in names:counts[type(node).__name__]+=1
    return counts


def equivalent(a,b):
    import sympy as sp
    if type(a) is bool or type(b) is bool or a is None or b is None:
        return a is b
    if isinstance(a,dict) or isinstance(b,dict):
        if not isinstance(a,dict) or not isinstance(b,dict) or len(a)!=len(b):return False
        return match_unordered(list(a.items()),list(b.items()))
    if isinstance(a,(list,tuple)) or isinstance(b,(list,tuple)):
        return (isinstance(a,(list,tuple)) and isinstance(b,(list,tuple)) and
                len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b)))
    if isinstance(a,sp.MatrixBase) or isinstance(b,sp.MatrixBase):
        return (isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and
                a.shape==b.shape and all(equivalent(x,y) for x,y in zip(a,b)))
    if isinstance(a,str) or isinstance(b,str):return a==b
    if a==b:return True
    if isinstance(a,sp.Equality) or isinstance(b,sp.Equality):
        return (isinstance(a,sp.Equality) and isinstance(b,sp.Equality) and
                equivalent(a.lhs,b.lhs) and equivalent(a.rhs,b.rhs))
    if isinstance(a,sp.Set) or isinstance(b,sp.Set):
        if not isinstance(a,sp.Set) or not isinstance(b,sp.Set):return False
        # Lambda-bound Dummy identities differ between independent processes.
        # Canonicalize only bound variables; free parameters remain significant.
        if a.as_dummy()==b.as_dummy():return True
        return sp.simplify(sp.SymmetricDifference(a,b))==sp.EmptySet
    difference=a-b
    if difference==0:return True
    # Exact checks only. Numerical substitution is never used to prove identity.
    if sp.cancel(difference)==0:return True
    return sp.simplify(difference)==0


def match_unordered(a,b):
    remaining=list(b)
    for x in a:
        for i,y in enumerate(remaining):
            if equivalent(x,y):
                remaining.pop(i);break
        else:return False
    return not remaining


def numeric_equal(a,b,digits):
    import sympy as sp
    if not (isinstance(a,sp.Basic) and isinstance(b,sp.Basic)):
        a,b=sp.sympify(a),sp.sympify(b)
    if a.free_symbols or b.free_symbols:return False
    tolerance=sp.Rational(1,10**max(1,digits-5))
    for actual,expected in zip(a.as_real_imag(),b.as_real_imag()):
        actual,expected=sp.N(actual,digits+30),sp.N(expected,digits+30)
        if not (actual.is_finite and expected.is_finite):return False
        # A tiny nonzero component must not silently become zero.
        error=abs(actual-expected) if expected==0 else abs((actual-expected)/expected)
        if error>tolerance:return False
    return True


def verify(case,actual,reference):
    """A check on this input/result, never a certificate for the whole API."""
    import sympy as sp
    a,b=decode(actual),decode(reference)
    if case.family=="assumptions" and (a is None or b is None):
        return dict(status="undecided_assumption",detail="An unknown fact is not a proved mathematical answer")
    ah,bh=held_counts(a),held_counts(b)
    if ah or bh:
        return dict(status="unevaluated",detail="Held work is excluded from speed rankings",
                    native_held=dict(ah),reference_held=dict(bh))
    if actual["kind"]!=reference["kind"]:
        return dict(status="result_structure_difference",detail="Return categories differ")
    if actual["kind"] in ("poly","basis"):
        if actual["domain"]!=reference["domain"]:
            return dict(status="domain_difference",detail=f"{actual['domain']} versus {reference['domain']}")
        if actual["kind"]=="basis" and actual["order"]!=reference["order"]:
            return dict(status="order_difference")
    def expressions(value):
        if isinstance(value,dict):return [e for pair in value.items() for item in pair for e in expressions(item)]
        if isinstance(value,(list,tuple,sp.MatrixBase)):return [e for item in value for e in expressions(item)]
        return [value] if isinstance(value,sp.Basic) else []
    if case.validation not in ("numeric","numeric_roots","array"):
        if any(e.has(sp.Float) for e in expressions(a)) and not any(e.has(sp.Float) for e in expressions(b)):
            return dict(status="exactness_difference",detail="Inexact output replaced an exact reference result")
    # Algebraic equality alone would wrongly reward a factor/cancel/simplify
    # implementation that simply returned its input without doing the work.
    if case.operation.startswith("factor_") and case.operation!="factor_modular":
        def factor_shape(expression):
            shape=[]
            for factor in sp.Mul.make_args(expression):
                base,power=factor.as_base_exp()
                if not (power.is_Integer and power>0):base,power=factor,sp.S.One
                symbols=sorted(base.free_symbols,key=str)
                if symbols:shape.append((sp.Poly(base,*symbols).total_degree(),int(power)))
            return sorted(shape)
        if factor_shape(a)!=factor_shape(b):
            return dict(status="factorization_quality_difference",detail="Factor degree/multiplicity structure differs")
    if case.family=="simplification":
        counts=dict(native=int(sp.count_ops(a)),reference=int(sp.count_ops(b)))
        if counts["native"]>counts["reference"]:
            return dict(status="simplification_quality_difference",operation_counts=counts,
                        detail="Equal expressions with greater output complexity do not qualify as a speed win")
    if case.operation=="cancel_sum":
        numerator,denominator=sp.fraction(a)
        symbols=sorted(a.free_symbols,key=str)
        try:
            pn,pd=sp.Poly(numerator,*symbols),sp.Poly(denominator,*symbols)
            reduced=sp.gcd(pn,pd).total_degree()==0
        except sp.PolynomialError:reduced=False
        if not reduced:return dict(status="cancellation_quality_difference")
    if case.operation=="apart":
        x=sp.Symbol("x")
        for term in sp.Add.make_args(a):
            numerator,denominator=sp.fraction(term)
            if denominator==1:continue
            factors=sp.factor_list(denominator,x)[1]
            if len(factors)!=1 or sp.degree(factors[0][0],x)!=1 or sp.degree(numerator,x)>0:
                return dict(status="partial_fraction_quality_difference")
    if case.family=="series":
        x=sp.Symbol("x");numerator,denominator=sp.fraction(sp.together(a.removeO()))
        try:
            finite_laurent=sp.Poly(numerator,x) is not None and len(sp.Poly(denominator,x).terms())==1
        except sp.PolynomialError:finite_laurent=False
        if a.getO()!=b.getO() or not finite_laurent:
            return dict(status="series_quality_difference")
        a,b=a.removeO(),b.removeO()
    if case.operation=="factor_modular":
        x=sp.Symbol("x")
        coefficient,factors=a
        if any(not sp.Poly(f,x,modulus=5).is_irreducible or k<=0 for f,k in factors):
            return dict(status="factorization_quality_difference")
        reconstructed=coefficient*sp.prod(f**k for f,k in factors)
        good=sp.Poly(reconstructed-(x**case.size-1),x,modulus=5).is_zero
        return dict(status="verified" if good else "observation_difference",scope="Exact modular identity and irreducibility of returned factors")
    if case.validation=="array":
        import numpy as np
        good=(a["shape"]==b["shape"] and
              np.allclose(a["values"],b["values"],rtol=1e-12,atol=1e-14,equal_nan=False))
        return dict(status="verified" if good else "numerical_difference",scope="Float64 arrays: rtol=1e-12, atol=1e-14")
    if case.validation=="text":
        return dict(status="verified" if a==b else "output_difference",
                    scope="Literal output equality; alternate formatting is not declared mathematically wrong")
    if case.validation=="numeric_roots":
        if len(a)!=len(b) or len(a)!=case.size:return dict(status="root_count_difference")
        remaining=list(b)
        for value in a:
            for i,target in enumerate(remaining):
                if numeric_equal(value,target,case.digits):remaining.pop(i);break
            else:return dict(status="numerical_difference",scope="Complete root count and componentwise relative precision")
        return dict(status="verified",scope="Complete root count and componentwise relative precision")
    if case.validation=="numeric":
        digits=case.digits or case.size
        # Evaluate a separate high-precision reference input, rather than
        # accepting agreement with a low-precision or underflowed answer.
        if case.family=="numerics":
            if case.operation=="constant":target=sp.exp(sp.pi)*sp.sqrt(2)+sp.log(3)
            elif case.operation=="tiny_sine":target=sp.sin(sp.Rational(1,10**1000))
            elif case.operation=="large_argument":target=sp.sin(sp.Integer(10)**100)
            else:target=sp.log(-2+sp.I)
            b=target.evalf(digits+40)
        elif case.operation=="nsolve":b=sp.nsolve(sp.cos(sp.Symbol('x'))-sp.Symbol('x'),1,prec=digits+40)
        good=numeric_equal(a,b,digits)
        return dict(status="verified" if good else "numerical_difference",
                    scope=f"Componentwise relative precision within 5 guard digits of {digits} requested digits")
    if case.validation=="antiderivative":
        x=sp.Symbol("x")
        # Equality of derivatives on their common analytic domain, not a proof
        # about arbitrary branches or distributional endpoint conventions.
        good=equivalent(sp.diff(a-b,x),sp.S.Zero)
        return dict(status="verified" if good else "observation_difference",scope="Equal derivative on common analytic domain")
    if case.validation=="ode":
        x=sp.Symbol("x");order=1 if case.operation=="linear_first" else 2
        def complete_solution(value):
            if not isinstance(value,sp.Equality) or value.lhs!=sp.Function("f")(x):return False
            rhs=value.rhs;constants=sorted(rhs.free_symbols-{x},key=str)
            if len(constants)!=order or rhs.subs(dict.fromkeys(constants,0))!=0:return False
            if any(sp.diff(rhs,c,d)!=0 for c in constants for d in constants):return False
            residual=sp.diff(rhs,x)-rhs if order==1 else sp.diff(rhs,x,2)+rhs
            if not equivalent(residual,sp.S.Zero):return False
            basis=[sp.diff(rhs,c) for c in constants]
            wronskian=sp.Matrix([[sp.diff(v,x,k).subs(x,0) for v in basis] for k in range(order)])
            return wronskian.det().is_zero is False
        good=complete_solution(a) and complete_solution(b)
        return dict(status="verified" if good else "observation_difference",scope="Complete independent basis for the specified homogeneous linear ODE")
    if case.validation=="unordered":good=match_unordered(a,b)
    elif case.validation=="nullspace":
        from performance_cases import integer_rows
        matrix=sp.Matrix(integer_rows(case.size,case.size+6,45))
        rank=matrix.rank();dimension=matrix.cols-rank
        good=(len(a)==len(b)==dimension and all(matrix*v==sp.zeros(matrix.rows,1) for v in a)
              and sp.Matrix.hstack(*a).rank()==dimension)
    else:good=equivalent(a,b)
    return dict(status="verified" if good else "observation_difference",
                scope="Exact comparison on the stated workload; broader behavior remains unaudited")
