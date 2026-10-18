"""Rational transfer functions, interconnections and state-space systems."""
import symox as S


def _tf(system):
    value=system.doit() if hasattr(system,'doit') else system
    if not isinstance(value,TransferFunction):raise TypeError('a SISO transfer function is required')
    return value


class TransferFunction:
    def __init__(self,num,den,var):
        self.num,self.den,self.var=map(S.sympify,(num,den,var))
        if not self.var.is_Symbol:raise TypeError('transfer-function variable must be a Symbol')
        if self.den==0:raise ValueError('denominator cannot be zero')
    @classmethod
    def from_rational_expression(cls,expr,var=None):
        expr=S.sympify(expr)
        if var is None:
            if len(expr.free_symbols)!=1:raise ValueError('specify the transfer-function variable')
            var=next(iter(expr.free_symbols))
        num,den=expr.as_numer_denom();return cls(num,den,var)
    @classmethod
    def from_coeff_lists(cls,num_list,den_list,var):
        var=S.sympify(var)
        polynomial=lambda values:S.Add(*(S.sympify(c)*var**(len(values)-i-1) for i,c in enumerate(values)))
        return cls(polynomial(num_list),polynomial(den_list),var)
    @property
    def args(self):return self.num,self.den,self.var
    def __repr__(self):return 'TransferFunction(%s, %s, %s)'%self.args
    def __eq__(self,other):return isinstance(other,TransferFunction) and self.args==other.args
    def to_expr(self):return self.num/self.den
    def doit(self,**kwargs):return self
    def simplify(self):return self.from_rational_expression(S.cancel(self.to_expr()),self.var)
    def expand(self):return TransferFunction(S.expand(self.num),S.expand(self.den),self.var)
    def evalf(self,*args,**kwargs):return TransferFunction(self.num.evalf(*args,**kwargs),self.den.evalf(*args,**kwargs),self.var)
    def subs(self,*args,**kwargs):return TransferFunction(self.num.subs(*args,**kwargs),self.den.subs(*args,**kwargs),self.var)
    def _compatible(self,other):
        other=_tf(other)
        if self.var!=other.var:raise ValueError('transfer-function variables must match')
        return other
    def __mul__(self,other):return Series(self,self._compatible(other))
    def __add__(self,other):return Parallel(self,self._compatible(other))
    def __neg__(self):return TransferFunction(-self.num,self.den,self.var)
    def __sub__(self,other):return self+-other
    def __truediv__(self,other):
        other=self._compatible(other);return Series(self,TransferFunction(other.den,other.num,self.var))
    def __pow__(self,power):
        power=S.sympify(power)
        if not power.is_Integer:raise ValueError('transfer-function power must be integral')
        return TransferFunction(self.num**power,self.den**power,self.var) if power>=0 else TransferFunction(self.den**(-power),self.num**(-power),self.var)
    def eval_frequency(self,other):return S.simplify(self.to_expr().subs(self.var,other))
    def dc_gain(self):return S.limit(self.to_expr(),self.var,0)
    def _roots(self,expression):
        roots=S.roots(expression,self.var)
        if isinstance(roots,dict) and sum(int(n) for n in roots.values())==int(S.degree(expression,self.var)):
            return [root for root,count in roots.items() for _ in range(int(count))]
        return S.Poly(expression,self.var).all_roots()
    def poles(self):return self._roots(self.den)
    def zeros(self):return [] if self.num==0 else self._roots(self.num)
    def is_stable(self):
        signs=[S.re(p).is_negative for p in self.poles()]
        return False if False in signs else True if all(s is True for s in signs) else None
    @property
    def is_proper(self):return S.degree(self.num,self.var)<=S.degree(self.den,self.var)
    @property
    def is_strictly_proper(self):return S.degree(self.num,self.var)<S.degree(self.den,self.var)
    @property
    def is_biproper(self):return S.degree(self.num,self.var)==S.degree(self.den,self.var)
    def rewrite(self,target,*args,**kwargs):
        if target is StateSpace:
            denominator=S.Poly(self.den,self.var).all_coeffs();numerator=S.Poly(self.num,self.var).all_coeffs();order=len(denominator)-1
            if len(numerator)>len(denominator):raise ValueError('state-space conversion requires a proper transfer function')
            lead=denominator[0];a=[v/lead for v in denominator];b=[S.Integer(0)]*(order+1-len(numerator))+[v/lead for v in numerator]
            A=S.zeros(order);B=S.zeros(order,1)
            if order:
                for j in range(order):A[0,j]=-a[j+1]
                for i in range(1,order):A[i,i-1]=1
                B[0,0]=1
            C=S.Matrix([[b[i+1]-b[0]*a[i+1] for i in range(order)]]) if order else S.zeros(1,0)
            return StateSpace(A,B,C,S.Matrix([[b[0]]]))
        return self


class PIDController(TransferFunction):
    def __init__(self,kp=S.Symbol('kp'),ki=S.Symbol('ki'),kd=S.Symbol('kd'),tf=0,var=S.Symbol('s')):
        self.kp,self.ki,self.kd,self.tf=map(S.sympify,(kp,ki,kd,tf));var=S.sympify(var)
        expression=self.kp+self.ki/var+self.kd*var/(1+self.tf*var)
        num,den=S.together(expression).as_numer_denom();super().__init__(num,den,var)


class _Connection(TransferFunction):
    def __init__(self,*args,evaluate=False):
        if not args:raise ValueError('at least one transfer function is required')
        self.systems=tuple(_tf(v) for v in args)
        if len({v.var for v in self.systems})!=1:raise ValueError('transfer-function variables must match')
        self.var=self.systems[0].var
    @property
    def args(self):return self.systems
    @property
    def num(self):return self.doit().num
    @property
    def den(self):return self.doit().den
    def __repr__(self):return type(self).__name__+'('+', '.join(map(str,self.systems))+')'
    def to_expr(self):return self.doit().to_expr()


class Series(_Connection):
    def doit(self,**kwargs):return TransferFunction(S.Mul(*(v.num for v in self.systems)),S.Mul(*(v.den for v in self.systems)),self.var)
class Parallel(_Connection):
    def doit(self,**kwargs):
        denominator=S.Mul(*(v.den for v in self.systems));numerator=S.Add(*(v.num*S.Mul(*(w.den for j,w in enumerate(self.systems) if j!=i)) for i,v in enumerate(self.systems)))
        return TransferFunction(numerator,denominator,self.var)


class Feedback(_Connection):
    def __init__(self,sys1,sys2=None,sign=-1):
        sys1=_tf(sys1);sys2=TransferFunction(1,1,sys1.var) if sys2 is None else _tf(sys2)
        if sign not in (-1,1):raise ValueError('feedback sign must be -1 or 1')
        super().__init__(sys1,sys2);self.sys1,self.sys2,self.sign=sys1,sys2,sign
    def doit(self,**kwargs):
        a,b=self.systems
        return TransferFunction(a.num*b.den,a.den*b.den-self.sign*a.num*b.num,self.var)
    @property
    def sensitivity(self):return 1/(1-self.sign*self.sys1.to_expr()*self.sys2.to_expr())


class TransferFunctionMatrix:
    def __init__(self,arg):
        self._entries=tuple(tuple(_tf(v) for v in row) for row in arg)
        if not self._entries or not self._entries[0] or any(len(r)!=len(self._entries[0]) for r in self._entries):raise ValueError('transfer-function matrix must be nonempty and rectangular')
        variables={v.var for row in self._entries for v in row}
        if len(variables)!=1:raise ValueError('matrix transfer-function variables must match')
        self.var=next(iter(variables));self.num_outputs,self.num_inputs=len(self._entries),len(self._entries[0]);self.shape=self.num_outputs,self.num_inputs
    @classmethod
    def from_Matrix(cls,matrix,var):return cls([[TransferFunction.from_rational_expression(matrix[i,j],var) for j in range(matrix.cols)] for i in range(matrix.rows)])
    def __getitem__(self,key):
        if isinstance(key,tuple):return self._entries[key[0]][key[1]]
        return self._entries[key]
    def to_expr(self):return S.Matrix([[v.to_expr() for v in row] for row in self._entries])
    def doit(self,**kwargs):return self
    def __neg__(self):return TransferFunctionMatrix([[-v for v in row] for row in self._entries])
    def __add__(self,other):return MIMOParallel(self,other)
    def __mul__(self,other):return MIMOSeries(other,self)
    def poles(self):return [[v.poles() for v in row] for row in self._entries]
    def zeros(self):return [[v.zeros() for v in row] for row in self._entries]
    def subs(self,*args,**kwargs):return TransferFunctionMatrix([[v.subs(*args,**kwargs) for v in row] for row in self._entries])
    def __repr__(self):return 'TransferFunctionMatrix(%s)'%(self._entries,)


class MIMOSeries:
    def __init__(self,*args,evaluate=False):
        self.args=tuple(v.doit() for v in args)
        if not self.args or any(not isinstance(v,TransferFunctionMatrix) for v in self.args):raise TypeError('MIMO series requires transfer-function matrices')
        self.var=self.args[0].var
    def doit(self,**kwargs):
        result=self.args[0].to_expr()
        for matrix in self.args[1:]:
            if matrix.var!=self.var:raise ValueError('transfer-function variables differ')
            result=matrix.to_expr()*result
        return TransferFunctionMatrix.from_Matrix(result,self.var)
    def to_expr(self):return self.doit().to_expr()


class MIMOParallel(MIMOSeries):
    def doit(self,**kwargs):
        result=self.args[0].to_expr()
        for matrix in self.args[1:]:result+=matrix.to_expr()
        return TransferFunctionMatrix.from_Matrix(result,self.var)


class MIMOFeedback:
    def __init__(self,sys1,sys2=None,sign=-1):
        self.sys1=sys1.doit();self.var=self.sys1.var
        self.sys2=sys2.doit() if sys2 is not None else TransferFunctionMatrix.from_Matrix(S.eye(self.sys1.num_inputs),self.var)
        if sign not in (-1,1):raise ValueError('feedback sign must be -1 or 1')
        self.sign=sign
    def doit(self,**kwargs):
        G,H=self.sys1.to_expr(),self.sys2.to_expr()
        return TransferFunctionMatrix.from_Matrix((S.eye(G.rows)-self.sign*G*H).inv()*G,self.var)
    def to_expr(self):return self.doit().to_expr()


class StateSpace:
    def __init__(self,A=None,B=None,C=None,D=None):
        self.state_matrix=S.Matrix(A) if A is not None else S.zeros(0)
        n=self.state_matrix.rows
        if self.state_matrix.cols!=n:raise ValueError('state matrix must be square')
        self.input_matrix=S.Matrix(B) if B is not None else S.zeros(n,0)
        self.output_matrix=S.Matrix(C) if C is not None else S.zeros(0,n)
        self.feedforward_matrix=S.Matrix(D) if D is not None else S.zeros(self.output_matrix.rows,self.input_matrix.cols)
        if self.input_matrix.rows!=n or self.output_matrix.cols!=n or self.feedforward_matrix.shape!=(self.output_matrix.rows,self.input_matrix.cols):raise ValueError('state-space matrix dimensions are inconsistent')
        self.num_states,self.num_inputs,self.num_outputs=n,self.input_matrix.cols,self.output_matrix.rows
    @property
    def A(self):return self.state_matrix
    @property
    def B(self):return self.input_matrix
    @property
    def C(self):return self.output_matrix
    @property
    def D(self):return self.feedforward_matrix
    def controllability_matrix(self):
        result=self.B
        for power in range(1,self.num_states):result=result.row_join(self.A**power*self.B)
        return result
    def observability_matrix(self):
        result=self.C
        for power in range(1,self.num_states):result=result.col_join(self.C*self.A**power)
        return result
    def is_controllable(self):return self.controllability_matrix().rank()==self.num_states
    def is_observable(self):return self.observability_matrix().rank()==self.num_states
    def rewrite(self,target,var=None):
        if target is TransferFunction:
            var=S.Symbol('s') if var is None else S.sympify(var)
            matrix=self.C*(var*S.eye(self.num_states)-self.A).inv()*self.B+self.D
            return [[TransferFunction.from_rational_expression(matrix[i,j],var) for j in range(matrix.cols)] for i in range(matrix.rows)]
        return self
    def to_transfer_function(self,var=None):return TransferFunctionMatrix(self.rewrite(TransferFunction,var))
    def append(self,other):return StateSpace(S.diag(self.A,other.A),S.diag(self.B,other.B),S.diag(self.C,other.C),S.diag(self.D,other.D))
    def __add__(self,other):
        if (self.num_inputs,self.num_outputs)!=(other.num_inputs,other.num_outputs):raise ValueError('parallel systems must have equal input/output dimensions')
        return StateSpace(S.diag(self.A,other.A),self.B.col_join(other.B),self.C.row_join(other.C),self.D+other.D)
    def __mul__(self,other):
        if self.num_inputs!=other.num_outputs:raise ValueError('series input/output dimensions differ')
        A=self.A.row_join(self.B*other.C).col_join(S.zeros(other.num_states,self.num_states).row_join(other.A))
        return StateSpace(A,(self.B*other.D).col_join(other.B),self.C.row_join(self.D*other.C),self.D*other.D)
    def dsolve(self,initial_conditions=None,input_vector=None,var=S.Symbol('t')):
        var=S.sympify(var);x0=S.zeros(self.num_states,1) if initial_conditions is None else S.Matrix(initial_conditions)
        value=(self.A*var).exp()*x0
        if input_vector is not None:
            tau=S.Dummy('tau',real=True);integrand=(self.A*(var-tau)).exp()*self.B*S.Matrix(input_vector).subs(var,tau)
            value+=integrand.applyfunc(lambda entry:S.integrate(entry,(tau,0,var)))
        return value


def gbt(tf,sample_per,alpha):
    tf=_tf(tf);T,alpha=map(S.sympify,(sample_per,alpha));z=S.Symbol('_z')
    if alpha.is_number and (alpha<0 or alpha>1):raise ValueError('alpha must lie in [0,1]')
    expression=S.cancel(tf.to_expr().subs(tf.var,(z-1)/(T*(alpha*z+1-alpha))))
    num,den=expression.as_numer_denom();a=S.Poly(den,z).all_coeffs();b=S.Poly(num,z).all_coeffs();lead=a[0]
    return [S.simplify(v/lead) for v in b],[S.simplify(v/lead) for v in a]

def bilinear(tf,sample_per):return gbt(tf,sample_per,S.Rational(1,2))
def forward_diff(tf,sample_per):return gbt(tf,sample_per,0)
def backward_diff(tf,sample_per):return gbt(tf,sample_per,1)


__all__=['TransferFunction','PIDController','Series','MIMOSeries','Parallel','MIMOParallel','Feedback','MIMOFeedback','TransferFunctionMatrix','StateSpace','gbt','bilinear','forward_diff','backward_diff']

# Imported lazily to avoid a cycle while the LTI classes are constructed.
def phase_margin(system):
    from .margins import phase_margin as margin
    return margin(system)
def gain_margin(system):
    from .margins import gain_margin as margin
    return margin(system)
