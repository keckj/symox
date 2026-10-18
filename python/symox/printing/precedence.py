"""Expression precedence compatible with SymPy's printer interface."""

PRECEDENCE = {
    "Lambda": 1, "Xor": 10, "Or": 20, "And": 30, "Relational": 35,
    "BitwiseOr": 36, "BitwiseXor": 37, "BitwiseAnd": 38,
    "Add": 40, "Mul": 50, "Pow": 60, "Func": 70, "Not": 100,
    "Atom": 1000,
}


def precedence(item):
    value = getattr(item, "precedence", None)
    if value is not None:
        return value
    name = getattr(item, "class_name", type(item).__name__)
    name = {"Half": "Rational", "Zero": "Integer", "One": "Integer", "NegativeOne": "Integer"}.get(name, name)
    if name in globals().get('PRECEDENCE_FUNCTIONS',{}):return PRECEDENCE_FUNCTIONS[name](item)
    if name in globals().get('PRECEDENCE_VALUES',{}):return PRECEDENCE_VALUES[name]
    if name == "UnevaluatedExpr":
        return precedence(item.args[0]) - 0.5
    if name in ("Integer", "Float", "Rational", "Mul"):
        if getattr(item, "is_negative", None) is True:
            return PRECEDENCE["Add"]
        if name == "Mul":
            args = item.args
            if args and getattr(args[0], "is_negative", None) is True:
                return PRECEDENCE["Add"]
            return PRECEDENCE["Mul"]
        if name == "Rational":
            return PRECEDENCE["Mul"]
        return PRECEDENCE["Atom"]
    if getattr(item, "is_Relational", False):
        return PRECEDENCE["Relational"]
    return PRECEDENCE.get(name, PRECEDENCE["Func"] if getattr(item, "is_Function", False) else PRECEDENCE["Atom"])


def precedence_Integer(item):return PRECEDENCE['Add'] if item<0 else PRECEDENCE['Atom']
def precedence_Float(item):return PRECEDENCE['Add'] if item<0 else PRECEDENCE['Atom']
def precedence_Rational(item):return PRECEDENCE['Add'] if item.p<0 else PRECEDENCE['Mul']
def precedence_Mul(item):return PRECEDENCE['Add'] if item.could_extract_minus_sign() else PRECEDENCE['Mul']
def precedence_UnevaluatedExpr(item):return precedence(item.args[0])-.5
def precedence_PolyElement(item):
    if item.is_generator:return PRECEDENCE['Atom']
    if item.is_ground:return precedence(item.as_expr())
    return PRECEDENCE['Mul'] if item.is_term else PRECEDENCE['Add']
def precedence_FracElement(item):return precedence_PolyElement(item.numer) if item.denom==1 else PRECEDENCE['Mul']
PRECEDENCE_FUNCTIONS={name:globals()['precedence_'+name] for name in 'Integer Float Rational Mul UnevaluatedExpr PolyElement FracElement'.split()}
PRECEDENCE_VALUES={**PRECEDENCE,'Equivalent':10,'Implies':10,'NegativeInfinity':40,'MatAdd':40,'MatPow':60,
 'MatMul':50,'MatrixSolve':50,'TensAdd':40,'TensMul':50,'HadamardProduct':50,'KroneckerProduct':50,'HadamardPower':60,'Mod':50}
PRECEDENCE_TRADITIONAL={**PRECEDENCE_VALUES,**dict.fromkeys('Integral Sum Product Limit Derivative TensorProduct DotProduct'.split(),PRECEDENCE['Mul'])}
PRECEDENCE_TRADITIONAL.update(dict.fromkeys('Union Intersection Complement SymmetricDifference ProductSet'.split(),PRECEDENCE['Xor']))

def precedence_traditional(item):return PRECEDENCE_TRADITIONAL.get(getattr(item,'class_name',type(item).__name__),precedence(item))
