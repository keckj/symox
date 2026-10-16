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
