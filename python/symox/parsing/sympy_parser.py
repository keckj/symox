"""Python-style expression parsing on native Symox objects.

Token transformations adapt SymPy 1.14 algorithms (see LICENSE-SYMPY).
Like Python eval, parse_expr accepts trusted input, not untrusted expressions.
"""
from __future__ import annotations
from tokenize import generate_tokens, untokenize, TokenError, NUMBER, STRING, NAME, OP, ENDMARKER, ERRORTOKEN, NEWLINE
from keyword import iskeyword
import ast
import unicodedata
from io import StringIO
import builtins
import types
from typing import Any, Callable
from functools import reduce
from functools import lru_cache
from symox.assumptions.ask import AssumptionKeys
from symox import Basic
from symox import Symbol
from symox import Function

def func_name(value):
    return value.__name__
from symox import Max, Min
null = ''
TOKEN = tuple[int, str]
DICT = dict[str, Any]
TRANS = Callable[[list[TOKEN], DICT, DICT], list[TOKEN]]

def _token_splittable(token_name: str) -> bool:
    """Predicate for whether a token name can be split into multiple tokens."""
    if '_' in token_name:
        return False
    try:
        return not unicodedata.lookup('GREEK SMALL LETTER ' + token_name)
    except KeyError:
        return len(token_name) > 1

def _token_callable(token: TOKEN, local_dict: DICT, global_dict: DICT, nextToken=None):
    """Predicate for whether a token name represents a callable function."""
    func = local_dict.get(token[1])
    if not func:
        func = global_dict.get(token[1])
    return callable(func) and (not isinstance(func, Basic) or getattr(func, "is_Lambda", False))

def _add_factorial_tokens(name: str, result: list[TOKEN]) -> list[TOKEN]:
    if result == [] or result[-1][1] == '(':
        raise TokenError()
    beginning = [(NAME, name), (OP, '(')]
    end = [(OP, ')')]
    diff = 0
    length = len(result)
    for (index, token) in enumerate(result[::-1]):
        (toknum, tokval) = token
        i = length - index - 1
        if tokval == ')':
            diff += 1
        elif tokval == '(':
            diff -= 1
        if diff == 0:
            if i - 1 >= 0 and result[i - 1][0] == NAME:
                return result[:i - 1] + beginning + result[i - 1:] + end
            else:
                return result[:i] + beginning + result[i:] + end
    return result

class ParenthesisGroup(list[TOKEN]):
    """List of tokens representing an expression in parentheses."""
    pass

class AppliedFunction:
    """A group of tokens representing a function and its arguments."""

    def __init__(self, function: TOKEN, args: ParenthesisGroup, exponent=None):
        if exponent is None:
            exponent = []
        self.function = function
        self.args = args
        self.exponent = exponent
        self.items = ['function', 'args', 'exponent']

    def expand(self) -> list[TOKEN]:
        """Return a list of tokens representing the function"""
        return [self.function, *self.args]

    def __getitem__(self, index):
        return getattr(self, self.items[index])

    def __repr__(self):
        return 'AppliedFunction(%s, %s, %s)' % (self.function, self.args, self.exponent)

def _flatten(result: list[TOKEN | AppliedFunction]):
    result2: list[TOKEN] = []
    for tok in result:
        if isinstance(tok, AppliedFunction):
            result2.extend(tok.expand())
        else:
            result2.append(tok)
    return result2

def _group_parentheses(recursor: TRANS):

    def _inner(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
        """Group tokens between parentheses with ParenthesisGroup."""
        result: list[TOKEN | ParenthesisGroup] = []
        stacks: list[ParenthesisGroup] = []
        stacklevel = 0
        for token in tokens:
            if token[0] == OP:
                if token[1] == '(':
                    stacks.append(ParenthesisGroup([]))
                    stacklevel += 1
                elif token[1] == ')':
                    stacks[-1].append(token)
                    stack = stacks.pop()
                    if len(stacks) > 0:
                        stacks[-1].extend(stack)
                    else:
                        inner = stack[1:-1]
                        inner = recursor(inner, local_dict, global_dict)
                        parenGroup = [stack[0]] + inner + [stack[-1]]
                        result.append(ParenthesisGroup(parenGroup))
                    stacklevel -= 1
                    continue
            if stacklevel:
                stacks[-1].append(token)
            else:
                result.append(token)
        if stacklevel:
            raise TokenError('Mismatched parentheses')
        return result
    return _inner

def _apply_functions(tokens: list[TOKEN | ParenthesisGroup], local_dict: DICT, global_dict: DICT):
    """Convert a NAME token + ParenthesisGroup into an AppliedFunction."""
    result: list[TOKEN | AppliedFunction] = []
    symbol = None
    for tok in tokens:
        if isinstance(tok, ParenthesisGroup):
            if symbol and _token_callable(symbol, local_dict, global_dict):
                result[-1] = AppliedFunction(symbol, tok)
                symbol = None
            else:
                result.extend(tok)
        elif tok[0] == NAME:
            symbol = tok
            result.append(tok)
        else:
            symbol = None
            result.append(tok)
    return result

def _implicit_multiplication(tokens: list[TOKEN | AppliedFunction], local_dict: DICT, global_dict: DICT):
    """Implicitly adds '*' tokens."""
    result: list[TOKEN | AppliedFunction] = []
    skip = False
    for (tok, nextTok) in zip(tokens, tokens[1:]):
        result.append(tok)
        if skip:
            skip = False
            continue
        if tok[0] == OP and tok[1] == '.' and (nextTok[0] == NAME):
            skip = True
            continue
        if isinstance(tok, AppliedFunction):
            if isinstance(nextTok, AppliedFunction):
                result.append((OP, '*'))
            elif nextTok == (OP, '('):
                if tok.function[1] == 'Function':
                    tok.function = (tok.function[0], 'Symbol')
                result.append((OP, '*'))
            elif nextTok[0] == NAME:
                result.append((OP, '*'))
        elif tok == (OP, ')'):
            if isinstance(nextTok, AppliedFunction):
                result.append((OP, '*'))
            elif nextTok[0] == NAME:
                result.append((OP, '*'))
            elif nextTok == (OP, '('):
                result.append((OP, '*'))
        elif tok[0] == NAME and (not _token_callable(tok, local_dict, global_dict)):
            if isinstance(nextTok, AppliedFunction) or (nextTok[0] == NAME and _token_callable(nextTok, local_dict, global_dict)):
                result.append((OP, '*'))
            elif nextTok == (OP, '('):
                result.append((OP, '*'))
            elif nextTok[0] == NAME:
                result.append((OP, '*'))
    if tokens:
        result.append(tokens[-1])
    return result

def _implicit_application(tokens: list[TOKEN | AppliedFunction], local_dict: DICT, global_dict: DICT):
    """Adds parentheses as needed after functions."""
    result: list[TOKEN | AppliedFunction] = []
    appendParen = 0
    skip = 0
    exponentSkip = False
    for (tok, nextTok) in zip(tokens, tokens[1:]):
        result.append(tok)
        if tok[0] == NAME and nextTok[0] not in [OP, ENDMARKER, NEWLINE]:
            if _token_callable(tok, local_dict, global_dict, nextTok):
                result.append((OP, '('))
                appendParen += 1
        elif tok[0] == NAME and nextTok[0] == OP and (nextTok[1] == '**'):
            if _token_callable(tok, local_dict, global_dict):
                exponentSkip = True
        elif exponentSkip:
            if isinstance(tok, AppliedFunction) or (tok[0] == OP and tok[1] == '*'):
                if not (nextTok[0] == OP and nextTok[1] == '*'):
                    if not (nextTok[0] == OP and nextTok[1] == '('):
                        result.append((OP, '('))
                        appendParen += 1
                    exponentSkip = False
        elif appendParen:
            if nextTok[0] == OP and nextTok[1] in ('^', '**', '*'):
                skip = 1
                continue
            if skip:
                skip -= 1
                continue
            result.append((OP, ')'))
            appendParen -= 1
    if tokens:
        result.append(tokens[-1])
    if appendParen:
        result.extend([(OP, ')')] * appendParen)
    return result

def function_exponentiation(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Allows functions to be exponentiated, e.g. ``cos**2(x)``."""
    result: list[TOKEN] = []
    exponent: list[TOKEN] = []
    consuming_exponent = False
    level = 0
    for (tok, nextTok) in zip(tokens, tokens[1:]):
        if tok[0] == NAME and nextTok[0] == OP and (nextTok[1] == '**'):
            if _token_callable(tok, local_dict, global_dict):
                consuming_exponent = True
        elif consuming_exponent:
            if tok[0] == NAME and tok[1] == 'Function':
                tok = (NAME, 'Symbol')
            exponent.append(tok)
            if tok[0] == nextTok[0] == OP and tok[1] == ')' and (nextTok[1] == '('):
                consuming_exponent = False
            if tok[0] == nextTok[0] == OP and tok[1] == '*' and (nextTok[1] == '('):
                consuming_exponent = False
                del exponent[-1]
            continue
        elif exponent and (not consuming_exponent):
            if tok[0] == OP:
                if tok[1] == '(':
                    level += 1
                elif tok[1] == ')':
                    level -= 1
            if level == 0:
                result.append(tok)
                result.extend(exponent)
                exponent = []
                continue
        result.append(tok)
    if tokens:
        result.append(tokens[-1])
    if exponent:
        result.extend(exponent)
    return result

def split_symbols_custom(predicate: Callable[[str], bool]):
    """Creates a transformation that splits symbol names."""

    def _split_symbols(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
        result: list[TOKEN] = []
        split = False
        split_previous = False
        for tok in tokens:
            if split_previous:
                split_previous = False
                continue
            split_previous = False
            if tok[0] == NAME and tok[1] in ['Symbol', 'Function']:
                split = True
            elif split and tok[0] == NAME:
                symbol = tok[1][1:-1]
                if predicate(symbol):
                    tok_type = result[-2][1]
                    del result[-2:]
                    i = 0
                    while i < len(symbol):
                        char = symbol[i]
                        if char in local_dict or char in global_dict:
                            result.append((NAME, '%s' % char))
                        elif char.isdigit():
                            chars = [char]
                            for i in range(i + 1, len(symbol)):
                                if not symbol[i].isdigit():
                                    i -= 1
                                    break
                                chars.append(symbol[i])
                            char = ''.join(chars)
                            result.extend([(NAME, 'Number'), (OP, '('), (NAME, "'%s'" % char), (OP, ')')])
                        else:
                            use = tok_type if i == len(symbol) else 'Symbol'
                            result.extend([(NAME, use), (OP, '('), (NAME, "'%s'" % char), (OP, ')')])
                        i += 1
                    split = False
                    split_previous = True
                    continue
                else:
                    split = False
            result.append(tok)
        return result
    return _split_symbols
split_symbols = split_symbols_custom(_token_splittable)

def implicit_multiplication(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT) -> list[TOKEN]:
    """Makes the multiplication operator optional in most cases."""
    res1 = _group_parentheses(implicit_multiplication)(tokens, local_dict, global_dict)
    res2 = _apply_functions(res1, local_dict, global_dict)
    res3 = _implicit_multiplication(res2, local_dict, global_dict)
    result = _flatten(res3)
    return result

def implicit_application(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT) -> list[TOKEN]:
    """Makes parentheses optional in some cases for function calls."""
    res1 = _group_parentheses(implicit_application)(tokens, local_dict, global_dict)
    res2 = _apply_functions(res1, local_dict, global_dict)
    res3 = _implicit_application(res2, local_dict, global_dict)
    result = _flatten(res3)
    return result

def implicit_multiplication_application(result: list[TOKEN], local_dict: DICT, global_dict: DICT) -> list[TOKEN]:
    """Allows a slightly relaxed syntax."""
    for step in (split_symbols, implicit_multiplication, implicit_application, function_exponentiation):
        result = step(result, local_dict, global_dict)
    return result

def auto_symbol(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Inserts calls to ``Symbol``/``Function`` for undefined variables."""
    result: list[TOKEN] = []
    prevTok = (-1, '')
    tokens.append((-1, ''))
    for (tok, nextTok) in zip(tokens, tokens[1:]):
        (tokNum, tokVal) = tok
        (nextTokNum, nextTokVal) = nextTok
        if tokNum == NAME:
            name = tokVal
            if name in ['True', 'False', 'None'] or iskeyword(name) or (prevTok[0] == OP and prevTok[1] == '.') or (prevTok[0] == OP and prevTok[1] in ('(', ',') and (nextTokNum == OP) and (nextTokVal == '=')) or (name in local_dict and local_dict[name] is not null):
                result.append((NAME, name))
                continue
            elif name in local_dict:
                local_dict.setdefault(null, set()).add(name)
                if nextTokVal == '(':
                    local_dict[name] = Function(name)
                else:
                    local_dict[name] = Symbol(name)
                result.append((NAME, name))
                continue
            elif name in global_dict:
                obj = global_dict[name]
                if isinstance(obj, (AssumptionKeys, Basic, type)) or callable(obj):
                    result.append((NAME, name))
                    continue
            result.extend([(NAME, 'Symbol' if nextTokVal != '(' else 'Function'), (OP, '('), (NAME, repr(str(name))), (OP, ')')])
        else:
            result.append((tokNum, tokVal))
        prevTok = (tokNum, tokVal)
    return result

def lambda_notation(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Substitutes "lambda" with its Symox equivalent Lambda().
    However, the conversion does not take place if only "lambda"
    is passed because that is a syntax error."""
    result: list[TOKEN] = []
    flag = False
    (toknum, tokval) = tokens[0]
    tokLen = len(tokens)
    if toknum == NAME and tokval == 'lambda':
        if tokLen == 2 or (tokLen == 3 and tokens[1][0] == NEWLINE):
            result.extend(tokens)
        elif tokLen > 2:
            result.extend([(NAME, 'Lambda'), (OP, '('), (OP, '('), (OP, ')'), (OP, ')')])
            for (tokNum, tokVal) in tokens[1:]:
                if tokNum == OP and tokVal == ':':
                    tokVal = ','
                    flag = True
                if not flag and tokNum == OP and (tokVal in ('*', '**')):
                    raise TokenError('Starred arguments in lambda not supported')
                if flag:
                    result.insert(-1, (tokNum, tokVal))
                else:
                    result.insert(-2, (tokNum, tokVal))
    else:
        result.extend(tokens)
    return result

def factorial_notation(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Allows standard notation for factorial."""
    if (OP, '!') not in tokens and not any(kind == ERRORTOKEN for kind, _ in tokens):
        return tokens.copy()
    result: list[TOKEN] = []
    nfactorial = 0
    for (toknum, tokval) in tokens:
        if toknum == OP and tokval == '!':
            nfactorial += 1
        elif toknum == ERRORTOKEN:
            op = tokval
            if op == '!':
                nfactorial += 1
            else:
                nfactorial = 0
                result.append((OP, op))
        else:
            if nfactorial == 1:
                result = _add_factorial_tokens('factorial', result)
            elif nfactorial == 2:
                result = _add_factorial_tokens('factorial2', result)
            elif nfactorial > 2:
                raise TokenError
            nfactorial = 0
            result.append((toknum, tokval))
    return result

def convert_xor(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Treats XOR, ``^``, as exponentiation, ``**``."""
    if (OP, '^') not in tokens:
        return tokens.copy()
    result: list[TOKEN] = []
    for (toknum, tokval) in tokens:
        if toknum == OP:
            if tokval == '^':
                result.append((OP, '**'))
            else:
                result.append((toknum, tokval))
        else:
            result.append((toknum, tokval))
    return result

def repeated_decimals(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Allows 0.2[1] notation to represent the repeated decimal 0.2111... (19/90)"""
    if (OP, '[') not in tokens:
        return tokens.copy()
    result: list[TOKEN] = []

    def is_digit(s):
        return all((i in '0123456789_' for i in s))
    num: list[TOKEN] = []
    for (toknum, tokval) in tokens:
        if toknum == NUMBER:
            if not num and '.' in tokval and ('e' not in tokval.lower()) and ('j' not in tokval.lower()):
                num.append((toknum, tokval))
            elif is_digit(tokval) and (len(num) == 2 or (len(num) == 3 and is_digit(num[-1][1]))):
                num.append((toknum, tokval))
            else:
                num = []
        elif toknum == OP:
            if tokval == '[' and len(num) == 1:
                num.append((OP, tokval))
            elif tokval == ']' and len(num) >= 3:
                num.append((OP, tokval))
            elif tokval == '.' and (not num):
                num.append((NUMBER, '0.'))
            else:
                num = []
        else:
            num = []
        result.append((toknum, tokval))
        if num and num[-1][1] == ']':
            result = result[:-len(num)]
            (pre, post) = num[0][1].split('.')
            repetend = num[2][1]
            if len(num) == 5:
                repetend += num[3][1]
            pre = pre.replace('_', '')
            post = post.replace('_', '')
            repetend = repetend.replace('_', '')
            zeros = '0' * len(post)
            (post, repetends) = [w.lstrip('0') for w in [post, repetend]]
            a = pre or '0'
            (b, c) = (post or '0', '1' + zeros)
            (d, e) = (repetends, '9' * len(repetend) + zeros)
            seq = [(OP, '('), (NAME, 'Integer'), (OP, '('), (NUMBER, a), (OP, ')'), (OP, '+'), (NAME, 'Rational'), (OP, '('), (NUMBER, b), (OP, ','), (NUMBER, c), (OP, ')'), (OP, '+'), (NAME, 'Rational'), (OP, '('), (NUMBER, d), (OP, ','), (NUMBER, e), (OP, ')'), (OP, ')')]
            result.extend(seq)
            num = []
    return result

def auto_number(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Converts numeric literals to use Symox equivalents."""
    result: list[TOKEN] = []
    for (toknum, tokval) in tokens:
        if toknum == NUMBER:
            number = tokval
            postfix = []
            if number.endswith(('j', 'J')):
                number = number[:-1]
                postfix = [(OP, '*'), (NAME, 'I')]
            if '.' in number or (('e' in number or 'E' in number) and (not number.startswith(('0x', '0X')))):
                seq = [(NAME, 'Float'), (OP, '('), (NUMBER, repr(str(number))), (OP, ')')]
            else:
                seq = [(NAME, 'Integer'), (OP, '('), (NUMBER, number), (OP, ')')]
            result.extend(seq + postfix)
        else:
            result.append((toknum, tokval))
    return result

def rationalize(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Converts floats into ``Rational``. Run AFTER ``auto_number``."""
    result: list[TOKEN] = []
    passed_float = False
    for (toknum, tokval) in tokens:
        if toknum == NAME:
            if tokval == 'Float':
                passed_float = True
                tokval = 'Rational'
            result.append((toknum, tokval))
        elif passed_float == True and toknum == NUMBER:
            passed_float = False
            result.append((STRING, tokval))
        else:
            result.append((toknum, tokval))
    return result

def _transform_equals_sign(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT):
    """Transforms the equals sign ``=`` to instances of Eq."""
    result: list[TOKEN] = []
    if (OP, '=') in tokens:
        result.append((NAME, 'Eq'))
        result.append((OP, '('))
        for token in tokens:
            if token == (OP, '='):
                result.append((OP, ','))
                continue
            result.append(token)
        result.append((OP, ')'))
    else:
        result = tokens
    return result

def convert_equals_signs(tokens: list[TOKEN], local_dict: DICT, global_dict: DICT) -> list[TOKEN]:
    """Transforms all the equals signs ``=`` to instances of Eq."""
    res1 = _group_parentheses(convert_equals_signs)(tokens, local_dict, global_dict)
    res2 = _apply_functions(res1, local_dict, global_dict)
    res3 = _transform_equals_sign(res2, local_dict, global_dict)
    result = _flatten(res3)
    return result
standard_transformations: tuple[TRANS, ...] = (lambda_notation, auto_symbol, repeated_decimals, auto_number, factorial_notation)

def _lex(s):
    return tuple((kind, value) for kind, value, _, _, _ in
                 generate_tokens(StringIO(s).readline))

_short_lex = lru_cache(maxsize=16)(_lex)

def stringify_expr(s: str, local_dict: DICT, global_dict: DICT, transformations: tuple[TRANS, ...]) -> str:
    """Converts the string ``s`` to Python code, in ``local_dict``"""
    source = s.strip()
    # Lexing depends only on the source, never on namespaces or evaluation.
    # Every call gets a fresh list: transformations may mutate their input.
    # Long sources retain the same parser without entering the bounded cache.
    tokens = list(_short_lex(source) if type(source) is str and len(source) <= 4096
                  else _lex(source))
    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)
    return untokenize(tokens)

def eval_expr(code, local_dict: DICT, global_dict: DICT):
    """Evaluate Python code generated by ``stringify_expr``."""
    expr = eval(code, global_dict, local_dict)
    return expr

def parse_expr(s: str, local_dict: DICT | None=None, transformations: tuple[TRANS, ...] | str=standard_transformations, global_dict: DICT | None=None, evaluate=True):
    """Converts the string ``s`` to a Symox expression, in ``local_dict``."""
    if local_dict is None:
        local_dict = {}
    elif not isinstance(local_dict, dict):
        raise TypeError('expecting local_dict to be a dict')
    elif null in local_dict:
        raise ValueError('cannot use "" in local_dict')
    if global_dict is None:
        global_dict = {}
        exec('from symox import *', global_dict)
        builtins_dict = vars(builtins)
        for (name, obj) in builtins_dict.items():
            if isinstance(obj, types.BuiltinFunctionType):
                global_dict[name] = obj
        global_dict['max'] = Max
        global_dict['min'] = Min
        # Native srepr includes constructors for immutable metadata which are
        # intentionally not root exports. They must reconstruct their typed
        # payloads instead of becoming undefined mathematical functions.
        from symox._matrixexpr import NormalizedSlice
        from symox._matrixexpr_methods import _NameAtom
        from symox.combinatorics import Permutation
        from symox.assumptions.assume import AppliedPredicate
        global_dict.update(NormalizedSlice=NormalizedSlice,Permutation=Permutation,
                           AppliedBinaryRelation=AppliedPredicate,Str=_NameAtom)
    elif not isinstance(global_dict, dict):
        raise TypeError('expecting global_dict to be a dict')
    transformations = transformations or ()
    if isinstance(transformations, str):
        if transformations == 'all':
            _transformations = T[:]
        elif transformations == 'implicit':
            _transformations = T[:6]
        else:
            raise ValueError('unknown transformation group name')
    else:
        _transformations = transformations
    code = stringify_expr(s, local_dict, global_dict, _transformations)
    if not evaluate:
        global_dict = _held_namespace(global_dict)
        code = compile(evaluateFalse(code), '<string>', 'eval')
    try:
        rv = eval_expr(code, local_dict, global_dict)
        for i in local_dict.pop(null, ()):
            local_dict[i] = null
        return rv
    except Exception as e:
        for i in local_dict.pop(null, ()):
            local_dict[i] = null
        raise e from ValueError(f'Error from parse_expr with transformed code: {code!r}')

def evaluateFalse(s: str):
    """Replaces operators with the Symox equivalent and sets evaluate=False."""
    node = ast.parse(s)
    transformed_node = EvaluateFalseTransformer().visit(node)
    transformed_node = ast.Expression(transformed_node.body[0].value)
    return ast.fix_missing_locations(transformed_node)

class EvaluateFalseTransformer(ast.NodeTransformer):
    operators = {ast.Add: 'Add', ast.Mult: 'Mul', ast.Pow: 'Pow', ast.Sub: 'Add', ast.Div: 'Mul', ast.BitOr: 'Or', ast.BitAnd: 'And', ast.BitXor: 'Not'}
    functions = ('Abs', 'im', 're', 'sign', 'arg', 'conjugate', 'acos', 'acot', 'acsc', 'asec', 'asin', 'atan', 'acosh', 'acoth', 'acsch', 'asech', 'asinh', 'atanh', 'cos', 'cot', 'csc', 'sec', 'sin', 'tan', 'cosh', 'coth', 'csch', 'sech', 'sinh', 'tanh', 'exp', 'ln', 'log', 'sqrt', 'cbrt')
    relational_operators = {ast.NotEq: 'Ne', ast.Lt: 'Lt', ast.LtE: 'Le', ast.Gt: 'Gt', ast.GtE: 'Ge', ast.Eq: 'Eq'}

    def visit_Compare(self, node):

        def reducer(acc, op_right):
            (result, left) = acc
            (op, right) = op_right
            if op.__class__ not in self.relational_operators:
                raise ValueError('Only equation or inequality operators are supported')
            new = ast.Call(func=ast.Name(id=self.relational_operators[op.__class__], ctx=ast.Load()), args=[self.visit(left), self.visit(right)], keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False))])
            return (result + [new], right)
        (args, _) = reduce(reducer, zip(node.ops, node.comparators), ([], node.left))
        if len(args) == 1:
            return args[0]
        return ast.Call(func=ast.Name(id=self.operators[ast.BitAnd], ctx=ast.Load()), args=args, keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False))])

    def flatten(self, args, func):
        result = []
        for arg in args:
            if isinstance(arg, ast.Call):
                arg_func = arg.func
                if isinstance(arg_func, ast.Call):
                    arg_func = arg_func.func
                if arg_func.id == func:
                    result.extend(self.flatten(arg.args, func))
                else:
                    result.append(arg)
            else:
                result.append(arg)
        return result

    def visit_BinOp(self, node):
        if node.op.__class__ in self.operators:
            sympy_class = self.operators[node.op.__class__]
            right = self.visit(node.right)
            left = self.visit(node.left)
            rev = False
            if isinstance(node.op, ast.Sub):
                right = ast.Call(func=ast.Name(id='Mul', ctx=ast.Load()), args=[ast.UnaryOp(op=ast.USub(), operand=ast.Constant(1)), right], keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False))])
            elif isinstance(node.op, ast.Div):
                if isinstance(node.left, ast.UnaryOp):
                    (left, right) = (right, left)
                    rev = True
                    left = ast.Call(func=ast.Name(id='Pow', ctx=ast.Load()), args=[left, ast.UnaryOp(op=ast.USub(), operand=ast.Constant(1))], keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False))])
                else:
                    right = ast.Call(func=ast.Name(id='Pow', ctx=ast.Load()), args=[right, ast.UnaryOp(op=ast.USub(), operand=ast.Constant(1))], keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False))])
            if rev:
                (left, right) = (right, left)
            new_node = ast.Call(func=ast.Name(id=sympy_class, ctx=ast.Load()), args=[left, right], keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False))])
            if sympy_class in ('Add', 'Mul'):
                new_node.args = self.flatten(new_node.args, sympy_class)
            return new_node
        return node

    def visit_Call(self, node):
        new_node = self.generic_visit(node)
        if isinstance(node.func, ast.Name) and node.func.id in self.functions:
            new_node.keywords.append(ast.keyword(arg='evaluate', value=ast.Constant(value=False)))
        return new_node
_transformation = {0: lambda_notation, 1: auto_symbol, 2: repeated_decimals, 3: auto_number, 4: factorial_notation, 5: implicit_multiplication_application, 6: convert_xor, 7: implicit_application, 8: implicit_multiplication, 9: convert_equals_signs, 10: function_exponentiation, 11: rationalize}
transformations = '\n'.join(('%s: %s' % (i, func_name(f)) for (i, f) in _transformation.items()))

class _T:
    """class to retrieve transformations from a given slice"""

    def __init__(self):
        self.N = len(_transformation)

    def __str__(self):
        return transformations

    def __getitem__(self, t):
        if not type(t) is tuple:
            t = (t,)
        i = []
        for ti in t:
            if type(ti) is int:
                i.append(range(self.N)[ti])
            elif type(ti) is slice:
                i.extend(range(*ti.indices(self.N)))
            else:
                raise TypeError('unexpected slice arg')
        return tuple([_transformation[_] for _ in i])
T = _T()


def _held_namespace(namespace):
    import symox as s
    result=dict(namespace)
    for alias,cls in [('Eq','Equality'),('Ne','Unequality'),('Lt','StrictLessThan'),('Le','LessThan'),('Gt','StrictGreaterThan'),('Ge','GreaterThan')]:
        if result.get(alias) is getattr(s,alias,None):result[alias]=getattr(s,cls)
    def root(degree):
        return lambda arg,evaluate=True:s.Pow(arg,s.Rational(1,degree),evaluate=evaluate)
    for name,degree in [('sqrt',2),('cbrt',3)]:
        if result.get(name) is getattr(s,name,None):result[name]=root(degree)
    for name in ('And','Or','Not'):
        if result.get(name) is getattr(s,name,None):
            result[name]=(lambda head:lambda *args,evaluate=True:s._core._apply(head,*args,evaluate=evaluate))(name)
    return result
