"""Classical ciphers and exact educational public-key cryptosystems.

Integer algorithms use Python's exact modular arithmetic and symox number
theory; Hill ciphers, feedback sequences and polynomials use native symox
matrices and coefficient domains. No computation delegates to SymPy.
"""
from collections import Counter as _Counter
from collections.abc import Mapping as _Mapping
from functools import lru_cache as _lru_cache
from math import gcd as _gcd, isqrt as _isqrt, lcm as _lcm, prod as _prod
import operator as _operator
import random as _random
import re as _re
from string import ascii_uppercase as _UPPER, digits as _DIGITS, printable as _PRINTABLE
import warnings as _warnings
_NON_UPPER = _re.compile("[^A-Z]+")

__all__ = [
    'AZ','NonInvertibleCipherWarning','bg_private_key','bg_public_key',
    'bifid5_square','bifid6_square','bifid_square','check_and_join','cycle_list',
    'decipher_affine','decipher_atbash','decipher_bg','decipher_bifid',
    'decipher_bifid5','decipher_bifid6','decipher_elgamal','decipher_gm',
    'decipher_hill','decipher_kid_rsa','decipher_railfence','decipher_rot13',
    'decipher_rsa','decipher_shift','decipher_vigenere','decode_morse',
    'dh_private_key','dh_public_key','dh_shared_key','elgamal_private_key',
    'elgamal_public_key','encipher_affine','encipher_atbash','encipher_bg',
    'encipher_bifid','encipher_bifid5','encipher_bifid6','encipher_elgamal',
    'encipher_gm','encipher_hill','encipher_kid_rsa','encipher_railfence',
    'encipher_rot13','encipher_rsa','encipher_shift','encipher_substitution',
    'encipher_vigenere','encode_morse','gm_private_key','gm_public_key',
    'kid_rsa_private_key','kid_rsa_public_key','lfsr_autocorrelation',
    'lfsr_connection_polynomial','lfsr_sequence','padded_key','rsa_private_key',
    'rsa_public_key',
]

bifid5 = _UPPER.replace('J','')
bifid6 = _UPPER+_DIGITS
bifid10 = _PRINTABLE


def _integer(value):
    try:return _operator.index(value)
    except TypeError:raise TypeError('an exact integer is required') from None


def AZ(s=None):
    """Uppercase and retain A-Z; without an argument return that alphabet."""
    if not s:return _UPPER
    def clean(text):return _NON_UPPER.sub('',text.upper())
    return clean(s) if isinstance(s,str) else [clean(text) for text in s]


def check_and_join(phrase,symbols=None,filter=None):
    """Join strings, optionally validating or filtering an alphabet."""
    text = ''.join(phrase)
    if symbols is None:return text
    allowed = set(''.join(symbols))
    invalid = set(text)-allowed
    if invalid and not filter:
        raise ValueError('characters in phrase but not symbols: "%s"' % ''.join(sorted(invalid)))
    return text.translate(dict.fromkeys(map(ord,invalid))) if invalid else text


def _alphabet(symbols):
    symbols = ''.join(symbols)
    if not symbols or len(set(symbols)) != len(symbols):
        raise ValueError('an alphabet must be nonempty and contain distinct characters')
    return symbols


def padded_key(key,symbols):
    """Remove repeated key characters and append unused alphabet characters."""
    syms = ''.join(dict.fromkeys(symbols))
    if len(syms) != len(symbols):
        extra = ''.join(sorted({i for i in symbols if symbols.count(i) > 1}))
        raise ValueError('duplicate characters in symbols: %s' % extra)
    extra = set(key) - set(syms)
    if extra:
        raise ValueError('characters in key but not symbols: %s' % ''.join(sorted(extra)))
    return ''.join(dict.fromkeys(key + syms))


def _prepare(message,key='',symbols=None,default=None):
    if not symbols and default is None:
        alphabet = _UPPER
        message = AZ(message) if message else ''
        key = AZ(key) if key else ''
    else:
        alphabet = _alphabet(symbols or default)
    return (check_and_join(message,alphabet,filter=True),
            check_and_join(key,alphabet,filter=True),alphabet)


def cycle_list(k,n):
    """Rotate range(n) by the signed offset k."""
    k,n = _integer(k),_integer(n)
    if n <= 0:raise ValueError('cycle length must be positive')
    return [(i+k)%n for i in range(n)]


def encipher_shift(msg,key,symbols=None):
    """Apply an additive alphabet shift."""
    text,_,alphabet = _prepare(msg,symbols=symbols)
    offset = _integer(key)%len(alphabet)
    return text.translate(str.maketrans(alphabet,alphabet[offset:]+alphabet[:offset]))


def decipher_shift(msg,key,symbols=None):
    return encipher_shift(msg,-_integer(key),symbols)


def encipher_rot13(msg,symbols=None):return encipher_shift(msg,13,symbols)
def decipher_rot13(msg,symbols=None):return decipher_shift(msg,13,symbols)


def encipher_affine(msg,key,symbols=None,_inverse=False):
    """Apply a*x+b to alphabet positions, optionally using its modular inverse."""
    text,_,alphabet = _prepare(msg,symbols=symbols)
    a,b = map(_integer,key)
    modulus = len(alphabet)
    if _gcd(a,modulus) != 1:raise ValueError('the affine multiplier must be coprime to the alphabet length')
    if _inverse:a,b = pow(a,-1,modulus),-b*pow(a,-1,modulus)
    transformed = ''.join(alphabet[(a*i+b)%modulus] for i in range(modulus))
    return text.translate(str.maketrans(alphabet,transformed))


def decipher_affine(msg,key,symbols=None):return encipher_affine(msg,key,symbols,_inverse=True)
def encipher_atbash(msg,symbols=None):return encipher_affine(msg,(25,25),symbols)
def decipher_atbash(msg,symbols=None):return decipher_affine(msg,(25,25),symbols)


def encipher_substitution(msg,old,new=None):
    """Replace characters or dictionary-defined substrings without mutating a key."""
    if old is None:return msg.translate(str.maketrans('','',new or ''))
    if not isinstance(old,_Mapping):
        if new is None or len(old) != len(new):raise ValueError('old and new alphabets must have equal length')
        return msg.translate(str.maketrans(old,new))
    characters,substrings = {},{}
    for source,target in old.items():
        if not isinstance(source,str) or not source or not isinstance(target,str):
            raise ValueError('substitution keys must be nonempty strings and values must be strings')
        if len(source) == len(target) == 1:characters[source] = target
        else:substrings[source] = target
    if substrings:
        pattern = '|'.join(_re.escape(key) for key in sorted(substrings,key=len,reverse=True))
        msg = _re.sub(pattern,lambda match:substrings[match.group()],msg)
    return msg.translate(str.maketrans(characters))


def _vigenere(msg,key,symbols,inverse):
    text,key,alphabet = _prepare(msg,key,symbols)
    if not key:raise ValueError('Vigenere requires a nonempty key in the alphabet')
    positions = {c:i for i,c in enumerate(alphabet)}
    n,period = len(alphabet),len(key)
    out = [''] * len(text)
    # one translation table per key letter, applied to each residue class
    for j,c in enumerate(key):
        shift = inverse*positions[c]
        table = str.maketrans(alphabet,''.join(alphabet[(i+shift)%n] for i in range(n)))
        out[j::period] = text[j::period].translate(table)
    return ''.join(out)


def encipher_vigenere(msg,key,symbols=None):return _vigenere(msg,key,symbols,1)
def decipher_vigenere(msg,key,symbols=None):return _vigenere(msg,key,symbols,-1)


def _hill_matrix(key):
    from .. import Matrix
    if not isinstance(key,Matrix) and not getattr(key,'_symox_sparse',False):
        raise TypeError('a Hill key must be a native matrix')
    if key.rows != key.cols or key.rows < 1:raise ValueError('a Hill key must be nonempty and square')
    for value in key:
        _integer(value)
    return key


def _hill_transform(text,key,alphabet,padding):
    positions = {c:i for i,c in enumerate(alphabet)}
    values = [positions[c] for c in text]
    size,modulus = key.rows,len(alphabet)
    values.extend([padding]*(-len(values)%size))
    coefficients = [[_integer(key[i,j]) for j in range(size)] for i in range(size)]
    result = []
    for offset in range(0,len(values),size):
        for row in coefficients:
            result.append(alphabet[sum(a*b for a,b in zip(row,values[offset:offset+size]))%modulus])
    return ''.join(result)


def encipher_hill(msg,key,symbols=None,pad='Q'):
    """Multiply column blocks by a native integral matrix modulo the alphabet size."""
    key = _hill_matrix(key)
    if not isinstance(pad,str) or len(pad) != 1:raise ValueError('padding must be one character')
    text,padding,alphabet = _prepare(msg,pad,symbols)
    if len(text)%key.rows and len(padding) != 1:raise ValueError('padding must belong to the alphabet')
    return _hill_transform(text,key,alphabet,alphabet.index(padding) if padding else 0)


def decipher_hill(msg,key,symbols=None):
    key = _hill_matrix(key)
    text,_,alphabet = _prepare(msg,symbols=symbols)
    inverse = key.inv_mod(len(alphabet))
    return _hill_transform(text,inverse,alphabet,0)


def _bifid_setup(msg,key,symbols):
    text,key,alphabet = _prepare(msg,key,symbols,default=bifid10)
    width = _isqrt(len(alphabet))
    if width*width != len(alphabet):raise ValueError('Bifid requires a square alphabet length')
    return text,padded_key(key,alphabet),width


def encipher_bifid(msg,key,symbols=None):
    """Fractionate square-table coordinates, joining all rows before columns."""
    text,table,width = _bifid_setup(msg,key,symbols)
    where = {c:i for i,c in enumerate(table)}
    indices = [where[c] for c in text]
    coordinates = [i//width for i in indices]+[i%width for i in indices]
    return ''.join(table[width*a+b] for a,b in zip(coordinates[::2],coordinates[1::2]))


def decipher_bifid(msg,key,symbols=None):
    text,table,width = _bifid_setup(msg,key,symbols)
    where = {c:i for i,c in enumerate(table)}
    coordinates = [coordinate for c in text for coordinate in divmod(where[c],width)]
    return ''.join(table[width*a+b] for a,b in zip(coordinates[:len(text)],coordinates[len(text):]))


def bifid_square(key):
    """Build a native symbolic matrix from distinct square-table characters."""
    from .. import Matrix,Symbol
    alphabet = ''.join(dict.fromkeys(''.join(key)))
    width = _isqrt(len(alphabet))
    if width*width != len(alphabet):raise ValueError('Bifid requires a square alphabet length')
    return Matrix(width,width,[Symbol(c) for c in alphabet])


def _bifid_key(key,alphabet):
    return padded_key(check_and_join((key or '').upper(),alphabet,filter=True),alphabet)


def encipher_bifid5(msg,key):return encipher_bifid(msg.upper(),'',_bifid_key(key,bifid5))
def decipher_bifid5(msg,key):return decipher_bifid(msg.upper(),'',_bifid_key(key,bifid5))
def encipher_bifid6(msg,key):return encipher_bifid(msg.upper(),'',_bifid_key(key,bifid6))
def decipher_bifid6(msg,key):return decipher_bifid(msg.upper(),'',_bifid_key(key,bifid6))
def bifid5_square(key=None):return bifid_square(_bifid_key(key,bifid5))
def bifid6_square(key=None):return bifid_square(_bifid_key(key,bifid6))


def _rail_positions(length,rails):
    rails = _integer(rails)
    if rails < 1:raise ValueError('the number of rails must be positive')
    if rails == 1:return [0]*length
    period = 2*(rails-1)
    return [min(i%period,period-i%period) for i in range(length)]


def encipher_railfence(message,rails):
    positions = _rail_positions(len(message),rails)
    rows = [[] for i in range(min(_integer(rails),len(message)))]
    for char,row in zip(message,positions):rows[row].append(char)
    return ''.join(''.join(row) for row in rows)


def decipher_railfence(ciphertext,rails):
    positions = _rail_positions(len(ciphertext),rails)
    counts = _Counter(positions)
    starts,offset = {},0
    for row in sorted(counts):starts[row],offset = offset,offset+counts[row]
    text = []
    for row in positions:
        text.append(ciphertext[starts[row]])
        starts[row] += 1
    return ''.join(text)


class NonInvertibleCipherWarning(RuntimeWarning):
    """RSA's repeated prime factors restrict the invertible message domain."""
    def __init__(self,msg):
        self.fullMessage = msg
        super().__init__(msg)
    def __str__(self):return '\n\t'+self.fullMessage
    def warn(self,stacklevel=3):_warnings.warn(self,stacklevel=stacklevel)


@_lru_cache(maxsize=128)
def _rsa_modulus(factors):
    from .. import factorint,isprime
    primes = _Counter()
    for value in factors:
        if value < 2:raise ValueError('RSA modulus factors must be greater than one')
        if isprime(value):primes[value] += 1
        else:primes.update({int(p):int(e) for p,e in factorint(value).items()})
    modulus = _prod(factors)
    euler = _prod((p-1)*p**(power-1) for p,power in primes.items())
    carmichael = 1
    for p,power in primes.items():
        term = 2**(power-2) if p == 2 and power >= 3 else (p-1)*p**(power-1)
        carmichael = _lcm(carmichael,term)
    return modulus,euler,carmichael,any(power > 1 for power in primes.values())


def _rsa_key(args,private,totient='Euler',index=None,multipower=None):
    if len(args) < 2:return False
    if totient not in ('Euler','Carmichael'):raise ValueError("totient must be 'Euler' or 'Carmichael'")
    if index is not None:
        index = _integer(index)
        if totient != 'Carmichael':raise ValueError('index requires the Carmichael totient')
    factors = tuple(_integer(value) for value in args[:-1])
    exponent = _integer(args[-1])
    modulus,euler,carmichael,repeated = _rsa_modulus(factors)
    if repeated and not multipower:
        NonInvertibleCipherWarning(
            'Repeated prime factors restrict invertible RSA messages to units modulo %s; '
            'use multipower=True to acknowledge this domain.' % modulus).warn(stacklevel=4)
    period = euler if totient == 'Euler' else carmichael
    if _gcd(exponent,period) != 1:return False
    if private:exponent = pow(exponent,-1,period)
    if index is not None:exponent = exponent%period+index*period
    return modulus,exponent


def rsa_public_key(*args,**kwargs):
    """Generate ordinary, multiprime or prime-power RSA public modulus/exponent."""
    return _rsa_key(args,False,**kwargs)


def rsa_private_key(*args,**kwargs):
    """Invert an RSA exponent using the Euler or Carmichael totient."""
    return _rsa_key(args,True,**kwargs)


def _rsa_power(message,key,factors):
    modulus,exponent = map(_integer,key)
    message = _integer(message)
    if modulus < 2:raise ValueError('RSA modulus must be greater than one')
    if factors:
        factors = tuple(map(_integer,factors))
        if all(p > 1 for p in factors) and _prod(factors) == modulus and all(_gcd(p,q) == 1 for j,p in enumerate(factors) for q in factors[j+1:]):
            # Pairwise-coprime CRT works for prime powers as well as primes.
            result = 0
            for p in factors:
                partial = modulus//p
                result += pow(message,exponent,p)*partial*pow(partial,-1,p)
            return result%modulus
    return pow(message,exponent,modulus)


def encipher_rsa(i,key,factors=None):return _rsa_power(i,key,factors)
def decipher_rsa(i,key,factors=None):return _rsa_power(i,key,factors)


def _kid_keys(a,b,A,B):
    a,b,A,B = map(_integer,(a,b,A,B))
    product = a*b-1
    if product == 0:raise ValueError('Kid RSA requires a*b != 1')
    public,private = A*product+a,B*product+b
    modulus = (public*private-1)//product
    if modulus < 2:raise ValueError('Kid RSA modulus must be greater than one')
    return modulus,public,private


def kid_rsa_public_key(a,b,A,B):
    n,e,d = _kid_keys(a,b,A,B)
    return n,e


def kid_rsa_private_key(a,b,A,B):
    n,e,d = _kid_keys(a,b,A,B)
    return n,d


def encipher_kid_rsa(msg,key):
    modulus,multiplier = map(_integer,key)
    if modulus < 2:raise ValueError('Kid RSA modulus must be greater than one')
    return _integer(msg)*multiplier%modulus


def decipher_kid_rsa(msg,key):return encipher_kid_rsa(msg,key)


_random_source = _random.Random()


def _draw(seed=None):
    if seed is None:return _random_source.randrange
    if isinstance(seed,(list,tuple)) or getattr(seed,'is_Tuple',False):
        values = iter(tuple(seed.args if getattr(seed,'is_Tuple',False) else seed))
        def sample(start,stop=None):
            start,stop = (0,start) if stop is None else (start,stop)
            if stop <= start:raise ValueError('empty random range')
            for value in values:
                value = _integer(value)
                if start <= value < stop:return value
            raise ValueError('seed sequence exhausted before an admissible value was found')
        return sample
    try:seed = _integer(seed)
    except TypeError:raise ValueError('seed must be an integer or sequence of integers') from None
    return _random.Random(seed).randrange


@_lru_cache(maxsize=32)
def _prime_group(digit):
    from .. import nextprime,primitive_root
    if digit < 1:raise ValueError('digit must be a positive bit exponent')
    p = int(nextprime(1<<digit))
    return p,int(primitive_root(p))


def elgamal_private_key(digit=10,seed=None):
    """Generate p=nextprime(2**digit), a primitive generator and a secret exponent."""
    p,g = _prime_group(_integer(digit))
    return p,g,_draw(seed)(2,p)


def elgamal_public_key(key):
    p,g,secret = map(_integer,key)
    if p <= 2 or not 0 < g < p:raise ValueError('invalid multiplicative group parameters')
    return p,g,pow(g,secret,p)


def encipher_elgamal(i,key,seed=None):
    p,g,public = map(_integer,key)
    i = _integer(i)
    if p <= 2 or not 0 < g < p or not 0 < public < p:raise ValueError('invalid ElGamal public key')
    if not 0 <= i < p:raise ValueError('message must be between zero and p-1')
    nonce = _draw(seed)(2,p)
    return pow(g,nonce,p),i*pow(public,nonce,p)%p


def decipher_elgamal(msg,key):
    p,_,secret = map(_integer,key)
    first,second = map(_integer,msg)
    if p <= 2:raise ValueError('invalid ElGamal modulus')
    return second*pow(first,-secret,p)%p


def dh_private_key(digit=10,seed=None):return elgamal_private_key(digit,seed)
def dh_public_key(key):return elgamal_public_key(key)


def dh_shared_key(key,b):
    """Combine a received public key with the other participant's exponent b."""
    p,g,public = map(_integer,key)
    b = _integer(b)
    if not 1 < b < p:raise ValueError('private exponent must be greater than one and below p')
    return pow(public,b,p)


def _prime_pair(p,q,blum=False):
    from .. import isprime
    p,q = _integer(p),_integer(q)
    if p == q or not isprime(p) or not isprime(q):raise ValueError('two distinct primes are required')
    if blum:
        if p%4 != 3 or q%4 != 3:raise ValueError('Blum primes must be congruent to 3 modulo 4')
    elif p == 2 or q == 2:raise ValueError('two odd primes are required')
    return p,q


def _legendre(value,prime):
    result = pow(value,(prime-1)//2,prime)
    return -1 if result == prime-1 else result


def gm_private_key(p,q,a=None):return _prime_pair(p,q)


def gm_public_key(p,q,a=None,seed=None):
    """Select a common quadratic nonresidue for a Goldwasser-Micali public key."""
    p,q = _prime_pair(p,q)
    modulus = p*q
    if a is not None:
        a = _integer(a)
        return (a,modulus) if _legendre(a,p) == _legendre(a,q) == -1 else False
    sample = _draw(seed)
    while True:
        candidate = sample(modulus)
        if _legendre(candidate,p) == _legendre(candidate,q) == -1:return candidate,modulus


def _message_bits(message):
    message = _integer(message)
    if message < 0:raise ValueError('message must be a nonnegative integer')
    return [int(bit) for bit in bin(message)[2:]] if message else []


def encipher_gm(i,key,seed=None):
    nonresidue,modulus = map(_integer,key)
    if modulus < 2 or _gcd(nonresidue,modulus) != 1:raise ValueError('invalid Goldwasser-Micali public key')
    sample = _draw(seed)
    encrypted = []
    for bit in _message_bits(i):
        while True:
            nonce = sample(modulus)
            if _gcd(nonce,modulus) == 1:break
        encrypted.append((nonce*nonce*(nonresidue if bit else 1))%modulus)
    return encrypted


def decipher_gm(message,key):
    p,q = _prime_pair(*key)
    value = 0
    for ciphertext in message:
        ciphertext = _integer(ciphertext)
        if not 0 < ciphertext < p*q or _gcd(ciphertext,p*q) != 1:raise ValueError('ciphertext components must be units modulo p*q')
        bit = not (_legendre(ciphertext,p) == _legendre(ciphertext,q) == 1)
        value = 2*value+bit
    return value


def bg_private_key(p,q):return _prime_pair(p,q,blum=True)
def bg_public_key(p,q):return _prod(bg_private_key(p,q))


def encipher_bg(i,key,seed=None):
    """XOR message bits with iterated-square parity, retaining the terminal state."""
    bits = _message_bits(i)
    modulus = _integer(key)
    if modulus <= 3:raise ValueError('Blum-Goldwasser modulus must be greater than three')
    nonce = _draw(seed)(2,modulus)
    state = nonce*nonce%modulus
    encrypted = []
    for bit in bits:
        encrypted.append(bit^(state&1))
        state = state*state%modulus
    return encrypted,state


def decipher_bg(message,key):
    p,q = bg_private_key(*key)
    encrypted,terminal = message
    terminal = _integer(terminal)
    if not 0 <= terminal < p*q:raise ValueError('terminal state is outside the residue ring')
    bits = [_integer(bit) for bit in encrypted]
    if any(bit not in (0,1) for bit in bits):raise ValueError('ciphertext must contain bits')
    # Invert each squaring on the quadratic residues. Keep these exponents
    # integral and reduce them modulo p-1/q-1; long messages never use floats.
    a = pow(terminal,pow((p+1)//4,len(bits),p-1),p)
    b = pow(terminal,pow((q+1)//4,len(bits),q-1),q)
    state = (a*q*pow(q,-1,p)+b*p*pow(p,-1,q))%(p*q)
    result = 0
    for bit in bits:
        result = 2*result+(bit^(state&1))
        state = state*state%(p*q)
    return result


_MORSE_LETTERS = '.- -... -.-. -.. . ..-. --. .... .. .--- -.- .-.. -- -. --- .--. --.- .-. ... - ..- ...- .-- -..- -.-- --..'.split()
_MORSE_DIGITS = '----- .---- ..--- ...-- ....- ..... -.... --... ---.. ----.'.split()
char_morse = dict(zip(_UPPER,_MORSE_LETTERS))
char_morse.update(zip(_DIGITS,_MORSE_DIGITS))
char_morse.update({'.':'.-.-.-',',':'--..--',':':'---...',';':'-.-.-.','?':'..--..','-':'-....-',
                  '_':'..--.-','(':'-.--.',')':'-.--.-',"'":'.----.','=':'-...-','+':'.-.-.',
                  '/':'-..-.','@':'.--.-.','$':'...-..-','!':'-.-.--'})
morse_char = {code:char for char,code in char_morse.items()}


def encode_morse(msg,sep='|',mapping=None):
    """Encode mapped characters, using a doubled separator between words."""
    if not isinstance(sep,str) or not sep:raise ValueError('Morse separator must be nonempty')
    mapping = char_morse if mapping is None else mapping
    if sep in mapping:raise ValueError('Morse separator conflicts with a mapped character')
    # Filtering precedes splitting: unsupported punctuation disappears while
    # original word boundaries and the trailing whitespace marker are retained.
    filtered = ''.join(c for c in msg if c in mapping or c.isspace())
    words = [sep.join(mapping[c] for c in word) for word in filtered.split()]
    return (2*sep).join(words)+(2*sep if msg and msg[-1].isspace() else '')


def decode_morse(msg,sep='|',mapping=None):
    if not isinstance(sep,str) or not sep:raise ValueError('Morse separator must be nonempty')
    mapping = morse_char if mapping is None else mapping
    text = msg
    while text.startswith(2*sep):text = text[len(2*sep):]
    while text.endswith(2*sep):text = text[:-len(2*sep)]
    if not text:return ''
    return ' '.join(''.join(mapping[code] for code in word.split(sep)) for word in text.split(2*sep))


def _feedback_field(values):
    from .. import GF,isprime
    if not values:raise ValueError('a finite-field sequence must be nonempty')
    first = values[0]
    parent = getattr(first,'parent',None)
    field = parent() if callable(parent) else None
    if field is not None:
        if not getattr(field,'is_FF',False):raise TypeError('feedback elements must belong to a prime field')
        modulus = int(field.characteristic())
    else:
        modulus_method = getattr(first,'modulus',None)
        if not callable(modulus_method):raise TypeError('feedback elements must belong to a prime field')
        modulus = _integer(modulus_method())
        field = GF(modulus)
    # The reference also runs over the ring Z/nZ of a composite GF(n).
    for value in values:
        parent = getattr(value,'parent',None)
        if callable(parent) and parent() != field:raise ValueError('all feedback elements must belong to the same field')
        method = getattr(value,'modulus',None)
        if callable(method) and int(method()) != modulus:raise ValueError('all feedback elements must have the same modulus')
    return field,modulus


def lfsr_sequence(key,fill,n):
    """Generate feedback terms over the key's native prime field."""
    if not isinstance(key,list) or not isinstance(fill,list):raise TypeError('key and fill must be lists')
    if not key or len(key) != len(fill):raise ValueError('key and fill must have the same positive length')
    n = _integer(n)
    if n < 0:raise ValueError('sequence length must be nonnegative')
    field,modulus = _feedback_field(key+fill)
    state = [int(field(value))%modulus for value in fill]
    coefficients = [int(value)%modulus for value in key]
    result = []
    for i in range(n):
        result.append(field(state[0]))
        following = sum(c*v for c,v in zip(coefficients,state))%modulus
        state = state[1:]+[following]
    return result


def lfsr_autocorrelation(L,P,k):
    """Return exact cyclic binary autocorrelation at an arbitrary signed lag."""
    from .. import Rational
    if not isinstance(L,list):raise TypeError('the sequence must be a list')
    P,k = _integer(P),_integer(k)
    if not 0 < P <= len(L):raise ValueError('period must lie between one and the sequence length')
    k %= P
    bits = [_integer(value) for value in L[:P]]
    return Rational(sum(1 if (bits[i]+bits[(i+k)%P])%2 == 0 else -1 for i in range(P)),P)


def lfsr_connection_polynomial(s):
    """Recover the minimal connection polynomial with Berlekamp-Massey."""
    from .. import Symbol,Integer,Add
    if not s:return Integer(1)
    field,modulus = _feedback_field(s)
    samples = [int(value)%modulus for value in s]
    connection,backup = [1],[1]
    length,gap,last = 0,1,1
    for position,sample in enumerate(samples):
        discrepancy = (sample+sum(connection[j]*samples[position-j] for j in range(1,min(length+1,len(connection)))))%modulus
        if discrepancy == 0:
            gap += 1
            continue
        previous = connection[:]
        scale = discrepancy*pow(last,-1,modulus)%modulus
        connection += [0]*max(0,len(backup)+gap-len(connection))
        for j,coefficient in enumerate(backup):connection[j+gap] = (connection[j+gap]-scale*coefficient)%modulus
        if 2*length <= position:
            length,backup,last,gap = position+1-length,previous,discrepancy,1
        else:gap += 1
    x = Symbol('x')
    return Add(*(Integer(coefficient)*x**j for j,coefficient in enumerate(connection) if coefficient))
