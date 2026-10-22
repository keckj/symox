"""Known vectors, differential checks and algebraic invariants for crypto."""
from copy import deepcopy
from math import gcd
from random import Random
import importlib
import subprocess
import sys

import pytest
import sympy as reference
import sympy.crypto.crypto as ref
import symox as s
import symox.crypto.crypto as crypto


CLASSICAL_CASES = [
    ('AZ',()),('AZ',('à bientôt, Straße! 123',)),('AZ',(['Ab c','123','déf'],)),
    ('check_and_join',(['ab','cd'],)),('check_and_join',('a!b c','abc',True)),
    ('padded_key',('BANANA','ABCDEFGHIJKLMNOPQRSTUVWXYZ')),
    ('encipher_substitution',('go navy','go navy','GO NAVY')),
    ('encipher_substitution',('abc',{'a':'x','bc':'yz'})),
    ('encipher_substitution',('abab',{'ab':'x','a':'y'})),
    ('encipher_substitution',('abc',None,'b')),
    ('encipher_hill',('ACT',lambda lib:lib.Matrix([[6,24,1],[13,16,10],[20,17,15]]))),
    ('encipher_hill',('ABA',lambda lib:lib.Matrix([[1,2],[3,5]]))),
    ('encipher_hill',('ABA',lambda lib:lib.Matrix([[1,2],[3,5]]),None,'Z')),
    ('decipher_hill',('CFA',lambda lib:lib.Matrix([[1,2],[3,5]]))),
    ('bifid5_square',('KEYWORD',)),('bifid6_square',('B0',)),('bifid_square',('abcd',)),
    ('encipher_bifid',('abc','b','abcd')),('decipher_bifid',('bdb','b','abcd')),
    ('encipher_bifid',('hello World!\n','key')),('decipher_bifid',('bdb','b','abcd')),
]


def _convert(value):
    if isinstance(value,s.Matrix):return value.to_sympy()
    if isinstance(value,s.Expr):return value.to_sympy()
    return value


@pytest.mark.parametrize('name,args',CLASSICAL_CASES)
def test_classical_differential(name,args):
    expected_args = [a(reference) if callable(a) else deepcopy(a) for a in args]
    actual_args = [a(s) if callable(a) else deepcopy(a) for a in args]
    assert _convert(getattr(crypto,name)(*actual_args)) == getattr(ref,name)(*expected_args)


@pytest.mark.parametrize('size',[1,2,5,26])
@pytest.mark.parametrize('shift',[-100,-1,0,1,31,10**20])
def test_cycles(size,shift):
    assert crypto.cycle_list(shift,size) == ref.cycle_list(shift,size)


@pytest.mark.parametrize('alphabet',[None,'abcxyz','0123456789','ABCD','αβγδ'])
@pytest.mark.parametrize('shift',[-53,-1,0,1,29])
def test_shift_differential_and_inverse(alphabet,shift):
    message = 'Hello, world! AZ 01234' if alphabet is None else alphabet*3
    ciphertext = crypto.encipher_shift(message,shift,alphabet)
    assert ciphertext == ref.encipher_shift(message,shift,alphabet)
    plain = crypto.decipher_shift(ciphertext,shift,alphabet)
    assert plain == (crypto.AZ(message) if alphabet is None else message)
    for name in ('encipher_rot13','decipher_rot13'):
        assert getattr(crypto,name)(message,alphabet) == getattr(ref,name)(message,alphabet)


@pytest.mark.parametrize('alphabet',['ABCDEFGHIJKLMNOPQRSTUVWXYZ','ABCD','0123456789','αβγδεζη'])
@pytest.mark.parametrize('a,b',[(1,0),(-1,1),(3,16),(7,-99)])
def test_affine_reference_and_inverse(alphabet,a,b):
    if gcd(a,len(alphabet)) != 1:
        with pytest.raises(ValueError):crypto.encipher_affine(alphabet,(a,b),alphabet)
        return
    text = alphabet*2
    ciphertext = crypto.encipher_affine(text,(a,b),alphabet)
    assert ciphertext == ref.encipher_affine(text,(a,b),alphabet)
    assert crypto.decipher_affine(ciphertext,(a,b),alphabet) == text


@pytest.mark.parametrize('name',['encipher_atbash','decipher_atbash'])
@pytest.mark.parametrize('alphabet',[None,'ABCD','abcdef'])
def test_atbash_reference(name,alphabet):
    msg = 'ABC xyz!' if alphabet is None else alphabet*2
    assert getattr(crypto,name)(msg,alphabet) == getattr(ref,name)(msg,alphabet)


@pytest.mark.parametrize('alphabet',[None,'ABCD','0123456789','αβγδε'])
@pytest.mark.parametrize('key_index',[1,2,3])
def test_vigenere_reference_and_inverse(alphabet,key_index):
    symbols = alphabet or crypto.AZ()
    key = symbols[:key_index]
    message = symbols*2
    encrypted = crypto.encipher_vigenere(message,key,alphabet)
    assert encrypted == ref.encipher_vigenere(message,key,alphabet)
    assert crypto.decipher_vigenere(encrypted,key,alphabet) == message
    assert crypto.decipher_vigenere(message,key,alphabet) == ref.decipher_vigenere(message,key,alphabet)


@pytest.mark.parametrize('size',[2,3])
@pytest.mark.parametrize('message',['ACT','PAY MORE MONEY','ABCDEF','S','ABCDEFGHIJKLMNOPQRSTUVWXYZ'])
def test_hill_vectors_and_inverse(size,message):
    rows = [[1,2],[3,5]] if size == 2 else [[6,24,1],[13,16,10],[20,17,15]]
    matrix = s.Matrix(rows)
    ciphertext = crypto.encipher_hill(message,matrix)
    assert ciphertext == ref.encipher_hill(message,reference.Matrix(rows))
    result = crypto.decipher_hill(ciphertext,matrix)
    clean = crypto.AZ(message)
    assert result == clean+'Q'*(-len(clean)%size)
    if size == 3 and message == 'ACT':assert ciphertext == 'POH'


@pytest.mark.parametrize('dimension',[5,6])
@pytest.mark.parametrize('key',['','KEYWORD','B0','A repeat key','SECRET'])
@pytest.mark.parametrize('message',['ATTACK AT DAWN','meet me at 12','JIGSAW 0123456789'])
def test_bifid_reference_and_inverse(dimension,key,message):
    encipher = getattr(crypto,'encipher_bifid%d' % dimension)
    decipher = getattr(crypto,'decipher_bifid%d' % dimension)
    ciphertext = encipher(message,key)
    assert ciphertext == getattr(ref,'encipher_bifid%d' % dimension)(message,key)
    alphabet = getattr(crypto,'bifid%d' % dimension)
    expected = ''.join(c for c in message.upper() if c in alphabet)
    assert decipher(ciphertext,key) == expected
    assert getattr(crypto,'bifid%d_square' % dimension)(key).to_sympy() == getattr(ref,'bifid%d_square' % dimension)(key)


@pytest.mark.parametrize('message',['WE ARE DISCOVERED. FLEE AT ONCE','hello world','αβ\n🙂🙂','x',''])
@pytest.mark.parametrize('rails',[1,2,3,4,20])
def test_railfence_differential_and_inverse(message,rails):
    ciphertext = crypto.encipher_railfence(message,rails)
    assert ciphertext == ref.encipher_railfence(message,rails)
    assert crypto.decipher_railfence(ciphertext,rails) == message


@pytest.mark.parametrize('message',['ABC','SMS ','SMS\n','',' ','12345','67890','!@#$%^&*()_+','HELLO  WORLD!','lowercase'])
@pytest.mark.parametrize('separator',['|','`','``',' / '])
def test_morse_differential(message,separator):
    encoded = crypto.encode_morse(message,separator)
    assert encoded == ref.encode_morse(message,separator)
    if encoded and encoded.strip(separator):
        assert crypto.decode_morse(encoded,separator) == ref.decode_morse(encoded,separator)


def test_morse_custom_mapping_is_unchanged():
    mapping = {'X':'.','Y':'-'}
    before = dict(mapping)
    encoded = crypto.encode_morse('XX YY ',mapping=mapping)
    assert encoded == '.|.||-|-||'
    assert mapping == before
    assert crypto.decode_morse(encoded,mapping={'.':'X','-':'Y'}) == 'XX YY'
    assert crypto.decode_morse('') == ''
    with pytest.raises(KeyError):crypto.decode_morse('.....----')


@pytest.mark.parametrize('factors,e',[
    ((61,53),17),((23,29),5),((3,5,7,11,13),7),((53,41,43,47),97),
    ((2,3),1),((5,),3),((5,5,7),7),((8,8),5),((15,77),13),
])
@pytest.mark.parametrize('totient',['Euler','Carmichael'])
def test_rsa_key_reference(factors,e,totient):
    for name in ('rsa_public_key','rsa_private_key'):
        actual = getattr(crypto,name)(*factors,e,totient=totient,multipower=True)
        expected = getattr(ref,name)(*factors,e,totient=totient,multipower=True)
        assert actual == expected


@pytest.mark.parametrize('index',[0,1,2,100])
def test_rsa_index_and_warning(index):
    for name in ('rsa_public_key','rsa_private_key'):
        assert getattr(crypto,name)(61,53,17,totient='Carmichael',index=index) == getattr(ref,name)(61,53,17,totient='Carmichael',index=index)
        with pytest.warns(crypto.NonInvertibleCipherWarning):
            assert getattr(crypto,name)(5,5,7,7)
    assert crypto.rsa_public_key(3,5,30) is False
    assert crypto.rsa_private_key(3,5,30) is False
    assert crypto.rsa_public_key(5) is False
    assert crypto.rsa_private_key() is False
    warning = crypto.NonInvertibleCipherWarning('message')
    assert warning.fullMessage == 'message' and str(warning) == '\n\tmessage'


@pytest.mark.parametrize('factors,e',[((61,53),17),((3,5,7,11),7),((5,5,7),7)])
def test_rsa_complete_residue_domain(factors,e):
    public = crypto.rsa_public_key(*factors,e,multipower=True)
    private = crypto.rsa_private_key(*factors,e,multipower=True)
    for message in range(public[0]):
        if len(set(factors)) != len(factors) and gcd(message,public[0]) != 1:continue
        encrypted = crypto.encipher_rsa(message,public)
        assert crypto.decipher_rsa(encrypted,private) == message
        assert encrypted == pow(message,e,public[0])
        assert crypto.encipher_rsa(message,public,factors=list(factors)) == encrypted
        assert crypto.decipher_rsa(encrypted,private,factors=list(factors)) == message


def test_rsa_large_integers_and_crt_fallback():
    p,q,e = (1<<127)-1,(1<<521)-1,65537
    public = crypto.rsa_public_key(p,q,e)
    private = crypto.rsa_private_key(p,q,e)
    message = (1<<500)+987654321
    encrypted = crypto.encipher_rsa(message,public,[p,q])
    assert encrypted == pow(message,e,p*q)
    assert crypto.decipher_rsa(encrypted,private,[p,q]) == message
    for invalid in ([p,p],[3,5],[],[0,p*q]):
        assert crypto.decipher_rsa(encrypted,private,invalid) == message
    key = (45,7)
    for message in range(45):
        assert crypto.encipher_rsa(message,key,[9,5]) == pow(message,7,45)


@pytest.mark.parametrize('parameters',[(1,2,1,1),(1,2,2,1),(1,2,1,2),(3,5,7,11)])
def test_kid_rsa_all_residues(parameters):
    public = crypto.kid_rsa_public_key(*parameters)
    private = crypto.kid_rsa_private_key(*parameters)
    assert public == ref.kid_rsa_public_key(*parameters)
    assert private == ref.kid_rsa_private_key(*parameters)
    for message in range(public[0]):
        assert crypto.decipher_kid_rsa(crypto.encipher_kid_rsa(message,public),private) == message


@pytest.mark.parametrize('bits',[2,5,10,32,100])
@pytest.mark.parametrize('seed',[0,1,3,42])
def test_elgamal_and_dh_seeded_reference(bits,seed):
    # SymPy can consume its global RNG on a cold primitive-root factor cache.
    # Group derivation must precede the seeded reference draw.
    reference.primitive_root(reference.nextprime(1<<bits))
    private = crypto.elgamal_private_key(bits,seed)
    assert private == ref.elgamal_private_key(bits,seed)
    assert private == crypto.dh_private_key(bits,seed) == ref.dh_private_key(bits,seed)
    public = crypto.elgamal_public_key(private)
    assert public == ref.elgamal_public_key(private) == crypto.dh_public_key(private)
    p,g,a = private
    assert s.isprime(p) and s.is_primitive_root(g,p)
    for message in (0,1,p//2,p-1):
        encrypted = crypto.encipher_elgamal(message,public,seed)
        assert encrypted == ref.encipher_elgamal(message,public,seed)
        assert crypto.decipher_elgamal(encrypted,private) == message
    other = crypto.dh_private_key(bits,seed+100)
    assert crypto.dh_shared_key(public,other[2]) == crypto.dh_shared_key(crypto.dh_public_key(other),a)
    assert crypto.dh_shared_key(public,other[2]) == ref.dh_shared_key(public,other[2])


@pytest.mark.parametrize('seed',[[0,1,5],[5],[100,0,1,7]])
def test_seed_sequences_are_copied_and_match_reference(seed):
    original = list(seed)
    assert crypto.elgamal_private_key(5,seed) == ref.elgamal_private_key(5,seed)
    assert seed == original
    assert crypto.encipher_bg(42,77,seed) == ref.encipher_bg(42,77,seed)
    assert seed == original
    with pytest.raises(ValueError):crypto.elgamal_private_key(5,[0,1,1000])


@pytest.mark.parametrize('p,q',[(3,5),(17,19),(67,79),(131,137)])
@pytest.mark.parametrize('seed',[0,1,7,51])
def test_goldwasser_micali_reference_and_inverse(p,q,seed):
    private = crypto.gm_private_key(p,q)
    public = crypto.gm_public_key(p,q,seed=seed)
    assert private == ref.gm_private_key(p,q)
    assert public == ref.gm_public_key(p,q,seed=seed)
    a,n = public
    assert pow(a,(p-1)//2,p) == p-1 and pow(a,(q-1)//2,q) == q-1
    assert crypto.gm_public_key(p,q,a) == public
    assert crypto.gm_public_key(p,q,1) is False
    for message in (0,1,2,12345,2**200+3):
        encrypted = crypto.encipher_gm(message,public,seed)
        assert encrypted == ref.encipher_gm(message,public,seed)
        assert crypto.decipher_gm(encrypted,private) == message


@pytest.mark.parametrize('p,q',[(3,7),(7,11),(23,31),(67,79)])
@pytest.mark.parametrize('seed',[0,1,7,51])
def test_blum_goldwasser_reference_and_inverse(p,q,seed):
    private = crypto.bg_private_key(p,q)
    public = crypto.bg_public_key(p,q)
    assert private == ref.bg_private_key(p,q)
    assert public == ref.bg_public_key(p,q)
    for message in (0,1,2,12345,2**31+3):
        encrypted = crypto.encipher_bg(message,public,seed)
        assert encrypted == ref.encipher_bg(message,public,seed)
        assert crypto.decipher_bg(encrypted,private) == message


def test_blum_goldwasser_long_message_uses_exact_exponents():
    private = crypto.bg_private_key(67,79)
    message = (1<<8192)+(1<<4096)+12345
    encrypted = crypto.encipher_bg(message,crypto.bg_public_key(*private),seed=17)
    assert len(encrypted[0]) == message.bit_length()
    assert crypto.decipher_bg(encrypted,private) == message


@pytest.mark.parametrize('prime',[2,3,5,7,11])
@pytest.mark.parametrize('order',[1,2,3,5])
def test_lfsr_reference_and_recovered_recurrence(prime,order):
    random = Random(100*prime+order)
    coefficients = [random.randrange(prime) for i in range(order)]
    initial = [random.randrange(prime) for i in range(order)]
    if not any(initial):initial[0] = 1
    F,R = s.GF(prime),reference.FF(prime)
    key,fill = [F(v) for v in coefficients],[F(v) for v in initial]
    sequence = crypto.lfsr_sequence(key,fill,8*order)
    expected = ref.lfsr_sequence([R(v) for v in coefficients],[R(v) for v in initial],8*order)
    assert list(map(int,sequence)) == list(map(int,expected))
    assert list(map(int,key)) == coefficients and list(map(int,fill)) == initial
    polynomial = crypto.lfsr_connection_polynomial(sequence)
    # Check recurrence residuals independently of the recovery algorithm.
    poly = s.Poly(polynomial,s.Symbol('x'))
    degree = int(poly.degree())
    for position in range(max(order,degree),len(sequence)):
        assert sum(int(poly.nth(j))*int(sequence[position-j]) for j in range(degree+1))%prime == 0
    try:expected_poly = ref.lfsr_connection_polynomial(expected)
    except reference.polys.polyerrors.GeneratorsNeeded:
        assert polynomial == 1  # SymPy's constant-polynomial corner case.
    else:assert polynomial.to_sympy() == expected_poly


@pytest.mark.parametrize('prime',[2,3,5,7])
def test_lfsr_zero_and_native_result_types(prime):
    F = s.GF(prime)
    assert crypto.lfsr_connection_polynomial([F(0)]*20) == 1
    assert crypto.lfsr_sequence([F(1)],[F(2%prime)],0) == []
    assert crypto.lfsr_connection_polynomial([]) == 1
    result = crypto.lfsr_sequence([F(1)],[F(1)],5)
    assert all(v.parent() == F for v in result)


@pytest.mark.parametrize('lag',list(range(-15,16)))
def test_lfsr_cyclic_autocorrelation(lag):
    F = s.GF(2)
    seq = crypto.lfsr_sequence([F(1),F(0),F(1)],[F(1),F(0),F(0)],7)
    actual = crypto.lfsr_autocorrelation(seq,7,lag)
    assert isinstance(actual,s.Rational)
    assert actual == crypto.lfsr_autocorrelation(seq,7,lag%7)
    assert actual.to_sympy() == ref.lfsr_autocorrelation(list(map(int,seq)),7,lag%7)
    assert actual == (1 if lag%7 == 0 else s.Rational(-1,7))


@pytest.mark.parametrize('name',[
    'encipher_shift','decipher_shift','encipher_rot13','decipher_rot13',
    'encipher_affine','decipher_affine','encipher_atbash','decipher_atbash',
    'encipher_vigenere','decipher_vigenere','encipher_hill','decipher_hill',
    'encipher_bifid','decipher_bifid','encipher_bifid5','decipher_bifid5',
    'encipher_bifid6','decipher_bifid6'])
def test_empty_plaintext_is_empty(name):
    args = {
        'shift':(3,), 'affine':((3,5),), 'vigenere':('KEY',),
        'hill':(s.eye(2),), 'bifid':('KEY',),
    }
    arguments = next((v for stem,v in args.items() if stem in name),())
    assert getattr(crypto,name)('',*arguments) == ''


@pytest.mark.parametrize('call',[
    lambda:crypto.padded_key('ab','abba'),lambda:crypto.padded_key('Z','abc'),
    lambda:crypto.encipher_affine('abc',(2,1)),lambda:crypto.encipher_vigenere('abc','!'),
    lambda:crypto.encipher_hill('ABC',s.Matrix([[1,2,3]])),
    lambda:crypto.encipher_hill('A',s.eye(2),symbols='ABC'),
    lambda:crypto.encipher_bifid('abc','a','abcde'),lambda:crypto.bifid_square('abcde'),
    lambda:crypto.encipher_railfence('abc',0),lambda:crypto.decipher_railfence('abc',-1),
    lambda:crypto.rsa_public_key(3,5,7,totient='unknown'),
    lambda:crypto.rsa_private_key(3,5,7,index=1),lambda:crypto.rsa_public_key(0,5,7),
    lambda:crypto.elgamal_private_key(0),lambda:crypto.elgamal_private_key(5,seed='bad'),
    lambda:crypto.encipher_elgamal(-1,(37,2,3)),lambda:crypto.encipher_elgamal(37,(37,2,3)),
    lambda:crypto.dh_shared_key((37,2,3),1),lambda:crypto.dh_shared_key((37,2,3),37),
    lambda:crypto.gm_private_key(2,3),lambda:crypto.gm_private_key(3,3),
    lambda:crypto.gm_private_key(9,11),lambda:crypto.encipher_gm(-1,(2,15)),
    lambda:crypto.decipher_gm([0],(3,5)),lambda:crypto.bg_private_key(13,17),
    lambda:crypto.bg_private_key(7,7),lambda:crypto.encipher_bg(-1,77),
    lambda:crypto.decipher_bg(([2],1),(7,11)),lambda:crypto.encode_morse('ABC',sep=''),
    lambda:crypto.decode_morse('.-',sep=''),lambda:crypto.lfsr_sequence([],[],1),
    lambda:crypto.lfsr_sequence([s.GF(2)(1)],[s.GF(3)(1)],3),
    lambda:crypto.lfsr_autocorrelation([1,0],3,0),
])
def test_invalid_inputs_raise_explicitly(call):
    with pytest.raises((ValueError,TypeError)):call()


def test_package_exports_and_no_runtime_sympy():
    package = importlib.import_module('symox.crypto')
    assert len(crypto.__all__) == 57
    assert all(hasattr(crypto,n) and hasattr(package,n) for n in crypto.__all__)
    code = r'''
import sys
class BlockSymPy:
 def find_spec(self,name,*args):
  if name == 'sympy' or name.startswith('sympy.'):
   raise RuntimeError('runtime SymPy dependency: '+name)
sys.meta_path.insert(0,BlockSymPy())
import symox as s
from symox.crypto import crypto as c
assert c.encipher_shift('ABC',1) == 'BCD'
M = s.Matrix([[1,2],[3,5]])
assert c.decipher_hill(c.encipher_hill('ABCD',M),M) == 'ABCD'
pub,pri = c.rsa_public_key(61,53,17),c.rsa_private_key(61,53,17)
assert c.decipher_rsa(c.encipher_rsa(42,pub),pri) == 42
pri = c.elgamal_private_key(10,seed=3)
assert c.decipher_elgamal(c.encipher_elgamal(123,c.elgamal_public_key(pri),seed=7),pri) == 123
for keygen,public,enc,dec in [(c.gm_private_key,c.gm_public_key,c.encipher_gm,c.decipher_gm),
                             (c.bg_private_key,c.bg_public_key,c.encipher_bg,c.decipher_bg)]:
 key = keygen(67,79)
 ciphertext = enc(2**2048+123,public(*key),seed=3)
 assert dec(ciphertext,key) == 2**2048+123
F = s.GF(3)
seq = c.lfsr_sequence([F(1),F(1)],[F(0),F(1)],20)
assert c.lfsr_connection_polynomial(seq) == 1+2*s.Symbol('x')+2*s.Symbol('x')**2
assert 'sympy' not in sys.modules
'''
    result = subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
    assert result.returncode == 0,result.stderr


def test_seeded_crypto_does_not_change_python_random_state():
    import random
    state = random.getstate()
    for i in range(3):
        key = crypto.elgamal_private_key(100,seed=0)
        assert key[2] == Random(0).randrange(2,key[0])
        assert crypto.gm_public_key(17,19,seed=3) == crypto.gm_public_key(17,19,seed=3)
    assert random.getstate() == state
