# Cryptography

[README](../README.md)

The 57 public callables in SymPy 1.14.0's `crypto.crypto` are available under
`symox.crypto.crypto` and `symox.crypto`. This module provides classical ciphers,
finite-field feedback sequences and textbook public-key arithmetic. It does
not import SymPy at runtime.

```python
from symox.crypto import (
    encipher_vigenere, decipher_vigenere,
    rsa_public_key, rsa_private_key, encipher_rsa, decipher_rsa,
)

ciphertext = encipher_vigenere('ATTACK AT DAWN', 'LEMON')
assert ciphertext == 'LXFOPVEFRNHR'
assert decipher_vigenere(ciphertext, 'LEMON') == 'ATTACKATDAWN'

public = rsa_public_key(61, 53, 17)
private = rsa_private_key(61, 53, 17)
assert public == (3233, 17)
assert private == (3233, 2753)
assert decipher_rsa(encipher_rsa(42, public), private) == 42
```

## Classical ciphers and alphabets

Shift, ROT13, affine, Atbash, substitution, Vigenere, Hill, Bifid and rail-fence
ciphers provide encryption and inverse operations. With the default Latin
alphabet, shift/affine/Vigenere/Hill inputs are uppercased and filtered to A-Z.
An explicit alphabet preserves case and accepts its own characters, including
spaces or non-Latin symbols. Alphabets must contain distinct characters.
`AZ`, `check_and_join`, `cycle_list` and `padded_key` are also public.

`encipher_substitution` accepts character tables or substring dictionaries.
Dictionary keys are copied; longer replacements are applied before the
single-character translation table. The input mapping remains unchanged.

Hill uses native integral matrices and native modular matrix inversion.
Encryption pads incomplete column blocks with Q by default, or the supplied
single-character `pad`. Decryption retains padding and completes an incomplete
ciphertext block with alphabet position zero, following the reference API.
A decryptable key must be invertible modulo the alphabet size.

Bifid supports arbitrary square alphabets, plus the 5 by 5 and 6 by 6 helpers
and native symbolic square matrices. The 5 by 5 alphabet excludes J; J is
filtered out rather than substituted with I. The 6 by 6 alphabet appends 0-9.
The general default is Python's 100-character `string.printable` alphabet.
Rail-fence operations retain the original message characters and whitespace.

Empty plaintext returns an empty ciphertext. Unlike SymPy 1.14's historical
normalization corner case, an empty string is not expanded to the Latin
alphabet by cipher preparation. Vigenere rejects an empty filtered key.
Standalone `AZ()` and `AZ('')` still return the reference Latin alphabet.
Atbash retains the reference affine key `(25,25)`, including for custom
alphabets; it is not redefined as reversal of every supplied alphabet.

## RSA and Kid RSA

RSA accepts individual primes, several primes, repeated primes or composite
factors. Composite arguments are factored with symox; prime inputs avoid that
factorization. Key generation supports `totient='Euler'` and
`totient='Carmichael'`. `index` selects an exponent representative and requires
Carmichael's convention. Non-coprime exponents return `False`.

Repeated prime factors issue `NonInvertibleCipherWarning` unless
`multipower=True`. In that case encryption/decryption is guaranteed on units
modulo the modulus; it need not be bijective on every residue. The warning
class exposes `fullMessage` and `warn` as in the reference interface.

Encryption and decryption use exact modular powers over arbitrary-size Python
integers. Optional pairwise-coprime `factors` whose product is the modulus
select CRT recombination; unsuitable factor lists fall back to the ordinary
modular power. Pairwise-coprime prime-power factors are supported as well.
Kid RSA implements the elementary key construction and modular multiplication.
These integer interfaces do not add message padding, byte encoding or a key
container format.

## ElGamal, Diffie-Hellman and probabilistic bit ciphers

ElGamal and Diffie-Hellman generate `p=nextprime(2**digit)` and a primitive root
using symox number theory. The parameter is a bit exponent, so p has
`digit+1` bits. Keys are integer triples. ElGamal messages lie in `0 <= m < p`.
`dh_shared_key(received_public_key, local_private_exponent)` produces the shared
integer from the other participant's public triple.

Goldwasser-Micali validates two distinct odd primes, chooses a common quadratic
nonresidue and encrypts the message bits using independent random squares.
An explicitly supplied invalid nonresidue returns `False`. Blum-Goldwasser
requires distinct primes congruent to 3 modulo 4, generates one parity bit per
squaring and returns the encrypted bits plus the final squaring state. Its
inverse uses exact integral exponents and CRT; long messages never pass
through floating-point exponentiation. Integer zero has an empty bit list.
These integer APIs do not retain leading zero bits or a fixed payload width.

The `seed` parameter accepts an integer or a copied sequence of candidate
integers. Out-of-range sequence entries are skipped; exhausting a sequence
raises `ValueError`. Integer seeds use a local `random.Random` instance and
are reproducible without changing Python's global random state or depending
on number-theory cache activity. This differs from SymPy's shared generator
side effects on a cold group-parameter calculation. The random generators
and raw protocols serve the same educational scope as SymPy's crypto module;
production padding, authentication and secure key lifecycle APIs are outside
this module.

## Morse and feedback sequences

Morse encoding supports the reference letters, digits and punctuation, custom
mappings and a nonempty separator. A doubled separator marks word boundaries
and trailing whitespace. Unsupported characters are filtered; text is not
automatically uppercased. Decoding an unknown code raises `KeyError`. Custom
mappings are not mutated, and decoding an empty string returns an empty string.

`lfsr_sequence` operates on native prime-field elements and leaves the key and
fill lists unchanged. Lists must have the same positive length and consistent
prime characteristic. Extension fields and composite-modulus rings are not
accepted. Returned elements retain the native field parent.

`lfsr_connection_polynomial` uses Berlekamp-Massey to recover a connection
polynomial with native symbolic coefficients, reduced to `0,...,p-1`.
The empty and all-zero sequences return one. `lfsr_autocorrelation` returns a
native exact rational for the supplied period, normalizing arbitrary positive
or negative lags modulo that period.

## Validation and coverage

The regression tests cover known vectors, differential results against SymPy,
seeded randomness, native result types, exhaustive RSA residue domains,
feedback recurrence residuals, long bit messages, invalid inputs and execution
with SymPy imports blocked. The 51 upstream SymPy crypto test groups have also
been exercised with the corresponding symox functions and native constructors.

The README percentage counts importable API names. It does not certify every
signature edge case or imply that arbitrary SymPy internal crypto helpers are
implemented. Invalid sizes, alphabet/key shapes and arithmetic domains raise
explicit errors rather than relying on disabled assertions or incidental
index errors from the reference implementation.
