# CC140/211 exact-rounding diagnosis (no source changes)

There is a genuine pre-existing native precision defect. It is not explained by presentation precision or a correct alternative rounding convention. The first faulty complex primitive is reciprocal, before later RREF differences propagate.

The frozen inexact_linear package, main and every received stage remain unchanged. The diagnostic runs against the private DEBUG `acd1ed8ff7d63c38c88bb2a56bd9957adf06ee014d0dc7f04bad252d3d7e3801` built from immutable r1b plus the frozen inexact dispatch patch. Every operation is the established DomainElement scalar operation; no new mathematical kernel was measured.

## Independent certificate

`audit.py` extracts direct native `.as_expr()` and reference `._mpc_` mantissa/exponent payloads, and separately verifies they equal both domains' `to_sympy` payloads. Every reference input is constructed from exactly the same dyadic numerators and denominators. There is no string-to-number conversion and no tolerance.

For each scalar operation in the original RREF loop, the oracle calculates the complex rational result with Python Fraction. It finds the binary exponent by integer bit lengths and an exact comparison, scales the absolute value to p bits, divides the two integers, then applies nearest/ties-to-even using the exact remainder. The returned dyad is thus the correctly rounded component. Every component of every reciprocal, multiply and subtract is checked against that operation's own exact stored inputs (even after the two trajectories diverge).

`complex-rounding-audit-v2.json` records all 20 scalar events per precision, operands, exact complex rationals, independently rounded components, actual dyads and error in ulps. `classification.json` proves the native source algorithm reproduces the bad reciprocal exactly.

## First errors

- **p140**, `pivot1.reciprocal`: inputs are still identical to the reference. The imaginary component should be `-940993504466598865261755653150685443041553 * 2**(-140)`. Native returns `-470496752233299432630877826575342721520777 * 2**(-139)`: error `-2**(-140)`, exactly **-1 ulp**. The reference is correctly rounded for both components.
- **p211**, `pivot0.reciprocal`: initial inputs are identical to the reference. The real component should be `1307450619184860710792925795146412643264894257215424150334213737 * 2**(-208)`. Native returns `2614901238369721421585851590292825286529788514430848300668427473 * 2**(-209)`: error `-2**(-209)`, exactly **-1 ulp**. The reference is correctly rounded for both components.

At p211 there is also a later native reciprocal real error of +1 ulp against its own local exact operands. At p140 the first reciprocal error listed above is the only native primitive in this witness failing correctly rounded local arithmetic. The full trace records later propagated differences separately.

## Source cause

`poly/domains.rs::DomainElement::binary` routes plain inexact division through `ops::div` and then `domain.convert`. Literal complex reciprocal is simplified by `ops.rs::pure_complex_power` / `gaussian_reciprocal` as:

1. `norm = r.mul(r).add(i.mul(i))` at the Float inputs' p bits;
2. `norm.inv()` at p bits;
3. evaluate `(r-i*I) * inverse_norm` and convert to the domain.

This prematurely rounds the norm and its reciprocal, then rounds the component products. Exact integer simulation of these p-bit steps produces exactly the incorrect native outputs in both witnesses. Merely removing eval_mp's p+4 presentation-style double rounding would not fix these errors: direct rounding exact reciprocal through p+4 does not reproduce the defective components. The primary defect is expression-level complex reciprocal simplification being used as field arithmetic.

The existing general `MpComplex::div` uses exact products but p+10 intermediate sums/divisions (matching mpmath's algorithm), which is better but is not itself a proof of correctly rounded complex division for all dyades. A correction should certify component rounding and should preserve existing fallback for representations it cannot certify, rather than change global symbolic reciprocal presentation.

## Reference also has a separate defect

On the reference's p140 trajectory, `pivot1.eliminate[0,1].product` has an imaginary component **+1 ulp** from its own correctly rounded exact product. Native's corresponding local product is correctly rounded. `mpmath.libmp.mpc_mul` reproduces the erroneous reference output directly; `libmp.from_rational` of the exact product agrees with the independent integer oracle. All operands and component numerators are in the trace. Therefore copying every reference bit is not a valid accuracy repair.

The reference p211 trajectory has no incorrectly rounded primitive in this witness. The p140 multiplication issue does not excuse the earlier native reciprocal defect and will not be adopted as a target convention.

## Proposed isolated repair scope

Route stock plain inexact DomainElement finite literal complex division to a native operation on the stored components at domain p, with independently certified quotient rounding. Keep domain selected precision authoritative; leave composite field/parameter algorithms, symbolic expression simplification, protocol wrappers, poles and existing fallbacks intact. Audit multiplication/subtraction once the early reciprocal error is removed, since they then see different legitimate inputs. Additional primitives should change only when another genuine native accuracy defect is proved.

No speed gain is claimed for repairing an old wrong value. No source patch has been made at this report stage.
