# Algorithmic parity requirement

Symox must aim for mathematical capabilities at least equivalent to the
reference SymPy version, including accepted domains, solution completeness,
precision and error behavior. Speed must not be obtained by reducing those
capabilities. This is a user requirement, not an optional benchmark goal.

API-name coverage and selected performance examples do not establish behavioral
parity. Preserve exact arithmetic, precision and established fallback paths.
Distinguish an unsupported or undecidable computation from a proven absence of
solutions. Track remaining gaps explicitly before claiming equivalence.

See [the parity contract and current gaps](docs/algorithmic-parity.md).
