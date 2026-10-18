"""Notes on Meijer-G integration adapters.

The public integration engine performs native symbolic integration. The
meijerint entry points report only results that the native engine resolved;
no imported SymPy formula engine is used. Unevaluated cases remain explicit.
"""
