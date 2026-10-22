"""Shared equilibrium solve of the native Arch and Cable classes."""
import symox as S


def _solve(equations, unknowns):
    solutions = S.solve(equations, unknowns, dict=True)
    if not solutions:raise ValueError('boundary conditions or equilibrium have no solution')
    solution = solutions[0]
    if any(v not in solution or solution[v].has(*unknowns) for v in unknowns):raise ValueError('insufficient independent boundary conditions')
    return solution
