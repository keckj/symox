"""Timings for the crypto / stats / physics scenes.

    python bench_modules.py LIB CASE  -> one JSON line (worker, fresh process)
    python bench_modules.py [CASE...] -> both libs, best of 3 fresh processes, writes cases4.json

Each worker runs the case setup, then a smaller call of the same kind (warm-up),
then times the case once. Results are compared for mathematical equality in the
parent (SymPy simplify of the difference when the printed forms differ).
"""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "data"
PY = sys.executable
TIMEOUT = 120

PRELUDE = """
sp = __import__(LIB)
from importlib import import_module
crypto = import_module(LIB + '.crypto.crypto')
stats = import_module(LIB + '.stats')
x, y, z, t, r, s, p, mu = sp.symbols('x y z t r s p mu')
sigma, k, theta, a, b, lam = sp.symbols('sigma k theta a b lam', positive=True)
"""

L = "{L}"
# name: (setup, warm-up, timed expression)
CASES = {
    # crypto & number theory
    "f7": ("", "sp.factorint(2**64 + 1)", "sp.factorint(2**128 + 1)"),
    "berlekamp": (f"from {L}.crypto.crypto import lfsr_sequence, lfsr_connection_polynomial; F = sp.FF(2); "
                  "key = [F(v) for v in (1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1)]; fill = [F(1)]*11; seq = lfsr_sequence(key, fill, 200)",
                  "lfsr_connection_polynomial(seq[:24])", "lfsr_connection_polynomial(seq)"),
    "hill": (f"from {L}.crypto.crypto import encipher_hill; msg = 'MEETMEATTHEUSUALPLACE' * 100; "
             "K = sp.Matrix([[6, 24, 1], [13, 16, 10], [20, 17, 15]])",
             "encipher_hill('HELLOWORLD', K)", "encipher_hill(msg, K)"),
    "rsa": (f"from {L}.crypto.crypto import rsa_private_key, rsa_public_key, encipher_rsa, decipher_rsa; "
            "P = sp.nextprime(10**60); Q = sp.nextprime(10**61); E = 65537",
            "decipher_rsa(encipher_rsa(5, rsa_public_key(61, 53, 17)), rsa_private_key(61, 53, 17))",
            "[decipher_rsa(encipher_rsa(m, rsa_public_key(P, Q, E)), rsa_private_key(P, Q, E)) for m in range(1, 200)][-1]"),
    "vigenere": (f"from {L}.crypto.crypto import encipher_vigenere; msg = 'ATTACKATDAWN' * 2000",
                 "encipher_vigenere('HELLO', 'KEY')", "encipher_vigenere(msg, 'LEMONADE')"),
    "lfsr": (f"from {L}.crypto.crypto import lfsr_sequence; F = sp.FF(2); key = [F(v) for v in (1, 0, 0, 1, 0, 1, 1)]; fill = [F(1)]*7",
             "lfsr_sequence(key, fill, 10)", "sum(int(v) for v in lfsr_sequence(key, fill, 3000))"),
    # stats
    "normal_moment": ("X = stats.Normal('X', mu, sigma)", "stats.E(stats.Normal('W', 0, 1)**2)", "sp.simplify(stats.E(X**4))"),
    "normal_var": ("X = stats.Normal('X', mu, sigma)", "stats.E(stats.Normal('W', 0, 1)**2)", "sp.expand(stats.variance(X**2))"),
    "binom_var": ("B = stats.Binomial('B', 12, p)", "stats.E(stats.Binomial('W', 2, sp.Rational(1, 2)))", "sp.factor(stats.variance(B))"),
    "binom_skew": ("B = stats.Binomial('B', 10, p)", "stats.E(stats.Binomial('W', 2, sp.Rational(1, 2)))", "sp.factor(stats.skewness(B))"),
    "dice": ("D1, D2, D3 = stats.Die('D1', 6), stats.Die('D2', 6), stats.Die('D3', 6)", "stats.E(stats.Die('W', 4))",
             "stats.P(D1 + D2 + D3 > 12)"),
    "dice_E": ("D1, D2 = stats.Die('D1', 6), stats.Die('D2', 6)", "stats.E(stats.Die('W', 4))", "stats.E((D1 + D2)**3)"),
    "exp_P": ("X = stats.Exponential('X', lam)", "stats.E(stats.Normal('W', 0, 1)**2)", "stats.P(X > 2)"),
    "gamma_mgf": ("G = stats.Gamma('G', k, theta)", "stats.E(stats.Normal('W', 0, 1)**2)",
                  "sp.simplify(stats.moment_generating_function(G)(t))"),
    "sum_normals": ("X = stats.Normal('X', 0, 1); Y = stats.Normal('Y', 0, 1)", "stats.E(stats.Normal('W', 0, 1)**2)",
                    "sp.simplify(stats.density(X + Y)(z))"),
    "beta_E": ("Bt = stats.Beta('Bt', a, b)", "stats.E(stats.Normal('W', 0, 1)**2)", "sp.simplify(stats.E(Bt))"),
    "normal_entropy": ("X = stats.Normal('X', mu, sigma)", "stats.E(stats.Normal('W', 0, 1)**2)", "sp.simplify(stats.entropy(X))"),
    "poisson_var": ("Pn = stats.Poisson('Pn', lam)", "stats.E(stats.Binomial('W', 2, sp.Rational(1, 2)))", "sp.simplify(stats.variance(Pn))"),
    "uniform_E": ("U = stats.Uniform('U', a, a + b)", "stats.E(stats.Normal('W', 0, 1)**2)", "sp.factor(stats.E(U**3))"),
    # physics
    "beam": (f"from {L}.physics.continuum_mechanics.beam import Beam; E_, I_ = sp.symbols('E I'); R1, R2 = sp.symbols('R1 R2'); "
             "bm = Beam(10, E_, I_); bm.apply_load(R1, 0, -1); bm.apply_load(R2, 10, -1); bm.apply_load(-8, 4, -1); "
             "bm.apply_load(-3, 0, 0, end=10); bm.bc_deflection = [(0, 0), (10, 0)]; "
             "w = Beam(2, E_, I_); w.apply_load(R1, 0, -1); w.apply_load(R2, 2, -1); w.apply_load(-1, 1, -1); w.bc_deflection = [(0, 0), (2, 0)]",
             "(w.solve_for_reaction_loads(R1, R2), w.deflection())",
             "(bm.solve_for_reaction_loads(R1, R2), sp.factor(bm.deflection()*E_*I_))[1]"),
    "hydrogen": (f"from {L}.physics.hydrogen import R_nl", "sp.integrate(R_nl(1, 0, r, 1)**2 * r**2, (r, 0, sp.oo))",
                 "sp.integrate(R_nl(5, 2, r, 1)**2 * r**2, (r, 0, sp.oo))"),
    "pauli": (f"from {L}.physics.paulialgebra import Pauli, evaluate_pauli_product; s1, s2, s3 = Pauli(1), Pauli(2), Pauli(3)",
              "evaluate_pauli_product(s1*s2)", "sp.factor(evaluate_pauli_product(sp.expand((s1*x + s2*y + s3*z)**6)))"),
    "lagrange2": (f"from {L}.physics.mechanics import dynamicsymbols, ReferenceFrame, Point, Particle, LagrangesMethod, Lagrangian; "
                  "q1, q2 = dynamicsymbols('q1 q2'); m_, l_, g_ = sp.symbols('m l g'); N = ReferenceFrame('N'); "
                  "A = N.orientnew('A', 'Axis', [q1, N.z]); B = N.orientnew('B', 'Axis', [q2, N.z]); O = Point('O'); O.set_vel(N, 0); "
                  "P1 = O.locatenew('P1', l_*A.x); P1.v2pt_theory(O, N, A); P2 = P1.locatenew('P2', l_*B.x); P2.v2pt_theory(P1, N, B); "
                  "Pa1 = Particle('Pa1', P1, m_); Pa2 = Particle('Pa2', P2, m_); "
                  "Pa1.potential_energy = -m_*g_*l_*sp.cos(q1); Pa2.potential_energy = -m_*g_*l_*(sp.cos(q1) + sp.cos(q2))",
                  "sp.simplify(sp.cos(q1)**2 + sp.sin(q1)**2)",
                  "sp.simplify(LagrangesMethod(Lagrangian(N, Pa1, Pa2), [q1, q2]).form_lagranges_equations())"),
    "tf": (f"from {L}.physics.control.lti import TransferFunction, Series, Feedback; "
           "G1 = TransferFunction(s + 3, s**3 + 4*s**2 + 5*s + 2, s); G2 = TransferFunction(s**2 + 1, s**2 + 2*s + 7, s)",
           "Feedback(TransferFunction(1, s, s), TransferFunction(1, 1, s)).doit(cancel=True)",
           "Feedback(Series(G1, G2), TransferFunction(1, 1, s)).doit(cancel=True, expand=True).to_expr()"),
    "grover": (f"from {L}.physics.quantum.qubit import IntQubit; from {L}.physics.quantum.grover import apply_grover; "
               "f = lambda q: q == IntQubit(5, nqubits=4); f2 = lambda q: q == IntQubit(1, nqubits=2)",
               "apply_grover(f2, 2)", "apply_grover(f, 4)"),
    "wigner9j": (f"from {L}.physics.wigner import wigner_9j", "wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0)", "wigner_9j(4, 4, 4, 4, 4, 4, 4, 4, 4)"),
    "cg": (f"from {L}.physics.wigner import clebsch_gordan", "clebsch_gordan(1, 1, 2, 0, 0, 0)", "clebsch_gordan(15, 12, 10, 3, -5, -2)"),
}


def worker(lib, case):
    setup, warm, expr = CASES[case]
    ns = {"LIB": lib}
    exec(PRELUDE.replace("LIB", repr(lib)), ns)
    exec(setup.replace("{L}", lib), ns)
    eval(warm, ns)
    def expired(signum, frame):
        raise TimeoutError('timed operation exceeded 120 seconds')
    signal.signal(signal.SIGALRM, expired)
    t0 = time.perf_counter()
    signal.setitimer(signal.ITIMER_REAL, TIMEOUT)
    try:
        r = eval(expr, ns)
    except TimeoutError:
        return {"timeout": TIMEOUT}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    dt = time.perf_counter() - t0
    return {"time": dt, "result": str(r)}


def run(lib, case):
    try:
        out = subprocess.run([PY, str(Path(__file__).resolve()), lib, case], capture_output=True, text=True,
                             timeout=TIMEOUT+30, env=os.environ.copy())
    except subprocess.TimeoutExpired:
        return {"error": "worker exceeded process budget; timed operation bound unavailable"}
    if out.returncode:
        return {"error": (out.stderr.strip().splitlines() or ["?"])[-1][:300]}
    return json.loads(out.stdout)


def best(lib, case, repeat=3):
    best_r = None
    for _ in range(repeat):
        r = run(lib, case)
        if "time" not in r:
            return r
        if best_r is None or r["time"] < best_r["time"]:
            best_r = r
        if r["time"] > 5:
            break
    return best_r


def same(a, b):
    if a == b:
        return True
    import sympy
    try:
        d = sympy.simplify(sympy.sympify(a) - sympy.sympify(b))
        return d == 0
    except Exception:
        return False


def main():
    out = DATA / "cases4.json"
    results = json.loads(out.read_text()) if out.exists() else {}
    for case in sys.argv[1:] or list(CASES):
        a, b = best("symox", case), best("sympy", case)
        eq = "result" in a and "result" in b and same(a["result"], b["result"])
        results[case] = {"symox": a, "sympy": b, "same": eq}
        fa = f"{a['time']*1e3:9.2f}ms" if "time" in a else list(a)[0].upper()
        fb = f"{b['time']*1e3:9.2f}ms" if "time" in b else list(b)[0].upper()
        ratio = f"x{b['time']/a['time']:.0f}" if "time" in a and "time" in b else ""
        print(f"{case:15s} {fa:>11s} {fb:>11s} {ratio:>8s} {'SAME' if eq else 'DIFF'}", flush=True)
        if not eq:
            print("    sx:", str(a.get("result", a))[:200])
            print("    sp:", str(b.get("result", b))[:200])
        out.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] in ("symox", "sympy"):
        print(json.dumps(worker(sys.argv[1], sys.argv[2])))
    else:
        main()
