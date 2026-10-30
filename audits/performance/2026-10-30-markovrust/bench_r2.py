import sys, importlib, time, random
lib = importlib.import_module(sys.argv[1]); st = importlib.import_module(sys.argv[1]+'.stats')
a, b, p, q = lib.symbols('a b p q', positive=True)
rng = random.Random(7)
def rnd(n):
    rows = []
    for i in range(n):
        w = [rng.random() for _ in range(n)]; s = sum(w); rows.append([lib.Float(x/s) for x in w])
    return lib.Matrix(rows)
cases = {
 'float6 stationary': lambda: st.DiscreteMarkovChain('X', list(range(6)), F6).stationary_distribution(),
 'float6 fundamental': lambda: st.DiscreteMarkovChain('X', list(range(6)), F6).fundamental_matrix(),
 'float4b fundamental': lambda: st.DiscreteMarkovChain('X', list(range(4)), FB).fundamental_matrix(),
 'sym3 stationary': lambda: st.DiscreteMarkovChain('X', [0,1,2], S3).stationary_distribution(),
 'ctmc sym3 limiting': lambda: st.ContinuousMarkovChain('C', [0,1,2], G3).limiting_distribution(),
 'ctmc float4 limiting': lambda: st.ContinuousMarkovChain('C', [0,1,2,3], GF4).limiting_distribution(),
 'ctmc P query': lambda: st.P(lib.Eq(C2(1), 1), lib.Eq(C2(0), 1)),
 'ctmc variance': lambda: st.variance(C2(lib.Rational(3,2)), lib.Eq(C2(0), 1)),
}
F6 = rnd(6)
FB = lib.Matrix([[0.25, 0.5, 0.25, 0], [0.125, 0.375, 0.25, 0.25], [0.5, 0, 0.25, 0.25], [0.0625, 0.4375, 0.25, 0.25]])
S3 = lib.Matrix([[1-p-q, p, q], [lib.Rational(1,3)]*3, [a, 0, 1-a]])
G3 = lib.Matrix([[-a, a, 0], [1, -2, 1], [0, b, -b]])
GF4 = lib.Matrix([[-1.5, 0.5, 0.5, 0.5], [0.25, -0.75, 0.25, 0.25], [0, 0.5, -1.0, 0.5], [0.375, 0.125, 0, -0.5]])
C2 = st.ContinuousMarkovChain('C', [0,1,2], lib.Matrix([[-1,1,0],[1,-1,0],[0,1,-1]]))
for k, f in cases.items():
    try: f()
    except Exception as e:
        print('%-22s ERR %s' % (k, type(e).__name__)); continue
    n = 0; t0 = time.perf_counter()
    while True:
        f(); n += 1
        if time.perf_counter() - t0 > 0.5 and n >= 3: break
    print('%-22s %10.3f ms' % (k, 1e3*(time.perf_counter()-t0)/n))
