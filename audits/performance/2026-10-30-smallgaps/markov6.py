import sys, importlib, random
lib = importlib.import_module(sys.argv[1]); st = importlib.import_module(sys.argv[1]+'.stats')
rng = random.Random(7)
rows = []
for i in range(6):
    w = [rng.random() for _ in range(6)]; s = sum(w); rows.append([lib.Float(x/s) for x in w])
F6 = lib.Matrix(rows)
print('det', lib.srepr((lib.eye(6) - F6).det()))
try:
    r = st.DiscreteMarkovChain('X', list(range(6)), F6).fundamental_matrix()
    print('fund', r[0, 0], r[5, 5])
except Exception as e:
    print('fund', type(e).__name__, e)
