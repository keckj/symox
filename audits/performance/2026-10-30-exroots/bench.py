"""exroots bench: python bench.py sympy|symox CASE_INDEX -> one json line.
First call in a fresh process, then the median of 5 repeated calls."""
import sys, time, json
if sys.argv[1] == 'sympy':
    from sympy import *
else:
    from symox import *
x = symbols('x')
CASES = [
    ('roots', '(3125/4 + 3125*sqrt(5)/4)*x**4 + (-125*sqrt(5)/4 - 125/2)*x**2 + 5*sqrt(5)/64 + 75/64'),
    ('roots', 'sqrt(3)*x**2 - 2*x + sqrt(2)'),
    ('roots', 'x**4 - sqrt(3)*x**2 - 1'),
    ('roots', 'sqrt(2)*x**4 + 3*x**2 + sqrt(5)'),
    ('simplify', '(2 + sqrt(5))/(50*(1 + sqrt(5)))'),
    ('simplify', '(1 + sqrt(3))/(2 - sqrt(3))'),
    ('simplify', '(2 + sqrt(5))/(50*(1 + sqrt(5))) + 1/(25*(1 + sqrt(5)))'),
    ('simplify', 'sqrt(-1/(25*(1 + sqrt(5))) + (2 + sqrt(5))/(50*(1 + sqrt(5))))'),
    ('integrate', '1/(x**3-2)'),
    ('integrate', '1/(x**3+x+1)'),
    ('integrate', '1/(x**3+2*x-1)'),
    ('integrate', '1/(x**5+1)'),
]
kind, src = CASES[int(sys.argv[2])]
f = sympify(src, locals={'x': x})
def call():
    if kind == 'roots':
        return roots(f, x, multiple=True)
    if kind == 'simplify':
        return simplify(f)
    return integrate(f, x)
t0 = time.perf_counter(); r = call(); first = time.perf_counter() - t0
reps = []
for _ in range(5):
    t0 = time.perf_counter(); call(); reps.append(time.perf_counter() - t0)
reps.sort()
print(json.dumps({'case': '%s(%s)' % (kind, src), 'first_ms': round(first * 1e3, 2),
                  'repeat_ms': round(reps[2] * 1e3, 2), 'len': len(str(r))}))
