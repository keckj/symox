from fractions import Fraction as F
from pathlib import Path
import json
namespace={'F':F}
script=(Path(__file__).parent/'audit.py').read_text()
exec(script[script.index('def quantum'):script.index('def trace')],namespace)
R=namespace['rounded'];Q=namespace['quantum'];encode=namespace.get('encode',lambda q:[q.numerator,q.denominator])
results=[]
for rec in json.loads((Path(__file__).parent/'complex-rounding-audit-v2.json').read_text()):
 p=rec['precision'];e=rec['first_divergence']['native'];r,i=[F(*q) for q in e['b']]
 norm_exact=r*r+i*i
 norm=R(R(r*r,p)+R(i*i,p),p);inv=R(1/norm,p)
 simulation=[R(r*inv,p),R(-i*inv,p)]
 assert simulation==[F(*q) for q in e['actual']]
 results.append(dict(precision=p,label=e['label'],norm_exact=encode(norm_exact),norm_rounded=encode(norm),
                     inverse_of_rounded_norm=encode(inv),source_algorithm_simulation=[encode(q) for q in simulation],
                     actual_matches_p_norm_reciprocal_simulation=True,
                     exact_component_errors=[encode(F(*a)-F(*q)) for a,q in zip(e['actual'],e['rounded'])],
                     p_plus_4_double_round_alone_explains=[R(R(F(*q),p+4),p)==F(*a) for q,a in zip(e['exact'],e['actual'])]))
print(json.dumps(results,indent=2))
