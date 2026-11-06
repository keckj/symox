| case | SymPy ms | before ms / same | after ms / same |
|---|---:|---:|---:|
| RootSum cubic id | 1.9 | 2.6 / yes | 0.4 / yes |
| RootSum cubic Lambda 1/(x+1) | 3.6 | 3.7 / yes | 0.6 / yes |
| RootSum cubic Lambda x**2 | 1.9 | 2.9 / yes | 0.4 / yes |
| RootSum cubic Lambda x/(x**2+2) | 4.5 | 4.1 / yes | 0.6 / yes |
| RootSum cubic log | 0.6 | 0.6 / yes | 0.4 / yes |
| RootSum cubic y*log(z-y) | 0.7 | 0.6 / yes | 0.4 / yes |
| RootSum cubic 2+3*log | 0.7 | 0.6 / yes | 0.4 / yes |
| RootSum cubic exp | 0.6 | 0.6 / yes | 0.3 / yes |
| RootSum cubic Function f | 0.6 | ERR / no | 0.4 / yes |
| RootSum constant Lambda | 0.2 | 0.4 / yes | 0.3 / yes |
| RootSum quartic a param 1/(x-a) | ERR | ERR / yes | ERR / yes |
| RootSum quintic sum x**3 | 2.9 | 4.5 / yes | 0.5 / yes |
| RootSum quintic 1/x | 3.8 | 4.3 / yes | 0.4 / yes |
| RootSum product quad*cubic log | 1.4 | 0.8 / yes | 0.5 / yes |
| RootSum quadratic=True log | 0.6 | 0.7 / yes | 0.4 / yes |
| RootSum quadratic=False log | 0.5 | 0.6 / yes | 0.3 / yes |
| RootSum linear factor log | 1.2 | 0.8 / yes | 0.5 / yes |
| RootSum repeated factor log | 0.7 | 0.7 / yes | 0.5 / yes |
| RootSum auto=False | 0.6 | 0.6 / yes | 0.3 / yes |
| RootSum scaled 4*x**3+2 | 1.0 | 0.6 / yes | 0.4 / yes |
| RootSum integer basis 8 | 0.6 | 0.7 / yes | 0.4 / yes |
| RootSum integer basis rational | 3.7 | 3.6 / yes | 0.6 / yes |
| RootSum symbolic coeff a | ERR | ERR / yes | ERR / yes |
| RootSum symbolic coeff rational | ERR | ERR / yes | ERR / yes |
| RootSum param p | ERR | ERR / yes | ERR / yes |
| RootSum param rational p | ERR | ERR / yes | ERR / yes |
| RootSum sqrt2 coeff log | 1.0 | 0.6 / yes | 0.4 / yes |
| RootSum sqrt2 coeff rational | 2.2 | 2.8 / yes | 0.5 / yes |
| RootSum x gen keyword | 0.6 | 0.4 / no | 0.2 / no |
| RootSum fraction coeff | 0.6 | 0.6 / yes | 0.4 / yes |
| RootSum quadratic=True param | ERR | ERR / yes | ERR / yes |
| RootSum multivariate error | ERR | ERR / yes | ERR / yes |
| RootSum bad func error | ERR | ERR / yes | ERR / yes |
| RootSum args | 0.7 | 0.6 / yes | 0.4 / yes |
| RootSum hash eq | 1.3 | 1.2 / yes | 0.7 / yes |
| RootSum doit cubic | 2.5 | 2.7 / yes | 2.5 / yes |
| RootSum doit quartic | 1.0 | 3.7 / yes | 3.3 / yes |
| RootSum doit quintic stays | 1.8 | 1.2 / yes | 1.0 / yes |
| RootSum evalf cubic | 2.2 | 1.3 / yes | 1.0 / yes |
| RootSum evalf quintic | 6.4 | 1.6 / yes | 1.4 / yes |
| RootSum evalf symbolic | ERR | ERR / yes | ERR / yes |
| RootSum diff | 8.0 | 8.0 / yes | 0.9 / yes |
| RootSum subs | 0.7 | 1.7 / yes | 1.3 / yes |
| integrate 1/(x**3+x+1) real=False | 15.2 | 1.7 / yes | 1.5 / yes |
| integrate 1/(x**5+x+1) | 78.7 | 11.0 / yes | 11.0 / yes |
| integrate x**2/(x**5-x+1) | 81.9 | 9.8 / yes | 9.6 / yes |
| integrate 1/(x**3+p*x+1) | 661.0 | 32.0 / yes | 35.2 / yes |
| integrate 1/(x**2+sqrt(2)*x+p) | 142.1 | 13.8 / yes | 14.6 / yes |
| integrate 1/(x**3+sqrt(2)*x+1) | ERR | 62.1 / no | 61.2 / no |
| integrate 1/(x**2+sqrt(2)*p) | 115.3 | 12.6 / yes | 12.4 / yes |
| integrate x/(x**4+a) | 39.3 | 7.3 / yes | 4.9 / yes |
| integrate 1/(x**4+a*x+1) | 494.1 | 54.2 / yes | 51.8 / yes |
| ratint 1/(x**3+a) real=False | 29.2 | 3.7 / yes | 3.2 / yes |
| logpart extension=True param | 97.3 | 7.9 / yes | 8.0 / yes |
| logpart QQ<s2> 1/(x**3+s2*x+1) | 245.3 | 26.9 / yes | 27.1 / yes |
| logpart QQ<s2>[p] 1/(x**2+s2*p) | 22.1 | ERR / no | ERR / no |
| logpart QQ<s2>[p] x/(x**2+s2*x+p) | 53.4 | ERR / no | ERR / no |
| logpart QQ<s2>[p] 1/(x**2+p) | 13.2 | ERR / no | ERR / no |
| logpart QQ<s2>[p] 1/(x**3+s2*p) | 23.4 | ERR / no | ERR / no |
| logpart QQ<s2>[p] 1/((x-p)(x+s2)) | 68.4 | ERR / no | ERR / no |
| logpart QQ<s2>(t) 1/(x**2+s2*t) | 27.7 | ERR / no | ERR / no |
| logpart QQ<s2>(t) x/(x**2+s2*x+t) | 63.1 | ERR / no | ERR / no |
| logpart QQ<s2>(t) 1/(x**3+t*x+s2) | 309.1 | ERR / no | ERR / no |
| logpart QQ<s2>(t) 1/((x-t)(x+s2)) | 75.0 | ERR / no | ERR / no |
| logpart QQ<s2>(t) (x+1)/(x**2+t*x+s2) | 181.7 | ERR / no | ERR / no |
| apart full 1/(x**3+x+1) | 3.1 | 2.0 / yes | 1.7 / yes |
| apart full doit 1/(x**3+x+1) | 4.9 | 2.4 / yes | 2.1 / yes |
| apart full x/(x**4+1) | 3.1 | 1.9 / yes | 1.6 / yes |
| apart full 1/(x**2+1)**2 | 9.2 | 3.0 / yes | 2.5 / yes |
| apart full 1/((x-1)(x**3+2)) | 4.5 | 2.4 / yes | 2.1 / yes |
| apart full 1/(x**5-x+1) | 3.9 | 2.0 / yes | 1.8 / yes |
| apart full a/(x**3+a) | 8.1 | 3.5 / yes | 3.0 / yes |
| apart full 1/(x**2+sqrt(2)) | 9.8 | 2.8 / yes | 2.6 / yes |
| EX sqrt2*(1+sqrt2) | 0.0 | 0.4 / no | 0.4 / no |
| EX (x+1)**2 | 0.0 | 0.2 / yes | 0.2 / yes |
