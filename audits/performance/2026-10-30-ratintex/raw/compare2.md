| case | SymPy ms | after ms / same | after2 ms / same |
|---|---:|---:|---:|
| logpart QQ<s2> 1/(x**2+s2) | 15.8 | 1.5 / yes | 1.6 / yes |
| logpart QQ<s2> 1/((x-s2)(x+1)) | 42.3 | 3.6 / yes | 3.9 / yes |
| logpart QQ<s2> x/(x**3+s2) | 11.5 | 1.1 / yes | 1.2 / yes |
| logpart QQ<s2> 1/(x**3-s2*x+1) | 233.7 | 22.3 / yes | 24.4 / yes |
| logpart QQ<s2> (x+s2)/(x**4+s2*x**2+1) | 356.5 | 30.1 / yes | 32.5 / yes |
| logpart QQ<s3> 1/(x**2-s3*x+1) | 11.7 | 1.0 / yes | 1.1 / yes |
| logpart QQ<s2> extension=True | 14.1 | 1.4 / yes | 1.5 / yes |
| logpart QQ<s2>/QQ<s3> 1/(x**2+s3) | 46.4 | missing | 1.7 / yes |
| logpart QQ<s2>/QQ<s3> (x+s2)/(x**2+s3*x+1) | 104.5 | missing | 8.7 / yes |
| logpart QQ<s2>/QQ<2**(1/3)> x/(x**3+2**(1/3)) | 73.9 | missing | 2.0 / yes |
| ratint tuple QQ<s2> | 91.2 | 8.0 / yes | 9.0 / yes |
| log_to_atan QQ/EX doctest | 0.9 | 0.4 / yes | 0.4 / yes |
| log_to_atan EX non-constant | 0.5 | 0.3 / yes | 0.3 / yes |
| log_to_atan EX deg2/deg1 | 50.6 | 4.6 / yes | 4.8 / yes |
| integrate 1/(x**4+2) | 97.7 | 11.1 / yes | 12.6 / yes |
| integrate 1/(x**2+sqrt(2)) | 97.8 | 0.2 / yes | 0.2 / yes |
| integrate 1/(x**3+sqrt(2)) | 87.6 | 10.4 / yes | 12.0 / yes |
| integrate 1/(x**2+p) | 40.6 | 7.0 / yes | 8.0 / yes |
| integrate 1/(x**5+1) | 766.5 | 32.3 / yes | 36.8 / yes |
| integrate x/(x**5+1) | 896.5 | 39.0 / yes | 44.1 / yes |
| integrate 1/(x**3-2) | 44.4 | 8.0 / yes | 9.1 / yes |
| integrate 1/(x**3+x+1) | 2189.8 | 322.8 / no | 336.2 / yes |
| integrate 1/(x**3+2*x-1) | 2122.6 | 291.5 / no | 311.7 / yes |
| integrate 1/((x**2+1)*(x**3+2)) | 1880.8 | 239.8 / yes | 252.1 / yes |
| ratint real=False x/(x**3-a*x+1) | 403.8 | 37.0 / yes | 41.4 / yes |
| ratint real=False 1/(x**4+x+a) | 420.2 | 37.4 / yes | 42.2 / yes |
| ratint real=False 1/(x**5+1) | 15.4 | 1.6 / yes | 1.8 / yes |
| ratint real=False 1/(x**2+2) | 9.1 | 1.5 / yes | 1.6 / yes |
| roots x**8-2*x**4+x**2-3 | 122.7 | 20.1 / no | 21.9 / no |
| roots (x**4-2*x**2-3)**2 | 2.2 | 0.1 / yes | 0.2 / yes |
| roots octic in x**2+1 | 118.1 | 21.4 / no | 23.0 / no |
| roots case 20 x**4+2*x**2+3 | 5.1 | 0.6 / yes | 0.7 / yes |
| roots x**6+x**2+1 | 25.6 | 11.0 / no | 12.2 / no |
| roots x**4+x**2+1 | 4.9 | 0.3 / yes | 0.3 / yes |
| roots x**4-x**2+3 | 5.7 | 0.5 / yes | 0.6 / yes |
| roots x**8+x**4+3 | 14.8 | 6.2 / yes | 6.7 / yes |
| roots x**4+a*x**2+1 | 10.5 | 0.4 / yes | 0.5 / yes |
| roots x**6-3*x**3+1 | 5.6 | 0.4 / yes | 0.4 / yes |
