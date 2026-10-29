| case | first before (ms) | first after (ms) | repeated before (ms) | repeated after (ms) |
|---|---|---|---|---|
| convert_array_to_matrix | 0.44 | 0.44 | 0.407 | 0.414 |
| Trace(A*B).diff(A) | 0.50 | 0.48 | 0.394 | 0.369 |
| sqrt(x.T*y).diff(x) | 0.27 | 0.26 | 0.187 | 0.183 |
| convert_indexed_to_array | 0.65 | 0.64 | 0.436 | 0.446 |
| (x.T*x).diff(x) | 0.18 | 0.18 | 0.139 | 0.137 |
| determinant.det(A) | n/a | 0.04 | n/a | 0.001 |
| (A*A).diff(A) | 0.63 | 1.50 | 0.554 | 0.543 |
| A.diff(A) | 1.48 | 1.46 | 0.214 | 0.210 |
| list(Flatten(20x30 array)) | n/a | 0.59 | n/a | 0.540 |
| convert_matrix_to_array | 1.21 | 0.35 | 0.267 | 0.267 |
