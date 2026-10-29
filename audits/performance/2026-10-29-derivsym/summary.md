| mode | case | before (ms) | after (ms) | before result | after result |
|---|---|---|---|---|---|
| fresh | import | 231.688 | 219.651 | - | - |
| fresh | diff_x_n | 0.167 | 0.264 | Derivative | Piecewise |
| fresh | diff_xsin_n | 0.117 | 0.320 | Derivative | Sum |
| fresh | diff_fg_n | 0.178 | 0.355 | Derivative | Sum |
| fresh | held_diff_x | 0.206 | 0.084 | Derivative | Derivative |
| fresh | held_subs_3 | 0.070 | 0.063 | Derivative | Derivative |
| fresh | leibniz_subs_doit | 0.507 | 1.573 | Mul | Add |
| fresh | rewrite_exp | 0.059 | 0.092 | Derivative | Derivative |
| fresh | sum_doit_derivs | 0.461 | 0.331 | Add | cos |
| fresh | eq_parity | 0.091 | 0.031 | Equality | BooleanFalse |
| fresh | diff_literal_5 | 0.260 | 0.295 | Mul | Mul |
| fresh | diff_literal_mixed | 0.222 | 0.219 | Mul | Mul |
| fresh | diff_sin_n | 0.058 | 0.057 | Derivative | Derivative |
| repeat | import | 230.037 | 221.550 | - | - |
| repeat | diff_x_n | 0.021 | 0.062 | Derivative | Piecewise |
| repeat | diff_xsin_n | 0.035 | 0.165 | Derivative | Sum |
| repeat | diff_fg_n | 0.041 | 0.078 | Derivative | Sum |
| repeat | held_diff_x | 0.029 | 0.029 | Derivative | Derivative |
| repeat | held_subs_3 | 0.021 | 0.020 | Derivative | Derivative |
| repeat | leibniz_subs_doit | 0.311 | 1.378 | Mul | Add |
| repeat | rewrite_exp | 0.018 | 0.029 | Derivative | Derivative |
| repeat | sum_doit_derivs | 0.307 | 0.252 | Add | cos |
| repeat | eq_parity | 0.045 | 0.005 | Equality | BooleanFalse |
| repeat | diff_literal_5 | 0.190 | 0.179 | Mul | Mul |
| repeat | diff_literal_mixed | 0.173 | 0.164 | Mul | Mul |
| repeat | diff_sin_n | 0.031 | 0.030 | Derivative | Derivative |
