| case | before first (s) | after first (s) | ratio | before repeat (s) | after repeat (s) | ratio | result |
|---|---|---|---|---|---|---|---|
| res_witness_200_150 | 15.36 | 0.2238 | - | - | 0.2166 | - | before refused (ValueError: pseudo-division resource limit: PRS storage); after ok |
| res_60_45 | 0.08558 | 0.00458 | 18.68x | 0.08517 | 0.003988 | 21.36x | identical |
| res_dense_40_30 | 0.004086 | 0.001346 | 3.04x | 0.00356 | 0.0008867 | 4.01x | identical |
| res_small_zz | 0.0004972 | 0.0004074 | 1.22x | 0.0001816 | 0.0001521 | 1.19x | identical |
| res_small_qq | 0.000503 | 0.0004583 | 1.10x | 0.0001768 | 0.0001624 | 1.09x | identical |
| res_prs_small | 0.0004798 | 0.0004463 | 1.08x | 0.0001681 | 0.0001645 | 1.02x | identical |
| res_bivariate | 0.0006477 | 0.0005863 | 1.10x | 0.0002731 | 0.0002918 | 0.94x | identical |
| poly_res_25_20 | 0.003842 | 0.001133 | 3.39x | 0.00344 | 0.0007332 | 4.69x | identical |
| disc_dense_25 | 0.001438 | 0.0008217 | 1.75x | 0.001052 | 0.0004951 | 2.12x | identical |
| disc_small | 0.0004829 | 0.0003796 | 1.27x | 0.0001703 | 0.0001493 | 1.14x | identical |
| cyclotomic_20001 | 1.406e-05 | 0.06309 | - | - | 0.06401 | - | before refused (ValueError: polynomial structure resource limit: cyclotomic ); after ok |
| cyclotomic_15015 | 0.824 | 0.06237 | 13.21x | 0.8224 | 0.06278 | 13.10x | identical |
| cyclotomic_19999 | 0.2037 | 0.08047 | 2.53x | 0.2059 | 0.08059 | 2.55x | identical |
| factor_x1155_m1 | 0.08885 | 0.07471 | 1.19x | 0.08534 | 0.07471 | 1.14x | identical |
| decompogen_4000 | 0.06439 | 0.06442 | 1.00x | 0.06402 | 0.06503 | 0.98x | identical |
| decompogen_4100 | 0.0004154 | 0.06696 | - | - | 0.06636 | - | before refused (NotImplementedError: decompogen: functional decomposition is); after ok |
