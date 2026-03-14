# CEPII Gravity PPML With Migrant Stock

## Fit Statistics

| Metric | Value |
|---|---:|
| Flow unit | annualized_per_year |
| Prediction target | E[flow|X] = exp(Xb) |
| Observations in panel | 207,998 |
| Observations used in estimation | 207,998 |
| Zero-flow observations included | 62,221 |
| Parameters | 17 |
| Converged | True |
| Iterations | 9 |
| McFadden pseudo-R-squared | 0.775158 |
| RMSE (flow) | 3863.566083 |
| Deviance | 164885528.462 |
| Poisson pseudo-log-likelihood | -82657645.294 |
| Null pseudo-log-likelihood | -367625935.264 |
| AIC | 165315324.588 |
| BIC | 165315498.758 |
| Mean observed flow | 454.390272 |
| Mean fitted flow | 454.390272 |

## Coefficients

| Term | Estimate | Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |
|---|---:|---:|---:|---:|---:|---:|
| const | -6.248136 | 0.472794 | -13.215356 | 0 | -7.174794 | -5.321477 |
| ln_origin_population_total | 0.290805 | 0.011679 | 24.899014 | 0 | 0.267914 | 0.313696 |
| ln_destination_population_total | 0.360243 | 0.013120 | 27.458346 | 0 | 0.334529 | 0.385957 |
| ln_origin_gdp_pc_ppp_constant | -0.017801 | 0.018814 | -0.946204 | 0.344045 | -0.054675 | 0.019072 |
| ln_destination_gdp_pc_ppp_constant | 0.157140 | 0.021208 | 7.409558 | 1.26565e-13 | 0.115574 | 0.198706 |
| ln_distance_km | -0.268243 | 0.024918 | -10.764943 | 0 | -0.317082 | -0.219404 |
| contig | 0.460906 | 0.062006 | 7.433215 | 1.06137e-13 | 0.339376 | 0.582435 |
| comlang_off | 0.291355 | 0.043543 | 6.691225 | 2.21312e-11 | 0.206013 | 0.376698 |
| colony | 0.200855 | 0.101371 | 1.981395 | 0.047547 | 0.002172 | 0.399539 |
| col45 | 0.297357 | 0.098497 | 3.018955 | 0.00253648 | 0.104307 | 0.490406 |
| smctry | -0.050646 | 0.078745 | -0.643166 | 0.520116 | -0.204983 | 0.103691 |
| ln_migrant_stock_both_sexes_plus1 | 0.399159 | 0.011135 | 35.847789 | 0 | 0.377335 | 0.420983 |
| period_1995 | -0.094441 | 0.070332 | -1.342780 | 0.179343 | -0.232290 | 0.043408 |
| period_2000 | -0.033741 | 0.068157 | -0.495049 | 0.620566 | -0.167326 | 0.099844 |
| period_2005 | -0.024940 | 0.063944 | -0.390030 | 0.696514 | -0.150268 | 0.100388 |
| period_2010 | -0.044877 | 0.069379 | -0.646834 | 0.51774 | -0.180857 | 0.091103 |
| period_2015 | -0.121602 | 0.061597 | -1.974148 | 0.048365 | -0.242331 | -0.000874 |

## Average Annual Flow Comparison Over Full Sample

Fitted values below are PPML fitted means.

| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| Spain (ESP) | 363,742.138 | 300,948.603 | 135,494.855 | 184,813.888 | 228,247.283 | 116,134.715 |
| United States (USA) | 2,285,537.328 | 1,711,236.713 | 943,916.082 | 441,121.056 | 1,341,621.246 | 1,270,115.657 |
| Italy (ITA) | 342,668.390 | 471,164.991 | 191,491.188 | 272,545.723 | 151,177.202 | 198,619.267 |
| France (FRA) | 350,714.070 | 655,718.465 | 267,801.926 | 325,580.281 | 82,912.145 | 330,138.184 |
| Germany (DEU) | 602,520.377 | 516,165.067 | 349,678.215 | 453,996.774 | 252,842.163 | 62,168.293 |

## Spain By Period

Observed and fitted values below are annualized flows. Fitted values are PPML fitted means. Balance is `inflow - outflow`.

| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| 1990-1995 | 170,883.848 | 157,624.569 | 45,122.692 | 172,081.556 | 125,761.156 | -14,456.987 |
| 1995-2000 | 216,706.606 | 182,185.149 | 50,693.896 | 160,131.497 | 166,012.710 | 22,053.651 |
| 2000-2005 | 575,094.972 | 244,687.205 | 69,236.416 | 178,653.008 | 505,858.556 | 66,034.197 |
| 2005-2010 | 571,463.868 | 374,914.509 | 122,850.870 | 191,213.092 | 448,612.998 | 183,701.417 |
| 2010-2015 | 213,632.720 | 440,385.178 | 303,163.170 | 202,392.684 | -89,530.450 | 237,992.494 |
| 2015-2020 | 434,670.814 | 405,895.009 | 221,902.086 | 204,411.493 | 212,768.728 | 201,483.516 |