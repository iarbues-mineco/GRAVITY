# CEPII Gravity OLS With Migrant Stock

## Fit Statistics

| Metric | Value |
|---|---:|
| Flow unit | annualized_per_year |
| Level prediction (median) | `exp(fitted_ln_flow)` |
| Level prediction (mean, smearing) | `exp(fitted_ln_flow) * smearing_factor` |
| Observations in panel | 207,998 |
| Observations used in estimation | 145,777 |
| Zero-flow observations excluded | 62,221 |
| Parameters | 17 |
| R-squared | 0.624530 |
| Adjusted R-squared | 0.624489 |
| RMSE (ln flow) | 2.526008 |
| Residual std. error | 2.526155 |
| Log-likelihood | -341931.436 |
| AIC | 683896.873 |
| BIC | 684065.000 |
| Smearing factor | 39.110108 |

## Coefficients

| Term | Estimate | Std. Error | t-stat | p-value | 95% CI Low | 95% CI High |
|---|---:|---:|---:|---:|---:|---:|
| const | -19.447124 | 0.148017 | -131.384095 | 0 | -19.737233 | -19.157016 |
| ln_origin_population_total | 0.406137 | 0.003388 | 119.875100 | 0 | 0.399497 | 0.412777 |
| ln_destination_population_total | 0.639867 | 0.003403 | 188.035214 | 0 | 0.633198 | 0.646537 |
| ln_origin_gdp_pc_ppp_constant | 0.421951 | 0.005733 | 73.596606 | 0 | 0.410714 | 0.433188 |
| ln_destination_gdp_pc_ppp_constant | 0.757438 | 0.006144 | 123.290593 | 0 | 0.745397 | 0.769479 |
| ln_distance_km | -1.043921 | 0.009200 | -113.465910 | 0 | -1.061953 | -1.025889 |
| contig | 0.152045 | 0.050994 | 2.981612 | 0.00286735 | 0.052098 | 0.251992 |
| comlang_off | 1.128466 | 0.018993 | 59.416370 | 0 | 1.091241 | 1.165691 |
| colony | 0.228676 | 0.088880 | 2.572855 | 0.0100863 | 0.054474 | 0.402879 |
| col45 | 1.164103 | 0.110698 | 10.516008 | 0 | 0.947139 | 1.381068 |
| smctry | 0.661877 | 0.066843 | 9.901930 | 0 | 0.530867 | 0.792888 |
| ln_migrant_stock_both_sexes_plus1 | 0.587321 | 0.002571 | 228.480417 | 0 | 0.582283 | 0.592360 |
| period_1995 | 0.006891 | 0.023562 | 0.292473 | 0.769925 | -0.039290 | 0.053072 |
| period_2000 | -0.119460 | 0.023404 | -5.104232 | 3.32141e-07 | -0.165331 | -0.073589 |
| period_2005 | -0.316986 | 0.023433 | -13.527263 | 0 | -0.362914 | -0.271058 |
| period_2010 | -0.499633 | 0.023305 | -21.438696 | 0 | -0.545310 | -0.453955 |
| period_2015 | -0.592778 | 0.023377 | -25.357216 | 0 | -0.638597 | -0.546960 |

## Average Annual Flow Comparison Over Full Sample

These averages are computed as total period flows over all periods divided by the total sample length in years.
Fitted values below use the median prediction `exp(Xb)`.

| Country | Observed Inflow | Fitted Inflow Median | Observed Outflow | Fitted Outflow Median | Observed Balance | Fitted Balance Median |
|---|---:|---:|---:|---:|---:|---:|
| Spain (ESP) | 363,742.138 | 487,338.902 | 135,494.855 | 369,908.963 | 228,247.283 | 117,429.939 |
| United States (USA) | 2,285,537.328 | 6,323,986.928 | 943,916.082 | 863,276.533 | 1,341,621.246 | 5,460,710.395 |
| Italy (ITA) | 342,668.390 | 1,130,642.415 | 191,491.188 | 1,016,215.338 | 151,177.202 | 114,427.077 |
| France (FRA) | 350,714.070 | 4,399,492.070 | 267,801.926 | 1,224,943.288 | 82,912.145 | 3,174,548.782 |
| Germany (DEU) | 602,520.377 | 2,333,685.675 | 349,678.215 | 2,091,152.657 | 252,842.163 | 242,533.018 |

## Spain By Period

Observed and fitted values below are annualized flows. Fitted values use the median prediction `exp(Xb)`. Balance is `inflow - outflow`.

| Period | Observed Inflow | Fitted Inflow Median | Observed Outflow | Fitted Outflow Median | Observed Balance | Fitted Balance Median |
|---|---:|---:|---:|---:|---:|---:|
| 1990-1995 | 170,883.848 | 251,582.020 | 45,122.692 | 357,315.272 | 125,761.156 | -105,733.252 |
| 1995-2000 | 216,706.606 | 327,116.019 | 50,693.896 | 383,814.280 | 166,012.710 | -56,698.261 |
| 2000-2005 | 575,094.972 | 428,228.326 | 69,236.416 | 408,932.365 | 505,858.556 | 19,295.961 |
| 2005-2010 | 571,463.868 | 624,882.316 | 122,850.870 | 376,010.205 | 448,612.998 | 248,872.111 |
| 2010-2015 | 213,632.720 | 671,331.808 | 303,163.170 | 336,553.287 | -89,530.450 | 334,778.521 |
| 2015-2020 | 434,670.814 | 620,892.924 | 221,902.086 | 356,828.367 | 212,768.728 | 264,064.556 |