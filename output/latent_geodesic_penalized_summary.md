# Penalized Latent Geodesic PPML

## Fit Statistics

| Metric | Value |
|---|---:|
| Flow unit | annualized_per_year |
| Prediction target | E[flow|X] = exp(Xb) |
| Countries | 191 |
| Free country coordinates | 189 |
| Anchors | USA, CHN |
| Observations used in estimation | 207,998 |
| Zero-flow observations included | 62,221 |
| Parameters | 395 |
| Coordinate iterations | 58 |
| Reference pseudo-log-likelihood | -83345993.326 |
| Final pseudo-log-likelihood | -63456008.256 |
| Log-likelihood improvement | 19889985.070 |
| McFadden pseudo-R-squared | 0.827390 |
| RMSE (flow) | 3276.068860 |
| Mean observed flow | 454.390272 |
| Mean fitted flow | 454.390272 |
| Mean displacement (km, free countries) | 1663.846 |
| Median displacement (km, free countries) | 1377.771 |
| Max displacement (km, free countries) | 7276.340 |
| Quadratic penalty weight | 2.000000 |
| Quadratic penalty value | 1821411.710 |
| Penalized objective | -65277419.966 |
| Penalized objective improvement | 18068573.359 |

## Coefficients

| Term | Estimate | Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |
|---|---:|---:|---:|---:|---:|---:|
| const | -8.604227 | 0.400294 | -21.494792 | 0 | -9.388788 | -7.819666 |
| ln_origin_population_total | 0.418127 | 0.010743 | 38.919478 | 0 | 0.397070 | 0.439183 |
| ln_destination_population_total | 0.496171 | 0.010833 | 45.800481 | 0 | 0.474938 | 0.517404 |
| ln_origin_gdp_pc_ppp_constant | 0.161958 | 0.014757 | 10.974700 | 0 | 0.133034 | 0.190882 |
| ln_destination_gdp_pc_ppp_constant | 0.412018 | 0.020086 | 20.513035 | 0 | 0.372651 | 0.451386 |
| contig | 0.219403 | 0.039167 | 5.601655 | 2.12314e-08 | 0.142636 | 0.296169 |
| comlang_off | 0.660448 | 0.043549 | 15.165691 | 0 | 0.575094 | 0.745802 |
| colony | 0.300735 | 0.092289 | 3.258626 | 0.00111953 | 0.119852 | 0.481618 |
| col45 | 0.301767 | 0.093121 | 3.240605 | 0.00119276 | 0.119254 | 0.484280 |
| smctry | 0.149780 | 0.072086 | 2.077791 | 0.0377287 | 0.008494 | 0.291067 |
| ln_migrant_stock_both_sexes_plus1 | 0.265329 | 0.008492 | 31.244448 | 0 | 0.248685 | 0.281973 |
| period_1995 | -0.088824 | 0.054486 | -1.630210 | 0.103057 | -0.195614 | 0.017967 |
| period_2000 | -0.093807 | 0.053040 | -1.768591 | 0.0769621 | -0.197764 | 0.010150 |
| period_2005 | -0.157653 | 0.053100 | -2.968977 | 0.00298793 | -0.261727 | -0.053579 |
| period_2010 | -0.233257 | 0.056102 | -4.157725 | 3.21432e-05 | -0.343215 | -0.123299 |
| period_2015 | -0.372708 | 0.052185 | -7.142082 | 9.19265e-13 | -0.474989 | -0.270428 |
| ln_latent_distance_km | -0.910079 | 0.027080 | -33.606978 | 0 | -0.963155 | -0.857003 |

## Average Annual Flow Comparison Over Full Sample

Fitted values below are latent-distance PPML fitted means.

| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| Spain (ESP) | 363,742.138 | 290,960.770 | 135,494.855 | 187,635.519 | 228,247.283 | 103,325.251 |
| United States (USA) | 2,285,537.328 | 1,829,916.984 | 943,916.082 | 606,353.857 | 1,341,621.246 | 1,223,563.127 |
| Italy (ITA) | 342,668.390 | 413,568.300 | 191,491.188 | 222,189.858 | 151,177.202 | 191,378.442 |
| France (FRA) | 350,714.070 | 455,299.356 | 267,801.926 | 228,432.805 | 82,912.145 | 226,866.551 |
| Germany (DEU) | 602,520.377 | 501,333.905 | 349,678.215 | 386,647.152 | 252,842.163 | 114,686.753 |

## Spain By Period

Observed and fitted values below are annualized flows. Fitted values are latent-distance PPML fitted means. Balance is `inflow - outflow`.

| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| 1990-1995 | 170,883.848 | 173,055.246 | 45,122.692 | 169,123.024 | 125,761.156 | 3,932.222 |
| 1995-2000 | 216,706.606 | 206,286.782 | 50,693.896 | 165,440.640 | 166,012.710 | 40,846.143 |
| 2000-2005 | 575,094.972 | 264,272.558 | 69,236.416 | 186,884.630 | 505,858.556 | 77,387.928 |
| 2005-2010 | 571,463.868 | 359,470.101 | 122,850.870 | 198,412.164 | 448,612.998 | 161,057.936 |
| 2010-2015 | 213,632.720 | 392,138.026 | 303,163.170 | 206,813.360 | -89,530.450 | 185,324.665 |
| 2015-2020 | 434,670.814 | 350,541.908 | 221,902.086 | 199,139.297 | 212,768.728 | 151,402.611 |