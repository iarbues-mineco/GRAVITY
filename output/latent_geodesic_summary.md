# Latent Geodesic PPML

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
| Final pseudo-log-likelihood | -62276719.600 |
| Log-likelihood improvement | 21069273.726 |
| McFadden pseudo-R-squared | 0.830598 |
| RMSE (flow) | 3351.624113 |
| Mean observed flow | 454.390272 |
| Mean fitted flow | 454.390272 |
| Mean displacement (km, free countries) | 3796.984 |
| Median displacement (km, free countries) | 2803.286 |
| Max displacement (km, free countries) | 12909.301 |

## Coefficients

| Term | Estimate | Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |
|---|---:|---:|---:|---:|---:|---:|
| const | -8.329965 | 0.405744 | -20.530109 | 0 | -9.125209 | -7.534722 |
| ln_origin_population_total | 0.380741 | 0.010402 | 36.601708 | 0 | 0.360353 | 0.401130 |
| ln_destination_population_total | 0.461279 | 0.010576 | 43.616816 | 0 | 0.440551 | 0.482007 |
| ln_origin_gdp_pc_ppp_constant | 0.196553 | 0.015246 | 12.891789 | 0 | 0.166671 | 0.226435 |
| ln_destination_gdp_pc_ppp_constant | 0.408134 | 0.019671 | 20.748041 | 0 | 0.369580 | 0.446688 |
| contig | 0.297185 | 0.041591 | 7.145348 | 8.97726e-13 | 0.215667 | 0.378702 |
| comlang_off | 0.552273 | 0.044199 | 12.495021 | 0 | 0.465643 | 0.638902 |
| colony | 0.345212 | 0.095337 | 3.620963 | 0.000293508 | 0.158355 | 0.532069 |
| col45 | 0.530763 | 0.092718 | 5.724484 | 1.03749e-08 | 0.349039 | 0.712486 |
| smctry | 0.387820 | 0.076959 | 5.039329 | 4.67166e-07 | 0.236984 | 0.538656 |
| ln_migrant_stock_both_sexes_plus1 | 0.272779 | 0.008085 | 33.737290 | 0 | 0.256932 | 0.288626 |
| period_1995 | -0.088381 | 0.055696 | -1.586860 | 0.112544 | -0.197543 | 0.020780 |
| period_2000 | -0.101849 | 0.054284 | -1.876206 | 0.060627 | -0.208244 | 0.004547 |
| period_2005 | -0.172171 | 0.054261 | -3.173038 | 0.00150853 | -0.278520 | -0.065822 |
| period_2010 | -0.251779 | 0.056529 | -4.453976 | 8.42945e-06 | -0.362573 | -0.140984 |
| period_2015 | -0.390439 | 0.054141 | -7.211540 | 5.53335e-13 | -0.496553 | -0.284325 |
| ln_latent_distance_km | -0.821494 | 0.022186 | -37.027877 | 0 | -0.864977 | -0.778011 |

## Average Annual Flow Comparison Over Full Sample

Fitted values below are latent-distance PPML fitted means.

| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| Spain (ESP) | 363,742.138 | 285,186.592 | 135,494.855 | 179,675.917 | 228,247.283 | 105,510.676 |
| United States (USA) | 2,285,537.328 | 1,851,406.213 | 943,916.082 | 627,223.530 | 1,341,621.246 | 1,224,182.683 |
| Italy (ITA) | 342,668.390 | 393,913.178 | 191,491.188 | 222,661.876 | 151,177.202 | 171,251.302 |
| France (FRA) | 350,714.070 | 417,556.432 | 267,801.926 | 205,089.528 | 82,912.145 | 212,466.903 |
| Germany (DEU) | 602,520.377 | 453,713.087 | 349,678.215 | 362,436.201 | 252,842.163 | 91,276.886 |

## Spain By Period

Observed and fitted values below are annualized flows. Fitted values are latent-distance PPML fitted means. Balance is `inflow - outflow`.

| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| 1990-1995 | 170,883.848 | 164,322.688 | 45,122.692 | 161,851.497 | 125,761.156 | 2,471.191 |
| 1995-2000 | 216,706.606 | 199,960.395 | 50,693.896 | 158,302.252 | 166,012.710 | 41,658.142 |
| 2000-2005 | 575,094.972 | 257,085.579 | 69,236.416 | 179,006.668 | 505,858.556 | 78,078.911 |
| 2005-2010 | 571,463.868 | 353,821.268 | 122,850.870 | 189,730.329 | 448,612.998 | 164,090.939 |
| 2010-2015 | 213,632.720 | 387,620.489 | 303,163.170 | 197,656.868 | -89,530.450 | 189,963.620 |
| 2015-2020 | 434,670.814 | 348,309.135 | 221,902.086 | 191,507.884 | 212,768.728 | 156,801.251 |