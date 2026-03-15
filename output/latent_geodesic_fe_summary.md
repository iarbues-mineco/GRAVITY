# Latent Geodesic PPML With Origin and Destination Fixed Effects

## Fit Statistics

| Metric | Value |
|---|---:|
| Flow unit | annualized_per_year |
| Prediction target | E[flow|X] = exp(alpha_i + delta_j + Xb) |
| Countries | 191 |
| Free country coordinates | 189 |
| Anchors | USA, CHN |
| Origin fixed effects | 191 |
| Destination fixed effects | 191 |
| Observations used in estimation | 207,998 |
| Zero-flow observations included | 62,221 |
| Parameters | 776 |
| Coordinate iterations | 57 |
| Slope iterations in last step | 0 |
| FE iterations in last step | 80 |
| Reference pseudo-log-likelihood | -1912288662.337 |
| Final pseudo-log-likelihood | -1777267390.750 |
| Log-likelihood improvement | 135021271.586 |
| McFadden pseudo-R-squared | -3.834445 |
| RMSE (flow) | 5795.430116 |
| Mean observed flow | 454.390272 |
| Mean fitted flow | 0.000308 |
| Mean displacement (km, free countries) | 2473.780 |
| Median displacement (km, free countries) | 1924.084 |
| Max displacement (km, free countries) | 13345.953 |

## Slope Coefficients

Standard errors below are conditional on the estimated country fixed effects.

| Term | Estimate | Conditional Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |
|---|---:|---:|---:|---:|---:|---:|
| ln_origin_population_total | 0.380741 | 74533.697328 | 0.000005 | 0.999996 | -146082.981655 | 146083.743138 |
| ln_destination_population_total | 0.461279 | 43685.964352 | 0.000011 | 0.999992 | -85622.455481 | 85623.378040 |
| ln_origin_gdp_pc_ppp_constant | 0.196553 | 27202.593406 | 0.000007 | 0.999994 | -53315.906810 | 53316.299916 |
| ln_destination_gdp_pc_ppp_constant | 0.408134 | 54468.976627 | 0.000007 | 0.999994 | -106756.824330 | 106757.640598 |
| contig | 0.297185 | 221104.855815 | 0.000001 | 0.999999 | -433357.257020 | 433357.851389 |
| comlang_off | 0.552273 | 143587.787844 | 0.000004 | 0.999997 | -281426.340521 | 281427.445066 |
| colony | 0.345212 | 494514.030031 | 0.000001 | 0.999999 | -969229.343499 | 969230.033922 |
| col45 | 0.530763 | 1316733.286864 | 0.000000 | 1 | -2580749.288736 | 2580750.350261 |
| smctry | 0.387820 | 416957.663963 | 0.000001 | 0.999999 | -817221.616625 | 817222.392266 |
| ln_migrant_stock_both_sexes_plus1 | 0.272779 | 17588.787369 | 0.000016 | 0.999988 | -34473.116996 | 34473.662554 |
| period_1995 | -0.088381 | 160432.888744 | -0.000001 | 1 | -314442.772255 | 314442.595493 |
| period_2000 | -0.101849 | 160979.663010 | -0.000001 | 0.999999 | -315514.443591 | 315514.239893 |
| period_2005 | -0.172171 | 154397.844408 | -0.000001 | 0.999999 | -302614.386501 | 302614.042159 |
| period_2010 | -0.251779 | 157170.040138 | -0.000002 | 0.999999 | -308047.869898 | 308047.366340 |
| period_2015 | -0.390439 | 149994.056784 | -0.000003 | 0.999998 | -293983.339630 | 293982.558752 |
| ln_latent_distance_km | -0.821494 | 25398.829641 | -0.000032 | 0.999974 | -49781.612840 | 49779.969852 |

Full country fixed effects are written to the technical CSV outputs, not listed here.

## Average Annual Flow Comparison Over Full Sample

Fitted values below are latent-distance PPML fitted means.

| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| Spain (ESP) | 363,742.138 | 0.108 | 135,494.855 | 0.000 | 228,247.283 | 0.108 |
| United States (USA) | 2,285,537.328 | 3.128 | 943,916.082 | 2.325 | 1,341,621.246 | 0.802 |
| Italy (ITA) | 342,668.390 | 0.080 | 191,491.188 | 0.000 | 151,177.202 | 0.080 |
| France (FRA) | 350,714.070 | 0.054 | 267,801.926 | 0.001 | 82,912.145 | 0.054 |
| Germany (DEU) | 602,520.377 | 0.359 | 349,678.215 | 0.001 | 252,842.163 | 0.358 |

## Spain By Period

Observed and fitted values below are annualized flows. Fitted values are latent-distance PPML fitted means. Balance is `inflow - outflow`.

| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| 1990-1995 | 170,883.848 | 0.086 | 45,122.692 | 0.000 | 125,761.156 | 0.086 |
| 1995-2000 | 216,706.606 | 0.089 | 50,693.896 | 0.000 | 166,012.710 | 0.089 |
| 2000-2005 | 575,094.972 | 0.105 | 69,236.416 | 0.000 | 505,858.556 | 0.105 |
| 2005-2010 | 571,463.868 | 0.123 | 122,850.870 | 0.000 | 448,612.998 | 0.123 |
| 2010-2015 | 213,632.720 | 0.127 | 303,163.170 | 0.000 | -89,530.450 | 0.127 |
| 2015-2020 | 434,670.814 | 0.118 | 221,902.086 | 0.000 | 212,768.728 | 0.118 |