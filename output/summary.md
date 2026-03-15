# CEPII Stock PPML With Origin and Destination FE

## Fit Statistics

| Metric | Value |
|---|---:|
| Flow unit | annualized_per_year |
| Prediction target | E[flow|X] = exp(alpha_i + delta_j + Xb) |
| Observations in panel | 207,998 |
| Observations used in estimation | 207,998 |
| Zero-flow observations included | 62,221 |
| Countries | 191 |
| Origin fixed effects | 191 |
| Destination fixed effects | 191 |
| Parameters | 397 |
| Slope iterations | 200 |
| FE iterations in last slope step | 35 |
| Converged | False |
| Log-likelihood | -52958261.880 |
| Null log-likelihood | -367625935.264 |
| McFadden pseudo-R-squared | 0.855945 |
| RMSE (flow) | 3080.539561 |
| Mean observed flow | 454.390272 |
| Mean fitted flow | 454.390272 |

## Slope Coefficients

Standard errors below are conditional on the estimated origin and destination fixed effects.

| Term | Estimate | Conditional Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |
|---|---:|---:|---:|---:|---:|---:|
| ln_origin_population_total | 0.940810 | 0.008184 | 114.951043 | 0 | 0.924769 | 0.956852 |
| ln_destination_population_total | -0.068303 | 0.008527 | -8.010317 | 1.11022e-15 | -0.085015 | -0.051590 |
| ln_origin_gdp_pc_ppp_constant | -0.140947 | 0.011978 | -11.767579 | 0 | -0.164422 | -0.117471 |
| ln_destination_gdp_pc_ppp_constant | 0.184936 | 0.014002 | 13.207516 | 0 | 0.157492 | 0.212380 |
| ln_distance_km | -0.577750 | 0.021932 | -26.343300 | 0 | -0.620735 | -0.534765 |
| contig | 0.213168 | 0.045033 | 4.733625 | 2.20545e-06 | 0.124906 | 0.301431 |
| comlang_off | 0.062716 | 0.037506 | 1.672152 | 0.0944943 | -0.010795 | 0.136228 |
| colony | 0.517879 | 0.093643 | 5.530364 | 3.19568e-08 | 0.334342 | 0.701416 |
| col45 | 0.652239 | 0.092054 | 7.085401 | 1.38645e-12 | 0.471817 | 0.832662 |
| smctry | 0.259010 | 0.062971 | 4.113193 | 3.90224e-05 | 0.135590 | 0.382430 |
| ln_migrant_stock_both_sexes_plus1 | 0.348020 | 0.006545 | 53.175734 | 0 | 0.335192 | 0.360847 |
| period_1995 | -0.118357 | 0.055871 | -2.118394 | 0.0341417 | -0.227862 | -0.008852 |
| period_2000 | -0.063205 | 0.051941 | -1.216879 | 0.22365 | -0.165007 | 0.038596 |
| period_2005 | -0.056741 | 0.051440 | -1.103060 | 0.270001 | -0.157562 | 0.044079 |
| period_2010 | -0.080659 | 0.056352 | -1.431347 | 0.152331 | -0.191106 | 0.029788 |
| period_2015 | -0.175063 | 0.050118 | -3.493009 | 0.00047761 | -0.273293 | -0.076833 |

Full country fixed effects are written to the technical CSV outputs, not listed here.

## Average Annual Flow Comparison Over Full Sample

Fitted values below are PPML fitted means.

| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| Spain (ESP) | 363,742.138 | 363,742.138 | 135,494.855 | 135,494.855 | 228,247.283 | 228,247.283 |
| United States (USA) | 2,285,537.328 | 2,285,537.328 | 943,916.082 | 943,916.082 | 1,341,621.246 | 1,341,621.246 |
| Italy (ITA) | 342,668.390 | 342,668.390 | 191,491.188 | 191,491.188 | 151,177.202 | 151,177.202 |
| France (FRA) | 350,714.070 | 350,714.070 | 267,801.926 | 267,801.926 | 82,912.145 | 82,912.145 |
| Germany (DEU) | 602,520.377 | 602,520.377 | 349,678.215 | 349,678.215 | 252,842.163 | 252,842.163 |

## Spain By Period

Observed and fitted values below are annualized flows. Fitted values are PPML fitted means. Balance is `inflow - outflow`.

| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |
|---|---:|---:|---:|---:|---:|---:|
| 1990-1995 | 170,883.848 | 221,978.068 | 45,122.692 | 128,614.736 | 125,761.156 | 93,363.332 |
| 1995-2000 | 216,706.606 | 243,095.696 | 50,693.896 | 115,354.255 | 166,012.710 | 127,741.442 |
| 2000-2005 | 575,094.972 | 313,976.136 | 69,236.416 | 128,310.218 | 505,858.556 | 185,665.919 |
| 2005-2010 | 571,463.868 | 440,453.514 | 122,850.870 | 139,968.492 | 448,612.998 | 300,485.022 |
| 2010-2015 | 213,632.720 | 498,913.695 | 303,163.170 | 151,764.332 | -89,530.450 | 347,149.364 |
| 2015-2020 | 434,670.814 | 464,035.718 | 221,902.086 | 148,957.098 | 212,768.728 | 315,078.621 |