# World Migration Gravity Data

Repository scaffold for the first step of a global gravity model of migration flows: reproducible data acquisition.

## Baseline scope

This repo is set up for a standard bilateral country-to-country gravity model with:

- Dependent variable: annualized bilateral migration flows derived from Abel-Cohen 5-year estimates, 1990-1995 through 2015-2020.
- Country covariates: annual World Bank indicators from 1990-2024.
- Bilateral controls: CEPII GeoDist variables such as distance, contiguity, common language, and colonial links.
- Network control: UN DESA migrant stock by origin-destination, 1990-2024.

## Variables included in the acquisition design

### Core outcome

- `flow`: annualized bilateral migration flow, obtained by dividing the 5-year period total by period length.
- `flow_total_period`: original 5-year bilateral flow estimate for the period.

### Country covariates

- `population_total`: World Bank `SP.POP.TOTL`, 1990-2024.
- `gdp_pc_ppp_constant`: World Bank `NY.GDP.PCAP.PP.KD`, 1990-2024.
- `unemployment_total_pct`: World Bank `SL.UEM.TOTL.ZS`, 1990-2024.
- `gini_index`: World Bank `SI.POV.GINI`, 1990-2024.

### Bilateral controls

- `distance_km`: CEPII `distw` population-weighted bilateral distance.
- `contig`
- `comlang_off`
- `comlang_ethno`
- `colony`
- `comcol`
- `curcol`
- `col45`
- `smctry`

### Network control

- `migrant_stock_both_sexes`: UN DESA bilateral migrant stock at mid-year.
- `migrant_stock_male`
- `migrant_stock_female`
- `ln_migrant_stock_both_sexes_plus1`: `log(1 + migrant_stock_both_sexes)`.

## Quick start

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -e .
python -m gravity_world.cli download --source world_bank_wdi
python -m gravity_world.cli build-country-reference
python -m gravity_world.cli normalize-abel-cohen
python -m gravity_world.cli assemble-country-year
python -m gravity_world.cli assemble-minimal-panel
python -m gravity_world.cli normalize-cepii
python -m gravity_world.cli assemble-cepii-panel
```

If you want the migrant-stock extension too, place the UN DESA `Destination and origin` workbook in `data/raw/un_desa_stock/` and run:

```bash
python -m gravity_world.cli normalize-un-desa-stock
python -m gravity_world.cli assemble-cepii-stock-panel
```

## Country harmonization

Build the canonical reference with:

```bash
python -m gravity_world.cli build-country-reference
```

This writes:

- `data/processed/reference/country_reference.csv`
- `data/processed/reference/source_country_overrides.csv`

## Abel-Cohen flow normalization

Normalize the raw Abel-Cohen flow file with:

```bash
python -m gravity_world.cli normalize-abel-cohen
```

This writes:

- `data/processed/flows/bilateral_flows_5y_abel_cohen.csv`
- `data/processed/flows/abel_cohen_unmatched_origins.csv`
- `data/processed/flows/abel_cohen_unmatched_destinations.csv`

## Country-year covariates

Build the World Bank country-year table with:

```bash
python -m gravity_world.cli assemble-country-year
```

This writes:

- `data/processed/country_year_covariates.csv`

## Minimal panel and model

Build the first estimation-ready panel with:

```bash
python -m gravity_world.cli assemble-minimal-panel
```

Estimate the baseline mass-only model with:

```bash
python -m gravity_world.cli estimate-minimal-model
```

Main outputs:

- `data/processed/panels/bilateral_panel_minimal.csv`
- `data/processed/models/minimal_gravity_ols_summary.txt`
- `data/processed/models/minimal_gravity_ols_spain_totals.csv`

## CEPII dyadic controls

Place or download `dist_cepii.zip` into `data/raw/cepii/`, then run:

```bash
python -m gravity_world.cli normalize-cepii
```

This writes:

- `data/processed/dyadic/cepii_geodist_controls.csv`
- `data/processed/dyadic/cepii_unmatched_origins.csv`
- `data/processed/dyadic/cepii_unmatched_destinations.csv`

## UN DESA migrant stock

Place the UN DESA `Destination and origin` workbook in `data/raw/un_desa_stock/`, then run:

```bash
python -m gravity_world.cli normalize-un-desa-stock
```

If you keep several UN DESA workbooks in that folder, the normalizer will automatically prefer the one whose filename contains both `destination` and `origin`.

This writes:

- `data/processed/stocks/bilateral_migrant_stock_un_desa.csv`
- `data/processed/stocks/un_desa_stock_unmatched_origins.csv`
- `data/processed/stocks/un_desa_stock_unmatched_destinations.csv`

## CEPII panel and models

Build the CEPII-extended panel with:

```bash
python -m gravity_world.cli assemble-cepii-panel
```

Estimate the extended model with:

```bash
python -m gravity_world.cli estimate-cepii-model
```

Build the CEPII plus migrant stock panel with:

```bash
python -m gravity_world.cli assemble-cepii-stock-panel
```

Estimate the stock-augmented model with:

```bash
python -m gravity_world.cli estimate-cepii-stock-model
```

Estimate the stock model plus origin and destination unemployment with:

```bash
python -m gravity_world.cli estimate-cepii-stock-unemployment-model
```

Estimate the stock model with PPML on the full nonnegative flow sample with:

```bash
python -m gravity_world.cli estimate-cepii-stock-ppml-model
```

Estimate the experimental latent geodesic PPML with the default anchors (`ESP FRA`) with:

```bash
python -m gravity_world.cli estimate-latent-geodesic-model
```

To change the fixed countries, pass ISO3 codes after `--anchors`. For example, to fix the United States and China and let Spain move:

```bash
python -m gravity_world.cli estimate-latent-geodesic-model --anchors USA CHN
```

To add a quadratic penalty that keeps countries closer to their real-world coordinates, use:

```bash
python -m gravity_world.cli estimate-latent-geodesic-model --anchors USA CHN --penalty-weight 1.0
```

The penalty is applied to squared chord displacement from the World Bank reference map, scaled so values around `0.25` to `2.0` are a reasonable first range to try.

To render country-surface maps from the latent coordinates, use:

```bash
python -m gravity_world.cli plot-latent-country-maps --model-prefix latent_geodesic_ppml
```

This writes two SVGs: one that keeps countries in place and draws centroid arrows to their displaced positions, and one that shifts each country shape by its latent displacement vector.

Main outputs:

- `data/processed/panels/bilateral_panel_cepii.csv`
- `data/processed/panels/bilateral_panel_cepii_stock.csv`
- `data/processed/models/cepii_gravity_ols_summary.txt`
- `data/processed/models/cepii_stock_gravity_ols_summary.txt`
- `output/summary.md`

The CEPII model adds bilateral distance and standard dyadic indicators to the mass-only benchmark. The stock model adds lagged bilateral migrant stock at the period start year. The root `output/` folder is the curated user-facing location for the current model summary.
