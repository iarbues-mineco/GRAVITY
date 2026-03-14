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
python -m gravity_world.cli estimate-minimal-model
python -m gravity_world.cli normalize-cepii
python -m gravity_world.cli assemble-cepii-panel
python -m gravity_world.cli estimate-cepii-model
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

## CEPII panel and model

Build the CEPII-extended panel with:

```bash
python -m gravity_world.cli assemble-cepii-panel
```

Estimate the extended model with:

```bash
python -m gravity_world.cli estimate-cepii-model
```

Main outputs:

- `data/processed/panels/bilateral_panel_cepii.csv`
- `data/processed/models/cepii_gravity_ols_summary.txt`
- `data/processed/models/cepii_gravity_ols_spain_totals.csv`

The CEPII model adds bilateral distance and standard dyadic indicators to the mass-only benchmark.