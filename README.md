# World Migration Gravity Data

Repository scaffold for the first step of a global gravity model of migration flows: reproducible data acquisition.

## Baseline scope

This repo is set up for a standard bilateral country-to-country gravity model with:

- Dependent variable: 5-year bilateral migration flows from Abel-Cohen, 1990-1995 through 2015-2020.
- Country covariates: annual World Bank indicators from 1990-2024.
- Bilateral controls: CEPII gravity variables such as distance, contiguity, common language, and colonial links.
- Network control: UN DESA migrant stock by origin-destination, 1990-2024.

## Variables included in the acquisition design

### Core outcome

- `flow_ijt`: bilateral migration flows, five-year periods, 1990-1995 to 2015-2020.

### Country covariates

- `population_total`: World Bank `SP.POP.TOTL`, 1990-2024.
- `gdp_pc_ppp_constant`: World Bank `NY.GDP.PCAP.PP.KD`, 1990-2024.
- `unemployment_total_pct`: World Bank `SL.UEM.TOTL.ZS`, 1990-2024.
- `gini_index`: World Bank `SI.POV.GINI`, 1990-2024.

### Bilateral controls

- `distance_km`
- `contiguity`
- `common_language`
- `colonial_tie`
- `ever_same_country`

### Network control

- `migrant_stock_ijt`: UN DESA destination-origin stock matrix, 1990-2024.

## Sources

- Abel-Cohen bilateral flows: figshare direct file download.
- World Bank WDI: official V2 API.
- UN DESA migrant stock: official landing page, with a manual or override hook because file links can change between revisions.
- CEPII gravity data: official landing page, with a manual or override hook because the public download path is not stable enough to hard-code.

## Repository layout

```text
config/
  harmonization/
data/
  raw/
  interim/
  processed/
  logs/
scripts/
src/gravity_world/
```

## Quick start

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -e .
python -m gravity_world.cli inventory
python -m gravity_world.cli download --source world_bank_wdi
python -m gravity_world.cli build-country-reference
python -m gravity_world.cli normalize-abel-cohen
python -m gravity_world.cli assemble-country-year
python -m gravity_world.cli assemble-minimal-panel
```

If UN DESA or CEPII direct links are known later, add them in `config/source_overrides.json` using the template in `config/source_overrides.example.json`.

## Country harmonization

The first normalization layer lives in `config/harmonization/`:

- `manual_countries.csv`: canonical manual additions and historical entities that are not reliably covered by live APIs.
- `source_overrides.csv`: source-specific code or name remaps for legacy labels and non-standard reporting.

After downloading World Bank metadata, build the canonical reference with:

```bash
python -m gravity_world.cli build-country-reference
```

This writes:

- `data/processed/reference/country_reference.csv`
- `data/processed/reference/source_country_overrides.csv`

## Abel-Cohen flow normalization

After downloading the raw Abel-Cohen flow file, normalize it with:

```bash
python -m gravity_world.cli normalize-abel-cohen
```

Optional:

```bash
python -m gravity_world.cli normalize-abel-cohen --preferred-flow da_pb_closed
```

This writes:

- `data/processed/flows/bilateral_flows_5y_abel_cohen.csv`
- `data/processed/flows/abel_cohen_unmatched_origins.csv`
- `data/processed/flows/abel_cohen_unmatched_destinations.csv`

The normalized flow table retains all available published flow-estimation variants and exposes one selected series as the baseline `flow` column.

## Minimal panel assembly

Build the first estimation-ready panel with:

```bash
python -m gravity_world.cli assemble-minimal-panel
```

This writes:

- `data/processed/panels/bilateral_panel_minimal.csv`
- `data/processed/panels/bilateral_panel_minimal_dropped_rows.csv`

The minimal panel uses:

- bilateral Abel-Cohen flow
- start-year origin population and GDP per capita
- start-year destination population and GDP per capita
- optional unemployment and Gini columns when available
- no distance and no migrant stock yet