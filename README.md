# asthma-lactylation-RCC2

Analysis code for the asthma lactylation / RCC2 project.

## Repository layout

- `methods/`: analysis scripts organized in a Methods-section order
- `scripts/`: helper CLIs (e.g., table generation)

## Environment

- R: 4.3.1
- Python: 3.x

### R (recommended)

This repository includes a `renv.lock` for dependency tracking.

```r
install.packages("renv")
renv::restore()
```

If you cannot use `renv`, see `R_sessionInfo.txt` for a reference environment.

### Python

```bash
python -m pip install -r requirements.txt
```

## Quick start: generate ML summary tables

The main reproducible artifact in this repo is the ML performance summary table.

```bash
python scripts/make_ml_summary_tables.py \
  --benchmark_csv path/to/benchmark_113.csv \
  --lmbd_csv path/to/lmbd.csv \
  --single_pathway_csv path/to/single_pathway.csv \
  --out_dir ml_summary_tables
```

Outputs:
- `ml_summary_tables/full_benchmark_ranked.csv` / `.xlsx`
- `ml_summary_tables/main_table_x.csv` / `.xlsx` / `.md`

## Reproducibility notes

- Many R scripts intentionally keep the original hard-coded `setwd("D:\\...")` patterns.
  Run scripts from the intended working directory or adjust paths locally.
- `methods/03_lactylation/01_pathway_scoring/lactylation_pathway_scoring.R` calls `source("./1.R")`.
  The helper file `1.R` is not included in this repository and must be provided locally to run the script as-is.
