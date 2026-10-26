# Climate-driven shifts in maize and cassava suitability in Uganda

R code accompanying the article

> Wangalwa, R. *et al.* **Projected climate-driven shifts in maize and cassava suitability in Uganda: implications for adaptation planning and food security.** *Discover Sustainability* (submitted).

The pipeline models current agro-climatic suitability for maize and cassava with an ensemble of machine-learning species distribution models (biomod2: XGBoost, Random Forest, MaxNet), projects it to 16 future climate scenarios (4 SSPs × 4 periods, 2021–2100), and quantifies suitable area and its change at national, regional and sub-regional level.

## Workflow

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 0 | `R/00_config.R` | Crop selection, paths and all analysis parameters | n/a |
| 1 | `R/01_extract_env_at_occurrences.R` | Summary of bioclimatic variables at occurrence points | `tables/*env_variable_stats*.csv` |
| 2 | `R/02_fit_ensemble_model.R` | SRE pseudo-absences (≥ 5 km from presences), XGBOOST/RF/MAXNET with 2-fold CV, EMmean/EMwmean ensembles, projection to current climate | `baseline/<crop>_ensemble_EMwmean_current.tif` |
| 3 | `R/03_evaluate_model.R` | AUC/TSS tables, variable importance, response curves, AUC bar chart and ROC curves | `tables/`, `figures/` |
| 4 | `R/04_project_future.R` | Ensemble projection for every future scenario folder | `projections/SSPxxx_YYYY-YYYY.tif` |
| 5 | `R/05_baseline_suitability_summary.R` | Area per suitability class (country, region, sub-region); donut, treemap, map, alluvial and lollipop plots | `tables/`, `figures/` |
| 6 | `R/06_future_suitability_area.R` | Five-class reclassification of future maps; national and sub-regional suitable area; line plot and heatmap | `tables/`, `figures/` |
| 7 | `R/07_change_detection.R` | Gain / loss / persistence of optimal suitability (≥ 0.6) vs. baseline, with areas | `change/`, `tables/`, `figures/` |
| 8 | `R/08_plot_future_maps.R` | 4 × 4 panel map of future suitability | `figures/` |

Shared helper functions are in `R/utils.R`. All outputs are written to `outputs/<crop>/`.

### Key settings (from `R/00_config.R`)

| Setting | Value |
|---|---|
| Projection | WGS 84 / UTM 36N (EPSG:32636); climate rasters aggregated by a factor of 2 |
| Pseudo-absences | SRE strategy, quantile 0.025, 2 repetitions, as many as presences; those within 5 km of a presence removed |
| Algorithms | XGBOOST, RF, MAXNET |
| Cross-validation | k-fold, k = 2, balanced on presences; prevalence 0.5; seed 123 |
| Ensemble | EMmean and EMwmean (proportional weights); models with validation AUC ≥ 0.6; seed 42 |
| Suitability classes | < 0.2 Not, 0.2–0.4 Marginally, 0.4–0.6 Moderately, 0.6–0.8 Highly, 0.8–1.0 Very highly suitable |
| Optimal suitability (change analysis) | ≥ 0.6 |

## Getting started

### 1. Requirements

- R ≥ 4.1
- R packages: `terra`, `sf`, `biomod2` (≥ 4.2), `xgboost`, `randomForest`, `maxnet`, `pROC`, `dplyr`, `tidyr`, `here`, `ggplot2`, `patchwork`, `scales`, `treemapify`, `ggalluvial`

```r
source("install_packages.R")
```

For exact reproducibility, use [renv](https://rstudio.github.io/renv/): `renv::init()` then `renv::snapshot()` to record package versions in `renv.lock`.

### 2. Data

**Data availability:** the input data will be made publicly available upon acceptance of the article. Until then they are not included in this repository. Place them under `data/` as described in [`data/README.md`](data/README.md). Paths can be changed in `R/00_config.R`.

### 3. Run

From the repository root (or open the folder as an RStudio project):

```bash
Rscript run_all.R              # both crops, all steps
Rscript run_all.R maize        # one crop
Rscript run_all.R cassava 5 8  # one crop, steps 5 to 8
```

Or run a single script for one crop:

```bash
CROP=cassava Rscript R/02_fit_ensemble_model.R
```

```r
Sys.setenv(CROP = "cassava"); source("R/02_fit_ensemble_model.R")
```

Step 2 fits the models and takes the longest. Steps 3 to 8 reuse its saved outputs.

## Repository structure

```
.
├── R/
│   ├── 00_config.R                      # paths and parameters
│   ├── utils.R                          # shared helpers
│   ├── 01_extract_env_at_occurrences.R
│   ├── 02_fit_ensemble_model.R
│   ├── 03_evaluate_model.R
│   ├── 04_project_future.R
│   ├── 05_baseline_suitability_summary.R
│   ├── 06_future_suitability_area.R
│   ├── 07_change_detection.R
│   └── 08_plot_future_maps.R
├── data/                # input data (not tracked, see data/README.md)
├── outputs/             # results (not tracked)
├── docs/CODE_REVIEW_NOTES.md
├── run_all.R
├── install_packages.R
├── CITATION.cff
└── LICENSE
```

## Citation

If you use this code, please cite the article above (see `CITATION.cff`).

## License

Code released under the MIT License (see `LICENSE`). Input data are subject to the licences of their original providers.

## Contact

Rapheal Wangalwa: open an issue in this repository.
