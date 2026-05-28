# paper_2025_joint_model_importance

Research code for the 2025 paper project on **Joint Model Importance** based on
Rashomon sets.

## What this project does

Rather than interpreting a single "best" model, this project studies the
**Rashomon set** — the set of all models whose predictive performance is within
`epsilon` of the best model — and investigates how these near-equally good
models differ in their behaviour and feature importance, and how they relate
to one another.

In concrete terms:

1. Trained models from an upstream AutoML run are loaded
   (`data/run_models_merged.rds`, `data/results_vic.RData`,
   `data/design.RData`).
2. For each task, the Rashomon set is determined from the performance scores
   (`get_RS()` in `init/functions.R`).
3. For every model in the Rashomon set, predictions are computed on a
   canonical validation split (2/3), parallelised via
   [`batchtools`](https://mlr-org.github.io/batchtools/) on a socket cluster
   (30 CPUs).
4. A pairwise distance matrix between models is built from these predictions
   (Euclidean / Manhattan).
5. The relationships between models are visualised and analysed using:
   - hierarchical clustering (`hclust`),
   - classical (metric) multidimensional scaling (`cmdscale`),
   - partitioning around medoids (`cluster::pam`).

## Datasets (tasks)

Defined in `init/assets/tasks.R`:

| Key | Task                                    | Type            |
|-----|-----------------------------------------|-----------------|
| gc  | German Credit                           | Classification  |
| cs  | COMPAS                                  | Classification  |
| bs  | Bike Sharing                            | Regression      |
| st  | Synthetic Task (10,000 observations)    | Regression      |

## Learners

Defined in `init/assets/learners.R` — both regression and classification
variants, with tuning spaces provided via `mlr3tuning` / `paradox`:

- `xgb`    — XGBoost
- `tree`   — rpart decision tree
- `nnet`   — single-hidden-layer neural network (`nnet`)
- `glmnet` — elastic net (`glmnet`)
- `svm`    — support vector machine (`e1071`); kernel is part of the tuning space

## Project structure

```
.
├── main.R                  # main script: load, build Rashomon sets, compute,
│                           # cluster, and visualise
├── init.R                  # sources everything under init/
├── init/
│   ├── source.R            # library imports + subdirectory loading
│   ├── functions.R         # get_performance_and_SVMkernel(), get_RS()
│   ├── zzz_settings.R      # task/learner lists, loads design + VIC
│   ├── assets/
│   │   ├── tasks.R         # mlr3 tasks (incl. synthetic task)
│   │   ├── learners.R      # mlr3 learners with tuning spaces
│   │   └── datasplits.R
│   └── batchtools/
│       └── registry.R      # batchtools registry setup
├── data/                   # pre-computed inputs (AutoML results, VIC, design)
├── renv/                   # renv project library
└── renv.lock               # pinned package versions
```

## Prerequisites

- R ≥ 4.3.0 (required for reproducible synthetic data via
  `RNGversion("4.3.0")`)
- Restore the pinned package versions with `renv`:
  ```r
  renv::restore()
  ```
- External model path: `main.R` reads individual sample models from
  `/media/external/rashomon/datafiles/<task>/<learner>/...` and writes the
  batchtools registry to
  `/media/external/ewaldf/paper_2025_joint_model_importance`. These paths must
  be adapted for local execution.

## Running the analysis

```r
source("main.R")
```

`main.R` is written as an interactive analysis script — cluster setup, job
definition, submission, and post-processing run step-by-step from top to
bottom. The central parameters at the top of the script are:

- `RS_epsilon = 0.01` — tolerance for the Rashomon set (1%)
- `distance_metrics = c("euclidean", "manhattan")`
