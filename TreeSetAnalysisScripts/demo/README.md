# TreeSet Demo

This directory contains a lightweight demo for the
`TreeSetAnalysisScripts` workflow. The demo is intentionally much smaller than
the full manuscript pipeline: it uses 3 example posterior trees, generates
time-sliced regression datasets, runs the main socio-ecological regression
model on those demo datasets, and produces a compact summary figure.

The demo is intended as a qualitative illustration of the workflow rather than
as a full reproduction of the manuscript analyses.

## Before running the demo

A standard `git clone` does **not** download the large GitHub release assets.
Before running the demo, make sure the following additional files have been
downloaded from the release page and unpacked into `TreeSetAnalysisScripts/`:

- `final_env_data.tgz`
  Unpack to `input_data/final_env_data.csv`
- `final_env_datapoints.tgz`
  Unpack to `input_data/final_env_datapoints.csv`
- `processed_threat_data_frame.tgz`
  Unpack to `input_data/processed_threat_data_frame.csv` if that file is not
  already present in your clone
- `islands_transformed.tgz`
  Unpack to `input_data/islands_transformed.r`
- `langa.tgz`
  Unpack to `input_data/langa/`

The 3 example trees used by the demo are already included in this repository at
`demo/all_trees_demo/`, so the demo does not require downloading the full
posterior treeset release asset.

## Run order

Run the following commands from `TreeSetAnalysisScripts/`:

```bash
Rscript code/demo/01_prepare_demo_datasets.R
Rscript code/demo/02_run_demo_regression.R
Rscript code/demo/03_plot_demo_results.R
```

Typical runtime on a normal desktop computer is approximately:

- step 1: `2-5 minutes`
- step 2: `5-15 minutes`
- step 3: `<1 minute`

## Folder contents

- `all_trees_demo/`
Three example `Global_all` trees used as the demo starting point.
- `datasets_and_trees_demo/`
Demo time-sliced datasets and matching backbone trees created by
`01_prepare_demo_datasets.R`.
- `cl_graph_demo/`
Spatial graph objects generated for the demo regression runs.
- `regression_results_demo/`
Per-run regression outputs created by `02_run_demo_regression.R`.
- `regression_slopes_demo/`
Per-run fitted response curve outputs created by
`02_run_demo_regression.R`.
- `demo_summary/`
Compact CSV summaries of demo runs and coefficient outputs.
- `demo_figures/`
Demo figure outputs created by `03_plot_demo_results.R`.

## Expected outputs

After running all 3 demo steps, the main outputs should include:

- `demo/datasets_and_trees_demo/*.csv`
- `demo/datasets_and_trees_demo/*.tre`
- `demo/regression_results_demo/*.csv`
- `demo/regression_slopes_demo/*.csv`
- `demo/demo_summary/demo_regression_run_summary.csv`
- `demo/demo_summary/demo_coefficients_summary.csv`
- `demo/demo_summary/demo_session_info.txt`
- `demo/demo_figures/demo_regression_summary.png`
- `demo/demo_figures/demo_regression_summary.pdf`

## Interpretation

- The demo uses 3 trees rather than the full posterior sample.
- The demo figure is intended to show the structure of the workflow and the
qualitative behavior of fitted effects.
- Full manuscript reproduction requires the larger release assets and the
canonical analysis pipeline described in the main README.
