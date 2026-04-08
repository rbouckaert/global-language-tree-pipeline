# TreeSetAnalysisScripts

This directory contains the R-based downstream analysis workflow used for the manuscript **"Global language diversification is linked to socio-ecology and threat status"**. It starts from a posterior sample of global language phylogenies and associated spatial/cultural/threat datasets, and produces the diversification analyses, EDGE calculations, summary tables, and manuscript figures.

This directory does **not** build the global language trees themselves. Tree construction and XML-generation belong to the separate `monos/` and BEAST workflow. `TreeSetAnalysisScripts/` covers the downstream analyses run on the posterior treeset.

Unless otherwise noted below, code in this directory was authored by David Redding and Artur Trebski.

## Repository contents

- `code/analysis/`: core analytical scripts (`001`-`014`) plus helper scripts.
- `code/figures/`: scripts for main and supplementary figures.
- `code/tables/`: scripts for supplementary tables.
- `input_data/`: input tables, shapefiles, and some derived inputs used by the analysis scripts.
- `outputs/`: intermediate and working outputs.

## Script guide

- `001_cut_trees_parallel.R`: processes global language phylogenies, subdivides them by continent and taxonomy, and summarizes tree-level statistics for each subtree. Outputs: subtrees and summary statistics files.
- `002_create_enviromental_variables.R`: generates environmental predictor datasets used downstream in the socio-ecological analyses. Outputs: cleaned environmental inputs and derived spatial summaries.
- `003_simulate_trees.R`: simulates birth-death phylogenetic trees using rates estimated from empirical trees and calculates lineage-through-time statistics. Outputs: simulated trees and LTT summaries.
- `004a_time_slice_datasets.R`: creates time-sliced datasets and extracts clade-level values from trees at specified time points, preparing data for downstream regressions. Outputs: time-sliced dataframes and matching trees.
- `004b_tree_slicing_three_points.R`: slices phylogenetic trees at three defined time points and summarizes clade counts and properties for each slice. Outputs: clade count tables and time-slice summaries.
- `005ac_regression_analyses.R`: performs spatio-temporal regression analyses on language diversification rates, including model summaries and fitted response outputs. Implements INLA-based spatial models following Martino and Rue and incorporates a phylogenetic covariance component adapted from Russell Dinnage.
- `005b_regression_regions.R`: runs regression analyses excluding specific regions to assess regional effects on diversification rates. Shares the INLA formulation and phylogenetic component credits noted above.
- `006_regression_summary.R`: visualizes regression analysis results, generating plots and summary statistics for interpretation. Outputs: regression plots and statistics files.
- `007_prepare_threat_dataframe.R`: prepares and merges threat data with language datasets, incorporating geographical and threat-status information. Outputs: threat-annotated dataframes.
- `008_EDGE_calculate.R`: calculates EDGE scores and phylogenetic diversity loss metrics for languages. Sources `EDGE2.R`, contributed by Ricky Gumbs and used here with attribution following Gumbs et al. (2022). Outputs: EDGE scores and diversity loss tables.
- `009_EDGE_Summarise_scores.R`: summarizes and compares EDGE scores across parameter sets, producing summary tables and visualizations. Outputs: EDGE summary statistics and comparison plots.
- `010_LTT_data.R`: generates lineage-through-time data and plots for each subtree, including continent-wise analyses. Outputs: LTT data files and plots.
- `011_language_analysis_bisse.R`: conducts BiSSE/HiSSE trait-dependent diversification analyses on phylogenetic trees and summarizes those results. Adapted from the BiSSE toolkit assembled by David Redding, incorporating routines originally authored by JM Beaulieu.
- `012_simulation_trees_summary_statsb.R`: calculates summary statistics for simulated tree sets, including diversity and structure metrics. Outputs: statistics tables for simulated trees.
- `013_plot_and_summarise_summary_stats.R`: plots and summarizes statistics from simulated trees, providing visual interpretation of simulation results. Outputs: summary plots and tables.
- `014_BISSE_summary.R`: summarizes BiSSE/HiSSE results, generating tables and plots for trait-dependent diversification metrics.

## Credits and sourcing

- `EDGE2.R` functionality is provided via the implementation by Ricky Gumbs, used here with attribution.
- Threat data inputs are derived from Bromham et al. (2022) *Nature Ecology & Evolution* doi:10.1038/s41559-021-01604-y.

## 1. System requirements

### Operating systems tested

- Tested locally on `macOS 26.3.1` (`Darwin 25.3.0`).
- The scripts are expected to run on Linux with a comparable R installation and geospatial stack, but Linux was not re-validated in this draft pass.

### Software requirements

- `R 4.5.2`
- A standard command-line environment capable of running `Rscript`

### Key R package dependencies

Core analysis scripts use the following packages:

- `ape 5.8-1`
- `picante 1.8.2`
- `phytools 2.5-2`
- `data.table 1.18.2.1`
- `here 1.0.2`
- `caper 1.0.4`
- `RPANDA 2.4`
- `diversitree 1.7-23`
- `paleotree 3.1-168`
- `phangorn 0.10-1`
- `sf 0.10.1`
- `sp 1.8-93`
- `raster 3.6-32`
- `terra 1.4-2`
- `INLA 2.0.0`
- `spdep 4.0.2`
- `progress 1.2.0`
- `tidyverse 0.6.3`
- `ggplot2 2.3`
- `cowplot 1.1-3`
- `ggpubr 0.5.7`
- `ggridges 5.2.0`

Additional figure/table scripts also call supporting packages such as `RColorBrewer`, `reshape2`, `forcats`, `readxl`, `countrycode`, `psych`, `viridis`, `rnaturalearth`, `rnaturalearthdata`, `magick`, `patchwork`, `ggh4x`, `ggnewscale`, `ggbeeswarm`, and `lm.beta`.

### Legacy/archived R packages

Some scripts still call older spatial or phylogenetic packages that may require archived installation or a legacy R library:

- `rgdal`
- `maptools`
- `rgeos`
- `apTreeshape`
- `hisse`
- `bazar`

For exact historical reproduction of the original workflow, it may be easier to use the maintained analysis environment associated with the manuscript submission rather than a clean modern R installation.

### Data requirements

The full workflow expects several large inputs, including:

- the posterior treeset `input_data/global-language-trees-6636-taxa.trees`
- environmental raster directories such as `input_data/temp_data_10000BC/`, `input_data/gcrop/`, `input_data/baseline/zip/`, and `input_data/2020_walking_only_friction_surface/`
- downstream derived inputs created by earlier scripts

If a lightweight clone omits the largest files, those assets must be added back before full pipeline can be run.

### Additional large input sources

Some large environmental and spatial input folders are not included in a
standard clone. In addition to the GitHub release assets described below, the
following Dropbox folders can be used to restore major inputs required by
`002_create_enviromental_variables.R`:

- `input_data/2020_motorized_friction_surface/`
[Dropbox link](https://www.dropbox.com/scl/fo/4x49mscdyh3whdoqbbdjj/AB0tLi5TsBtLI7YXPpPhnYU/input_data/2020_motorized_friction_surface?rlkey=oir9k3ka3qvsk2i7i7st860w6&e=1&st=ismcmfxj&subfolder_nav_tracking=1&dl=0)
- `input_data/2020_walking_only_friction_surface/`
[Dropbox link](https://www.dropbox.com/scl/fo/4x49mscdyh3whdoqbbdjj/AGhC3UAOGSDgFtmJ4MM72G4/input_data/2020_walking_only_friction_surface?rlkey=oir9k3ka3qvsk2i7i7st860w6&subfolder_nav_tracking=1&st=k5uzap08&dl=0)
- `input_data/baseline/`
[Dropbox link](https://www.dropbox.com/scl/fo/4x49mscdyh3whdoqbbdjj/AJo8yKDJ-2jSI9Y_H37bHSw/input_data/baseline?rlkey=oir9k3ka3qvsk2i7i7st860w6&subfolder_nav_tracking=1&st=8vp12gqa&dl=0)
- `input_data/gcrop/`
[Dropbox link](https://www.dropbox.com/scl/fo/4x49mscdyh3whdoqbbdjj/ACQncQHd42nqa1F9OlIxayc/input_data/gcrop?rlkey=oir9k3ka3qvsk2i7i7st860w6&subfolder_nav_tracking=1&st=v6442v54&dl=0)
- `input_data/langa/`
[Dropbox link](https://www.dropbox.com/scl/fo/4x49mscdyh3whdoqbbdjj/AE2LYNjkgLkMuC1MtMqwnlY/input_data/langa?rlkey=oir9k3ka3qvsk2i7i7st860w6&subfolder_nav_tracking=1&st=s972iap9&dl=0)
- `input_data/temp_data_10000BC/`
[Dropbox link](https://www.dropbox.com/scl/fo/4x49mscdyh3whdoqbbdjj/AOWHfiF6KXC-JF1MjQkrm0I/input_data/temp_data_10000BC?rlkey=oir9k3ka3qvsk2i7i7st860w6&subfolder_nav_tracking=1&st=zt3896n6&dl=0)

### Hardware requirements

- No non-standard hardware is required.
- A multi-core desktop or workstation is strongly recommended for full reruns.
- At least `16 GB RAM` is recommended for the heavier tree, raster, and regression steps.

## 2. Installation guide

### Installation instructions

1. Clone the repository and enter this directory:

```bash
git clone https://github.com/rbouckaert/global-language-tree-pipeline.git
cd global-language-tree-pipeline/TreeSetAnalysisScripts
```

1. Install R package dependencies. A minimal installation example is:

```r
install.packages(c(
  "pacman","here","data.table","ape","picante","moments","phytools","doParallel",
  "foreach","caper","phyloTop","RPANDA","raster","sf","exactextractr","imputeTS",
  "dismo","e1071","phangorn","mice","paleotree","diversitree","sp","nlme","ade4",
  "terra","rworldmap","spdep","progress","tidyverse","ggplot2","cowplot","ggpubr",
  "gridExtra","RColorBrewer","ggridges","forcats","reshape2","ggthemes",
  "rnaturalearth","rnaturalearthdata","psych","viridis","readxl","coda",
  "stringr","magrittr","maps","MASS"
))
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable")
)
```

1. If exact reproduction is required, install any archived/legacy packages needed by the specific scripts you intend to run.
2. Download the required large data assets from the GitHub Releases page and unpack them into the expected folders. A standard `git clone` does **not** download GitHub release assets automatically.

GitHub Releases:
[https://github.com/rbouckaert/global-language-tree-pipeline/releases](https://github.com/rbouckaert/global-language-tree-pipeline/releases)

### Release assets required after cloning

At minimum, the following release assets should be downloaded and unpacked into `TreeSetAnalysisScripts/`:

- `global-language-trees-6636-taxa.tgz`
Unpack so the file is available at `input_data/global-language-trees-6636-taxa.trees`
- `all_tree.tgz`
Unpack so the directory is available at `all_trees/`
- `final_env_data.tgz`
Unpack so the file is available at `input_data/final_env_data.csv`
- `final_env_datapoints.tgz`
Unpack so the file is available at `input_data/final_env_datapoints.csv`
- `processed_threat_data_frame.tgz`
Unpack so the file is available at `input_data/processed_threat_data_frame.csv`
- `islands_transformed.tgz`
Unpack so the file is available at `input_data/islands_transformed.r`
- `langa.tgz`
Unpack so the directory is available at `input_data/langa/`

If you want to start from precomputed time-sliced datasets rather than regenerating them from `004a`, the following release assets can also be unpacked into `datasets_and_trees/`:

- `3500_csv.tgz`
- `3500_tre.tgz`
- `4250_csv.tgz`
- `4250_tre.tgz`
- `5000_csv.tgz`
- `5000_tre.tgz`

When those six archives are unpacked, the paired `.csv` and `.tre` files should sit directly inside `datasets_and_trees/`.

### Typical install time

- On a normal desktop computer with an existing R toolchain: approximately `15-45 minutes`.
- Installation can take longer if geospatial packages or archived packages need to be compiled from source.

## 3. Demo

Because the full manuscript workflow is computationally heavy, we provide a  
lightweight demo that exercises the main TreeSet socio-ecological regression  
workflow on 3 example trees and 3 time points.

### Demo inputs

- `demo/all_trees_demo/`
- supporting input tables in `input_data/`

The demo also depends on several large supporting files that are distributed through GitHub release assets rather than through a standard clone. Before running the demo, make sure the following paths exist:

- `input_data/final_env_data.csv`
from `final_env_data.tgz`
- `input_data/final_env_datapoints.csv`
from `final_env_datapoints.tgz`
- `input_data/processed_threat_data_frame.csv`
from `processed_threat_data_frame.tgz` if it is not already present in your clone
- `input_data/islands_transformed.r`
from `islands_transformed.tgz`
- `input_data/langa/`
from `langa.tgz`

### Demo commands

Run the following commands from `TreeSetAnalysisScripts/`:

```bash
Rscript code/demo/01_prepare_demo_datasets.R
Rscript code/demo/02_run_demo_regression.R
Rscript code/demo/03_plot_demo_results.R
```

### What the demo does

- reads 3 example posterior trees from `demo/all_trees_demo/`
- creates demo time-sliced datasets at `3500`, `4250`, and `5000` years before present
- runs the main demo regression workflow on those generated datasets
- writes a compact demo summary figure

### Expected demo output

The demo writes outputs under:

- `demo/datasets_and_trees_demo/`
- `demo/regression_results_demo/`
- `demo/regression_slopes_demo/`
- `demo/demo_summary/`
- `demo/demo_figures/`

Expected files include:

- `demo/datasets_and_trees_demo/*.csv`
- `demo/datasets_and_trees_demo/*.tre`
- `demo/regression_results_demo/*.csv`
- `demo/regression_slopes_demo/*.csv`
- `demo/demo_summary/demo_regression_run_summary.csv`
- `demo/demo_summary/demo_coefficients_summary.csv`
- `demo/demo_summary/demo_session_info.txt`
- `demo/demo_figures/demo_regression_summary.png`
- `demo/demo_figures/demo_regression_summary.pdf`

See [demo/README.md](/Users/arturtrebski/Documents/Coding_Projects/global-language-tree-pipeline/TreeSetAnalysisScripts/demo/README.md) for the demo folder layout and output descriptions.

### Expected demo runtime

- On a normal desktop computer: approximately `10-20 minutes` total.
- Typical breakdown:
  - dataset preparation: `2-5 minutes`
  - demo regressions: `5-15 minutes`
  - plotting: `<1 minute`

## 4. Instructions for use

### Working directory

Run all scripts from the `TreeSetAnalysisScripts/` directory unless a script explicitly states otherwise.

### Canonical analysis order

The full downstream analysis is organized in the following order:

1. `code/analysis/001_cut_trees_parallel.R`
2. `code/analysis/002_create_enviromental_variables.R`
3. `code/analysis/003_simulate_trees.R`
4. `code/analysis/004a_time_slice_datasets.R`
5. `code/analysis/004b_tree_slicing_three_points.R`
6. `code/analysis/005ac_regression_analyses.R` or `code/analysis/005a_regression_analyses.R`
7. `code/analysis/005b_regression_regions.R`
8. `code/analysis/006_regression_summary.R`
9. `code/analysis/007_prepare_threat_dataframe.R`
10. `code/analysis/008_EDGE_calculate.R`
11. `code/analysis/009_EDGE_Summarise_scores.R`
12. `code/analysis/010_LTT_data.R`
13. `code/analysis/011_language_analysis_bisse.R`
14. `code/analysis/012_simulation_trees_summary_statsb.R`
15. `code/analysis/013_plot_and_summarise_summary_stats.R`
16. `code/analysis/014_BISSE_summary.R`

After the analysis scripts complete, run the figure and table scripts in:

- `code/figures/`
- `code/tables/`

### Typical usage

Examples:

```bash
Rscript code/analysis/001_cut_trees_parallel.R
Rscript code/analysis/005ac_regression_analyses.R
Rscript code/figures/Fig_4.R
Rscript code/tables/Table_S1.R
```

### Major output locations

- `summary/`: per-tree summary statistics from tree cutting
- `all_trees/`: cut tree files generated from the posterior treeset
- `simulationtrees/`: simulated tree sets
- `datasets_and_trees/`: time-sliced datasets and backbone trees
- `demo/`: lightweight demo inputs, outputs, and figures
- `regression_results_*`: regression model outputs
- `EDGE_scores/` and `outputs/EDGEscores/`: EDGE metrics and summaries
- `BISSE_runs/`: BiSSE/HiSSE outputs
- `outputs/`: the main script-written location for tables, summaries, and intermediate products
- `outputs/figures/`: the main script-written location for figure exports

## 5. Reproduction instructions

### Reproducing the manuscript analyses

To reproduce the quantitative results reported in the manuscript:

1. Assemble the complete manuscript data bundle, including the posterior treeset and all large raster/spatial inputs.
2. Run the analysis scripts in canonical order from `001` to `014`.
3. Run the figure and table scripts.
4. Compare generated outputs with the key CSV summaries and figure exports written to `outputs/`.

For most users, step 1 requires downloading the large GitHub release assets and unpacking them into the locations listed above, because those files are not included in a standard clone of the repository.

### Key reproduction targets

Examples of manuscript-linked outputs include:

- `outputs/All_Trees_Summary.csv`
- `outputs/in_text_data/DR_rates_across_tips.csv`
- `outputs/EDGEscores/ALL_EDGE.csv`
- `outputs/EDGEscores/NEW_EDGE_top_100.csv`
- `outputs/bisse_results1.csv`
- `outputs/DR_test_summary.csv`
- `outputs/figures/Fig_1/`
- `outputs/figures/Fig_2/`
- `outputs/figures/Fig_3/`
- `outputs/figures/Fig_4/`
- `outputs/figures/Fig_S1/`
- `outputs/figures/Fig_S2/`
- `outputs/figures/Fig_S3/`
- `outputs/figures/Fig_S4/`
- `outputs/figures/Fig_S5/`
- `outputs/figures/Fig_S6/`
- `outputs/figures/Fig_S7/`

### Expected full-workflow runtime

- A full rerun is substantially heavier than the demo and should be treated as a long-running analysis.
- On a normal multi-core desktop or workstation, expect the complete workflow to take from `many hours to multiple days`, depending on whether all steps are rerun and whether intermediate files are already present.

## Citation

If you use this code or data, please cite the associated manuscript and repository release used for the review or analysis.