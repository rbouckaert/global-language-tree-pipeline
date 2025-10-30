## Code to analyse language diversification rates for Bouckeart et al.

### Scripts

Code authored by David Redding and Artur Trebski unless noted otherwise below.


-   **001_cut_trees.R**: Processes global language phylogenies, subdivides them by continent and taxonomy, and summarizes tree-level statistics for each subtree. Outputs: subtrees and summary statistics files.
-   **002_prepare_datasets.R**: Generates input datasets of socio-environmental variables for each subtree, merging relevant metadata. Outputs: cleaned datasets for further analyses.
-   **003_simulate_trees.R**: Simulates birth-death phylogenetic trees using rates estimated from empirical trees, and calculates lineage-through-time (LTT) data. Outputs: simulated trees and LTT statistics.
-   **004a_time_slice_datasets.R**: Creates time-sliced datasets and extracts values from trees at specified time points, preparing data for downstream analyses. Outputs: time-sliced dataframes.
-   **004b_tree_slicing_three_points.R**: Slices phylogenetic trees at three defined time points and summarizes clade counts and properties for each slice. Outputs: clade count tables and time-slice summaries.
-   **005ac_regression_analyses.R**: Performs spatio-temporal regression analyses on language diversification rates, including cross-validation and model selection. Outputs: regression results and model summaries. Implements INLA-based spatial models following Martino & Rue (Martino S, Rue H. Implementing approximate Bayesian inference using Integrated Nested Laplace Approximation: A manual for the inla program. Department of Mathematical Sciences, NTNU, Norway. 2009) and incorporates a phylogenetic covariance component adapted from Russell Dinnage.
-   **005b_regression_regions.R**: Runs regression analyses excluding specific regions to assess regional effects on diversification rates. Outputs: regional exclusion regression summaries and comparison plots. Shares the INLA formulation and phylogenetic component credits noted above.
-   **006_regression_summary.R**: Visualizes regression analysis results, generating plots and summary statistics for interpretation. Outputs: regression plots and statistics files.
-   **007_prepare_threat_dataframe.R**: Prepares and merges threat data with language datasets, incorporating geographical and threat status information. Outputs: threat-annotated dataframes.
-   **008_EDGE_calculate.R**: Calculates EDGE (Evolutionarily Distinct and Globally Endangered) scores and phylogenetic diversity loss metrics for languages. Sources the EDGE2 function contributed by Ricky Gumbs (see Gumbs et al. 2022). Outputs: EDGE scores and diversity loss tables.
-   **009_EDGE_Summarise_scores.R**: Summarizes and compares EDGE scores across parameter sets, producing summary tables and visualizations. Outputs: EDGE summary statistics and comparison plots.
-   **010_LTT_data.R**: Generates lineage-through-time (LTT) data and plots for each subtree, including continent-wise analyses. Outputs: LTT data files and plots.
-   **011_language_analysis_bisse.R**: Conducts BiSSE/HiSSE trait-dependent diversification analyses on phylogenetic trees, processing and summarizing results. Outputs: BiSSE/HiSSE analysis summaries. Adapted from the BiSSE toolkit assembled by David Redding, incorporating routines originally authored by JM Beaulieu.
-   **012_simulation_trees_summary_statsb.R**: Calculates summary statistics for simulated tree sets, including diversity and structure metrics. Outputs: statistics tables for simulated trees.
-   **013_plot_and_summarise_summary_stats.R**: Plots and summarizes statistics from simulated trees, providing visual interpretation of simulation results. Outputs: summary plots and tables.
-   **014_BISSE_summary.R**: Summarizes BiSSE/HiSSE results, generating tables and plots for trait-dependent diversification metrics. Outputs: BiSSE/HiSSE summary plots and tables.

### Credits and sourcing

-   EDGE2 functionality is provided via the implementation by Ricky Gumbs, used here with attribution (Gumbs et al. 2022).
-   Threat data (`res5.csv` / `res5.r`) are derived from Bromham et al. (2021) *Nat Ecol Evol* doi:10.1038/s41559-021-01604-y.

### Folders:

-   code

    - analysis

    - figures

    - tables
   
-   input_data
