
## Core Input Data (input_data/ directory)

| FileName                                                                      | Format    | Description                                                                                           | OriginalSource | used_by/output_by       |
|-------------------------------------------------------------------------------|-----------|-------------------------------------------------------------------------------------------------------|----------------|-------------------------|
| top27families.csv                                                             | CSV       | Top 27 language families, including tips (languages) associated with each family.                     |                | 001                     |
| Lexicongraphic_familes.csv                                                    | CSV       | Information about 9 lexicographic language families, data sources and glottocode mapping to taxa.     |                | 001                     |
| languoid.csv                                                                  | CSV       | Glottolog languoid data.                                                                              |                | 002                     |
| languages_and_dialects_geoNEW.csv                                             | CSV       | Additional language/dialect metadata [20,000 glottocodes].                                            |                | 002                     |
| urbanspatial-hist-urban-pop-3700bc-ad2000-xlsx.csv                            | CSV       | Historical urban population and location data from 3700 BC to 2000 AD.                                |                | 002                     |
| 2020_walking_only_friction_surface/<br>2020_walking_only_friction_surface.geotiff | Raster    | Global friction surface (walking) for 2020.                                                           |                | 002                     |
| temp_data_10000BC/                                                            | Directory | Raster files for temperature.                                                                         |                | 002                     |
| gcrop/                                                                        | Directory | Raster files for cropland.                                                                            |                | 002                     |
| baseline/zip/                                                                 | Directory | Files matching the pattern "popd": These files contain population density data.                       |                | 002                     |
| langa/langa.shp                                                               | Shapefile | Language spatial data as polygons.                                                                    |                | 002, 004, 007           |
| languages_and_dialects_geoFINALupdates6.csv                                   | CSV       | Language metadata, including family, macroarea, geographical coordinates, glottocodes, threat status. |                | 002, 004, 007, 009, 010 |
| groups_data1.csv                                                              | CSV       | Group names for simulations.                                                                          |                | 003                     |
| global-language-trees-6636-taxa.trees                                         | Nexus     | Language phylogenies in Nexus format. 1000 posterior trees with 6636 tips each.                       |                | 003, 004_tree_slicing   |
| islands_transformed.R                                                         | RData     | Island polygons, used to determine which languages are on islands.                                    |                | 004, 007                |
| threat_data1.csv                                      | CSV       | Threat data derived from Bromham et al. 2022 “Global predictors of language endangerment and the future of linguistic diversity” (_Nat Ecol Evol_, DOI: 10.1038/s41559-021-01604-y, Supplementary Data 4, CC BY 4.0).                                                                          |                | 007                     |
| bisse_data_in.csv                                                             | CSV       | BiSSE input data.                                                                                     |                | 011                     |
| tips_by_continent.csv                                                         | CSV       | Tips by continent.                                                                                    |                | 011                     |
| bisse_param_names.r                                                           | RData     | Parameter names for BiSSE.                                                                            |                | 014                     |

## Helper Scripts (code/ directory)

| FileName                        | Format   | Description                     | used_by/output_by |
|---------------------------------|----------|---------------------------------|-------------------|
| EDGE2.R                         | R script | EDGE score calculation function | 008               |
| AUX2_traitDependent_functions.R | R script | Trait-dependent functions       | 011               |

## Intermediate Data (Generated by previous scripts, used as input by later scripts)

| FileName | Format | Description | Used by / Output by |
|---|---|---|---|
| all_trees/ | Directory | Cut/processed phylogenetic trees | 003, 004, 008, 010, 011 |
| input_data/processed_threat_data_frame.csv | CSV | Processed threat data | 007 (output), 004, 008 (input) |
| input_data/final_env_data.csv | CSV | Environmental data for polygons | 002 (output), 004 (input) |
| input_data/final_env_datapoints.csv | CSV | Environmental data for points | 002 (output), 004 (input) |
| simulationtrees/ | Directory | Simulated trees.<br>Contains files matching pattern "extant" (used by 014). | 003 (output), 004a, 012 (input) |
| datasets_and_trees/ | Directory | Language/environmental data and trees for regression. Global at 3 time points (3500, 4250, 5000 Years).<br>Contains files matching *_2.csv<br>Contains files matching *_2.tre | 004a (output), 005a_regression_analyses, 005b_regression_regions (input) |
| regression_results_no_alt/ | Directory | Regression results.<br>Contains files matching regressions25_*.csv<br>Contains files matching slopes25_*.csv | 005a (output), 006 (input) |
| regression_results_regions/<br>results_*_regressions25.csv | CSV | Regression results for macroregion exclusion | 005b |
| EDGE_scores/ | Directory | EDGE score files.<br>Contains files matching *_EDGESCORES2.csv | 008 (output), 009 (input) |
| EDGE_scores/[tree_name]_EDGESCORES2.csv | CSV | EDGE scores for each tree | 008 [output] |
| EDGE_scores/[tree_name]_PDLOST2.csv | CSV | PD loss calculations for each tree | 008 [output] |
| edscores/ | Directory | Additional/intermediate EDGE scores | 001 (output), 009 (input) |
| summary_sim/ | Directory | Per-tree summary statistics for simulations | 012 (output), 013(input) |
| BISSE_runs/ | Directory | BiSSE analysis results.<br>Contains files matching pattern "ISSE_runs3.csv" | 011 (output), 014 (input) |
| summary/*.csv | CSV | Summary statistics for each cut tree | 001 [output] |
| outputs/glotto_languages_cites_states3points.csv | CSV | Main output for point environmental data | 002 |
| outputs/glotto_languages_cites_states3.csv | CSV | Main output for polygon environmental data | 002 |
| outputs/mean_frictionpoints.csv | CSV | Friction data for points | 002 |
| outputs/mean_frictionpoly.csv | CSV | Friction data for polygons | 002 |


## Output Files Summary (Scripts 000\_\* to 016\_)

| FileName                                            | Format | Description                                            | CreatedInScript |
|-----------------------------------------------------|--------|--------------------------------------------------------|-----------------|
| outputs/All_Trees_Summary.csv                       | CSV    | Combined summary of all trees                          | 001             |
| outputs/in_text_data/DR_rates_across_tips.csv       | CSV    | Summary of DR (diversification rate) rates across tips | 001             |
| outputs/new_trees_sliced/\*\_clade_counts.csv       | CSV    | Number of clades at each time slice per tree           | 004b            |
| outputs/Fig_2/not_zero_table4_linear.csv            | CSV    | Table of covariate proportions above/below zero        | 006             |
| outputs/EDGEscores/ALL_EDGE.csv                     | CSV    | All EDGE scores summary                                | 009             |
| outputs/EDGEscores/ALL_ED_wTREES.csv                | CSV    | EDGE scores with numbered trees                        | 009             |
| outputs/EDGEscores/NEW_EDGE_top_100.csv             | CSV    | Top 100 for 226 languages                              | 009             |
| outputs/EDGEscores/ALL_EDGE_top_100.csv             | CSV    | Top 100 summary with classification                    | 009             |
| outputs/Fig_1/LTT_Plot1_Data.csv                    | CSV    | Data for Lineage Through Time plots                    | 010             |
| outputs/summary_statsMEAN\_<date>.csv               | CSV    | Mean summary statistics, date-stamped                  | 013             |
| outputs/simulated_summary_statsALL_HPDI\_<date>.csv | CSV    | Simulated summary statistics with HPDI, date-stamped   | 013             |
| outputs/bisse_results1.csv                          | CSV    | Combined BiSSE run results                             | 014             |
| outputs/DR_test_summary.csv                         | CSV    | DR test summary statistics                             | 014             |
| outputs/AIC_test_summary.csv                        | CSV    | AIC summary statistics                                 | 014             |
| outputs/bisse_parameter_differences.csv             | CSV    | Parameter differences between states                   | 014             |
| outputs/bisse_parameter_directions.csv              | CSV    | Proportion of runs with higher State 1 parameters      | 014             |
| outputs/bisse_null_AIC_proportions.csv              | CSV    | Proportion of runs where BiSSE outperforms null model  | 014             |
