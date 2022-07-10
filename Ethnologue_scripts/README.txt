
Code for anaylsis of diversification of global human languages

All coded in R language.

Author: David Redding (unless stated)

Date: 05/07/2022

calculate_EDGE2_scores.R - R script that reads newick files from the data folder 'all_trees' and applies the EDGE2 function (EDGE2.R) using threat data from 'res5.csv'.

calculate_summary_statistics. R - R script that reads newick files from the data folder 'all_trees' and applies a series of functions to obtain mean phylogenetic diversity, evolutionary distinctness and other tree measurements to produce a summary table.

EDGE2.R - R script written by Ricky Gumbs to calculate the evolutionary distinctness of a species on a tree transformed by threat status (Gumbs et al. https://www.biorxiv.org/content/10.1101/2022.05.17.492313v1).

language_analysis_bisse7.R - R script that takes trait data (all_tips_by_year6ALL.csv) and trees in newick format and calculates diversification rate parameters for different trait states (BISSE analyses adapted from script by JM Beaulieu)

language_analysis_regression.R - R script that performs a phylogenetic and spatial regression in INLA (Martino S, Rue H. Implementing approximate Bayesian inference using Integrated Nested Laplace Approximation: A manual for the inla program. Department of Mathematical Sciences, NTNU, Norway. 2009.)using matching geographical data and time-sliced language phylogenies (in the folder 'datasets_and_trees'. Phylogenetic component adapted from code by Russell Dinnage.

res5.csv - Threat data for each global language from Bromham, L. et al. Global predictors of language endangerment and the future of linguistic diversity. Nat Ecol Evol (2021) doi:10.1038/s41559-021-01604-y.

Folders:
all_trees - 901 posterior global language trees in newick format. 
Available by donloading and uncompressing https://github.com/rbouckaert/global-language-tree-pipeline/releases/download/v1.0.0/all_trees.tgz

simulation_trees - 901 matching simulated trees with same crown height and approximate birth-death parameters - can be swapped in for 'all_trees' to create expectation on random trees.
Available by donloading and uncompressing https://github.com/rbouckaert/global-language-tree-pipeline/releases/download/v1.0.0/simulation_trees.tgz

data_sets_and_trees - 901 backbone trees from 3 time points alongside matching geographical summaries for the subclades at those tips.
Available by donloading and uncompressing from https://github.com/rbouckaert/global-language-tree-pipeline/releases/download/v1.0.0/ 3500_csv.tgz, 3500_tre.tgz, 4250_csv.tgz, 4250_tre.tgz, 5000_csv.tgz and 5000_tre.tgz
