This repository contains the functions and supplementary materials developed for my Master's thesis in Bioinformatics and Biostatistics at the Open University of Catalonia and the University of Barcelona. The files contain:

sample_processing.Rmd: a file containing the DADA2 processing scripts for the synthetic and environmental sample communities.

k_optimization.Rmd: contains the optimization scripts for the various evaluated algorithms referenced in the thesis.

k_validation.Rmd: contains the validation scripts for the optimized algorithms in a real-world context, along with their evaluation. ANOVA and PERMANOVA results are also included.

Project Functions.R: It contains all the functions necessary for the development of the scripts described in the .Rmd files.

ASV_table.csv: ASV count table for the synthetic community

ASV_table_env.csv: ASV count table for the environmental community


Taxonomia_100.csv: taxonomy file after the search performed in Blast to complete NAs

"all_plot_list_k(x).rds" files: rds format files containing all the plots evaluated during the optimization of each algorithm.

Supplementary material: contains the evaluation metrics for the optimized algorithms, as well as the performance evaluation metrics for the k2 and k4 filters. These are described in the final thesis.