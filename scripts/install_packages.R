# R vesrion
#Minimum: R ≥ 4.0.0
#Recommended: R ≥ 4.2.0 (stable and compatible with major MR packages)
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("ieugwasr")
install.packages("dplyr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
devtools::install_github("MRCIEU/MRInstruments")
install.packages(c('TwoSampleMR','gwasvcf', 'gwasglue', 'MRInstruments'), repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))
install.packages("tidyverse")
install.packages("ggplot2")
