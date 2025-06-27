# Install required CRAN packages
cran_packages <- c(
  "remotes",
  "ieugwasr",
  "dplyr",
  "tidyverse",
  "ggplot2"
)
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# Install Bioconductor package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("VariantAnnotation", ask = FALSE, update = FALSE)

# Install GitHub packages
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/MRInstruments")
remotes::install_github("mrcieu/gwasglue")
