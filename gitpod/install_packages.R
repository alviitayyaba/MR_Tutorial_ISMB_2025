install.packages(c(
  "devtools", "tidyverse", "ggplot2", "ggforce", "dplyr"
), repos = "https://cloud.r-project.org")

# BiocManager needed for gwasvcf
install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "gwasvcf"
))

devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
devtools::install_github("MRCIEU/ieugwasr")
devtools::install_github("MRCIEU/gwasglue2")

# Install genetics.binaRies for plink binary
devtools::install_github("explodecomputer/genetics.binaRies")

# Download PLINK binary
genetics.binaRies::get_plink_binary()
