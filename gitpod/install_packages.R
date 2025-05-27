install.packages(c(
  "TwoSampleMR",
  "MRInstruments",
  "tidyverse",
  "ieugwasr",
  "gwasvcf",
  "gwasglue2",
  "ggplot2",
  "ggforce",
  "dplyr"
), repos = "https://cloud.r-project.org")

# Install from GitHub if needed:
if (!requireNamespace("genetics.binaRies", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("MRCIEU/genetics.binaRies")
}
