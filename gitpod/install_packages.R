install.packages(c(
  "remotes",           # Needed for install_github
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

# Install PLINK helper package
if (!requireNamespace("genetics.binaRies", quietly = TRUE)) {
  remotes::install_github("MRCIEU/genetics.binaRies")
}

# Optional: Download PLINK binary if needed
genetics.binaRies::get_plink_binary()
