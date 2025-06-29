FROM rocker/tidyverse:4.2.3

# Install system dependencies if needed
RUN apt-get update && apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN Rscript -e "install.packages(c('dplyr', 'R.utils', 'ggplot2'), repos='https://cloud.r-project.org')" \
    -e "if (!requireNamespace('devtools', quietly=TRUE)) install.packages('devtools')" \
    -e "devtools::install_github('MRCIEU/MRInstruments')" \
    -e "devtools::install_github('MRCIEU/ieugwasr')" \
    -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes')" \
    -e "remotes::install_github('MRCIEU/TwoSampleMR')"