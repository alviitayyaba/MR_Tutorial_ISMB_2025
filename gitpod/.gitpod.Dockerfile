# Start from rocker R 4.5.1 base image (check if available)
FROM ghcr.io/rocker-org/devcontainer/r-ver:4.2

# Install system dependencies (if needed)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages used locally (optional)
RUN R -e "install.packages(c('tidyverse', 'radian', 'ieugwasr', 'ggplot2'), repos='https://cran.rstudio.com/')"

# Install RStudio Server
RUN apt-get update && apt-get install -y gdebi-core \
    && wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.06.2-561-amd64.deb \
    && gdebi -n rstudio-server-2023.06.2-561-amd64.deb \
    && rm rstudio-server-2023.06.2-561-amd64.deb \
    && rm -rf /var/lib/apt/lists/*

# Set working directory to Gitpod workspace folder
WORKDIR /workspace

# Expose RStudio port
EXPOSE 8787