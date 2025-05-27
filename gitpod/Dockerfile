FROM rocker/rstudio:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libv8-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    plink1.9 \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R
