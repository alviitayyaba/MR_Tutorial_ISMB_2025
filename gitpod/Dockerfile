FROM rocker/rstudio:4.3.2

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
    libboost-all-dev \
    zlib1g-dev \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Copy R package installation script
COPY gitpod/install_packages.R /tmp/install_packages.R

# Install R packages
RUN Rscript /tmp/install_packages.R
