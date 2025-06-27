FROM rocker/rstudio:4.2.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libzstd-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && apt-get clean

# Install devtools, remotes, BiocManager
RUN R -e "install.packages(c('devtools', 'remotes', 'BiocManager'), repos='https://cloud.r-project.org')"
