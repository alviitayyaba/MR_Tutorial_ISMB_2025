FROM rocker/tidyverse:latest

# Install system dependencies and RStudio Server
RUN apt-get update && apt-get install -y \
    gdebi-core \
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
    wget \
    plink1.9 && \
    wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2023.06.1-524-amd64.deb && \
    gdebi -n rstudio-server-2023.06.1-524-amd64.deb && \
    rm rstudio-server-2023.06.1-524-amd64.deb

# Add R packages
COPY gitpod/install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R

# Allow login without password (default user is `gitpod`)
RUN echo "auth-none=1" >> /etc/rstudio/rserver.conf
