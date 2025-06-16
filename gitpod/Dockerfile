FROM rocker/tidyverse:latest

# Install RStudio Server
RUN apt-get update && \
    apt-get install -y gdebi-core && \
    wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2023.06.1-524-amd64.deb && \
    gdebi -n rstudio-server-2023.06.1-524-amd64.deb && \
    rm rstudio-server-2023.06.1-524-amd64.deb

# Install dependencies (as before)
COPY gitpod/install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R
