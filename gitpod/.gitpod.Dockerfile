FROM gitpod/workspace-base:latest

USER root

# Install R base
RUN apt-get update && apt-get install -y \
    dirmngr \
    gnupg \
    software-properties-common \
    wget \
    gdebi-core \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libgit2-dev

# Add CRAN repo and install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor > /etc/apt/trusted.gpg.d/cran.gpg \
 && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
 && apt-get update \
 && apt-get install -y r-base

# Install RStudio Server
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.06.1-524-amd64.deb && \
    gdebi -n rstudio-server-2023.06.1-524-amd64.deb && \
    rm rstudio-server-2023.06.1-524-amd64.deb

# Create a default user (Gitpod uses "gitpod")
RUN echo "gitpod:gitpod" | chpasswd

EXPOSE 8787
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize=0"]

USER gitpod
