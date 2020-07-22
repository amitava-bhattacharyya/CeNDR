# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest
##FROM rocker/r-ver:3.2.2

## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libxml2-dev \
libcairo2-dev \
libsqlite3-dev \
libmariadbd-dev \
libpq-dev \
libssh2-1-dev \
unixodbc-dev \
libcurl4-openssl-dev \
libssl-dev

RUN R -e "install.packages('lme4', repos=c('http://lme4.r-forge.r-project.org/repos', getOption('repos')[['CRAN']]))"
RUN R -e "install.packages('dplyr', dependencies=TRUE)"

# Create directories
RUN mkdir -p /home/data /home/output /home/script

WORKDIR /home/data

COPY ./H2_script.R /home/script/H2_script.R

# Run the script
CMD Rscript /home/script/H2_script.R

