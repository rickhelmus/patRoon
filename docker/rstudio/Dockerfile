FROM rocker/rstudio:4.0.0-ubuntu18.04

ENV SETUPDIR=/usr/local/setup

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget libv8-3.14-dev openbabel \
        libxml2-dev pngquant openjdk-11-jdk libmagick++-dev pandoc git pngquant texinfo libfribidi-dev \
        zlib1g-dev libxml2-dev libnetcdf-dev libglpk-dev tzdata libnetcdf-dev netcdf-bin && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://abibuilder.informatik.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/2.5.0/OpenMS-2.5.0-Debian-Linux-x86_64.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-2.5.0-Debian-Linux-x86_64.deb && \
    rm -rf $SETUPDIR && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    R CMD javareconf

USER rstudio

COPY --from=patroon:latest --chown=rstudio:rstudio /home/patRoon/deps /home/rstudio/deps
COPY --from=patroon:latest --chown=rstudio:rstudio /home/patRoon/.Rprofile /home/rstudio

RUN sed -i -e 's/\/home\/patRoon/\/home\/rstudio/g' ~/.Rprofile

ENV OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0 R_MAX_NUM_DLLS=150 R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

USER root

COPY --from=patroon:latest /usr/local/lib/R/site-library/. /usr/local/lib/R/site-library/
