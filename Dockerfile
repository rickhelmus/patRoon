FROM rocker/r-ver:3.4.4

ENV SETUPDIR=/usr/local/setup

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget \
        libxml2-dev libnetcdf-dev netcdf-bin pngquant openjdk-8-jdk libmagick++-dev pandoc \
        r-cran-checkmate r-cran-data.table r-cran-withr r-cran-digest r-cran-xml r-cran-xml2 r-cran-dbi \
        r-cran-rsqlite r-cran-rjava r-cran-dplyr r-cran-rcolorbrewer \
        r-cran-htmlwidgets r-cran-shiny r-cran-knitr r-cran-r.utils \
        r-cran-ggplot2 r-cran-jsonlite r-cran-igraph r-cran-hmisc \
        r-cran-robustbase r-cran-testthat r-bioc-biocinstaller \
        r-bioc-biocparallel r-bioc-affy r-bioc-biocgenerics r-bioc-biobase \
        r-bioc-rbgl r-bioc-s4vectors r-bioc-biocparallel r-bioc-multtest && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://github.com/OpenMS/OpenMS/releases/download/Release2.3.0/OpenMS-2.3.0-Linux.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-2.3.0-Linux.deb && \
    rm -rf $SETUPDIR && \
    useradd -ms /bin/bash patRoon && \
    addgroup patRoon staff && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER patRoon
WORKDIR /home/patRoon

COPY ./docker/install_deps.R ./DESCRIPTION ./

RUN wget http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.3-CL.jar && \
    wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0/sirius-4.0-linux64-headless.zip && \
    unzip sirius-4.0-linux64-headless.zip && rm sirius-4.0-linux64-headless.zip && \
    Rscript install_deps.R

ADD --chown=patRoon . patRoon

ENV OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0
ARG FAIL_TESTS=1

WORKDIR patRoon
RUN Rscript docker/run_tests.R || [ $FAIL_TESTS -eq 0 ] && cat ~/olog/* && \
    cp ~/patRoon/junit.xml ~/ && rm -rf ~/patRoon /tmp/Rtmp* ~/install_deps.R

