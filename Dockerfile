FROM rocker/r-apt:bionic

LABEL maintainer="Rick Helmus <r.helmus@uva.nl>" \
    org.label-schema.name="patRoon" \
    org.label-schema.description="Docker image used for patRoon CI" \
    org.label-schema.schema-version="1.0" \
    org.label-schema.vcs-url="https://github.com/rickhelmus/patRoon" \
    org.label-schema.vendor="patRoon" \
    docker.cmd.test="docker run -t patroonorg/patroon /bin/bash -c 'cd patRoon; Rscript docker/run_tests.R'"

ENV SETUPDIR=/usr/local/setup

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget openbabel libv8-3.14-dev \
        libxml2-dev libnetcdf-dev netcdf-bin pngquant openjdk-11-jdk libmagick++-dev pandoc git pngquant texinfo \
        r-cran-checkmate r-cran-data.table r-cran-withr r-cran-digest r-cran-xml r-cran-xml2 r-cran-dbi \
        r-cran-rsqlite r-cran-rjava r-cran-dplyr r-cran-rcolorbrewer \
        r-cran-htmlwidgets r-cran-shiny r-cran-knitr r-cran-r.utils \
        r-cran-ggplot2 r-cran-jsonlite r-cran-igraph r-cran-hmisc \
        r-cran-robustbase r-cran-testthat r-cran-biocmanager \
        r-bioc-biocparallel r-bioc-affy r-bioc-biocgenerics r-bioc-biobase \
        r-bioc-rbgl r-bioc-s4vectors r-bioc-biocparallel r-bioc-multtest && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://abibuilder.informatik.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/2.5.0/OpenMS-2.5.0-Debian-Linux-x86_64.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-2.5.0-Debian-Linux-x86_64.deb && \
    rm -rf $SETUPDIR && \
    useradd -ms /bin/bash patRoon && \
    addgroup patRoon staff && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    R CMD javareconf

USER patRoon
WORKDIR /home/patRoon

COPY ./docker/install_deps.R ./DESCRIPTION ./

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN wget http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-CL.jar && \
    wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.4.17/sirius-4.4.17-linux64-headless.zip && \
    unzip sirius-4.4.17-linux64-headless.zip && rm sirius-4.4.17-linux64-headless.zip && \
    echo 'options(patRoon.path.metFragCL = "~/MetFrag2.4.5-CL.jar")' >> .Rprofile && \
    echo 'options(patRoon.path.SIRIUS = "~/sirius-linux64-headless-4.4.17/bin")' >> .Rprofile && \
    echo 'options(patRoon.progress.opts = list(style = 1))' >> .Rprofile && \
    Rscript install_deps.R && rm -f ~/install_deps.R ~/DESCRIPTION

ADD --chown=patRoon . patRoon

RUN Rscript -e 'devtools::install(pkg = "patRoon", upgrade = FALSE)'

ENV OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0 R_MAX_NUM_DLLS=150

CMD ["R"]
