FROM rocker/r-ver

WORKDIR /usr/local/setup

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget \
        libxml2-dev libnetcdf-dev netcdf-bin pngquant openjdk-8-jdk libmagick++-dev pandoc \
        r-cran-checkmate r-cran-data.table r-cran-withr r-cran-digest r-cran-xml r-cran-xml2 r-cran-dbi \
        r-cran-rsqlite r-cran-rjava \
        r-cran-htmlwidgets r-cran-shiny r-cran-knitr r-cran-r.utils \
        r-cran-ggplot2 r-cran-jsonlite r-cran-igraph r-cran-hmisc \
        r-cran-robustbase r-cran-testthat \
        r-bioc-biocparallel r-bioc-affy r-bioc-biocgenerics && \
    wget https://abibuilder.informatik.uni-tuebingen.de/archive/openms/OpenMSInstaller/nightly/OpenMS-2.2.0-Linux.deb && \
    apt-get install -y ./OpenMS-2.2.0-Linux.deb && \
    useradd -ms /bin/bash patRoon && \
    addgroup patRoon staff

USER patRoon
WORKDIR /home/patRoon

COPY ./install_deps.R .

RUN wget http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.3-CL.jar && \
    wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0/sirius-4.0-linux64-headless.zip && \
    unzip sirius-4.0-linux64-headless.zip && \
    Rscript install_deps.R

ADD --chown=patRoon . patRoon

ENV OPENMS_DATA_PATH=/usr/share/OpenMS

WORKDIR patRoon
RUN Rscript -e 'options(patRoon.path.metFragCL = "~/MetFrag2.4.3-CL.jar", patRoon.path.SIRIUS = "~/sirius-linux64-headless-4.0/bin"); devtools::install(); devtools::test()'
