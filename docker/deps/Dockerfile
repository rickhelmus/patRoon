# syntax=docker/dockerfile:experimental

FROM rocker/rstudio:4.2

ENV SETUPDIR=/usr/local/setup

RUN --mount=type=cache,id=cache_apt,target=/var/cache/apt --mount=type=cache,id=lib_apt,target=/var/lib/apt apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget libnode-dev openbabel \
        libxml2-dev pngquant openjdk-11-jdk libmagick++-dev pandoc git pngquant texinfo libfribidi-dev ghostscript \
        zlib1g-dev libxml2-dev libnetcdf-dev libglpk-dev tzdata libnetcdf-dev netcdf-bin ccache && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://github.com/OpenMS/OpenMS/releases/download/Release2.7.0/OpenMS-2.7.0-Debian-Linux-x86_64.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-2.7.0-Debian-Linux-x86_64.deb && \
    rm -rf $SETUPDIR && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    echo 'PATH=/home/rstudio/bin:${PATH}' >> /usr/local/lib/R/etc/Renviron && R CMD javareconf

USER rstudio
WORKDIR /home/rstudio

COPY --chown=rstudio docker/install_deps.R DESCRIPTION .
COPY --chown=rstudio docker/rstudio/ccache.conf /home/rstudio/.ccache/
COPY --chown=rstudio docker/rstudio/Makevars /home/rstudio/.R/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0 R_MAX_NUM_DLLS=150

RUN  mkdir deps && cd deps && \
    wget --progress=dot:mega https://github.com/ipb-halle/MetFragRelaunched/releases/download/v2.4.8/MetFragCommandLine-2.4.8.jar && \
    wget --progress=dot:mega https://github.com/boecker-lab/sirius/releases/download/v5.5.7/sirius-5.5.7-linux64-headless.zip && \
    wget --progress=dot:mega https://zenodo.org/record/6503754/files/PubChemLite_exposomics_20220429.csv && \
    wget --progress=dot:mega https://zenodo.org/record/3472781/files/CompTox_07March19_WWMetaData.csv && \
    unzip sirius-5.5.7-linux64-headless.zip && rm sirius-5.5.7-linux64-headless.zip && \
    git clone --depth 1 https://bitbucket.org/djoumbou/biotransformer.git && cd biotransformer && \
    wget --progress=dot:mega https://github.com/rickhelmus/patRoonDeps/raw/master/ext/biotransformer-3.0.0.jar && cd ~ && \
    echo 'options(patRoon.path.MetFragCL = "/home/rstudio/deps/MetFragCommandLine-2.4.8.jar")' >> .Rprofile && \
    echo 'options(patRoon.path.SIRIUS = "/home/rstudio/deps/sirius/bin")' >> .Rprofile && \
    echo 'options(patRoon.path.MetFragCompTox = "/home/rstudio/deps/CompTox_07March19_WWMetaData.csv")' >> .Rprofile && \
    echo 'options(patRoon.path.MetFragPubChemLite = "/home/rstudio/deps/PubChemLite_exposomics_20220429.csv")' >> .Rprofile && \
    echo 'options(patRoon.path.BioTransformer = "/home/rstudio/deps/biotransformer/biotransformer-3.0.0.jar")' >> .Rprofile && \
    echo 'options(java.parameters = "-Djava.awt.headless=true")' >> .Rprofile

RUN --mount=type=cache,id=ccache,target=/home/rstudio/ccache,uid=1000,gid=1000 Rscript install_deps.R && rm -f install_deps.R DESCRIPTION

RUN rm -rf /home/rstudio/.R/Makevars

CMD ["R"]
