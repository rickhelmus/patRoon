FROM rocker/r-rspm:18.04

LABEL maintainer="Rick Helmus <r.helmus@uva.nl>" \
    org.label-schema.name="patRoon" \
    org.label-schema.description="Docker image used for patRoon CI" \
    org.label-schema.schema-version="1.0" \
    org.label-schema.vcs-url="https://github.com/rickhelmus/patRoon" \
    org.label-schema.vendor="patRoon" \
    docker.cmd.test="docker run -t patroonorg/patroon /bin/bash -c 'cd patRoon; Rscript docker/run_tests.R'"

ENV SETUPDIR=/usr/local/setup

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget libv8-3.14-dev openbabel \
        libxml2-dev pngquant openjdk-11-jdk libmagick++-dev pandoc git pngquant texinfo libfribidi-dev \
        zlib1g-dev libxml2-dev libnetcdf-dev libglpk-dev tzdata libnetcdf-dev netcdf-bin && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://abibuilder.informatik.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/2.5.0/OpenMS-2.5.0-Debian-Linux-x86_64.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-2.5.0-Debian-Linux-x86_64.deb && \
    rm -rf $SETUPDIR && \
    useradd -ms /bin/bash patRoon && \
    addgroup patRoon docker && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    R CMD javareconf

USER patRoon
WORKDIR /home/patRoon

COPY ./docker/install_deps.R ./DESCRIPTION ./

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN mkdir deps && cd deps && \
    wget -q http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-CL.jar && \
    wget -q https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.4.29/sirius-4.4.29-linux64-headless.zip && \
    wget -q https://zenodo.org/record/3611238/files/PubChemLite_14Jan2020_tier1.csv && \
    wget -q ftp://newftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/Chemistry_Dashboard/MetFrag_metadata_files/CompTox_17March2019_SelectMetaData.csv && \
    unzip sirius-4.4.29-linux64-headless.zip && rm sirius-4.4.29-linux64-headless.zip && cd ~ && \
    echo 'options(patRoon.path.MetFragCL = "/home/patRoon/deps/MetFrag2.4.5-CL.jar")' >> .Rprofile && \
    echo 'options(patRoon.path.SIRIUS = "/home/patRoon/deps/sirius-linux64-headless-4.4.29/bin")' >> .Rprofile && \
    echo 'options(patRoon.path.MetFragCompTox = "/home/patRoon/deps/CompTox_17March2019_SelectMetaData.csv")' >> .Rprofile && \
    echo 'options(patRoon.path.MetFragPubChemLite = "/home/patRoon/deps/PubChemLite_14Jan2020_tier1.csv")' >> .Rprofile && \
    echo 'options(patRoon.progress.opts = list(style = 1))' >> .Rprofile && \
    Rscript install_deps.R && rm -f ~/install_deps.R ~/DESCRIPTION

ADD --chown=patRoon . patRoon

RUN Rscript -e 'devtools::install(pkg = "patRoon", upgrade = FALSE)'

ENV OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0 R_MAX_NUM_DLLS=150

CMD ["R"]
