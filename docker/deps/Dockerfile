# syntax=docker/dockerfile:experimental

FROM rocker/rstudio:4.5

ENV SETUPDIR=/usr/local/setup

# NOTE: we remove system libs in /usr/lib from OpenMS (libc etc, all except libOpenMS*) that shouldn't be installed --> they conflict with RStudio

RUN --mount=type=cache,id=cache_apt,target=/var/cache/apt --mount=type=cache,id=lib_apt,target=/var/lib/apt apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget libnode-dev openbabel \
        libxml2-dev openjdk-17-jdk libmagick++-dev pandoc git texinfo libfribidi-dev ghostscript \
        zlib1g-dev libxml2-dev libnetcdf-dev libglpk-dev tzdata libnetcdf-dev netcdf-bin ccache bzip2 \
        python3-venv python3-pip python3-dev libqt6network6t64 libsvm3 libxerces-c3.2t64 libboost-regex1.74.0 \
        libsnappy1v5 && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://github.com/OpenMS/OpenMS/releases/download/release%2F3.5.0/OpenMS-3.5.0-Debian-Linux-x86_64.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-3.5.0-Debian-Linux-x86_64.deb && \
    rm -rf $SETUPDIR && \
    find /usr/lib -maxdepth 1 -type f ! -name 'libO*' ! -name 'libth*' ! -name 'libC*' -print0 | xargs -0 rm && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    echo 'PATH=/home/rstudio/bin:${PATH}' >> /usr/local/lib/R/etc/Renviron && R CMD javareconf

USER rstudio
WORKDIR /home/rstudio

COPY --chown=rstudio docker/install_deps.R DESCRIPTION .
COPY --chown=rstudio docker/rstudio/ccache.conf /home/rstudio/.ccache/
COPY --chown=rstudio docker/rstudio/Makevars /home/rstudio/.R/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0 R_MAX_NUM_DLLS=250

RUN  echo 'options(java.parameters = "-Djava.awt.headless=true")' >> .Rprofile && \
     echo 'options(pkg.build_extra_flags = FALSE)' >> .Rprofile

ARG GIT_BRANCH=master
ENV GIT_BRANCH=${GIT_BRANCH}
RUN --mount=type=cache,id=ccache,target=/home/rstudio/ccache,uid=1000,gid=1000 Rscript install_deps.R && rm -f install_deps.R DESCRIPTION

RUN Rscript -e 'remove.packages("patRoon")' && rm -rf /home/rstudio/.R/Makevars

CMD ["R"]
