# syntax=docker/dockerfile:experimental

FROM rocker/rstudio:4.3

ENV SETUPDIR=/usr/local/setup

# NOTE: we remove system libs from OpenMS (libc etc) that shouldn't be installed --> they conflict with RStudio

RUN --mount=type=cache,id=cache_apt,target=/var/cache/apt --mount=type=cache,id=lib_apt,target=/var/lib/apt apt-get update -y && \
    apt-get install -y --no-install-recommends libssl-dev libssh2-1-dev wget libnode-dev openbabel \
        libxml2-dev openjdk-11-jdk libmagick++-dev pandoc git pngquant texinfo libfribidi-dev ghostscript \
        zlib1g-dev libxml2-dev libnetcdf-dev libglpk-dev tzdata libnetcdf-dev netcdf-bin ccache bzip2 && \
    mkdir -p $SETUPDIR && \
    wget -P $SETUPDIR https://github.com/OpenMS/OpenMS/releases/download/Release3.0.0/OpenMS-3.0.0-Debian-Linux-x86_64.deb && \
    apt-get install -y --no-install-recommends $SETUPDIR/OpenMS-3.0.0-Debian-Linux-x86_64.deb && \
    rm -rf $SETUPDIR && \
    rm /usr/lib/libc.so.6 /usr/lib/libm.so.6 /usr/lib/libpthread.so.0 /usr/lib/libstdc++.so.6 /usr/lib/libglib-2.0.so.0 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    echo 'PATH=/home/rstudio/bin:${PATH}' >> /usr/local/lib/R/etc/Renviron && R CMD javareconf

USER rstudio
WORKDIR /home/rstudio

COPY --chown=rstudio docker/install_deps.R DESCRIPTION .
COPY --chown=rstudio docker/rstudio/ccache.conf /home/rstudio/.ccache/
COPY --chown=rstudio docker/rstudio/Makevars /home/rstudio/.R/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true OPENMS_DATA_PATH=/usr/share/OpenMS _R_CHECK_FORCE_SUGGESTS_=0 R_MAX_NUM_DLLS=250

RUN  echo 'options(java.parameters = "-Djava.awt.headless=true")' >> .Rprofile

ARG GIT_BRANCH=master
ENV GIT_BRANCH=${GIT_BRANCH}
RUN --mount=type=cache,id=ccache,target=/home/rstudio/ccache,uid=1000,gid=1000 Rscript install_deps.R && rm -f install_deps.R DESCRIPTION

RUN rm -rf /home/rstudio/.R/Makevars

CMD ["R"]
