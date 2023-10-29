# syntax=docker/dockerfile:experimental

FROM uva-hva.gitlab.host:4567/r.helmus/patroon/patroondeps

ADD --chown=rstudio . patRoon
RUN --mount=type=cache,id=ccache,target=/home/rstudio/ccache,uid=1000,gid=1000 Rscript -e 'remotes::install_local("patRoon", upgrade = FALSE)'
