FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:02ab-main

SHELL ["/usr/bin/env", "bash", "-c"]


# Allow --mount=cache to do its job
# https://github.com/moby/buildkit/blob/86c33b66e176a6fc74b88d6f46798d3ec18e2e73/frontend/dockerfile/docs/syntax.md#run---mounttypecache
RUN rm /etc/apt/apt.conf.d/docker-clean
RUN echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache


# Generic installation dependencies
# wget - obvious
# software-properties-common - `add-apt-repository`
# dirmngr - GPG key manager, loads from `/etc/apt/trusted.gpg.d/`
RUN apt-get update && apt-get install --yes --no-install-recommends \
  wget \
  software-properties-common \
  dirmngr


#
# R
#

# >>> Install R
# https://cloud.r-project.org/bin/linux/debian/
# https://github.com/rocker-org/rocker-versioned2/blob/f3325b2cf88d8899ddcb2f0945aa9f87ad150cd7/scripts/install_R_ppa.sh
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/debian $(lsb_release --codename --short)-cran40/"

RUN apt-get update && apt-get install --yes \
    r-base \
    r-base-dev \
    locales

RUN apt-mark hold r-base r-base-dev

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen
RUN locale-gen en_US.utf8
RUN /usr/sbin/update-locale LANG="en_US.UTF-8"

# >>> R packages
RUN apt-get update
RUN apt install --yes \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(version = "3.15")'
RUN Rscript -e 'BiocManager::install(c( \
  "purrr", \
  "dplyr", \
  "tibble", \
  "readr", \
  "readxl", \
  "stringr", \
  "vctrs", \
  "clusterProfiler", \
  "DOSE", \
  "ggridges", \
  "enrichplot", \
  "ggplot2", \
  "msigdbr", \
  "pathview" \
  ), update=FALSE)'

RUN Rscript -e 'BiocManager::install(c( \
  "org.Hs.eg.db", \
  "org.Mm.eg.db", \
  "org.Rn.eg.db", \
  "org.Dm.eg.db", \
  "org.At.tair.db", \
  "org.Sc.sgd.db", \
  "org.Dr.eg.db", \
  "org.Ce.eg.db", \
  "org.Bt.eg.db", \
  "org.Ss.eg.db", \
  "org.Gg.eg.db", \
  "org.Mmu.eg.db", \
  "org.Cf.eg.db", \
  "org.EcK12.eg.db", \
  "org.Xl.eg.db", \
  "org.Ag.eg.db", \
  "org.Pt.eg.db", \
  "org.EcSakai.eg.db", \
  "org.Mxanthus.db" \
  ), update=FALSE)'

RUN python3 -m pip install --upgrade latch imagesize jinja2

# todo: automatically build frontend in Dockerfile

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY ./report/build/index.html /root/template.html
COPY ./go_pathway.r /root/go_pathway.r
COPY wf /root/wf

ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
