FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="83c68a2161fc26211943185ce48005234d9ea6f6f85ca72c3d696c30ab82c0af"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.5.1/bio/fastqc/environment.yaml
#   prefix: /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.5.1/bio/fastqc/environment.yaml /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.5.1/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.21
RUN mkdir -p /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.5.1/bio/samtools/index/environment.yaml /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yml
#   prefix: /conda-envs/c3af500606597b948946bf39685e5a98
#   name: deeptools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - deeptools=4.1.0
RUN mkdir -p /conda-envs/c3af500606597b948946bf39685e5a98
COPY workflow/envs/deeptools.yml /conda-envs/c3af500606597b948946bf39685e5a98/environment.yaml

# Conda environment:
#   source: workflow/envs/deseq2.yml
#   prefix: /conda-envs/d5cd58d499e48dfcf2ddba99fbbc70ed
#   name: deseq2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - r-tidyverse=2.0.0
#     - r-rcolorbrewer=1.1_3
#     - r-pheatmap=1.0.12
#     - r-ggrepel=0.9.6
#     - r-cowplot=1.1.3
#     - bioconductor-deseq2=1.42.0
#     - bioconductor-genomicfeatures=1.54.1
#     - bioconductor-tximport=1.30
#     - r-reshape2=1.4.4
#     - bioconductor-limma=3.58.1
#     - conda-forge::curl=8.20
#   #  - conda-forge::icu=75.1
RUN mkdir -p /conda-envs/d5cd58d499e48dfcf2ddba99fbbc70ed
COPY workflow/envs/deseq2.yml /conda-envs/d5cd58d499e48dfcf2ddba99fbbc70ed/environment.yaml

# Conda environment:
#   source: workflow/envs/mapping.yml
#   prefix: /conda-envs/8bda2bb25358eea446fa54587cd1f747
#   name: mapping
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - star=2.7.11a
#     - pigz=2.8
#     - multiqc=1.19
#     - curl
RUN mkdir -p /conda-envs/8bda2bb25358eea446fa54587cd1f747
COPY workflow/envs/mapping.yml /conda-envs/8bda2bb25358eea446fa54587cd1f747/environment.yaml

# Conda environment:
#   source: workflow/envs/rtracklayer.yml
#   prefix: /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582
#   name: rtracklayer
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - conda-forge::r-tidyverse=2.0.0
#     - bioconductor-rtracklayer=1.62.0
RUN mkdir -p /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582
COPY workflow/envs/rtracklayer.yml /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582/environment.yaml

# Conda environment:
#   source: workflow/envs/te.yml
#   prefix: /conda-envs/c4b565a9ee40a6be82be171208e594cd
#   name: tetranscripts
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - tetranscripts=2.2.3
#     - pip
#     - pip:
#         - TElocal==1.1.3
RUN mkdir -p /conda-envs/c4b565a9ee40a6be82be171208e594cd
COPY workflow/envs/te.yml /conda-envs/c4b565a9ee40a6be82be171208e594cd/environment.yaml

# Conda environment:
#   source: workflow/envs/trim_galore.yml
#   prefix: /conda-envs/4df970466a873ab6bad20311e5c763da
#   name: trim_galore
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     # Pinned to this exact build: build h9ee0642_0 was linked against a newer
#     # glibc than older hosts (e.g. Ubuntu 20.04, glibc 2.31) provide, causing
#     # `GLIBC_2.32/2.33/2.34 not found` errors. hf1b6044 builds fix this (glibc <=2.16).
#     # >=2.6.0 needed: 2.2.0 silently ignored -a2/--adapter2 unless -a was also given.
#     - trim-galore=2.6.0=hf1b6044_1
RUN mkdir -p /conda-envs/4df970466a873ab6bad20311e5c763da
COPY workflow/envs/trim_galore.yml /conda-envs/4df970466a873ab6bad20311e5c763da/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b --file /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b/environment.yaml && \
    conda env create --prefix /conda-envs/c3af500606597b948946bf39685e5a98 --file /conda-envs/c3af500606597b948946bf39685e5a98/environment.yaml && \
    conda env create --prefix /conda-envs/d5cd58d499e48dfcf2ddba99fbbc70ed --file /conda-envs/d5cd58d499e48dfcf2ddba99fbbc70ed/environment.yaml && \
    conda env create --prefix /conda-envs/8bda2bb25358eea446fa54587cd1f747 --file /conda-envs/8bda2bb25358eea446fa54587cd1f747/environment.yaml && \
    conda env create --prefix /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582 --file /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582/environment.yaml && \
    conda env create --prefix /conda-envs/c4b565a9ee40a6be82be171208e594cd --file /conda-envs/c4b565a9ee40a6be82be171208e594cd/environment.yaml && \
    conda env create --prefix /conda-envs/4df970466a873ab6bad20311e5c763da --file /conda-envs/4df970466a873ab6bad20311e5c763da/environment.yaml && \
    conda clean --all -y