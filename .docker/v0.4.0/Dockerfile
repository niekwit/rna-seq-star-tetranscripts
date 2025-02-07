FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="5f257ff4258181a416226253990f4afbb50651a383c1ad7c663b82e8589578a2"

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
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.5.1/bio/trim_galore/pe/environment.yaml
#   prefix: /conda-envs/975a4a78cc02737b050b54b2645b9ecc
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/975a4a78cc02737b050b54b2645b9ecc
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.5.1/bio/trim_galore/pe/environment.yaml /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yml
#   prefix: /conda-envs/cb3f51030520fe1abe74523b3e058594
#   name: deeptools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - deeptools=3.5.4
RUN mkdir -p /conda-envs/cb3f51030520fe1abe74523b3e058594
COPY workflow/envs/deeptools.yml /conda-envs/cb3f51030520fe1abe74523b3e058594/environment.yaml

# Conda environment:
#   source: workflow/envs/deseq2.yml
#   prefix: /conda-envs/31db22e8b168012308fe67eb4bf010cf
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
#   #  - conda-forge::icu=75.1
RUN mkdir -p /conda-envs/31db22e8b168012308fe67eb4bf010cf
COPY workflow/envs/deseq2.yml /conda-envs/31db22e8b168012308fe67eb4bf010cf/environment.yaml

# Conda environment:
#   source: workflow/envs/mapping.yml
#   prefix: /conda-envs/3157fa53c5396cf74c280fed3c76bb09
#   name: mapping
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - star=2.7.11a
#     - pigz=2.8
#     - multiqc=1.19
RUN mkdir -p /conda-envs/3157fa53c5396cf74c280fed3c76bb09
COPY workflow/envs/mapping.yml /conda-envs/3157fa53c5396cf74c280fed3c76bb09/environment.yaml

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
#   prefix: /conda-envs/8870bac9f447214976803485b3a9fe25
#   name: tetranscripts
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - tetranscripts=2.2.3
RUN mkdir -p /conda-envs/8870bac9f447214976803485b3a9fe25
COPY workflow/envs/te.yml /conda-envs/8870bac9f447214976803485b3a9fe25/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b --file /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b/environment.yaml && \
    conda env create --prefix /conda-envs/975a4a78cc02737b050b54b2645b9ecc --file /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml && \
    conda env create --prefix /conda-envs/cb3f51030520fe1abe74523b3e058594 --file /conda-envs/cb3f51030520fe1abe74523b3e058594/environment.yaml && \
    conda env create --prefix /conda-envs/31db22e8b168012308fe67eb4bf010cf --file /conda-envs/31db22e8b168012308fe67eb4bf010cf/environment.yaml && \
    conda env create --prefix /conda-envs/3157fa53c5396cf74c280fed3c76bb09 --file /conda-envs/3157fa53c5396cf74c280fed3c76bb09/environment.yaml && \
    conda env create --prefix /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582 --file /conda-envs/ad2411a11f8886f06eaf5ae0af7a2582/environment.yaml && \
    conda env create --prefix /conda-envs/8870bac9f447214976803485b3a9fe25 --file /conda-envs/8870bac9f447214976803485b3a9fe25/environment.yaml && \
    conda clean --all -y