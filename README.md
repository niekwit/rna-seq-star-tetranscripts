# Snakemake workflow: `rna-seq-star-tetranscripts`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/niekwit/rna-seq-star-tetranscripts/workflows/Tests/badge.svg?branch=main)](https://github.com/niekwit/rna-seq-star-tetranscripts/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for transposable element RNA-Seq using [TEtranscripts](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts).


## Usage

### Install Snakemake and Snakedeploy

Snakemake and Snakedeploy are best installed via the Mamba package manager (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via Mambaforge. For other options see here.

Given that Mamba is installed, run

```
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```
to install both Snakemake and Snakedeploy in an isolated environment. For all following commands ensure that this environment is activated via

```
conda activate snakemake
```
### Deploy workflow

Given that Snakemake and Snakedeploy are installed and available (first step), the workflow can be deployed as follows.

First, create an appropriate project working directory on your system and enter it:

```
mkdir -p path/to/project-workdir
cd path/to/project-workdir
```

In all following steps, we will assume that you are inside of that directory.
Second, run

```
snakedeploy deploy-workflow https://github.com/niekwit/rna-seq-star-tetranscripts . 

```
Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a Snakemake module, the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main Snakefile in the workflow subfolder.
Third, consider to put this directory under version control, e.g. by managing it via a (private) Github repository

### Configure workflow

**General configuration**

To configure this workflow, modify *config/config.yaml* according to your needs, following the explanations provided in the file.

**Sample sheet**

The default sample sheet is *config/samples.csv*. Each sample refers to an actual physical sample, and replicates (both biological and technical) may be specified as separate samples. For each sample, you will always have to specify a sample name in the *sample* columm. 

**Strandedness of library preparation protocol**

To get the correct data from TEcount output, you can provide information on the strandedness of the library preparation protocol used for an experiment. More information on this can be found [here](https://pypi.org/project/TEtranscripts/). 

### Run workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

Fow running the workflow while deploying any necessary software via conda (using the Mamba package manager by default), run Snakemake with

```
snakemake --cores all --use-conda 
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module that has been defined by the deployment in Deploy workflow step.

For further options, e.g. for cluster and cloud execution, see the [docs](https://snakemake.readthedocs.io/).

### Generate report

After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser via

```
snakemake --report report.zip
```

The resulting report.zip file can be passed on to collaborators, provided as a supplementary file in publications, or uploaded to a service like [Zenodo](https://zenodo.org/) in order to obtain a citable [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier).

## Citation

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

