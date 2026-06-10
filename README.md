[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10211476.svg)](https://doi.org/10.5281/zenodo.10211476)
[![Tests](https://github.com/niekwit/rna-seq-star-tetranscripts/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/rna-seq-star-tetranscripts/actions/workflows/main.yml)

# Snakemake workflow: `rna-seq-star-tetranscripts`

A Snakemake workflow for transposable element RNA-Seq using [TEtranscripts](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts).

## Usage

### Install Snakemake

Snakemake is best installed via [Conda](https://docs.conda.io/en/latest/)/[Mamba](https://github.com/mamba-org/mamba):

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake>=8.25.5
conda activate snakemake
```

### Deploy workflow

Clone the repository:

```bash
git clone https://github.com/niekwit/rna-seq-star-tetranscripts.git
cd rna-seq-star-tetranscripts
```

Place your raw reads (`.fastq.gz`) in the `reads/` directory. For single-end data, files should be named `{sample}.fastq.gz`. For paired-end data, files should follow the pattern `{sample}_R1_001.fastq.gz` / `{sample}_R2_001.fastq.gz`.

### Configure workflow

#### `config/config.yaml`

| Parameter                       | Description                                                                                                                     |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `genome`                        | Genome build (`mm39`, `hg38`, `hg19`, `mm38`, `T2T-CHM13v2.0`)                                                                  |
| `ensembl_genome_build`          | Ensembl release number (e.g. `113`)                                                                                             |
| `strand`                        | Library strandedness: `forward`, `reverse`, or `unstranded` (see [TEtranscripts docs](https://pypi.org/project/TEtranscripts/)) |
| `fdr_cutoff`                    | FDR threshold for volcano plots (default: `0.05`)                                                                               |
| `fc_cutoff`                     | log₂ fold change threshold for volcano plots (default: `0.5`)                                                                   |
| `spike_in.apply`                | Set to `True` to enable spike-in normalisation                                                                                  |
| `spike_in.name`                 | Pattern matching spike-in feature names (e.g. `ERCC-`)                                                                          |
| `spike_in.gtf`                  | Path to spike-in GTF file                                                                                                       |
| `spike_in.fasta`                | Path to spike-in FASTA file                                                                                                     |
| `mapping.outFilterMultimapNmax` | Max number of multi-mapping loci allowed (default: `100`)                                                                       |
| `mapping.winAnchorMultimapNmax` | Max multi-mapping loci for anchor windows (default: `100`)                                                                      |
| `mapping.extra_params`          | Any additional STAR arguments                                                                                                   |
| `split_bam`                     | Split BAMs by chromosome before TE counting — recommended for large datasets on a cluster                                       |
| `deeptools.binsize`             | Bin size (bp) for bigWig generation                                                                                             |
| `deeptools.normalisation`       | Normalisation method: `RPKM`, `FPKM`, or `TPM`                                                                                  |

#### `config/samples.csv`

A comma-separated file describing the samples:

```
sample,genotype,treatment,reference,batch
WT_1,WT,None,yes,1
KO_1,KO,None,no,1
```

| Column      | Description                                                                               |
| ----------- | ----------------------------------------------------------------------------------------- |
| `sample`    | Unique sample name matching the reads filename (without `.fastq.gz`)                      |
| `genotype`  | Genotype of the sample                                                                    |
| `treatment` | Treatment condition                                                                       |
| `reference` | Set to `yes` for control samples used as the reference in differential expression         |
| `batch`     | Sequencing batch. If more than one batch is present, a batch effect is modelled in DESeq2 |

### Run workflow

**Locally** (Conda environments resolved automatically):

```bash
snakemake --cores <N> --use-conda
```

**On an HPC with Slurm**:

```bash
snakemake --cores <N> --use-conda --slurm --default-resources slurm_account=<account> slurm_partition=<partition>
```

**With Apptainer** (using the pre-built Docker image):

```bash
snakemake --cores <N> --use-apptainer
```

To perform a dry run first:

```bash
snakemake -n
```

## Output

Results are written to the `results/` directory:

| Path                | Description                                                                    |
| ------------------- | ------------------------------------------------------------------------------ |
| `results/qc/`       | FastQC and MultiQC reports                                                     |
| `results/trimmed/`  | Adapter-trimmed reads                                                          |
| `results/mapped/`   | STAR-aligned BAM files and per-sample mapping stats                            |
| `results/te_count/` | TEcount count tables (genes + TEs)                                             |
| `results/bigwig/`   | RPKM-normalised bigWig tracks                                                  |
| `results/deseq2/`   | DESeq2 results CSVs for each comparison (genes and TEs separately)             |
| `results/plots/`    | Volcano plots, PCA, sample distance heatmap, mapping rates, and TE class plots |

## Citation

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).
