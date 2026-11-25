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

## Using 3′ mRNA-seq data (e.g. QuantSeq, 3′ Tag-Seq)

The workflow is designed for full-length RNA-seq. It will run on 3′ mRNA-seq reads, and gene-level differential expression is valid (DESeq2 uses raw counts and does not normalise by gene length), but TE quantification is much more limited.

### Limitations

- **Only 3′ ends are sampled.** Reads come from the few hundred bases upstream of a polyadenylation site. A TE is only counted if it lies near the 3′ end of a polyadenylated transcript, either in a host gene's 3′ UTR or at the polyA signal of an autonomous element (e.g. LTR 3′ ends, the L1 3′ UTR). Transcription across the body of a TE is largely invisible, so TE family counts partly reflect where elements sit relative to polyA sites.
- **Internal oligo-dT priming.** Oligo-dT also primes on genomic A-rich stretches in pre-mRNA and introns. SINEs (Alu, B1/B2) and LINE1 elements end in A-rich tails, so these families pick up spurious, highly multi-mapping reads. A change in host-gene transcription or intron retention between conditions can then appear as differential TE expression.
- **Multi-mapper assignment is biased.** The TEtranscripts EM algorithm assumes reads are spread along the full length of each feature. 3′ tags break that assumption, so the redistribution of multi-mapping reads across families and loci is biased.
- **Low depth and short reads.** Typical 3′ mRNA-seq runs (~5–10 M reads, single-end) give few uniquely mapping TE reads.

### Polyadenylated TE classes

Several TE classes produce polyadenylated transcripts of their own and therefore give genuine 3′ tags. These are the classes for which 3′ mRNA-seq carries the most real signal. DNA transposons are inactive in mammals and have no polyA signal of their own, so any 3′ tags in them come from host transcripts.

#### LINE1

Full-length L1 elements carry their own polyA signal, so polyadenylated L1 transcripts do produce 3′ tags, and family-level comparisons between conditions (e.g. L1 counts in WT vs KO) are reasonable to test. Keep in mind:

- Most genomic L1 copies are 5′-truncated but retain the 3′ end, so 3′ tags cannot distinguish full-length, potentially active elements from truncated copies co-transcribed with host genes.
- L1 subfamilies differ mainly in their 5′ UTR (in mouse, L1Md_A/Tf/Gf are defined by their 5′ monomers), so subfamily assignment from 3′ tags is unreliable.
- Internal priming on the genomic A-tail of an L1 copy within a host transcript gives a tag at exactly the same position as genuine L1 polyadenylation, so the two cannot be told apart.
- Young elements (human L1HS/L1PA2–3, mouse L1Md_Tf/Gf/A) are nearly identical at their 3′ ends, so locus-level calls will mostly come from older, more diverged copies.
- The L1 polyA signal is weak, and read-through transcripts are polyadenylated downstream (3′ transduction). These tags fall in unique flanking sequence and are not counted as L1 by TEtranscripts.

#### ERVs and other LTR retrotransposons

Endogenous retroviruses and related LTR elements (e.g. human HERV-K/HML-2 and HERV-H, mouse IAP, MusD/ETn and MERVL) are transcribed by RNA polymerase II from the 5′ LTR and polyadenylated at a signal in the 3′ LTR. Keep in mind:

- RepeatMasker annotates LTRs and internal regions as separate subfamilies. 3′ tags land in the 3′ LTR, so they are counted under the LTR subfamily (e.g. `LTR5_Hs`, `MT2_Mm`, `IAPLTR1_Mm`), not under the internal region (e.g. `HERVK-int`, `MERVL-int`, `IAPEz-int`). Look at the LTR subfamilies when interpreting ERV results.
- In young elements the 5′ and 3′ LTRs are identical, so each tag maps equally well to both LTRs of the same provirus and becomes a multi-mapper.
- Solo LTRs (left behind by recombination between the two LTRs) are far more common than full proviruses. Many act as promoters or polyA sites for host genes, so an LTR count can reflect a host transcript rather than ERV expression.

#### SVA (human and other hominids)

SVA elements are Pol II-transcribed, carry their own polyA signal at their 3′ end and, like L1, often read through into flanking sequence. They are young, GC-rich and highly similar to each other, so expect mostly multi-mapping reads and limited subfamily resolution.

#### SINEs

- Alu and B1 elements are transcribed by RNA polymerase III and are not normally polyadenylated, but their genome-encoded A-rich tails can still anneal to oligo-dT. Most of their 3′ tags come from copies embedded in the 3′ UTRs and introns of host genes.
- Some SINEs, including mouse B2, contain an AATAAA signal and their Pol III transcripts can be polyadenylated (Borodulina & Kramerov, 2008, *RNA*), so they can give genuine 3′ tags. Copies embedded in host transcripts are still the main confound.

Treat TE results from 3′ mRNA-seq as a screen. Validate key findings with an assay that targets a region away from the 3′ end: for L1, the 5′ UTR or ORF1 (RT-qPCR, ORF1p western blot or immunofluorescence); for ERVs, the internal region (RT-qPCR on `gag`/`env`, or Env/Gag protein). Full-length or long-read RNA-seq covers all classes.

### Suggested settings

- **Pre-trim the reads.** The built-in Trim Galore step removes standard adapters only. Remove poly(A) tails (and, for QuantSeq FWD, the first 12 nt from random priming, as recommended by Lexogen) before placing the reads in `reads/`, for example:

  ```bash
  cutadapt -u 12 -a "A{20}" -a AGATCGGAAGAGC -n 2 -q 20 -m 20 \
      -o reads/{sample}.fastq.gz raw/{sample}.fastq.gz
  ```

  Drop `-u 12` for kits that do not use random priming for second-strand synthesis.
- **`strand`**: set to match the kit. QuantSeq FWD is `forward`, while QuantSeq REV and other read-2-based protocols are usually `reverse`. Check this against the TEcount logs or the kit documentation.
- **`telocal.apply: True`** with **`telocal.extra_params: "--mode uniq"`**: locus-level counts from uniquely mapping reads only. These are easier to defend than family-level EM estimates, and let you check individual loci for internal priming (A-rich genomic sequence immediately downstream) or host-gene read-through.
- **`mapping.outFilterMultimapNmax` / `winAnchorMultimapNmax`**: keep the defaults (`100`) so that multi-mapping TE reads are retained for TEtranscripts.
- **Check host-gene confounding.** For differentially expressed TE families, check whether intronic or host-gene signal changes in the same direction. Prioritise intergenic loci or loci antisense to their host gene.

## Citation

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).
