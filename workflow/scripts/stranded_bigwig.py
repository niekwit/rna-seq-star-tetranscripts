"""Create plus- and minus-strand bigWig files with deepTools bamCoverage.

Plus-strand signal is positive and minus-strand signal is negative, so both
tracks can be overlaid in a genome browser.

bamCoverage --filterRNAstrand normalises each strand to the reads kept on that
strand, which puts the two strands on different scales. Instead, the scale
factor is computed here from all mapped alignments, using the same formula and
read total (sum of the BAM index mapped counts) that bamCoverage uses for the
unstranded bigWig, and passed with --normalizeUsing None --scaleFactor. Both
strand tracks are then on the same scale as each other and as the unstranded
bigWig.
"""

import shlex
import subprocess

import pysam

bam = snakemake.input["bam"]
binsize = int(snakemake.params["binsize"])
norm = snakemake.params["norm"]
extra = shlex.split(snakemake.params["extra"])
strand = snakemake.params["strand"]
threads = snakemake.threads
log = snakemake.log[0]

# --filterRNAstrand assumes a dUTP (reverse-stranded) library, where
# "forward" keeps reads from genes on the plus strand. For forward-stranded
# libraries the two values swap.
if strand == "reverse":
    rna_strand = {"plus": "forward", "minus": "reverse"}
elif strand == "forward":
    rna_strand = {"plus": "reverse", "minus": "forward"}
else:
    raise ValueError(
        f"Stranded bigWig files need a stranded library, but strand is '{strand}'"
    )

# Total mapped alignments, as counted by bamCoverage for normalisation
with pysam.AlignmentFile(bam) as fh:
    mapped = sum(s.mapped for s in fh.get_index_statistics())

if mapped == 0:
    raise ValueError(f"No mapped reads in {bam}")

# Same formulas as bamCoverage --normalizeUsing RPKM/CPM
if norm == "RPKM":
    scale = 1e6 / mapped * 1000 / binsize
elif norm == "CPM":
    scale = 1e6 / mapped
elif norm == "None":
    scale = 1.0
else:
    raise ValueError(
        f"Normalisation '{norm}' is not supported for stranded bigWig files "
        "(use RPKM, CPM or None)"
    )

with open(log, "w") as log_fh:
    log_fh.write(f"Mapped alignments: {mapped}\nScale factor ({norm}): {scale}\n")
    log_fh.flush()

    for direction, sign in (("plus", 1), ("minus", -1)):
        cmd = [
            "bamCoverage",
            "-p",
            str(threads),
            "-b",
            bam,
            "-o",
            snakemake.output[direction],
            "--binSize",
            str(binsize),
            "--filterRNAstrand",
            rna_strand[direction],
            "--normalizeUsing",
            "None",
            "--scaleFactor",
            str(sign * scale),
            *extra,
        ]
        log_fh.write(f"\n{shlex.join(cmd)}\n")
        log_fh.flush()
        subprocess.run(cmd, check=True, stdout=log_fh, stderr=log_fh)
