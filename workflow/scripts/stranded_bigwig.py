"""Create plus- and minus-strand bigWig files with deepTools bamCoverage.

Plus-strand signal is positive and minus-strand signal is negative, so both
tracks can be overlaid in a genome browser.

bamCoverage --filterRNAstrand normalises each strand to the reads kept on that
strand, which puts the two strands on different scales. Instead, the scale
factor that bamCoverage computed for the unstranded bigWig (same BAM, same
normalisation, same extra_params, so the same filters and chromosome
exclusions) is reused with --normalizeUsing None --scaleFactor. Both strand
tracks are then on the same scale as each other and as the unstranded bigWig.
"""

import math
import shlex
import subprocess

bam = snakemake.input["bam"]
binsize = str(snakemake.params["binsize"])
extra = shlex.split(snakemake.params["extra"])
strand = snakemake.params["strand"]
threads = str(snakemake.threads)
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

# Scale factor computed by bamCoverage for the unstranded bigWig
with open(snakemake.input["scale"]) as fh:
    scale = float(fh.read().strip())

if not math.isfinite(scale) or scale <= 0:
    raise ValueError(f"Invalid scale factor {scale} in {snakemake.input['scale']}")

with open(log, "w") as log_fh:
    log_fh.write(f"Scale factor from unstranded bamCoverage run: {scale}\n")
    log_fh.flush()

    for direction, sign in (("plus", 1), ("minus", -1)):
        # --binSize comes before extra_params so that an override there
        # matches the unstranded run; the normalisation and strand options
        # come last so that extra_params cannot override them.
        cmd = [
            "bamCoverage",
            "-p",
            threads,
            "-b",
            bam,
            "-o",
            snakemake.output[direction],
            "--binSize",
            binsize,
            *extra,
            "--filterRNAstrand",
            rna_strand[direction],
            "--normalizeUsing",
            "None",
            "--scaleFactor",
            str(sign * scale),
        ]
        log_fh.write(f"\n{shlex.join(cmd)}\n")
        log_fh.flush()
        subprocess.run(cmd, check=True, stdout=log_fh, stderr=log_fh)
