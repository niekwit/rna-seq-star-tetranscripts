## config.yaml
Use config.yaml to set the parameters of the analysis. 

The `time` setting in resources will only be relevant when running the pipeline on an HPC with the `--slurm` option.

## samples.csv

Columns:
`sample`: unique sample name (name of the reads file without `.fastq.gz`)
`genotype`: genotype of the sample
`treatment`: treatment condition of the sample 
`reference`: add `yes` if this sample should be set as a control in the differential transcript analysis
`batch`: batch for each sample (e.g. different sequencing run). If more than one batch is present, a batch effect will be modelled during differential transcript analysis with DESeq2


