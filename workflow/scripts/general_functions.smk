import pandas as pd
import datetime
import os
from scripts.resources import Resources
from scripts import general_functions as utils
from snakemake.utils import min_version, validate


def import_samples():
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"]
    
    # check if sample names match file names
    not_found = []
    for sample in SAMPLES:
        r1= f"reads/{sample}_R1_001.fastq.gz"
        r2= f"reads/{sample}_R2_001.fastq.gz"
        if not os.path.isfile(r1):
            not_found.append(r1)
        if not os.path.isfile(r2):
            not_found.append(r2)
    if len(not_found) != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"Following files not found:\n{not_found}")
        
    return SAMPLES    


def star_arguments(config):
    """
    Returns multimapping arguments for STAR set in config.yaml
    """
    ofmn = config["mapping"]["outFilterMultimapNmax"]
    wamn = config["mapping"]["winAnchorMultimapNmax"]
    extra_params = config["mapping"]["extra_params"]
    star_extra = f"--outFilterMultimapNmax {ofmn} --winAnchorMultimapNmax {wamn} {extra_params}"
    
    return star_extra


def comparisons():
    """
    Create pairwise comparison strings from samples.csv
    """
    sample_info = pd.read_csv("config/samples.csv")
    
    if len(sample_info["treatment"].unique()) == 1:
        sample_info["condition"] = sample_info["genotype"]
    else:
        sample_info["condition"] = sample_info[["genotype","treatment"]].agg('_'.join, axis=1)

    # Get reference conditions
    reference_conditions = sample_info[sample_info["reference"] == "yes"]["condition"].unique().tolist()
    assert len(reference_conditions) > 0, "No reference conditions found"
    
    # Get test conditions
    test_conditions = sample_info["condition"].unique().tolist()
    
    # Create strings for comparisons
    comparisons = []
    for test in test_conditions:
        for ref in reference_conditions:
            if test != ref:
                comparisons.append(f"{test}_vs_{ref}")
    
    return comparisons


def get_chromosomes(genome):
    """
    Return a list of chromosomes for the reference genome
    """
    if "hg" in genome or "T2T" in genome or genome == "test":
        return [i for i in range(1,23)] + ["X", "Y"]
    elif "mm" in genome:
        return [i for i in range(1,20)] + ["X", "Y"]
    else:
        raise ValueError(f"Genome {genome} not found")
    
    
def mapping_input(wildcards):
    """
    Return the input files for mapping based on whether
    spike-in should be applied or not.
    """
    base_idx = f"resources/{genome}_{resources.build}_index_star/"
    
    input_dict = {}
    input_dict["idx"] = base_idx
    
    idx_files = [
        "chrLength.txt", 
        "chrNameLength.txt", 
        "chrName.txt", 
        "chrStart.txt", 
        "exonGeTrInfo.tab", 
        "exonInfo.tab",
        "geneInfo.tab",
        "Genome",
        "genomeParameters.txt",
        "Log.out",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab"
    ]
    idx_files = [f"{base_idx}{i}" for i in idx_files]
    input_dict["idx_files"] = idx_files
    
    if config["spike_in"]["apply"]:
        val1 = "results/spike_in/{sample}_1.fq.gz"
        val2 = "results/spike_in/{sample}_2.fq.gz"
    else:
        val1 = "results/trimmed/{sample}_val_1.fq.gz",
        val2 = "results/trimmed/{sample}_val_2.fq.gz",
        
    input_dict["val1"] = val1
    input_dict["val2"] = val2
    
    return input_dict
