import pandas as pd
import os

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
        raise ValueError(f"ERROR: following files not found:\n{not_found}")
        
    return SAMPLES    


def star_arguments(config):
    """Returns multimapping arguments for STAR set in config.yaml
    """
    ofmn = config["mapping"]["outFilterMultimapNmax"]
    wamn = config["mapping"]["winAnchorMultimapNmax"]
    extra_params = config["mapping"]["extra_params"]
    star_extra = f"--outFilterMultimapNmax {ofmn} --winAnchorMultimapNmax {wamn} {extra_params}"
    
    return star_extra
 