import os

class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build
    """
    
    # create genome directory
    os.makedirs("resources/", exist_ok=True)
    
    def __init__(self, genome, build):
        self.genome = genome
        self.build = build
                
        # base URLs
        base_url_ens = f"https://ftp.ensembl.org/pub/release-{build}/"
        base_url_te = "https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/"
                
        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"
            elif genome == "hg38":
                name = "GRCh38"
            else:
                raise ValueError(f"Genome {genome} not found")
                
            # create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
            self.tegtf_url = f"{base_url_te}{name}_Ensembl_rmsk_TE.gtf.gz"
        
        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"
            elif genome == "mm39":
                name = "GRCm39"
            else:
                raise ValueError(f"Genome {genome} not found")
                
            # create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
            self.tegtf_url = f"{base_url_te}{name}_Ensembl_rmsk_TE.gtf.gz"
        
        elif genome == "test":
            # Download fasta for only one chromosome            
            self.fasta_url = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
            self.gtf_url = "https://github.com/niekwit/rna-seq-star-tetranscripts/raw/main/.test/Homo_sapiens.GRCh38.111_chr22gtf.gz"
            self.tegtf_url = "https://github.com/niekwit/rna-seq-star-tetranscripts/raw/main/.test/GRCh38_Ensembl_rmsk_TE_chr22.gtf.gz"
            
        else:
            raise ValueError(f"Genome {genome} not found/available")
                        
        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        self.tegtf = self._file_from_url(self.tegtf_url)
        
    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file
        """
        
        return f"resources/{os.path.basename(url).replace('.gz','')}"
    
    def _return_urls(self, name, build):
        """Returns URLs for fasta and GTF files
        """
        pass
  
        
        
            
            
    