###Downloaded test bam data and convert to fastq.gz:
#https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_PE/
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_PE/testData_control_rep1_PE.bam
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_PE/testData_control_rep2_PE.bam
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_PE/testData_treatment_rep1_PE.bam
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/test_data/testdata_PE/testData_treatment_rep2_PE.bam

# Convert bam to fastq
samtools fastq -N -F 4 testData_control_rep1_PE.bam -1 control_1_R1_001.fastq -2 control_1_R2_001.fastq
samtools fastq -N -F 4 testData_control_rep2_PE.bam -1 control_2_R1_001.fastq -2 control_2_R2_001.fastq
samtools fastq -N -F 4 testData_treatment_rep1_PE.bam -1 treatment_1_R1_001.fastq -2 treatment_1_R2_001.fastq
samtools fastq -N -F 4 testData_treatment_rep2_PE.bam -1 treatment_2_R1_001.fastq -2 treatment_2_R2_001.fastq

# Compress fastq files
for FILE in *.fastq; do pigz -f $FILE; done


### Prepare GTF files for Chr22 only
# This will speed up the test run

## TE GTF was downloaded from Hammell lab Dropbox
pigz -d GRCh38_Ensembl_rmsk_TE.gtf.gz

# Keep on Chr22 related entries
grep "^22" GRCh38_Ensembl_rmsk_TE.gtf > GRCh38_Ensembl_rmsk_TE_chr22.gtf  
pigz GRCh38_Ensembl_rmsk_TE_chr22.gtf

## Gene GTF was downloaded from Ensembl
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
pigz -d Homo_sapiens.GRCh38.111.gtf.gz

# Keep only Chr22 related entries
grep "^22" Homo_sapiens.GRCh38.111.gtf > Homo_sapiens.GRCh38.111_chr22gtf
pigz Homo_sapiens.GRCh38.111_chr22gtf