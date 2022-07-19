#!/bin/bash
#SBATCH -o /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR6.%A.out.txt
#SBATCH -e /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR_36.%A.err.txt
#SBATCH --ntasks=1
#SBATCH --time=1-12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8G

vcf="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/merged_snps_filtered.vcf"
transcriptome="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/background_ddTxi3_100622";
#Oh wait we need to filter these by the ones we include in the study.
#And I think also rename the chromosomes. Let's write a quick perl script to do this...
number=$(perl /home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_37_filterrename.pl $transcriptome $vcf);
echo "number $number";
module load R/4.2.0;
module load python-h5py/2.8.0;
sh /home/b.bssc1d/SNPhylo-20180901/snphylo.sh -v $vcf.included.vcf -P $vcf.snpphylo140722 -a 500000 -b;
#Note: to make the bootstrapping work...
#sh /home/b.bssc1d/SNPhylo-20180901/snphylo.sh -v $vcf.included.vcf -P $vcf.snpphylo140722 -t 8 -a 500000;
