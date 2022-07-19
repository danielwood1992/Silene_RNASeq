#!/bin/sh
#SBATCH -o /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR_35.%A.out.txt
#SBATCH -e /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR_35.%A.err.txt
#SBATCH --ntasks=12
#SBATCH --time=1-12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8G

module load parallel;
module load bcftools;
module load samtools;
#Looks to be a bit of a pain running this in parallel, so I will just do all the jobs in parallel instead...
genome="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta.trans_filtered.longesti";
bam_list="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/bam_list_incomplete.txt";
#bam_test="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/GR10.sorted.bam";

#So you need to do this, then maybe filter for QUAL later. 
parallel -N 1 -j 12 --delay 0.2 "samtools sort {1} -o {1}.sorted.bam && bcftools mpileup -Ou --gvcf 20 -d 2000 -f $genome {1}.sorted.bam | bcftools call -Ou -m --gvcf 20 | bcftools norm -m +any --fasta-ref $genome | bcftools +fill-tags -- -t all | bcftools plugin setGT - -- -t q -n . -i \"FMT/DP<6\" | bcftools view - -Ob -o {1}.tomerge.bcf && echo {1} Complete >> $bam_list.progress" :::: $bam_list;


