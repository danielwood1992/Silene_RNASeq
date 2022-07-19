#!/bin/bash --login
#SBATCH -o batch_LM_4.5_.%A.out.txt
#SBATCH -e batch_LM_4.5_.%A.err.txt
#SBATCH --ntasks=1
#SBATCH --time=1-12:00:00
#SBATCH --partition=compute

thing_list="/home/b.bssc1d/scripts/formal_RNAseq_scripts/samples_file.txt";
#thing_list="/home/b.bssc1d/scripts/formal_RNAseq_scripts/samples_file.txt_B";

ref="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta.trans_filtered.longesti";
out="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree";

module load samtools/1.10
module load hisat2/2.1.0-binary;
module load parallel;
parallel --colsep "\t" -j 1 " hisat2 -1 {3} -2 {4} -x $ref | samtools view -bSq 20 -o $out/{2}.hisat2.bam - && echo Done {2} >> $out/summary.txt" :::: $thing_list;
