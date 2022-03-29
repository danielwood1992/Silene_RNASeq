#!/bin/bash --login
#SBATCH --partition=highmem
#SBATCH --ntasks=3
#SBATCH --time=1-12:00
#SBATCH --mem=170G
#SBATCH -o error_messages/2.2_%A.out.txt
#SBATCH -e error_messages/2.2_%A.err.txt

#Returns transcripts with >100AA ORFs 
#based on /home/b.bssc1d/scripts/rna_seq_paper/2B_2.2_TransDecoder.sh

#Run this in the working directory...
dir="/scratch/b.bssc1d/Silene_RNA_seq/2_TranscriptomeCleanup/";
fasta="SU_RedSamp_k25mincov2.fasta";

mkdir $dir/Transdecoder_dir;
 ~/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $dir/$fasta -O $dir/Transdecoder_dir;
perl ~/scripts/formal_RNAseq_scripts/sub_2B_2.2_TransDecoder.pl $dir/$fasta $dir/Transdecoder_dir/longest_orfs.pep;
