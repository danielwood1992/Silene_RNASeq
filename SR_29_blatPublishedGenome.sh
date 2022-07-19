#!/bin/bash --login
#SBATCH --partition=compute
#SBATCH --ntasks=24
#SBATCH --time=0-12:00
#SBATCH -o %A.out.txt
#SBATCH -e %A.err.txt

#Based on SR_5
#So will use the contigs to blast against the published assembly. Keep the ones that blast against this. Hopefully that will recover more of the decent stuff.

#Downloaded the reference sequence from genbank
module load parallel;
query="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_29/SU_RedSamp_k25mincov2.fasta.trans_filtered";
list="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_29/list.txt";
prog="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_29/SR_29_Progress.txt";

#1 For published genome
outcode="SR29";
ref="/home/b.bssc1d/GCA_018983105.1_ASM1898310v1_genomic.fna";

#2 For the high quality unpublished genome
#outcode="SR29.B";
#ref="/scratch/b.bssc1d/Linkage_Mapping/TGS_GC_fmlrc.scaff_seqs.fa";

parallel -N 1 -j 24 "/home/b.bssc1d/blat/blat/blat $ref {1} -t=dnax -q=rnax {1}.${outcode}.psl && echo {1} Complete >> $prog" :::: $list;

#parallel -N 1 -j 1 "echo /home/b.bssc1d/blat/blat/blat $ref {1} -t=dnax -q=rnax {1}.SR29.psl && echo {1} Complete >> $prog" :::: $list;
