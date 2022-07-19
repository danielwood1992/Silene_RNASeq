#!/bin/bash --login
#SBATCH --partition=compute
#SBATCH --ntasks=24
#SBATCH --time=0-12:00
#SBATCH -o %A.out.txt
#SBATCH -e %A.err.txt

my_dir="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_29";

#Assumes separate blat dirs, with query files split into multiple files xaa xab xac etc...
#Run from in the directory where everuthing is stored
#Note: SR_29.B was from Owen's genome, don't use this.

#cat $my_dir/x*psl > $my_dir/sorted_merged_outfile.psl #Need to figure out how to delete the headers, but can't be bothered so do it manually
#truncate -s 0 $my_dir/sorted_merged_outfile.psl;

head -n 5 $my_dir/xaa.SR29.psl > $my_dir/sorted_merged_outfile.psl;
#head -n 5 $my_dir/xaa.SR29.B.psl > $my_dir/SR29B.sorted_merged_outfile.psl;

for file in $my_dir/x*psl; do tail -n +6 $file >> $my_dir/sorted_merged_outfile.psl; done
#for file in $my_dir/x*SR29.B*psl; do tail -n +6 $file >> $my_dir/SR29B.sorted_merged_outfile.psl; done

mkdir $my_dir/tmpdir_1 
/home/b.bssc1d/blat/pslSort dirs $my_dir/sorted_merged_outfile.psl $my_dir/tmpdir_1 .
/home/b.bssc1d/blat/pslReps $my_dir/sorted_merged_outfile.psl $my_dir/sorted_merged_outfile_best.psl $my_dir/sorted_merged_outfile_best.psr;

/home/b.bssc1d/blat/pslScore $my_dir/sorted_merged_outfile_best.psl > $my_dir/sorted_merged_outfile_best.scores;

#/home/b.bssc1d/blat/pslSort dirs $my_dir/SR29B.sorted_merged_outfile.psl $my_dir/tmpdir_1 .
#/home/b.bssc1d/blat/pslReps $my_dir/SR29B.sorted_merged_outfile.psl $my_dir/SR29B.sorted_merged_outfile_best.psl $my_dir/SR29B.sorted_merged_outfile_best.psr;

#/home/b.bssc1d/blat/pslScore $my_dir/SR29B.sorted_merged_outfile_best.psl > $my_dir/SR29B.sorted_merged_outfile_best.scores;

