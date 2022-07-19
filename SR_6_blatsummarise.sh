#Assumes separate blat dirs, with query files split into multiple files xaa xab xac etc...
#Run from in the directory where everuthing is stored

#my_dir="/home/dwood/blat_Slatifolia/RNAX_DNAX/SUniflora";
my_dir="/home/dwood/blat_oldsilenegenome";
cat my_dir/x*psl > $my_dir/sorted_merged_outfile.psl
mkdir $my_dir/tmpdir_1 
/home/dwood/bin/x86_64/pslSort dirs $my_dir/sorted_merged_outfile.psl $my_dir/tmpdir_1 .
/home/dwood/bin/x86_64/pslReps $my_dir/sorted_merged_outfile.psl $my_dir/sorted_merged_outfile_best.psl $my_dir/sorted_merged_outfile_best.psr;
perl ~/scripts/pslScore.pl $my_dir/sorted_merged_outfile_best.psl > $my_dir/sorted_merged_outfile_best.scores;

SU_RedSamp_k25mincov2.fasta.trans_filtered.cds
merged_outfiles.psl
#/home/dwood/blat_Slatifolia/RNAX_DNAX/SUniflora

