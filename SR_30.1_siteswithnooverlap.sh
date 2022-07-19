#OK so what we need to do here is...
#a) sort by scaffold and then b) sort by hit start position, right?
#Then go through them and then only keep hits where nothing is overlapping...
#Right so these are the BLAT hits to the published uniflora genome.
file="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_29/sorted_merged_outfile_best.psl";
#So these are the genes that are actually used
expr_names="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/background_ddTxi3_100622.clean";
#grep -v -f $expr_names $file > $file.expr; #This is now done. 30/06
file=$file.expr
module load bedtools;
#1 Add total length
perl /home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_30.0.pl $file;
#2 Rejig so the "start" position is the first one along the chromosome
perl /home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_30.1.pl $file.filt;
#3 Remove header | sort by start then end | put unique ones (?) into new file
tail -n +6 $file.filt.re | sort -k14,14 -k16,16n | uniq > $file.filt.re.sorted;
#4 prints out the sites which have only a single contig mapping to them...
perl /home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_30.2.pl $file.filt.re.sorted;
#For those regions with only a single contig mapping to them: how many contigs are shared? How many don't have a hit on a different scaffold?
perl /home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_30.3.pl $file.filt.re.sorted.not;

