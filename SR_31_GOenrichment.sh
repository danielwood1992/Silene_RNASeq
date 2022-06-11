#GO enrichment script...
#Temp directory - copy the lists across

dir="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO";
list_of_files="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/list_of_files.txt";
sub1="/home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR31_1.0_cleanfromR.sh";
sh $sub1 $list_of_files;

#Ok so now we have a nice file list.
#What do we need to do?

#1 Get sequence lengths
goseq="/home/b.bssc1d/trinityrnaseq-v2.10.0/Analysis/DifferentialExpression/run_GOseq.pl";
fasta="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta.trans_filtered.longesti";
sub2="/home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR31_2.0_lengths.sh";
sh $sub2 $fasta;
annotations="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/go_annotations.txt"; #recovered from previous runs on Hal

#so then...
background=$(head -n 1 $list_of_files).clean;
echo $background;
lengths=$fasta.lengths;
tail -n +2 $list_of_files > $list_of_files.todo.txt;
module load R/4.2.0;
while read gene_list; do perl $goseq --genes_single_factor ${gene_list}.clean --GO_assignments $annotations --lengths $lengths --background $background; done < $list_of_files.todo.txt; 

