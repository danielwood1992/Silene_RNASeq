#Using some trinity statistics to assess assembly quality.
#See //github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats 
#/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/MA1_TranscriptomeQuality/MA1_2_FilteredContigQC/SU_RedSamp_k25mincov2.fasta.trans_filtered.sub
trinity_exN50="/home/b.bssc1d/trinityrnaseq-v2.10.0/util/misc/contig_ExN50_statistic.pl";

#Full transcriptome
full="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta";
#Reduced transcriptome
red="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta.trans_filtered.sub";
#Kallisto expression matrix
kall_exp="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/MA1_TranscriptomeQuality/SR_27/kallisto.isoform.counts.matrix";

$trinity_exN50 $kall_exp $full > $full.stats;

#Get reduced kallisto matrix...
grep '>' $red | cut -c2- | cut -f1 -d' ' > $red.names;
perl /home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_27.pl $red.names $kall_exp; #output is $kall_exp.red
$trinity_exN50 $kall_exp.red $red > $red.stats;


