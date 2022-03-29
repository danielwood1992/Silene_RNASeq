#Based on: run_kallisto_trinityutils.sh

#transcripts="/home/dwood/Hal_HDD1/dwood/Corset_Kallisto/SU_RedSamp_k25mincov2.fasta.trans_filtered.clustonly";
#output_dir="/home/dwood/Hal_HDD1/dwood/Corset_Kallisto";
#gene_map="/home/dwood/Hal_HDD1/dwood/Corset_Kallisto/Corset-clusters.txt.mod";
#samples_file="/home/dwood/scripts/RNA_seq_scripts/samples_file.txt";

#06/05/20
#samples_file="/home/dwood/scripts/RNA_seq_scripts/samples_file.txt";
#transcripts="/home/dwood/Hal_HDD1/dwood/Corset_Kallisto/Not_Corset_Filtered/SU_RedSamp_k25mincov2.fasta.trans_filtered";
#output_dir="/home/dwood/Hal_HDD1/dwood/Corset_Kallisto/Not_Corset_Filtered";
#gene_map="/home/dwood/Hal_HDD1/dwood/Corset_Kallisto/Not_Corset_Filtered/Trinity-clusters.txt.mod";

#16/05/20
#samples_file="/home/dwood/scripts/RNA_seq_scripts/samples_file.txt";
#transcripts="/home/dwood/blat_Suniflora/SU_RedSamp_k25mincov2.fasta.trans_filtered.blat";
#output_dir="~/blat_Suniflora";

#18/05/20
#samples_file="/home/dwood/scripts/RNA_seq_scripts/samples_file.txt";
#transcripts="/home/dwood/blastdb_prot/Ref_protein_sets/Transdecoder_dir/blast_results/test_mapping/SU_RedSamp_k25mincov2.fasta.trans_filtered";
#output_dir="/home/dwood/blastdb_prot/Ref_protein_sets/Transdecoder_dir/blast_results/test_mapping";

#FOR ANALYSIS - DEFAULT 22/06/20
samples_file="/home/dwood/scripts/formal_RNAseq_scripts/samples_file.txt";
transcripts="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/SU_RedSamp_k25mincov2.fasta.trans_filtered";
output_dir="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq";
perl /home/dwood/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --seqType fq --transcripts $transcripts --samples_file $samples_file --est_method kallisto --output_dir $output_dir --thread_count 10 --prep_reference;

#TEMP1 - RAW ALIGNMENT
#samples_file="/home/dwood/scripts/formal_RNAseq_scripts/samples_file.txt";
#transcripts="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/SU_RedSamp_k25mincov2.fasta";
#output_dir="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/raw_transcriptome";
#perl /home/dwood/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --seqType fq --transcripts $transcripts --samples_file $samples_file --est_method kallisto --output_dir $output_dir --thread_count 20 --prep_reference;


#TEMP2 - FILTERED ALIGNMENT
#samples_file="/home/dwood/scripts/formal_RNAseq_scripts/samples_file.txt";
#transcripts="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/SU_RedSamp_k25mincov2.fasta.trans_filtered.txi";
#output_dir="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/filtered_transcriptome";
#perl /home/dwood/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --seqType fq --transcripts $transcripts --samples_file $samples_file --est_method kallisto --output_dir $output_dir --thread_count 20 --prep_reference;





