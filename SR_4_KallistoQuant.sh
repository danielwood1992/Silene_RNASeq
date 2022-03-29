#Based on: run_kallisto_trinityutils.sh
#FOR ANALYSIS - DEFAULT 22/06/20
samples_file="/home/dwood/scripts/formal_RNAseq_scripts/samples_file.txt";
transcripts="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/SU_RedSamp_k25mincov2.fasta.trans_filtered";
output_dir="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq";
perl /home/dwood/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --seqType fq --transcripts $transcripts --samples_file $samples_file --est_method kallisto --output_dir $output_dir --thread_count 10 --prep_reference;
