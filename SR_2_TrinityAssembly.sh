#Note - this s based off /home/dwood/scripts/RNA_seq_scripts/2B_2.1_submitTrinity_finalvused.sh
#run this with nohup &

file="/home/dwood/scripts/formal_RNAseq_scripts/samples_file.txt_reducedsamples";

~/trinityrnaseq-v2.10.0/Trinity --seqType fq --CPU 35 --max_memory 110G --samples_file $file --verbose --min_kmer_cov 2  --output ~/Hal_HDD1/dwood/trinity_15040_combined4samples 
