#Estimating sequence legnths
script="/home/dwood/trinityrnaseq-v2.10.0/util/abundance_estimates_to_matrix.pl";
#This is also best to do manually.
#In practice, take the trinity file and do the following in vim:
#:v/>/d
#:%s/>//g
#%s/ .*//g
#%s/\(.*\)/\1\t\t/g
#:%s/_i.*\t/\t/g
#This will give you the file you need.
gene_trans_map="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/Trinity_trans_map";
#Quant file - probably best to generate this yourself. 
#In practice it'll be going to the relevant directory with each of your BD11Z/ etc. files in it
#Then do ls $PWD/*/*tsv > quant_file.txt
#(As of 03/07/20 this is ~/Hal_HDD1/dwood/Formal_RNA_seq/, doing ls */*abundance.tsv > quant_file. 
#Have this as the quant file for the time being, but it could change.
quant_file="/home/dwood/Hal_HDD1/dwood/Formal_RNASeq/quant_file.txt";


perl $script --est_method kallisto --gene_trans_map $gene_trans_map --quant_files $quant_file --name_sample_by_basedir




