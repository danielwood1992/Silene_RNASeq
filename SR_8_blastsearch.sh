#For my annotation
#Do this for each of the xaa split peptide files from the longestiorf files....
dir="/home/dwood/blastdb_prot/Ref_protein_sets/Transdecoder_dir/blast_results/"
blastp -query $1 -db /home/dwood/Trinotate-Trinotate-v3.2.1/admin/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > $1.blast


