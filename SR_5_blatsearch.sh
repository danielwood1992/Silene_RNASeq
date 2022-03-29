#NOTE - the way to use this script:
#1. If not already done, use fasta_oneline.pl on your query.
#2. Split the query using split -l 200 $query or whatever, to get N files of equal size equal to the number of thereads you want to use

#i.e. for file in x*; do sh $this_script; done

ref="/home/dwood/blat_Suniflora/TGS_GC_fmlrc.scaff_seqs.fa";
dir="/home/dwood/blat_Slatifolia/RNAX_DNAX/SUniflora";
#query="/home/dwood/blat_Suniflora/SU_RedSamp_k25mincov2.fasta.trans_filtered";
query=$1;
~/bin/x86_64/blat $dir/$ref $query -t=dnax -q=rnax $query.psl

