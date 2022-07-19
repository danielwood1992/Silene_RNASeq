use POSIX;
use strict;
#Based on SR_9. Use this to find what taxa the various contigs likely come from...

#Get the taxon IDs and the numbers...
my ($line, @temp, %hash, $hash2);
#my $blast_results_list = "/scratch/b.bssc1d/040322_HalBackups/blastdb_prot/blast_results/blastp.outfmt6.longesti";
my $blast_results_list = "/scratch/b.bssc1d/040322_HalBackups/blastdb_prot/blast_results/blastp.outfmt6.longesti.contam";

`cut -f2 $blast_results_list | cut -f2 -d'_' | sort | uniq -c > $blast_results_list.tax`;
`perl -p -i -e 's/ +/ /g' $blast_results_list.tax`;
`perl -p -i -e 's/^ //g' $blast_results_list.tax`;

#Get all the embryophte taxon IDs...

open(IN, "<$blast_results_list.tax");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/ /, $line;
	$hash{$temp[1]} = $temp[0];
}


my $swissprot_list = "/home/b.bssc1d/scripts/formal_RNAseq_scripts/speclist_uniprot_2020_02.txt_mod";

my (%hash2);
open(IN, "<$swissprot_list");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	if (exists $hash{$temp[0]}){
			print "bork\t";
			$hash2{$temp[1]} = $temp[0];
	}
}

my $namedump = "/home/b.bssc1d/scripts/formal_RNAseq_scripts/270422_taxdump/fullnamelineage.dmp";
open(OUT, ">$blast_results_list.tax.sum");
open(IN, "<$namedump");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	if (exists ($hash2{$temp[0]})){
		print OUT "$hash{$hash2{$temp[0]}}\t$line\n";
		
	}
}

