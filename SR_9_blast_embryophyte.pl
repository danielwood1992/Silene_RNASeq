use POSIX;
use strict;
#based on blast_embryophyte_swissprot.pl

#my $blast_results_list = $ARGV[0];

#Need to run in appropriate dir: 
#/home/dwood/blastdb_prot/Ref_protein_sets/Transdecoder_dir/
#Need to have xaa.blast, the blast results, from the previous script.

#To get this list...
#See Bioinformatics_formal_methods_v2 om google docs.

#my $namedump = "/home/dwood/Hal_HDD1/dwood/nr_database_0705020/nt_database/fullnamelineage.dmp";

#my @files = grep(/x.*pc/, readdir(DIR));
#18_05_20 for swissprot



#Get all the embryophte taxon IDs...
my ($line, @temp, %hash);
open(IN, "<$namedump");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	if ($line =~ m/Embryophyta/){
		$hash{$temp[0]}= 0;
	}
}
print $temp[0]."\n";

#So this is then a list of all the embryophyte taxon IDS...
#So then we just need to update this to take into account the swissprot name IDs...
#Note - to download this file, see Bioinfomratics_Formal_Methods_v2, Filtering section. 
#my $swissprot_list = "/home/dwood/scripts/RNA_seq_scripts/speclist_uniprot_2020_02.txt_mod";
my (%hash2);
open(IN, "<$swissprot_list");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	if (exists $hash{$temp[1]}){
			print "bork\t";
			$hash2{$temp[0]} = "";
	}
}

opendir(DIR, ".");
my @files = grep(/x.*blast$/, readdir(DIR));

#This is...the thing, then...
my ($item, $taxid);
foreach $item (@files){
	open(IN, "<$item");
	open(OUT, ">$item.emb");
	while(!eof(IN)){
		$line = readline *IN;
		chomp $line;
		$taxid = (split/_/, (split/\t/, $line)[1])[1];
		if (exists $hash2{$taxid}){
			print OUT "$line\temb_yes\n";
		}else{
			print OUT "$line\temb_no\n";
		}
	}
	close IN;
}



