my $blastp = "/home/dwood/blastdb_prot/Ref_protein_sets/Transdecoder_dir/blast_results/blastp.outfmt6.emb";
my $transcriptome = "/home/dwood/blastdb_prot/Ref_protein_sets/Transdecoder_dir/blast_results/SU_RedSamp_k25mincov2.fasta.trans_filtered";


#my $transcriptome = $ARGV[0];
#my $blastp = $ARGV[1];
my ($line, $seq, $name, $clust, %length, %iname);
open(IN, "<$transcriptome");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($line =~ m/^>/){
		$seq = readline *IN;
		chomp $seq;
		$name = $line;
		$name =~ s/^>//g;
		$name =~ s/ .*//g;
		$clust = $name;
		$clust =~ s/_i.*//g;
		if ($length{$clust} < length($seq)){
			$length{$clust} = length($seq);
			$iname{$clust} = $name;
		}
	}else{
		die "something went horrily wrong\n";
	}
}
my (%hash, $item);
foreach $item (keys %iname){
	print $iname{$item}."\t";
	$hash{$iname{$item}} = "";
}

#So then filter the blast database for these guys (if they've got more than one protein...well,go for the one with the top bitscore...
my (%toprint);
open(IN, "<$blastp");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	$blastname = $temp[0];
	$blastname =~ s/\..*//g;
	if (exists ($hash{$blastname})){
		if (exists $toprint{$blastname}){
			if ($temp[11] > (split/\t/, $toprint{$blastname})[11]){
				$toprint{$blastname} = $line;
			}
		}
		else{
			$toprint{$blastname} = $line;
		}
	}
}
open(OUT, ">$blastp.longesti");
foreach $item (keys %toprint){
	print OUT $toprint{$item}."\n";
}

