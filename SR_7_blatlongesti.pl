#my $transcriptome = $ARGV[0];
#my $blat = $ARGV[1];
#
#}

my $dir = "/home/dwood/blat_Slatifolia/RNAX_DNAX/SUniflora";
my $transcriptome = "SU_RedSamp_k25mincov2.fasta.trans_filtered.cds";
my $blat = "merged_outfiles.psl";
;

my ($line, $seq, $name, $clust, %length, %iname);
open(IN, "<$dir/$transcriptome");
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
foreach $item (keys %iname){	print $iname{$item}."\t";
	$hash{$iname{$item}} = "";
}

#So then filter the blast database for these guys (if they've got more than one protein...well,go for the one with the top bitscore...
my (%toprint);
open(IN, "<$dir/$blat");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	$blastname = $temp[9];
	$blastname =~ s/\..*//g;
	if (exists ($hash{$blastname})){
		print "woof\t";
		if (exists $toprint{$blastname}){
			if ($temp[0] > (split/\t/, $toprint{$blastname})[0]){
				$toprint{$blastname} = $line;
			}
		}
		else{
			$toprint{$blastname} = $line;
		}
	}
}
open(OUT, ">$dir/$blat.longesti");
foreach $item (keys %toprint){
	print OUT $toprint{$item}."\n";

