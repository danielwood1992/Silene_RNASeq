#Gets an N50 from lengths.
my $fasta = $ARGV[0];
my ($name, $seq, @array, $length);
open(IN, "<$fasta");
while(!eof(IN)){
	$name = readline *IN;
	$seq = readline *IN;
	$length = length($seq);
	push @array, $length;
	$tot = $tot+$length;
	if ($name =~ m/^/){
	}else{
		die "wrong format";
	}
}
@array = sort { $a <=> $b }  @array;
#Half the bases:
my $half = $tot/2;
my $first = "T";
foreach $item (@array){
	$half_tot = $half_tot + $item;
	if ($first eq "T" && $half_tot > $half){
		print $item."\n";
		$first = "F";
	}
}
