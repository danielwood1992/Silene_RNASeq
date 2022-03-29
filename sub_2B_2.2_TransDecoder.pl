use POSIX;
use strict;

#This script will take out the pep files, or whatever, and use these to return the original transcriptome (so not just the cds, but everything...)
#Obviously this assumes a specific file format for both the peps and transcriptome...


my $transcriptome = $ARGV[0];
my $peps = $ARGV[1];

my $line;
my %hash;

open(IN, "<$peps");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($line =~ m/^>/){
		$line =~ s/\..*//g;
		$hash{$line} = "";
	}
}
close(IN);
#This will give you a hash with the transcript names for everything in longestorfs...
#Then need to extract the sequences from the original transcriptome...

my ($num, $n, $trans_filtered);

$num = keys %hash;
$n = 0;

$trans_filtered = $transcriptome.".trans_filtered";

my ($name, $seq, $namecheck);

open(IN, "<$transcriptome");
open(OUT, ">$trans_filtered");
while(!eof(IN)){
	$name = readline *IN;
	chomp $name;
	if ($name =~ m/^>/){
		$seq = readline *IN;	
		chomp $seq;
		$namecheck = $name;
		$namecheck =~ s/ .*//g;
		if (exists $hash{$namecheck}){
			$n++;
			print OUT "$name\n$seq\n";
		}
	}else{
		die "Something went horribly, horribly, horribly, horribly wrong.\n";
	}
}
print "$n / $num sequences retrieved succesfully, printed in $trans_filtered";	
