use POSIX;
use strict;

my ($list, $kall, $line, @temp, %hash);
$list = $ARGV[0];
$kall = $ARGV[1];

open(IN, "<$list");
while(!eof(IN)){$line = readline *IN; chomp $line; $hash{$line} = "";}

open(IN, "<$kall");
my $first = "T";
open(OUT, ">$kall.red");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($first eq "T"){
		$first = "F";
		print OUT $line."\n";
	}
	if (exists ($hash{(split/\t/, $line)[0]})){
		print OUT $line."\n";
	}
}










