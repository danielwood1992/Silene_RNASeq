use POSIX;
use strict;

#Rejigs the blat output so the start is always actually at the start...
my $file = $ARGV[0];
open(OUT, ">$file.filt");
open(IN, "<$file");
my ($line, @temp, $length);
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	$length = $temp[0]+$temp[1];
	if ($length > 500){
		print OUT $line."\n";
	}
}
