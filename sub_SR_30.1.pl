use POSIX;
use strict;

#Rejigs the blat output so the start is always actually at the start...
my $file = $ARGV[0];
open(OUT, ">$file.re");
open(IN, "<$file");
my ($line, @temp, $start, $end, $new);
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	$start = $temp[15];
	$end = $temp[16];
	if ($start > $end){
		$temp[15] = $end;
		$temp[16] = $start;
	}
	$new = join "\t", @temp;
	print OUT $new."\n";
}
