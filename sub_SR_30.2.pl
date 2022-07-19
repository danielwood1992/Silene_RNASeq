use POSIX;
use strict;

#So this script needs to...
#From the whatever file...find a way of finding those that don't overlap...

my $file = $ARGV[0];
open(OUT, ">$file.not");
open(IN, "<$file");
my $limit = 300;
my ($match, $start, $end, $scaff, $prev_scaff, $first1, $first2, $int_end, $line, @temp, $N, $prev);
$first1 = "T";
$first2 = "T";
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	$match = $temp[0]+$temp[1];
	$start = $temp[15];
	$end = $temp[16];
	$scaff = $temp[13];

	#So if this is the first ever one...
	if ($scaff ne $prev_scaff){
		#If this is the first ever scaffold, will set things anew.
		#If it isn't, but the scaffold has changed, then sort things out...
		#
		#So if the last hit doesn't have any overlaps, print it out...
		if ($N == 1){
			print OUT "$prev.\n";
		}
		$first2 = "T";
		$int_end = $end;
		$prev_scaff = $scaff;
		$scaff = $prev_scaff;
		$N = 0; 

	}
	#If it's still the same scaffold, do the thing...
	if ($match > $limit){
		print "b";
		#So how do we do this.
		#So if with the scaffold it's the first, you need to set up the end and the starting point I guess?
		if ($first2 eq "T"){
			$first2 = "F";
			$int_end = $end;
			$N = 1;
			$prev = $line;
		}else{
			#So if the start is < end of the other one, they overlap.
			if ($start < $int_end){
				print "c";
				$N++;
				#Need to check if the new end is bigger than the old end. If it is, extend the end.
				if ($end > $int_end){
						$int_end = $end;
						
				}
			#So if the start is greater than the end, it doesn't overlap wth the one in front.	
			}else{
				print "d";
				#i) If the previous one is the only thing in that cluster, nothing overlaps with it - print it out.	
				if ($N == 1){
					print OUT $prev."\n"; #So we need to get the previous line somehow...
				}			
				#ii) Start a new cluster...
				$N = 1;
				$prev = $line;
				$int_end = $end;

			}	
					
		}
	}else{
	}
}
if ($N == 1){
	print OUT "$prev.\n";
}
print "\n";
