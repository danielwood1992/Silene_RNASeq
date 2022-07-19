use POSIX;
use strict;

#Ok so I guess things to do...
#a) See what the percentage match length decrease on the next longest scaffold is...(unless this is what you'd expect where genes are fused together? But I think probably enough to say that it's unlikely to be two recently duplicated genes, no?
#
#my $psl_file = "/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_29/sorted_merged_outfile_best.psl";
#So I guess the problems are... 1) if it's 
my $psl_file = $ARGV[0];
open(OUT, ">$psl_file.SR30");
my $to_skip = 5;
my $i = 0;
my ($line, @temp, $length, $chrom, $query, %hash, $pc);
my $length_lim = "200";
my $pc_lim = "0.9";
open(IN, "<$psl_file");

while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($i < $to_skip){
		$i++;
	}else{
		@temp = split/\t/, $line;
		$length = $temp[0]+$temp[1];
		$chrom = $temp[13];
		$query = $temp[9];
		$pc = $temp[0]/($temp[0]+$temp[1]);
		if ($length > $length_lim){
			if (exists ($hash{$query}{$chrom})){
				if ($length > $hash{$query}{$chrom}[0]){

				}
			}else{
				$hash{$query}{$chrom}[0] = $length;
				$hash{$query}{$chrom}[1] = $pc;
			}
		}
	}

		
}
my ($item, $item2, @array, $ratio, $tot, $NAc); 	
foreach $item (keys %hash){
	$tot++;
#	print $item."\n";
	#So for each chromosome thing...
	if (scalar(keys $hash{$item}) == 1){
		print OUT "$item\tNA\tNA\tNA\n";
#		print "NA\n";

		$NAc++;	
	}else{
		foreach $item2 (keys $hash{$item}){
		#	print "$hash{$item}{$item2}[1] ";
		}
		@array = values $hash{$item};
##		print "$array[0][1] ";
		@array = sort { $a->[0] <=> $b->[0] } @array;
		$ratio = $array[-2][0]/$array[-1][0];
		print OUT "$item\t$ratio\t$array[-2][0]\t$array[-2][1]\t$array[-1][0]\t$array[-1][1]\n";	
	}
}
print "Of $tot contigs, $NAc don't have a hit on a different scaffold\n";
