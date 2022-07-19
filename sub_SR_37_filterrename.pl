use POSIX;
use strict;

my $transcriptome = $ARGV[0];
my $vcf = $ARGV[1];

my $first = "T";
my ($line, @temp, %hash);
my $i = "1";
open(IN, "<$transcriptome");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($first eq "T"){
		$first = "F";
	}else{
		@temp = split/"/, $line;
		#print $temp[3]."\n";
		$hash{$temp[3]} = "";
	}
}

#So now we need to open the vcf...
open(IN, "<$vcf");
open(OUT, ">$vcf.included.vcf");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($line =~ m/^##contig/){
		@temp = split/=/, $line;
		$temp[2] =~ s/_i.*//g;
#		print $temp[2]."\n";
		if (exists($hash{$temp[2]})){
#			print "woof\n";
			$line =~ s/$temp[2]/$i/g;
			$line =~ s/_i.*,/,/g;
			$hash{$temp[2]} = $i;
			$i++;
			print OUT $line."\n";
		}
	}elsif ($line =~ m/^TRINITY/){
		@temp = split/\t/, $line;
		$temp[0] =~ s/_i.*//g;
		if (exists($hash{$temp[0]})){
			$temp[0] = $hash{$temp[0]};
			$line = join("\t", @temp);
			print OUT $line."\n";
		}		
	}else{
		print OUT $line."\n";
	}	
}
print $i."\n";
