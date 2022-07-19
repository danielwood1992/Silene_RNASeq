use POSIX;
use strict;

my $vcf = $ARGV[0];

my ($line, @temp, %hash);
my $i = "1";

#So now we need to open the vcf...
open(IN, "<$vcf");
open(OUT, ">$vcf.mod.vcf");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	if ($line =~ m/^##contig/){
		@temp = split/=/, $line;
		$line =~ s/$temp[2]/$i,length/g;
		$line =~ s/_i.*,/,/g;
		my $name = $temp[2];
 		$name =~ s/,.*//g;
		$name =~ s/_i.*//g;
		$hash{$name} = $i;
#		print $name."\n";
		$i++;
		print OUT $line."\n";
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
