#!/usr/bin/perl
#Not actually used, but should probably work...

my $samples_file = "/home/b.bssc1d/scripts/formal_RNAseq_scripts/raw_reads_list.txt";
my $output_dir = "scratch/b.bssc1d/Silene_RNA_seq/1_CleanedReads/";
my $script = "/home/b.bssc1d/scripts/formal_RNAseq_scripts/sub_SR_1_CleanReads.sh";

my ($line, @temp, $read1, $read2, $base1, $base2);

open(IN, "<$samples_file");
open(OUT, ">SR_1_commands.txt");
while(!eof(IN)){
	$line = readline *IN;
	chomp $line;
	@temp = split/\t/, $line;
	$read1 = $temp[0];
	$read2 = $temp[1];
	$base1 = $read1;
	$base2 = $read2;
	$base1 =~ s/^.*\///g;
	$base2 =~ s/^.*\///g; 	

	print OUT "$script $read1 $read2 $output_dir woof $base1 $base2\n";	
}
#Uncomment this when you're ready to rumble. 
#`sbatch ~/scripts/24h_submission SR_1_commands.txt`;

