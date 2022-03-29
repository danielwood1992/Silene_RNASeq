#!/bin/bash --login
#SBATCH --partition=compute
#SBATCH --ntasks=3
#SBATCH --time=0-12:00
#SBATCH -o %A.out.txt
#SBATCH -e %A.err.txt

#where_we_are	both:_filtered.gz.trimmo.paired.gz

export PERL5LIB=~/perl5/lib/perl5/
module add FastQC/0.11.8
module add trimmomatic/0.39

read1=$1
read2=$2 
outputdir=$3
bgi_seqs="/home/b.bssc1d/F20FTSEUHT0112_SILghtE_primer_seqs.txt";
read1_base=$5
read2_base=$6

#NOTE - THIS SCRIPT DOESN'T WORK AT THE MINUTE, NEEDS TO GET THE TRIMMOIMATIC INPUT FILES AS THE OUTPUT FILES OF ILLUQC

fastqc $read1 -outdir=$outputdir
fastqc $read2 -outdir=$outputdir
#perl ~/bin/NGSQCToolkit_v2.3.3/QC/IlluQC.pl -pe $read1 $read2 $bgi_seqs A -l 70 -s 20 -t 2 -z g -o $outputdir
perl ~/bin/NGSQC_Toolkit_v2.3.3./QC/IlluQC.pl -pe $read1 $read2 1 2 $bgi_seqs -l 70 -s 20 -t 2 -z g -o $outputdir

java -jar $TRIMMOMATIC PE $outputdir/${read1_base}_filtered.gz $outputdir/${read2_base}_filtered.gz $outputdir/${read1_base}_filtered.gz.trimmo.paired.gz $outputdir/${read1_base}_filtered.gz.trimmo.unpaired.gz $outputdir/${read2_base}_filtered.gz.trimmo.paired.gz $outputdir/${read2_base}_filtered.gz.trimmo.unpaired.gz LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:70 -threads 3 
fastqc $outputdir/${read1_base}_filtered.gz.trimmo.paired.gz -outdir=$outputdir
fastqc $outputdir/${read2_base}_filtered.gz.trimmo.paired.gz -outdir=$outputdir
