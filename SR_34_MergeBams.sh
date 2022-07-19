#!/bin/sh
#SBATCH -o /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR_34.%A.out.txt
#SBATCH -e /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR_34.%A.err.txt
#SBATCH --ntasks=2
#SBATCH --time=1-12:00:00
#SBATCH --partition=compute
#$SBATCH --mem-per-cpu=9G


#sample_names="/home/b.bssc1d/scripts/formal_RNAseq_scripts/samples_file.txt2";
sample_names="/home/b.bssc1d/scripts/formal_RNAseq_scripts/samples_file.txt3";

dir="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree";
module load picard/2.20.2;
module load parallel;
module load samtools;

#The whole thing
#parallel -j 2 "java -Xmx8096m -jar $PICARD AddOrReplaceReadGroups INPUT=$dir/{1}C.hisat2.bam OUTPUT=$dir/{1}C.hisat2.bam.tomerge RGID={1} RGPL=\"illumina\" RGLB={1} RGSM={1} RGPU={1} && java -Xmx8096m -jar $PICARD AddOrReplaceReadGroups INPUT=$dir/{1}Z.hisat2.bam OUTPUT=$dir/{1}Z.hisat2.bam.tomerge RGID={1} RGPL=\"illumina\" RGLB={1} RGSM={1} RGPU={1} && samtools sort $dir/{1}C.hisat2.bam.tomerge && samtools sort $dir/{1}Z.hisat2.bam && samtools -f merge $dir/{1}.bam $dir/{1}C.hisat2.bam.tomerge $dir/{1}Z.hisat2.bam.tomerge && samtools sort $dir/{1}.bam -o $dir/{1}.sorted.bam && echo {1} Done >> $dir/SR_34_log.txt" :::: $sample_names;

#Second bit
parallel -j 2 "samtools merge -f $dir/{1}.bam $dir/{1}C.hisat2.bam.tomerge $dir/{1}Z.hisat2.bam.tomerge && echo {1} Done >> $dir/SR_34_log.txt" :::: $sample_names;
