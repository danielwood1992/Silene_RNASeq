#!/bin/bash
#SBATCH -o /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR6.%A.out.txt
#SBATCH -e /scratch/b.bssc1d/Silene_RNA_seq/joblogs/SR_36.%A.err.txt
#SBATCH --ntasks=12
#SBATCH --time=1-12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8G

module load parallel;
module load bcftools;
genome="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta.trans_filtered.longesti";
list="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/bcf_list.txt";
metalist="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/bcf_metalist.txt";
metalist2="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/bcf_metalist2.txt";


#Throws some errors but seems to have worked somehow...
#parallel -j 12 "bgzip {1} && bcftools index {1}.gz" :::: $list;
dir=$(dirname $list);
#parallel -j 4 "bcftools merge --gvcf $genome --file-list {1} -o {1}.merged.bcf.gz && bcftools index {1}.merged.bcf.gz" :::: $metalist;
##so this should hopefully give you your series of merged files.
#Then you need to filter for missing data...
#parallel -j 4 "bcftools filter -e 'AN<4' {1}.merged.bcf.gz -o {1}.merged.atleast2.bcf.gz && bcftools index {1}.merged.atleast2.bcf.gz" :::: $metalist;

file_string=$(cat $metalist2 | awk '{print $0}' ORS=" ");
echo $file_string;
#bcftools isec -c all -n=4 -p shared_snps $file_string;

#for file in $dir/shared_snps/*vcf; do bgzip $file; done;
#for file in $dir/shared_snps/*vcf.gz; do bcftools index $file; done;

#ls $dir/shared_snps/*vcf.gz > $dir/another_list.txt;

#bcftools merge --gvcf $genome --file-list $dir/another_list.txt -o -Obz $dir/merged_snps.bcf.gz;
#bcftools index $dir/merged_snps.bcf.gz;
bcftools view -m2 -M2 -v snps $dir/merged_snps.vcf -Obz -o  $dir/merged_snps_tmp.bcf.gz;
bcftools filter -i 'QUAL>20 & AC>1' $dir/merged_snps_tmp.bcf.gz -Ov -o $dir/merged_snps_filtered.vcf; 
