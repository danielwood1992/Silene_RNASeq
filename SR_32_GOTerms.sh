#KR3: stress related terms in resIresJ up
resIresJ="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/resIresJ_samedir_up_100622.csv.clean.GOseq.enriched";
echo "grep stress in resIresJ_upregulated";
grep -c "stress" $resIresJ;
echo "grep Stress in resIresJ_upregulated";
grep -c "Stress" $resIresJ;
cut -f1-9 $resIresJ > /scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/Table_SP1_110622.tsv;

#KR4: stress terms in resFresG
resFresG="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/resFresG_sig_overlap_samedir_100622.csv.clean.GOseq.enriched";
echo "grep stress in resFresG";
grep -c "stress" $resFresG;
echo "grep Stress in resFresG";
grep -c "Stress" $resFresG;
cut -f1-9 $resFresG > /scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/Table_SP2_110622.tsv;

#KR38: metal tolerance terms...
resBresC="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/resBresC_sig_overlap_samedir_100622.csv.clean.GOseq.enriched";
longesti="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SU_RedSamp_k25mincov2.fasta.trans_filtered.longesti";
echo "resBresC metal";
grep -c "metal" $resBresC;
echo "resBresC zinc";
grep -c "zinc" $resBresC;
grep "zinc ion transport" $resBresC | rev | cut -f1 | rev | sed 's/, /\n/g' | grep -f - $longesti -A 1 > $resBresC.zinciontransport.txt;
cut -f1-9 $resBresC > /scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/Table_SP4_040722.tsv;


#KR39 DPN metal tolerance terms
resKresLDPN="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/resKresL_DPN_sig_overlap_samedir_100622.csv.clean.GOseq.enriched";
echo "resKresLDPN metal";
grep -c "metal" $resKresLDPN;
grep "glutathione metabolic process" $resKresLDPN | rev | cut -f1 | rev | sed 's/, /\n/g' | grep -f - $longesti -A 1 > $resKresLDPN.glutathione.txt;
cut -f1-9 $resKresLDPN > /scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/Table_SP5_040722.tsv;

#KR? GO terms in plast vs. ancplast...

plast="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/resKresL_DPN_sig_overlap_samedir_ancplast_100622.csv.clean.GOseq.enriched";
noplast="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/resKresL_DPN_sig_overlap_samedir_noancplast_100622.csv.clean.GOseq.enriched";
dir="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO";
cut -f6 $plast | sort > $dir/plast.txt;
cut -f6 $noplast | sort > $dir/noplast.txt;
comm -23 $dir/plast.txt $dir/noplast.txt > $dir/plast_only.txt;
cp $dir/plast_only.txt $dir/plast_only.txt.tmp;
perl -pe 's/ /\n/g' $dir/plast_only.txt.tmp | sort | uniq -c | sort -nr > $dir/plast_only.txt.list;
comm -13 $dir/plast.txt $dir/noplast.txt > $dir/noplast_only.txt.tmp;
perl -pe 's/ /\n/g' $dir/noplast_only.txt.tmp | sort | uniq -c | sort -nr > $dir/noplast_only.txt.list;
cut -f1-9 $plast > /scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/Table_SP3_040722.tsv;
cut -f1-9 $noplast > /scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_31_GO/Table_SP9_040722.tsv;


perl -pe 's/ /\n/g' $dir/noplast_only.txt.tmp | sort | uniq -c | sort -nr > $dir/noplast_only.txt.list;
echo "stress ancplast";
grep -c "stress" $plast;
echo "stress noancplast";
grep -c "stress" $noplast;

echo "metal ancplast";
grep -c "metal" $plast;
echo "metal noancplast";
grep -c "metal" $noplast;

comm -12 $dir/plast.txt $dir/noplast.txt > $dir/plastandnoplast.txt.tmp;

