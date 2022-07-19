#Before running this script, you will need to load module R/4.2.0

library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(phytools)

col_coast = "#87CEFA"
col_mine = "orange"

col_PPBD = "#046C9A"
col_GRSA = "#F21A00"


snp_tree = read.tree(file="/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/merged_snps_filtered.vcf.snpphylo140722.bs.tree");
ggtree(snp_tree)+geom_tiplab()+geom_text(aes(label=node))+geom_treescale()
snp_tree$edge.length
#To find which branch length to halve...
plotTree(snp_tree)
edgelabels(snp_tree$eedge.length)
#snp_tree = reroot(snp_tree, 18, (0.03134)/2)
#snp_tree = reroot(snp_tree, 18, (0.0285)/2)
snp_tree = reroot(snp_tree, 18, 0.0285)

ggtree(snp_tree)+geom_tiplab()+geom_text(aes(label=node))
snp_tree2 = groupClade(snp_tree, .node=c(22, 20, 17, 15))

woof = ggtree(snp_tree2, aes(color=group), size = 1.5)+
#  scale_color_manual(values=c("black", col_mine, col_coast, col_coast, col_mine))+
  scale_color_manual(values=c("black", col_mine, col_coast, col_mine, col_coast))+

  theme(legend.position = "none")+
  geom_cladelabel(node=22, label = "T1", offset = 0.02, align = T, col = col_GRSA, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_cladelabel(node=20, label = "S1", offset = 0.02, align = T, col = col_GRSA, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_cladelabel(node=17, label = "T2", offset = 0.02, align = T, col = col_PPBD, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_cladelabel(node=15, label = "S2", offset = 0.02, align = T, col = col_PPBD, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_treescale()+
  
  xlim(0,0.5)
#  geom_text(aes(label=node))+
#  geom_tiplab()
woof
gridExtra::grid.arrange(flip(woof, 19, 14) %>% flip(22, 20), ncol=1)

pdf(file = "/scratch/b.bssc1d/Silene_RNA_seq/NEE_Revisions/SR_33_Tree/merged_snps_filtered.vcf.snpphylo140722.bs.tree.SR38.pdf")
woof
dev.off()


