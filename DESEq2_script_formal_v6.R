#0 LOADING IN THE DATA, FILTERING
remove(list=ls())
`%notin%` = Negate(`%in%`)
setwd(dir = "/Users/Danie/Documents/RNA_Seq")

###########################################
#1 FUNCTIONS
###########################################-----

#Fisher test for two overlapping sets (of genes)
my_fisher_test = function(overlap, set1, set2, background) {
  test_m = matrix(c(overlap, set1-overlap, set2-overlap, background-(set1+set2-overlap)), nrow = 2)
  x = fisher.test(test_m, alternative = "greater")
  return(x)
}

#Plot first 2 PCs (no centroids)
plotPCA12 <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  if (!all(intgroup1 %in% names(colData(object)))) {
    stop("the argument 'intgroup1' should specify columns of colData(dds)")
  }
  
  intgroup1.df <- as.data.frame(colData(object)[, intgroup1, 
                                                drop = FALSE])
  group1 <- if (length(intgroup1) > 1) {
    factor(apply(intgroup1.df, 1, paste, collapse = " : "))
  }else {
    colData(object)[[intgroup1]]
  }
  
  
  if (!all(intgroup2 %in% names(colData(object)))) {
    stop("the argument 'intgroup2' should specify columns of colData(dds)")
  }
  intgroup2.df <- as.data.frame(colData(object)[, intgroup2, 
                                                drop = FALSE])
  group2 <- if (length(intgroup2) > 1) {
    factor(apply(intgroup2.df, 1, paste, collapse = " : "))
  }else {
    colData(object)[[intgroup2]]
  }
  
  
  if (!all(intgroup3 %in% names(colData(object)))) {
    stop("the argument 'intgroup3' should specify columns of colData(dds)")
  }
  intgroup3.df <- as.data.frame(colData(object)[, intgroup3, 
                                                drop = FALSE])
  group3 <- if (length(intgroup3) > 1) {
    factor(apply(intgroup3.df, 1, paste, collapse = " : "))
  }else {
    colData(object)[[intgroup3]]
  }
  
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], 
                  group1 = group1, group2 = group2, group3 = group3, name = colnames(object))
  intgrop3d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group1 = group1, 
                          intgroup1.df, group2 = group2, intgroup2.df, name = colnames(object))
  
  
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1,2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 8, stroke=2, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                                                                                                                            100), "% variance")) + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

#Plot first 2 PCs (centroids)

plotPCA12_centroids <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  pca_thing = pca
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  if (!all(intgroup1 %in% names(colData(object)))) {
    stop("the argument 'intgroup1' should specify columns of colData(dds)")
  }
  
  intgroup1.df <- as.data.frame(colData(object)[, intgroup1, 
                                                drop = FALSE])
  group1 <- if (length(intgroup1) > 1) {
    factor(apply(intgroup1.df, 1, paste, collapse = " : "))
  }else {
    colData(object)[[intgroup1]]
  }
  
  
  if (!all(intgroup2 %in% names(colData(object)))) {
    stop("the argument 'intgroup2' should specify columns of colData(dds)")
  }
  intgroup2.df <- as.data.frame(colData(object)[, intgroup2, 
                                                drop = FALSE])
  group2 <- if (length(intgroup2) > 1) {
    factor(apply(intgroup2.df, 1, paste, collapse = " : "))
  }else {
    colData(object)[[intgroup2]]
  }
  
  
  if (!all(intgroup3 %in% names(colData(object)))) {
    stop("the argument 'intgroup3' should specify columns of colData(dds)")
  }
  intgroup3.df <- as.data.frame(colData(object)[, intgroup3, 
                                                drop = FALSE])
  group3 <- if (length(intgroup3) > 1) {
    factor(apply(intgroup3.df, 1, paste, collapse = " : "))
  }else {
    colData(object)[[intgroup3]]
  }
  
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], 
                  group1 = group1, group2 = group2, group3 = group3, name = colnames(object))
  intgrop3d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group1 = group1, 
                          intgroup1.df, group2 = group2, intgroup2.df, name = colnames(object))
  
  #Need to add the centroids somehow.
  
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[3,4]
    return(d)
  }
  return(pca_thing)
  ggplot(data = d, aes_string(x = "PC1", y = "PC2")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #  ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 10, stroke=4, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                                                                                                                             100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) + 
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

###################################
#2 Fisher Table setup
###################################
Fisher_Table = data.frame("Comparison", "All_pvalue (0 = <2.2e-16)", "All_oddsratio", "Samedir_pvalue", "Samedir_oddsratio", stringsAsFactors = FALSE)

##################################
#3 Colours
##################################


ctrl_green = "#00CC33"
tmp = col2rgb(ctrl_green)
ctrl_green = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.3)
col_coast = "#87CEFA"
col_mine = "orange"
col_PPBD = "#046C9A"
col_GRSA = "#F21A00"
col_control = "antiquewhite"
tmp = col2rgb("lightgreen")
col_zinc = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.3)
col_zinc2 = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha =1)
zinc_col = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255)
ctrl_col = col_control
ctrl_col_temp = col2rgb(ctrl_col)
ctrl_col_temp2 = rgb(ctrl_col_temp[1]/255, ctrl_col_temp[2]/255, ctrl_col_temp[3]/255)

tmp = col2rgb("yellow1")
gold_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.3)
gold_plot_light = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.1)
tmp = col2rgb("red")
red_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.3)
purple_col = col2rgb("plum")
purple_plot = rgb(purple_col[1]/255, purple_col[2]/255, purple_col[3]/255, alpha = 0.3)
purple_col2 = rgb(purple_col[1]/255, purple_col[2]/255, purple_col[3]/255)
another_plot_purple = purple_col2

ftmp = col2rgb(col_mine)
mine_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.2)
tmp = col2rgb(col_coast)
coast_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.2)
tmp = col2rgb(col_GRSA)
GRSA_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.6)
tmp = col2rgb(col_PPBD)
PPBD_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.6)

###########################################
#4 Load relevant libraries
###########################################

library(DESeq2) #1.26.0
library(tximport) #1.14.2
#library(readr)
library(ggplot2)
#library(gridExtra)
#library(UpSetR)
#library(ComplexHeatmap)
#library(SuperExactTest)
library(reshape2)
#library(ashr)
library(VennDiagram)
#library(tidyverse)
library(gghalves)
library(ggnewscale)
library(rcompanion)

##########################################
#5 Setup and filtering
##########################################

samples = read.table("/Users/Danie/Documents/RNA_Seq/samples_file.txt_forslueth.tsv", sep = "\t", header = TRUE)                     
files = file.path("/Users/Danie/Documents/RNA_Seq/Formal_RNASeq/", samples$Sample, "abundance.h5")
samples

#Adjust samples table to add "Pop_Cond" variable
samples$Pop_Cond = paste(samples$Pop, samples$Cond, sep="")
samples$Ecotype = c("C", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M","M", "C", "C", "C", "C", "C", "C")
samples$Geog = c("E", "E", "E", "E", "E", "E", "W", "W", "W", "W", "W", "W", "E", "E", "E", "E", "E", "E", "W", "W", "W", "W", "W", "W")
samples$Ecotype = as.factor(samples$Ecotype)
samples$Geog = as.factor(samples$Geog)
samples


#Read in Trinity gene-transcript map 
tx2gene = read.csv("/Users/Danie/Documents/RNA_Seq/Trinity-clusters.txt", sep = "\t", header = FALSE)
tx2gene = tx2gene[c("V2", "V1")]
head(tx2gene)

#Read in counts - summarise to gene level
tx.kallisto = tximport(files, type = "kallisto", txOut = FALSE, tx2gene=tx2gene, ignoreTxVersion = TRUE)
head(tx.kallisto)

#Convert tx object to DESeq object (design matrix can change later)
ddTxi = DESeqDataSetFromTximport(tx.kallisto, colData = samples, design = ~ Cond + Pop)
remove(tx.kallisto)
dim(ddTxi)
ddTxi_backup = ddTxi

#Filtering 1 - filter contigs by sequence criteria 
#Either a >50bp match to latifolia genome, and/or a >50bp match to an Embryophyte

#Load in data
blast_info = read.csv("/Users/Danie/Documents/RNA_Seq/blastp.outfmt6.emb.longesti", header = FALSE, sep = "\t")
blat_info = read.csv("/Users/Danie/Documents/RNA_Seq/merged_outfiles.psl_forlongesti.longesti", sep = "\t", header = FALSE)
GC_content = read.csv("/Users/Danie/Documents/RNA_Seq/SU_RedSamp_k25mincov2.fasta.trans_filtered.GClongesti", header=FALSE, sep="\t")
emb_yes = blast_info[blast_info$V13 == "emb_yes",]
emb_no = blast_info[blast_info$V13 == "emb_no",]

#Filter by blat match length
blat_info$match_length = blat_info$V1+blat_info$V2
blat_match_50plus = blat_info[blat_info$match_length > 50,]
blat_50plus_list = blat_match_50plus$V10

#Filter by embryophyte match length
emb_50plus = emb_yes[emb_yes$V4 > 50,]
dim(emb_50plus)
length(blat_50plus_list)

#GC used to get full list of names (and check GC content before/after)
GC_content_Filter2 = GC_content[GC_content$V1 %in% blat_50plus_list | GC_content$V1 %in% emb_50plus$V1,]
dim(GC_content_Filter2)

#Filtering
ddTxi = subset(ddTxi, rownames(ddTxi) %in% GC_content_Filter2$V1)

#Check GC contnent
par(mfrow=c(2,1))
hist(GC_content$V2, breaks = 100)
hist(GC_content_Filter2$V2, breaks = 100)

#Filtering 2 - remove genes with a total of <10 counts total
keep = rowSums(counts(ddTxi)) >= 10
ddTxi2 = ddTxi[keep,]
dim(ddTxi2)

rownames(ddTxi2)
write.csv(rownames(ddTxi2), "background_ddTxi2")

#Adding some more factors

ddTxi2$Cond = factor(ddTxi$Cond, levels = c("C", "Z"))
ddTxi2$Pop_Cond = factor(ddTxi$Pop_Cond, levels = c("BDZ", "BDC", "GRZ", "GRC", "PPZ", "PPC", "SAZ", "SAC"))

#Background fit - just use Pop_Cond
#design(ddTxi2) <- ~Cond+Ecotype+Geog
design(ddTxi2) <- ~Pop_Cond

#Fit the DESeq2 model
dds = DESeq(ddTxi2)

#Examine the overall patterns of variation
par(mfrow  = c(1,2))
boxplot(log2(counts(dds, normalized = FALSE)))
boxplot(log2(counts(dds, normalized = TRUE)))

####################################################
#6 Counts, PCAs in control/zinc
####################################################

vst_all = vst(dds)

object = vst_all
intgroup1 = "Cond"
intgroup2 = "Geog"
intgroup3 = "Ecotype"
ntop = 1000000

#a) PCAs

dds_C = dds[,dds$Cond == "C"]
vst_C = vst(dds_C)
dds_Z = dds[,dds$Cond == "Z"]
vst_Z = vst(dds_Z)

pdf()

#PCA all individuals, all conditions (with centoirds and arrows)
#Figure 2C (20/07/21)

#https://stackoverflow.com/questions/36208909/how-to-calculate-centroids-in-pca - based on this
PCA_all = plotPCA12_centroids(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-300, 200), my_ylims = c(-150, 150), my_background_col = "lightgray")
PCA_all

blork.x = as.data.frame(PCA_all$x)
blork.x$groups = samples$Pop_Cond
pca.centroids = aggregate(blork.x[,1:2], list(Type = blork.x$groups), mean)

PCA_all_plot = PCA_all+
  geom_segment(aes(x = pca.centroids[7,2], y = pca.centroids[7,3], xend = pca.centroids[8,2], yend = pca.centroids[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[7,2], y = pca.centroids[7,3], xend = pca.centroids[8,2], yend = pca.centroids[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[1,2], y = pca.centroids[1,3], xend = pca.centroids[2,2], yend = pca.centroids[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5,col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[1,2], y = pca.centroids[1,3], xend = pca.centroids[2,2], yend = pca.centroids[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot,  lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[5,2], y = pca.centroids[5,3], xend = pca.centroids[6,2], yend = pca.centroids[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[5,2], y = pca.centroids[5,3], xend = pca.centroids[6,2], yend = pca.centroids[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[3,2], y = pca.centroids[3,3], xend = pca.centroids[4,2], yend = pca.centroids[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[3,2], y = pca.centroids[3,3], xend = pca.centroids[4,2], yend = pca.centroids[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_point(data = pca_centroids_c, size = 4, shape = 21, fill = c(PPBD_plot, GRSA_plot, PPBD_plot, GRSA_plot), alpha = 1)

PCA_all_plot

pdf(file = "Fig2A_200721_PCA_all_centroids.pdf")
grid.draw(PCA_all_plot)
dev.off()

#PCA all individuals, control conditions only
#Figure 2C (20/07/21)

Fig2A = plotPCA12(vst_C, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_xlims = c(-75,90), my_ylims = c(-75,50), my_background_col = ctrl_green)
Fig2A
pdf(file = "Fig2A_200721_ControlPCA.pdf")
grid.draw(Fig2A)
dev.off()

#################################################
#7 Within condition comparisons
#################################################

ddTxi3 = ddTxi2
design(ddTxi3) <- ~Pop_Cond
dds2 = DESeq(ddTxi3)
resultsNames(dds2)
colnames(dds2) = paste(dds2$IND, dds2$Cond, sep="")
blork = counts(dds2, normalized=T)
head(blork)
write.csv(blork, "~/Normalised_Counts,tsv", sep = "\t")

resA = results(dds2, contrast=c("Pop_Cond", "SAC", "BDC"))
resA$Names = rownames(resA)
resB = results(dds2, contrast=c("Pop_Cond", "SAC", "GRC"))
resB$Names = rownames(resB)
resC = results(dds2, contrast=c("Pop_Cond", "BDC", "PPC"))
resC$Names = rownames(resC)
resD = results(dds2, contrast=c("Pop_Cond", "GRC", "PPC"))
resD$Names = rownames(resD)
resE = results(dds2, contrast=c("Pop_Cond", "SAZ", "BDZ"))
resE$Names = rownames(resE)
resF = results(dds2, contrast=c("Pop_Cond", "SAZ", "GRZ"))
resF$Names = rownames(resF)
resG = results(dds2, contrast=c("Pop_Cond", "BDZ", "PPZ"))
resG$Names = rownames(resG)
resH = results(dds2, contrast=c("Pop_Cond", "GRZ", "PPZ"))
resH$Names = rownames(resH)

#Remove NAs: removes p-values for genes where counts are too low to get reliable estimate, as specified in user guide.
resA_sig = resA[is.na(resA$padj) == FALSE,]
resB_sig = resB[is.na(resB$padj) == FALSE,]
resC_sig = resC[is.na(resC$padj) == FALSE,]
resD_sig = resD[is.na(resD$padj) == FALSE,]
resE_sig = resE[is.na(resE$padj) == FALSE,]
resF_sig = resF[is.na(resF$padj) == FALSE,]
resG_sig = resG[is.na(resG$padj) == FALSE,]
resH_sig = resH[is.na(resH$padj) == FALSE,]

#Significantly different genes: Adjusted p-value < 0.05. 
resA_sig = resA_sig[resA_sig$padj < 0.05,]
resB_sig = resB_sig[resB_sig$padj < 0.05,]
resC_sig = resC_sig[resC_sig$padj < 0.05,]
resD_sig = resD_sig[resD_sig$padj < 0.05,]
resE_sig = resE_sig[resE_sig$padj < 0.05,]
resF_sig = resF_sig[resF_sig$padj < 0.05,]
resG_sig = resG_sig[resG_sig$padj < 0.05,]
resH_sig = resH_sig[resH_sig$padj < 0.05,]

resA_sig$Names = rownames(resA_sig)
resB_sig$Names = rownames(resB_sig)
resC_sig$Names = rownames(resC_sig)
resD_sig$Names = rownames(resD_sig)
resE_sig$Names = rownames(resE_sig)
resF_sig$Names = rownames(resF_sig)
resG_sig$Names = rownames(resG_sig)
resH_sig$Names = rownames(resH_sig)

write.csv(resA_sig$Names, "resA_sig_210721")
write.csv(resB_sig$Names, "resB_sig_210721")
write.csv(resC_sig$Names, "resC_sig_210721")
write.csv(resD_sig$Names, "resD_sig_210721")
write.csv(resE_sig$Names, "resE_sig_210721")
write.csv(resF_sig$Names, "resF_sig_210721")
write.csv(resG_sig$Names, "resG_sig_210721")
write.csv(resH_sig$Names, "resH_sig_210721")

dim(resA_sig)
dim(resB_sig)
dim(resC_sig)
dim(resD_sig)
dim(resE_sig)
dim(resF_sig)
dim(resG_sig)
dim(resH_sig)



###############################################
#####8 between condition comparisons
###############################################

#I/J/K/L
#Includes individual as factor (woo clones)

ddTxi_ind = ddTxi2
ddTxi_ind$Ind_plant = as.factor(c(1, 2, 2, 1, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3))
design(ddTxi_ind) <- ~Pop + Pop:Ind_plant + Pop:Cond
dds3 = DESeq(ddTxi_ind)

resI = results(dds3, contrast = list("PopSA.CondZ"))
resI$Names = rownames(resI)
resJ = results(dds3, contrast = list("PopBD.CondZ"))
resJ$Names = rownames(resJ)
resK = results(dds3, contrast = list("PopGR.CondZ"))
resK$Names = rownames(resK)
resL = results(dds3, contrast = list("PopPP.CondZ"))
resL$Names = rownames(resL)

resI_sig = resJ[is.na(resI$padj) == FALSE,]
resJ_sig = resJ[is.na(resJ$padj) == FALSE,]
resK_sig = resK[is.na(resK$padj) == FALSE,]
resL_sig = resL[is.na(resL$padj) == FALSE,]

resI_sig = resI_sig[resI_sig$padj < 0.05,]
resJ_sig = resJ_sig[resJ_sig$padj < 0.05,]
resK_sig = resK_sig[resK_sig$padj < 0.05,]
resL_sig = resL_sig[resL_sig$padj < 0.05,]

resI_sig$Names = rownames(resI_sig)
resJ_sig$Names = rownames(resJ_sig)
resK_sig$Names = rownames(resK_sig)
resL_sig$Names = rownames(resL_sig)

dim(resI_sig)
dim(resJ_sig)
dim(resK_sig)
dim(resL_sig)

write.csv(resI_sig$Names, "resI_sig_210721.csv")
write.csv(resJ_sig$Names, "resJ_sig_210721.csv")
write.csv(resK_sig$Names, "resK_sig_210721.csv")
write.csv(resL_sig$Names, "resL_sig_210721.csv")

############################################
#9 Overlapping sets (includes venn diagrams)
############################################

#Overlapping gene sets of interest

#Bit of a mess in here at the

resBresC_sig_overlap = merge(as.data.frame(resB_sig), as.data.frame(resC_sig), by = "Names")
resFresG_sig_overlap = merge(as.data.frame(resF_sig), as.data.frame(resG_sig), by = "Names")
resAresE_sig_overlap = merge(as.data.frame(resA_sig), as.data.frame(resE_sig), by = "Names")
resDresH_sig_overlap = merge(as.data.frame(resD_sig), as.data.frame(resH_sig), by = "Names")
resAresD_sig_overlap = merge(as.data.frame(resA_sig), as.data.frame(resD_sig), by = "Names")

dim(resBresC_sig_overlap)
dim(resFresG_sig_overlap)
dim(resAresE_sig_overlap)
dim(resDresH_sig_overlap)


#resBresC overlaps
resBresC_sig_overlap$sign = resBresC_sig_overlap$log2FoldChange.x*resBresC_sig_overlap$log2FoldChange.y
resBresC_sig_overlap_samedir = resBresC_sig_overlap[resBresC_sig_overlap$sign > 0,]
dim(resBresC_sig_overlap_samedir)

write.csv(resBresC_sig_overlap_samedir$Names, "resBresC_sig_overlap_samedir_210721.csv")
resBresC_sig_overlap_samedir_UP = resBresC_sig_overlap_samedir[resBresC_sig_overlap_samedir$log2FoldChange.x > 0,]
resBresC_sig_overlap_samedir_DOWN = resBresC_sig_overlap_samedir[resBresC_sig_overlap_samedir$log2FoldChange.x < 0,]
dim(resBresC_sig_overlap_samedir_UP)
dim(resBresC_sig_overlap_samedir_DOWN)
write.csv(resBresC_sig_overlap_samedir_UP$Names, "resBresC_sig_overlap_samedir_up_210721.csv")
write.csv(resBresC_sig_overlap_samedir_DOWN$Names, "resBresC_sig_overlap_samedir_down_210721.csv")

resIresJ_sig_overlap = merge(as.data.frame(resI_sig), as.data.frame(resJ_sig), by = "Names")
dim(resIresJ_sig_overlap)[1]/dim(resI)[1] #% of transcriptome DE in sensitive populations between treatments

resKresL_sig_overlap = merge(as.data.frame(resK_sig), as.data.frame(resL_sig), by = "Names")
resKresL_notInotJ = resKresL_sig_overlap[resKresL_sig_overlap$Names %notin% resIresJ_sig_overlap$Names,]
dim(resKresL_notInotJ)

#Multiply FCs; those with positive value are in same direction. 
resIresJ_sig_overlap$Sign = resIresJ_sig_overlap$log2FoldChange.x*resIresJ_sig_overlap$log2FoldChange.y
resIresJ_sig_overlap_samedir = resIresJ_sig_overlap[resIresJ_sig_overlap$Sign > 0,]
resKresL_sig_overlap$Sign = resKresL_sig_overlap$log2FoldChange.x*resKresL_sig_overlap$log2FoldChange.y
resKresL_sig_overlap_samedir = resKresL_sig_overlap[resKresL_sig_overlap$Sign > 0,]

resIresJ_samedir_up = resIresJ_sig_overlap_samedir[resIresJ_sig_overlap_samedir$log2FoldChange.x > 0,]
resIresJ_samedir_down = resIresJ_sig_overlap_samedir[resIresJ_sig_overlap_samedir$log2FoldChange.x < 0,]
resKresL_samedir_up = resKresL_sig_overlap_samedir[resKresL_sig_overlap_samedir$log2FoldChange.x > 0,]
resKresL_samedir_down = resKresL_sig_overlap_samedir[resKresL_sig_overlap_samedir$log2FoldChange.x < 0,]
dim(resIresJ_samedir_up)
dim(resIresJ_samedir_down)
dim(resKresL_samedir_up)
dim(resKresL_samedir_down)

write.csv(resIresJ_samedir_up$Names, "resIresJ_samedir_up_210721.csv")
write.csv(resIresJ_samedir_down$Names, "resIresJ_samedir_down_210721.csv")
write.csv(resKresL_samedir_up$Names, "resKresL_samedir_up_210721.csv")
write.csv(resKresL_samedir_down$Names, "resKresL_samedir_down_210721.csv")
write.csv(resIresJ_sig_overlap_samedir, "resIresJ_sig_overlap_samedir_210721.csv")
write.csv(resKresL_sig_overlap_samedir, "resKresL_sig_overlap_samedir_210721.csv")

resIresJ_resKresL = merge(resKresL_sig_overlap_samedir, resIresJ_sig_overlap_samedir, by = "Names")
write.csv(resIresJ_resKresL$Names, "resIresJ_resKresL_210721.csv" )
resIresJ_resKresL_UP = resIresJ_resKresL[resIresJ_resKresL$log2FoldChange.x.x > 0,]
dim(resIresJ_resKresL_UP)
resIresJ_resKresL_DOWN = resIresJ_resKresL[resIresJ_resKresL$log2FoldChange.x.x < 0,]
dim(resIresJ_resKresL_DOWN)
write.csv(resIresJ_resKresL_UP$Names, "resIresJ_resKresL_UP_210721.csv")
write.csv(resIresJ_resKresL_DOWN$Names, "resIresJ_resKresL_DOWN_210721.csv")

resIresJ_resFresG = merge(resIresJ_sig_overlap, resFresG_sig_overlap, by = "Names")
dim(resFresG_sig_overlap)
dim(resIresJ_resFresG)[1]/dim(resFresG_sig_overlap)[1] #% of genes DE in zinc between both T and S that are DE in both S between conditions. 
resKresL_resFresG = merge(resKresL_sig_overlap, resFresG_sig_overlap, by = "Names")
dim(resKresL_resFresG)[1]/dim(resFresG_sig_overlap)[1] #% of genes DE in zinc between both T and S that are DE in both S between conditions. 

Supplementary_Figure_2A = venn.diagram(x= list(resKresL_sig_overlap_samedir$Names, resBresC_sig_overlap_samedir$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), rotation.degree=180, cex= 2) 
Supplementary_Figure_2A
pdf(file = "Supplementary_Figure_2A.pdf")
grid.draw(Supplementary_Figure_2A)
dev.off()
resBresC_resKresL = intersect(resBresC_sig_overlap_samedir$Names, resKresL_sig_overlap_samedir$Names)
length(resBresC_resKresL) #26...
write.csv(resBresC_resKresL, "resBresC_resKresL_150820.csv")

fish_BCKL = my_fisher_test(length(resBresC_resKresL), length(resBresC_sig_overlap_samedir$Names), length(resKresL_sig_overlap_samedir$Names), length(rownames(ddTxi_ind)))
table_tmp = data.frame("BCKL", NA, NA, fish_BCKL$p.value, fish_BCKL$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)
Fisher_Table

resIresJ_resFresG_notKnotL = resIresJ_resFresG[resIresJ_resFresG$Names %notin% resKresL_sig_overlap_samedir$Names,]
dim(resIresJ_resFresG_notKnotL)
venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = "VennDiagram_FG_KL_IJ_180820.png", output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, fill = c(col_mine, col_coast, col_zinc), width = 2000, height = 2000) 
Fig2F = venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = NULL, output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, fill = c(col_mine, col_coast, col_zinc), width = 2000, height = 2000) 
pdf(file = "Supplementary_Fig2F_1908.pdf")
grid.draw(Fig2F)
dev.off()

#Fisher tests
fish_IJ_samedir = my_fisher_test(length(resIresJ_sig_overlap_samedir$Names), length(resI_sig$Names), length(resJ_sig$Names), length(rownames(ddTxi_ind)))
fish_IJ = my_fisher_test(length(resIresJ_sig_overlap$Names), length(resI_sig$Names), length(resJ_sig$Names), length(rownames(ddTxi_ind)))
table_tmp = data.frame("IJ", fish_IJ$p.value, fish_IJ$estimate, fish_IJ_samedir$p.value, fish_IJ_samedir$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)
Fisher_Table

fish_KL_samedir = my_fisher_test(length(resKresL_sig_overlap_samedir$Names), length(resK_sig$Names), length(resL_sig$Names), length(rownames(ddTxi_ind)))
fish_KL = my_fisher_test(length(resKresL_sig_overlap$Names), length(resK_sig$Names), length(resL_sig$Names), length(rownames(ddTxi_ind)))
table_tmp = data.frame("KL", fish_KL$p.value, fish_KL$estimate, fish_KL_samedir$p.value, fish_KL_samedir$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)
Fisher_Table

dim(resIresJ_resKresL)
all_DGE = resIresJ_sig_overlap[resIresJ_sig_overlap$Names %in% resKresL_sig_overlap$Names,]
fish_IJKL = my_fisher_test(length(all_DGE$Names), length(resIresJ_sig_overlap$Names), length(resKresL_sig_overlap$Names), length(rownames(ddTxi_ind)))
fish_IJKL

table_tmp = data.frame("KL", fish_IJKL$p.value, fish_IJKL$estimate, fish_IJKL_samedir$p.value, fish_IJKL_samedir$estimate, stringsAsFactors = FALSE)
table_tmp
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)
Fisher_Table

#"The adaptive components of evolutuonary change (AP and CEC genes) make up only
resKresL_resBresC_overlap = resKresL_sig_overlap_samedir[resKresL_sig_overlap_samedir$Names %in% resBresC_sig_overlap_samedir$Names,]
(dim(resKresL_sig_overlap_samedir)[1]+dim(resBresC_sig_overlap_samedir)[1]-dim(resKresL_resBresC_overlap)[1])/dim(ddTxi)[1]
#% of the transcriptome"


#Within-comparison Venn diagrams
#Fig 1B
venn.diagram(x= list(resB_sig$Names, resC_sig$Names), category.names = c("", ""), filename = "VennDiagram_BvC_140820.png", output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_control), col = c(col_GRSA, col_PPBD), rotation.degree = 180, height = 2000, width = 2000, imagetype = "png") 
woof = venn.diagram(x= list(resK_sig$Names, resL_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_control), col = c(col_GRSA, col_PPBD), rotation.degree = 180, height = 2000, width = 2000, imagetype = "png") 
pdf(file = "VennDiagram_KL.pdf")
grid.draw(woof)
dev.off()
Fig2B = venn.diagram(x= list(resB_sig$Names, resC_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_control), col = c(col_GRSA, col_PPBD), rotation.degree = 180, height = 2000, width = 2000, imagetype = "png") 
library(grDevices)
pdf(file = "VennDiagram_BC.pdf")
grid.draw(Fig2B)
dev.off()
venn.diagram(x= list(resF_sig$Names, resG_sig$Names), category.names = c("", ""), filename = "VennDiagram_FvG_140820.png", output = TRUE, cat.pos = c(0,0), fill = c(col_zinc, col_zinc), col = c(col_GRSA, col_PPBD), rotation.degree = 180, width = 2000, height = 2000) 
Fig2D = venn.diagram(x= list(resF_sig$Names, resG_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), fill = c(col_zinc, col_zinc), col = c(col_GRSA, col_PPBD), rotation.degree = 180, width = 2000, height = 1000) 
pdf(file = "VennDiagram_FG.pdf")
grid.draw(Fig2D)
dev.off()

#resFresG
resFresG_sig_overlap$sign = resFresG_sig_overlap$log2FoldChange.x*resFresG_sig_overlap$log2FoldChange.y
resFresG_sig_overlap_samedir = resFresG_sig_overlap[resFresG_sig_overlap$sign > 0,]
dim(resFresG_sig_overlap_samedir)
write.csv(resFresG_sig_overlap_samedir$Names, "resFresG_sig_overlap_samedir_150820.csv")

resFresG_sig_overlap_samedir_UP = resFresG_sig_overlap_samedir[resFresG_sig_overlap_samedir$log2FoldChange.x > 0,]
resFresG_sig_overlap_samedir_DOWN = resFresG_sig_overlap_samedir[resFresG_sig_overlap_samedir$log2FoldChange.x < 0,]
write.csv(resFresG_sig_overlap_samedir_UP$Names, "resFresG_sig_overlap_samedir_up_150820.csv")
write.csv(resFresG_sig_overlap_samedir_DOWN$Names, "resFresG_sig_overlap_samedir_down_150820.csv")

venn.diagram(x= list(resI_sig$Names, resJ_sig$Names,resK_sig$Names, resL_sig$Names), category.names = c("S1", "S2", "T1", "T2"), filename = "VennDiagram_5_140820.png", output = TRUE, cat.pos = c(0,0,0,0), col = c(col_coast, col_coast, col_mine, col_mine)) 
Fig3A = venn.diagram(x= list(resI_sig$Names, resJ_sig$Names, resK_sig$Names, resL_sig$Names), category.names = c("S1", "S2", "T1", "T2"), filename = NULL, output = TRUE, cat.pos = c(0,0,0,0), col = c(col_coast, col_coast, col_mine, col_mine)) 
Fig3A
pdf(file = "VenDiagram_IJKL.pdf")
grid.draw(Fig3A)
dev.off()

#Others for supplement etc. 

fish_BC_samedir = my_fisher_test(length(resBresC_sig_overlap_samedir$Names), length(resB_sig$Names), length(resC_sig$Names), length(rownames(ddTxi2)))
fish_BC = my_fisher_test(length(resBresC_sig_overlap$Names), length(resB_sig$Names), length(resC_sig$Names), length(rownames(ddTxi2)))
table_tmp = data.frame("BC", fish_BC$p.value, fish_BC$estimate, fish_BC_samedir$p.value, fish_BC_samedir$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)

Fisher_Table



###############################################
#10 PCAs of resKresL [AP] / resBresC [CEC]
###############################################

vst_resKresL = vst_all[rownames(vst_all) %in% resKresL_sig_overlap$Names,]
PCA_resKresL = plotPCA12_centroids(vst_resKresL, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-80, 80),my_ylims = c(-40, 40), my_background_col = col_mine)
blork_resKresL.x = as.data.frame(PCA_resKresL$x)
blork_resKresL.x$groups = samples$Pop_Cond
pca.centroids_resKresL = aggregate(blork_resKresL.x[,1:2], list(Type = blork_resKresL.x$groups), mean)
PCA_resKresL = plotPCA12(vst_resKresL, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-80, 80),my_ylims = c(-40, 40), my_background_col = "lightgray")
PCA_resKresL_centroids_SupplementaryFigure = PCA_resKresL+
  geom_segment(aes(x = pca.centroids_resKresL[7,2], y = pca.centroids_resKresL[7,3], xend = pca.centroids_resKresL[8,2], yend = pca.centroids_resKresL[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[7,2], y = pca.centroids_resKresL[7,3], xend = pca.centroids_resKresL[8,2], yend = pca.centroids_resKresL[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[1,2], y = pca.centroids_resKresL[1,3], xend = pca.centroids_resKresL[2,2], yend = pca.centroids_resKresL[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5,col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[1,2], y = pca.centroids_resKresL[1,3], xend = pca.centroids_resKresL[2,2], yend = pca.centroids_resKresL[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot,  lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[5,2], y = pca.centroids_resKresL[5,3], xend = pca.centroids_resKresL[6,2], yend = pca.centroids_resKresL[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[5,2], y = pca.centroids_resKresL[5,3], xend = pca.centroids_resKresL[6,2], yend = pca.centroids_resKresL[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[3,2], y = pca.centroids_resKresL[3,3], xend = pca.centroids_resKresL[4,2], yend = pca.centroids_resKresL[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresL[3,2], y = pca.centroids_resKresL[3,3], xend = pca.centroids_resKresL[4,2], yend = pca.centroids_resKresL[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_point(data = pca_centroids_c, size = 4, shape = 21, fill = c(PPBD_plot, GRSA_plot, PPBD_plot, GRSA_plot), alpha = 1)
PCA_resKresL_centroids_SupplementaryFigure


vst_resBresC = vst_all[rownames(vst_all) %in% resBresC_sig_overlap_samedir$Names,]
plotPCA12(vst_resBresC, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-50, 50),my_ylims = c(-25, 25), my_background_col = gold_plot_light) 

###########################################################################################################



######################################################################################
######11 FOLD CHANGES/plots 
######################################################################################

#c) Pairwise correlation of estimated counts in the control/zinc

#Shrunken log fold changes...
lfcA = lfcShrink(dds2, contrast=c("Pop_Cond", "SAC", "BDC"), type="ashr")
lfcA$Names = rownames(lfcA)
lfcB = lfcShrink(dds2, contrast=c("Pop_Cond", "SAC", "GRC"), type = "ashr")
lfcB$Names = rownames(lfcB)
lfcC = lfcShrink(dds2, contrast=c("Pop_Cond", "BDC", "PPC"), type = "ashr")
lfcC$Names = rownames(lfcC)
lfcD = lfcShrink(dds2, contrast=c("Pop_Cond", "GRC", "PPC"), type = "ashr")
lfcD$Names = rownames(lfcD)
lfcE = lfcShrink(dds2, contrast=c("Pop_Cond", "SAZ", "BDZ"), type = "ashr")
lfcE$Names = rownames(lfcE)
lfcF = lfcShrink(dds2, contrast=c("Pop_Cond", "SAZ", "GRZ"), type = "ashr")
lfcF$Names = rownames(lfcF)
lfcG = lfcShrink(dds2, contrast=c("Pop_Cond", "BDZ", "PPZ"), type = "ashr")
lfcG$Names = rownames(lfcG)
lfcH = lfcShrink(dds2, contrast=c("Pop_Cond", "GRZ", "PPZ"), type = "ashr")
lfcH$Names = rownames(lfcH)
lfcI = lfcShrink(dds3, contrast=list("PopSA.CondZ"), type = "ashr")
lfcI$Names = rownames(lfcI)
lfcJ = lfcShrink(dds3, contrast=list("PopBD.CondZ"), type = "ashr")
lfcJ$Names = rownames(lfcJ)
lfcK = lfcShrink(dds3, contrast=list("PopGR.CondZ"), type = "ashr")
lfcK$Names = rownames(lfcK)
lfcL = lfcShrink(dds3, contrast=list("PopPP.CondZ"), type = "ashr")
lfcL$Names = rownames(lfcL)


###############
#11A - resBresC lfc plots

lfcA_resBresC_bothdir = lfcA[lfcA$Names %in% resBresC_sig_overlap$Names,]
lfcB_resBresC_bothdir = lfcB[lfcB$Names %in% resBresC_sig_overlap$Names,]
lfcC_resBresC_bothdir = lfcC[lfcC$Names %in% resBresC_sig_overlap$Names,]
lfcD_resBresC_bothdir = lfcD[lfcD$Names %in% resBresC_sig_overlap$Names,]
lfcE_resBresC_bothdir = lfcE[lfcE$Names %in% resBresC_sig_overlap$Names,]
lfcF_resBresC_bothdir = lfcF[lfcF$Names %in% resBresC_sig_overlap$Names,]
lfcG_resBresC_bothdir = lfcG[lfcG$Names %in% resBresC_sig_overlap$Names,]
lfcH_resBresC_bothdir = lfcH[lfcH$Names %in% resBresC_sig_overlap$Names,]

lfc_BvC_resBresC = merge(as.data.frame(lfcB_resBresC_bothdir), as.data.frame(lfcC_resBresC_bothdir), by = "Names")

Supplementary_Figure_1 = ggplot()+
  xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_fill_viridis_c()+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.5, data = lfc_BvC_resBresC)+
  ylab("T2 v S2 log2(Fold Change)")+
  xlab("T1 v S1 log2(Fold Change)")+
  theme(panel.border = element_rect(colour = "lightgray", fill = NA, size = 2), panel.background = element_blank())
Supplementary_Figure_1
pdf(file = "Supplementary_Figure_1_1.pdf")
grid.draw(Supplementary_Figure_1)
dev.off()

summary(lm(lfc_BvC_resBresC$log2FoldChange.y ~ lfc_BvC_resBresC$log2FoldChange.x))
#So not a very positive correlaiton...

#3B) LFC between treatments

lfcKvI = merge(as.data.frame(lfcK), as.data.frame(lfcI), by = "Names")
lfcIvJ = merge(as.data.frame(lfcI), as.data.frame(lfcJ), by = "Names")
lfcLvJ = merge(as.data.frame(lfcL), as.data.frame(lfcJ), by = "Names")
lfcKvL = merge(as.data.frame(lfcK), as.data.frame(lfcL), by = "Names")

#For each get lists for K, L and KvL genes. 
lfcIvJ_resKresL = lfcIvJ[lfcIvJ$Names %in% resKresL_sig_overlap$Names,]
lfcKvI_resKresL = lfcKvI[lfcKvI$Names %in% resKresL_sig_overlap$Names,]
lfcLvJ_resKresL = lfcLvJ[lfcLvJ$Names %in% resKresL_sig_overlap$Names,]
lfcKvL_resKresL = lfcKvL[lfcKvL$Names %in% resKresL_sig_overlap$Names,]

lfcIvJ_resIresJ = lfcIvJ[lfcIvJ$Names %in% resIresJ_sig_overlap$Names,]

resKresL_LFC_Plot = ggplot()+
  #xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcIvJ, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = lfcKvL_resKresL$log2FoldChange.x, y = lfcKvL_resKresL$log2FoldChange.y), binwidth=0.3)+
  scale_fill_viridis_c()+
  xlim(-10, 10)+
  ylim(-10, 10)+
  geom_abline()+
  
  ylab("T2 Control to Zinc - log2(Fold Change)")+
  xlab("T1 Control to Zinc - log2(Fold Change)")+
  #  theme(panel.background = element_rect(fill = coast_plot))
  theme(text= element_text(size = 20), panel.border = element_rect(colour = col_mine, fill = NA, size = 5), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))
resKresL_LFC_Plot

summary(lm(lfcKvL_resKresL$log2FoldChange.y ~ lfcKvL_resKresL$log2FoldChange.x))
#1.077, adjusted R^2 = 0.94
summary(lm(lfcKvI_resKresL$log2FoldChange.x ~ lfcKvI_resKresL$log2FoldChange.y))
#Slope = 0.886, adj R^2 = 0.857
summary(lm(lfcLvJ_resKresL$log2FoldChange.x ~ lfcLvJ_resKresL$log2FoldChange.y))
#Slope = 0.911, adj. R^2 = 0.828
summary(lm(lfcIvJ_resKresL$log2FoldChange.y ~ lfcIvJ_resKresL$log2FoldChange.x))
#Slope = 1.04, adj R^2 = 0.95

bin2dplot_lfC_resIresJ_sigoverlap = ggplot()+
  #xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcIvJ, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = lfcIvJ_resIresJ$log2FoldChange.x, y = lfcIvJ_resIresJ$log2FoldChange.y), binwidth=0.3)+
  scale_fill_viridis_c()+
  
  ylab("resI - S1)")+
  xlab("resJ - S2)")+
  xlim(-10, 10)+
  ylim(-10, 10)+
  geom_abline()+
  #  theme(panel.background = element_rect(fill = coast_plot))
  theme(panel.border = element_rect(colour = col_coast, fill = NA, size = 2), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))
bin2dplot_lfC_resIresJ_sigoverlap


###############################################
#########12 "resids"/|FC| values
###############################################

#called resids for historical reasons, can't be bothered to change now. 

######################
###12A - transcriptome wide resids

resids_countsA = abs(lfcA$log2FoldChange)
resids_countsB = abs(lfcB$log2FoldChange)
resids_countsC = abs(lfcC$log2FoldChange)
resids_countsD = abs(lfcD$log2FoldChange)
resids_countsE = abs(lfcE$log2FoldChange)
resids_countsF = abs(lfcF$log2FoldChange)
resids_countsG = abs(lfcG$log2FoldChange)
resids_countsH = abs(lfcH$log2FoldChange)

median(resids_countsA)
median(resids_countsB)
median(resids_countsC)
median(resids_countsD)
median(resids_countsE)
median(resids_countsF)
median(resids_countsG)
median(resids_countsH)

both_CZ = melt(data.frame(resids_countsA, resids_countsE, resids_countsB, resids_countsF, resids_countsC, resids_countsG, resids_countsD, resids_countsH), na.rm = TRUE)
both_CZ_reduced = melt(data.frame(resids_countsA, resids_countsE, resids_countsD, resids_countsH), na.rm = TRUE)
both_C = melt(data.frame(resids_countsA, resids_countsB, resids_countsC, resids_countsD), na.rm = TRUE)
both_Z = melt(data.frame(resids_countsE, resids_countsF, resids_countsG, resids_countsH), na.rm = TRUE)

both_CZ_wilcox = pairwise.wilcox.test(both_CZ$value, both_CZ$variable, p.adjust.method = "BH")
pairwiseMedianTest(value ~ variable, data=both_CZ) #Alternative stats test
write.csv(as.data.frame(both_CZ_wilcox$p.value), "both_CZ_wilcox.csv")

###################################
#Supplementary Figure - Transcriptome-wide |FC| in control
Fig5A_nocol = ggplot()+
  geom_boxplot(fill = ctrl_green, data = both_C, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15),  panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_nocol
pdf(file = "TranscriptomeWide_FC_Control.pdf")
grid.draw(Fig5A_nocol)
dev.off()

####################################################
#resBresC_sig_overlap_samedir; "resids" of these
####################################################
resids_countsA_resBresC = abs(lfcA[lfcA$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsB_resBresC = abs(lfcB[lfcB$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsC_resBresC = abs(lfcC[lfcC$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsD_resBresC = abs(lfcD[lfcD$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsE_resBresC = abs(lfcE[lfcE$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsF_resBresC = abs(lfcF[lfcF$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsG_resBresC = abs(lfcG[lfcG$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsH_resBresC = abs(lfcH[lfcH$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)

median(resids_countsA_resBresC)
median(resids_countsB_resBresC)
median(resids_countsC_resBresC)
median(resids_countsD_resBresC)
median(resids_countsE_resBresC)
median(resids_countsF_resBresC)
median(resids_countsG_resBresC)
median(resids_countsH_resBresC)

2^median(resids_countsA_resBresC)
2^median(resids_countsB_resBresC)
2^median(resids_countsC_resBresC)
2^median(resids_countsD_resBresC)
2^median(resids_countsE_resBresC)
2^median(resids_countsF_resBresC)
2^median(resids_countsG_resBresC)
2^median(resids_countsH_resBresC)

both_CZ_resBresC_all = melt(data.frame(resids_countsA_resBresC, resids_countsE_resBresC, resids_countsB_resBresC, resids_countsF_resBresC, resids_countsC_resBresC, resids_countsG_resBresC, resids_countsD_resBresC, resids_countsH_resBresC), na.rm = TRUE)
both_CZ_resBresC_wilcox = pairwise.wilcox.test(both_CZ_resBresC_all$value, both_CZ_resBresC_all$variable, p.adjust.method = "BH")
both_CZ_resBresC_wilcox
write.csv(as.data.frame(both_CZ_resBresC_wilcox$p.value), "both_CZ_resBresC_wilcox.csv")
pairwiseMedianTest(value ~ variable, data=both_CZ_resBresC_all) #Could also include this? 
both_CZ_resBresC = melt(data.frame(resids_countsA_resBresC, resids_countsB_resBresC, resids_countsC_resBresC, resids_countsD_resBresC), na.rm = TRUE)

boxplot_resBresC_FC_220721 = ggplot()+
  geom_boxplot(fill = ctrl_green, data = both_CZ_resBresC, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  # scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1))
boxplot_resBresC_FC_220721

both_CZ_resBresC_two = melt(data.frame(resids_countsA_resBresC,  resids_countsD_resBresC), na.rm = TRUE)

boxplot_resBresC_FC_220721 = ggplot()+
  geom_boxplot(fill = ctrl_green, data = both_CZ_resBresC_two, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  # scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15))+ylim(0,0.4)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1))
boxplot_resBresC_FC_220721

##########################
########## resKresL resids

#########Compare across conditions

resids_countsA_resKresL = abs(lfcA[lfcA$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsB_resKresL = abs(lfcB[lfcB$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsC_resKresL = abs(lfcC[lfcC$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsD_resKresL = abs(lfcD[lfcD$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsE_resKresL = abs(lfcE[lfcE$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsF_resKresL = abs(lfcF[lfcF$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsG_resKresL = abs(lfcG[lfcG$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsH_resKresL = abs(lfcH[lfcH$Names %in% resKresL_sig_overlap_samedir$Names,]$log2FoldChange)

median(resids_countsA_resKresL)
median(resids_countsB_resKresL)
median(resids_countsC_resKresL)
median(resids_countsD_resKresL)
median(resids_countsE_resKresL)
median(resids_countsF_resKresL)
median(resids_countsG_resKresL)
median(resids_countsH_resKresL)

2^median(resids_countsA_resKresL)
2^median(resids_countsB_resKresL)
2^median(resids_countsC_resKresL)
2^median(resids_countsD_resKresL)
2^median(resids_countsE_resKresL)
2^median(resids_countsF_resKresL)
2^median(resids_countsG_resKresL)
2^median(resids_countsH_resKresL)

resids_CZ_resKresL = melt(data.frame(resids_countsA_resKresL, resids_countsE_resKresL, resids_countsB_resKresL, resids_countsF_resKresL, resids_countsC_resKresL, resids_countsG_resKresL, resids_countsD_resKresL, resids_countsH_resKresL), na.rm = TRUE)
resids_CZ_resKresL_wilcox = pairwise.wilcox.test(resids_CZ_resKresL$value, resids_CZ_resKresL$variable, p.adjust.method = "BH")
resids_CZ_resKresL_wilcox
write.csv(as.data.frame(resids_CZ_resKresL_wilcox$p.value), "resids_CZ_resKresL_wilcox.csv")
pairwiseMedianTest(value ~ variable, data=resids_CZ_resKresL) #Could also include this? 

resids_CZ_resKresL_eco = melt(data.frame(resids_countsA_resKresL, resids_countsE_resKresL, resids_countsD_resKresL, resids_countsH_resKresL), na.rm  = TRUE)

Supplementary_Figure_5A = ggplot()+
  #  geom_rect(aes(xmin = 0.43, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  # geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #geom_half_violin(data = resids_CZ_resKresL_eco, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_coast, col_coast, col_mine, col_mine))+
  #geom_half_violin(data = resids_CZ_resKresL_eco, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(fill = c(ctrl_green, another_plot_purple, ctrl_green, another_plot_purple), data = resids_CZ_resKresL_eco, aes(x = variable, y=value), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("S1C v S2 (C)", "S1 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+ylim(0,0.825)+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2), text = element_text(size = 15))

Supplementary_Figure_5A
pdf(file = "Supplementary_Figure_5A.pdf")
grid.draw(Supplementary_Figure_5A)
dev.off()

#####Comparing resids between resKresL, resBresC and whole T in Z
tmp1 = as.data.frame(resids_countsH_resKresL)
names(tmp1) = "value"
tmp1$variable = "resKresL"
tmp2 = as.data.frame(resids_countsH_resBresC)
names(tmp2) = "value"
tmp2$variable = "resBresC"
tmp3 = as.data.frame(resids_countsH)
names(tmp3) = "value"
tmp3$variable = "T"
resids_Z_subset = rbind(tmp1, tmp2, tmp3)
resids_Z_subset_wilcox = pairwise.wilcox.test(resids_Z_subset$value, resids_Z_subset$variable, p.adjust.method = "BH")
resids_Z_subset_wilcox

##################################################################################################
#13 Plasticity - REV/RI classification, |FC| of plastic vs. nonplastic, parametric bootstrapping
##################################################################################################

#Need to calculate the Lo, Lp and La for various sets of genes. 
#Then include those where Lp and La both are above a certain cutoff...

#Previous authiors used TMM normalised counts measured in edgeR. 

######################################
#Simulating (a/b) ~ (c/b)

AC = runif(100) #Ancestral control
DC = runif(100) #Derived control
AS = runif(100) #Ancestral stressor
plot(log2(DC/AC), log2(AS/AC))
hist(log2(DC/AC)*log2(AS/AC))

AC = rnbinom(100, 2, 0.5) #Ancestral control
DC = rnbinom(100, 2, 0.5)  #Derived control
AS = rnbinom(100, 2, 0.5)  #Ancestral stressor

plot(AS-AC, DC-AS) 
plot(log2(DC/AC), log2(AS/AC))
hist(log2(DC/AC)*log2(AS/AC))
#Seems to be a bit of a positive correlation bias, zoiks. 

#######################################

norm_counts = as.data.frame(counts(dds, normalized=TRUE))
head(norm_counts)
colnames(norm_counts) = paste(dds$IND, dds$Cond, sep = "")
norm_counts = norm_counts[complete.cases(norm_counts),]
dim(norm_counts)

#Assigning Lo, Lp and La values. 
norm_counts$Lo = rowMeans(norm_counts[,c("BD11C", "BD12C", "BD1C", "SA6C", "SA7C", "SA8C")])
norm_counts$Lp = rowMeans(norm_counts[,c("BD11Z", "BD12Z", "BD1Z", "SA6Z", "SA7Z", "SA8Z")])
norm_counts$La = rowMeans(norm_counts[,c("GR10Z", "GR2Z", "GR8Z", "PP1Z", "PP4Z", "PP8Z")])
norm_counts$PC = norm_counts$Lp - norm_counts$Lo
norm_counts$GC = norm_counts$La - norm_counts$Lp
norm_counts$Lp_SD = rowSds(as.matrix(norm_counts[,c("BD11Z", "BD12Z", "BD1Z", "SA6Z", "SA7Z", "SA8Z")]))

cutoff = 0.2
#Get those above the cutoff...
norm_counts$Cutoff_GC = ifelse(abs(norm_counts$GC) > cutoff*norm_counts$Lo, "Y", "N")
norm_counts$Cutoff_PC = ifelse(abs(norm_counts$PC) > cutoff*norm_counts$Lo, "Y", "N")

#Is Lp or Lo closer to La? 
norm_counts$Distance = abs(norm_counts$Lo-norm_counts$La) - abs(norm_counts$Lp-norm_counts$La)
norm_counts$Distance = abs(norm_counts$Lo-norm_counts$La) - abs(norm_counts$Lp-norm_counts$La)
norm_counts$Closer = ifelse(norm_counts$Distance > 0, "Y", "N")

norm_counts_cut = norm_counts[norm_counts$Cutoff_GC == "Y" & norm_counts$Cutoff_PC == "Y",]
norm_counts_noPC = norm_counts[norm_counts$Cutoff_GC == "Y" & norm_counts$Cutoff_PC == "N",]
norm_counts_noGC = norm_counts[norm_counts$Cutoff_GC == "N" & norm_counts$Cutoff_PC == "Y",]

norm_counts_cut$DirStat = norm_counts_cut$GC*norm_counts_cut$PC
norm_counts_cut$SameDir = ifelse(norm_counts_cut$DirStat > 0, "Y", "N")

#Reinforcement: DirStat > 0 [changes in same dir]
#Reversion: DirStat < 0 [changes in oppo dir] & GC > 0.5PC [i.e. moving backwards more than 50% of PC change]
#Overshooting: DirStat < 0 [changes in oppo dir] & GC < 0.5PC [i.e. moving backwards less than 50% of PC change]
norm_counts_cut$Type = ifelse(norm_counts_cut$DirStat > 0, "RI", ifelse(abs(norm_counts_cut$GC) > 0.5*abs(norm_counts_cut$PC), "REV", "OVERSHOOT"))

#1) Whole transcriptome

100*table(norm_counts_cut$Type)[1]/nrow(norm_counts_cut) #% overshooting
100*table(norm_counts_cut$Type)[2]/nrow(norm_counts_cut) #% reversion
100*table(norm_counts_cut$Type)[3]/nrow(norm_counts_cut) #% reinforcement

head(norm_counts_cut)
T_table = table(norm_counts_cut$Type)
T_vec = c(0,0,0)
T_vec[2] = table(norm_counts_cut$Type)[1] #% overshooting
T_vec[1] = table(norm_counts_cut$Type)[2] #% reversion
T_vec[3] = table(norm_counts_cut$Type)[3] #% reinforcement
T_vec = data.frame(val = T_vec)
T_vec = cbind(T_vec, type = c("REV", "OVER", "RI"))
T_vec$type= factor(T_vec$type, levels = T_vec$type)

#|FC| values for whole transcriptome
lfcA_T_REV = lfcA[lfcA$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "REV",]),]
lfcE_T_REV = lfcE[lfcE$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "REV",]),]
lfcD_T_REV = lfcD[lfcD$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "REV",]),]
lfcH_T_REV = lfcH[lfcH$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "REV",]),]
resids_countsA_T_REV = abs(lfcA_T_REV$log2FoldChange)
resids_countsE_T_REV = abs(lfcE_T_REV$log2FoldChange)
resids_countsD_T_REV = abs(lfcD_T_REV$log2FoldChange)
resids_countsH_T_REV = abs(lfcH_T_REV$log2FoldChange)

lfcA_T_RI = lfcA[lfcA$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "RI",]),]
lfcE_T_RI = lfcE[lfcE$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "RI",]),]
lfcD_T_RI = lfcD[lfcD$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "RI",]),]
lfcH_T_RI = lfcH[lfcH$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "RI",]),]
resids_countsA_T_RI = abs(lfcA_T_RI$log2FoldChange)
resids_countsE_T_RI = abs(lfcE_T_RI$log2FoldChange)
resids_countsD_T_RI = abs(lfcD_T_RI$log2FoldChange)
resids_countsH_T_RI = abs(lfcH_T_RI$log2FoldChange)

lfcA_T_OVERSHOOT = lfcA[lfcA$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "OVERSHOOT",]),]
lfcE_T_OVERSHOOT = lfcE[lfcE$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "OVERSHOOT",]),]
lfcD_T_OVERSHOOT = lfcD[lfcD$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "OVERSHOOT",]),]
lfcH_T_OVERSHOOT = lfcH[lfcH$Names %in% rownames(norm_counts_cut[norm_counts_cut$Type == "OVERSHOOT",]),]
resids_countsA_T_OVERSHOOT = abs(lfcA_T_OVERSHOOT$log2FoldChange)
resids_countsE_T_OVERSHOOT = abs(lfcE_T_OVERSHOOT$log2FoldChange)
resids_countsD_T_OVERSHOOT = abs(lfcD_T_OVERSHOOT$log2FoldChange)
resids_countsH_T_OVERSHOOT = abs(lfcH_T_OVERSHOOT$log2FoldChange)

resids_countsH_T_RI = as.data.frame(resids_countsH_T_RI)
colnames(resids_countsH_T_RI) = "VALUE"
resids_countsH_T_RI$Names = "RI"
resids_countsH_T_REV = as.data.frame(resids_countsH_T_REV)
colnames(resids_countsH_T_REV) = "VALUE"
resids_countsH_T_REV$Names = "REV"
resids_countsH_T_OVERSHOOT = as.data.frame(resids_countsH_T_OVERSHOOT)
colnames(resids_countsH_T_OVERSHOOT) = "VALUE"
resids_countsH_T_OVERSHOOT$Names = "OVERSHOOT"

median(resids_countsH_T_RI$VALUE) #Reinforcement |FC|
median(resids_countsH_T_REV$VALUE) #Reversion |FC|
median(resids_countsH_T_OVERSHOOT$VALUE) #Overshooting |FC|


Plast_FC_T = rbind(resids_countsH_T_OVERSHOOT, resids_countsH_T_REV, resids_countsH_T_RI)

median(resids_countsH_T_RI$VALUE)
median(resids_countsH_T_REV$VALUE)
median(resids_countsH_T_OVERSHOOT$VALUE)

median(resids_countsE_T_RI)
median(resids_countsE_T_REV)
median(resids_countsE_T_OVERSHOOT)

Plast_FC_T_wilcox = pairwise.wilcox.test(Plast_FC_T$VALUE, Plast_FC_T$Names, p.adjust.method = "BH")
Plast_FC_T_wilcox
pairwiseMedianTest(value ~ variable, data=both_CZ) #Might be better to do this, really. 

#So...how do we compare these?

#2) reKresL_samedir genes (AP genes)
norm_counts_resKresL = norm_counts_cut[rownames(norm_counts_cut) %in% resKresL_sig_overlap_samedir$Names,]
norm_counts_resKresL = norm_counts_cut[rownames(norm_counts_cut) %in% resKresL_sig_overlap_samedir$Names,]
dim(norm_counts_resKresL)

100*table(norm_counts_resKresL$Type)[1]/nrow(norm_counts_resKresL) #% overshooting
100*table(norm_counts_resKresL$Type)[2]/nrow(norm_counts_resKresL) #% reversion
100*table(norm_counts_resKresL$Type)[3]/nrow(norm_counts_resKresL) #% reinforcement

library(cowplot)

KL_Table = table(norm_counts_resKresL$Type)
KL_vec = c(0,0,0)
KL_vec[2] = table(norm_counts_resKresL$Type)[1] #% overshooting
KL_vec[1] = table(norm_counts_resKresL$Type)[2] #% reversion
KL_vec[3] = table(norm_counts_resKresL$Type)[3] #% reinforcement
KL_vec = data.frame(val = KL_vec)
KL_vec = cbind(KL_vec, type = c("REV", "OVER", "RI"))
KL_vec$type= factor(KL_vec$type, levels = KL_vec$type)


binom.test(x = c((KL_Table[1]+KL_Table[3]), KL_Table[2]), p = (T_table[1]+T_table[3])/(T_table[1]+T_table[2]+T_table[3]))

#Are "adaptive" plastic changes overrepresented in this gene set?M 

#|FC| values for whole transcriptome
lfcA_KL_REV = lfcA[lfcA$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "REV",]),]
lfcE_KL_REV = lfcE[lfcE$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "REV",]),]
lfcD_KL_REV = lfcD[lfcD$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "REV",]),]
lfcH_KL_REV = lfcH[lfcH$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "REV",]),]
resids_countsA_KL_REV = abs(lfcA_KL_REV$log2FoldChange)
resids_countsE_KL_REV = abs(lfcE_KL_REV$log2FoldChange)
resids_countsD_KL_REV = abs(lfcD_KL_REV$log2FoldChange)
resids_countsH_KL_REV = abs(lfcH_KL_REV$log2FoldChange)

lfcA_KL_RI = lfcA[lfcA$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "RI",]),]
lfcE_KL_RI = lfcE[lfcE$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "RI",]),]
lfcD_KL_RI = lfcD[lfcD$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "RI",]),]
lfcH_KL_RI = lfcH[lfcH$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "RI",]),]
resids_countsA_KL_RI = abs(lfcA_KL_RI$log2FoldChange)
resids_countsE_KL_RI = abs(lfcE_KL_RI$log2FoldChange)
resids_countsD_KL_RI = abs(lfcD_KL_RI$log2FoldChange)
resids_countsH_KL_RI = abs(lfcH_KL_RI$log2FoldChange)

lfcA_KL_OVERSHOOT = lfcA[lfcA$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "OVERSHOOT",]),]
lfcE_KL_OVERSHOOT = lfcE[lfcE$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "OVERSHOOT",]),]
lfcD_KL_OVERSHOOT = lfcD[lfcD$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "OVERSHOOT",]),]
lfcH_KL_OVERSHOOT = lfcH[lfcH$Names %in% rownames(norm_counts_resKresL[norm_counts_resKresL$Type == "OVERSHOOT",]),]
resids_countsA_KL_OVERSHOOT = abs(lfcA_KL_OVERSHOOT$log2FoldChange)
resids_countsE_KL_OVERSHOOT = abs(lfcE_KL_OVERSHOOT$log2FoldChange)
resids_countsD_KL_OVERSHOOT = abs(lfcD_KL_OVERSHOOT$log2FoldChange)
resids_countsH_KL_OVERSHOOT = abs(lfcH_KL_OVERSHOOT$log2FoldChange)


norm_counts_noPC_resKresL = norm_counts_noPC[rownames(norm_counts_noPC) %in% resKresL_sig_overlap_samedir$Names,]
dim(norm_counts_noPC_resKresL)
lfcA_KL_NOANCPLAST = lfcA[lfcA$Names %in% rownames(norm_counts_noPC_resKresL),]
lfcE_KL_NOANCPLAST = lfcE[lfcE$Names %in% rownames(norm_counts_noPC_resKresL),]
lfcD_KL_NOANCPLAST = lfcD[lfcD$Names %in% rownames(norm_counts_noPC_resKresL),]
lfcH_KL_NOANCPLAST = lfcH[lfcH$Names %in% rownames(norm_counts_noPC_resKresL),]

norm_counts_noGC_resKresL = norm_counts_noGC[rownames(norm_counts_noGC) %in% resKresL_sig_overlap_samedir$Names,]
dim(norm_counts_noGC_resKresL)

resids_countsA_KL_NOANCPLAST = abs(lfcA_KL_NOANCPLAST$log2FoldChange)
resids_countsE_KL_NOANCPLAST = abs(lfcE_KL_NOANCPLAST$log2FoldChange)
resids_countsD_KL_NOANCPLAST = abs(lfcD_KL_NOANCPLAST$log2FoldChange)
resids_countsH_KL_NOANCPLAST = abs(lfcH_KL_NOANCPLAST$log2FoldChange)
resids_countsH_KL_NOANCPLAST = as.data.frame(resids_countsH_KL_NOANCPLAST)
colnames(resids_countsH_KL_NOANCPLAST) = "VALUE"
resids_countsH_KL_NOANCPLAST$Names = "NOANCPLAST"
dim(resids_countsH_KL_NOANCPLAST)
median(resids_countsH_KL_NOANCPLAST$VALUE)
#0.12. So the same as all the rest of them, then. 

resids_countsH_KL_RI = as.data.frame(resids_countsH_KL_RI)
colnames(resids_countsH_KL_RI) = "VALUE"
resids_countsH_KL_RI$Names = "RI"
resids_countsH_KL_REV = as.data.frame(resids_countsH_KL_REV)
colnames(resids_countsH_KL_REV) = "VALUE"
resids_countsH_KL_REV$Names = "REV"
resids_countsH_KL_OVERSHOOT = as.data.frame(resids_countsH_KL_OVERSHOOT)
colnames(resids_countsH_KL_OVERSHOOT) = "VALUE"
resids_countsH_KL_OVERSHOOT$Names = "OVERSHOOT"

Plast_FC_KL = rbind(resids_countsH_KL_OVERSHOOT, resids_countsH_KL_REV, resids_countsH_KL_RI)
Plast_FC_KL_wilcox = pairwise.wilcox.test(Plast_FC_KL$VALUE, Plast_FC_KL$Names, p.adjust.method = "BH")
Plast_FC_KL_wilcox
median(resids_countsH_KL_RI$VALUE)
median(resids_countsH_KL_REV$VALUE)
median(resids_countsH_KL_OVERSHOOT$VALUE)

temp_Plast_FC_KL = Plast_FC_KL
temp_Plast_FC_KL$Names = "PLAST"
PlastNoPlast_FC_KL = rbind(temp_Plast_FC_KL, resids_countsH_KL_NOANCPLAST)
PlastNoPlast_FC_KL_wilcox = pairwise.wilcox.test(PlastNoPlast_FC_KL$VALUE, PlastNoPlast_FC_KL$Names, p.adjust.method = "BH")
PlastNoPlast_FC_KL_wilcox

KL_PlastNoPlast_Boxplot = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_KL, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("P", "NP"))+
  xlab("")+ylab("|FC|")+ylim(0,0.625)+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2), text = element_text(size = 15))
KL_PlastNoPlast_Boxplot

median(resids_countsE_KL_RI)
median(resids_countsE_KL_REV)
median(resids_countsE_KL_OVERSHOOT)

#Ok...I think this probably needs some more thought. 

pairwiseMedianTest(value ~ variable, data=both_CZ) #Might be better to do this, really. 

#3) resBresC_samedir genes (CEC genes)

norm_counts_resBresC = norm_counts_cut[rownames(norm_counts_cut) %in% resBresC_sig_overlap_samedir$Names,]
norm_counts_resBresC = norm_counts_cut[rownames(norm_counts_cut) %in% resBresC_sig_overlap_samedir$Names,]

100*table(norm_counts_resBresC$Type)[1] #% overshooting
100*table(norm_counts_resBresC$Type)[2] #% reversion
100*table(norm_counts_resBresC$Type)[3] #% reinforcement


BC_Table = table(norm_counts_resBresC$Type)
BC_Table

BC_vec = c(0,0,0)
BC_vec[2] = table(norm_counts_resBresC$Type)[1] #% overshooting
BC_vec[1] = table(norm_counts_resBresC$Type)[2] #% reversion
BC_vec[3] = table(norm_counts_resBresC$Type)[3] #% reinforcement
BC_vec = data.frame(val = BC_vec)
BC_vec = cbind(BC_vec, type = c("REV", "OVER", "RI"))
BC_vec$type= factor(BC_vec$type, levels = BC_vec$type)

binom.test(x = c((BC_Table[1]+BC_Table[3]), BC_Table[2]), p = (KL_Table[1]+KL_Table[3])/(KL_Table[1]+KL_Table[2]+KL_Table[3]))
binom.test(x = c((BC_Table[3]), BC_Table[2]+BC_Table[1]), p = (KL_Table[3])/(KL_Table[1]+KL_Table[2]+KL_Table[3]))
binom.test(x = c((BC_Table[3]), BC_Table[2]+BC_Table[1]), p = (T_table[3])/(T_table[1]+T_table[2]+T_table[3]))

#|FC| values for whole transcriptome
lfcA_BC_REV = lfcA[lfcA$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "REV",]),]
lfcE_BC_REV = lfcE[lfcE$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "REV",]),]
lfcD_BC_REV = lfcD[lfcD$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "REV",]),]
lfcH_BC_REV = lfcH[lfcH$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "REV",]),]
resids_countsA_BC_REV = abs(lfcA_BC_REV$log2FoldChange)
resids_countsE_BC_REV = abs(lfcE_BC_REV$log2FoldChange)
resids_countsD_BC_REV = abs(lfcD_BC_REV$log2FoldChange)
resids_countsH_BC_REV = abs(lfcH_BC_REV$log2FoldChange)

lfcA_BC_RI = lfcA[lfcA$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "RI",]),]
lfcE_BC_RI = lfcE[lfcE$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "RI",]),]
lfcD_BC_RI = lfcD[lfcD$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "RI",]),]
lfcH_BC_RI = lfcH[lfcH$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "RI",]),]
resids_countsA_BC_RI = abs(lfcA_BC_RI$log2FoldChange)
resids_countsE_BC_RI = abs(lfcE_BC_RI$log2FoldChange)
resids_countsD_BC_RI = abs(lfcD_BC_RI$log2FoldChange)
resids_countsH_BC_RI = abs(lfcH_BC_RI$log2FoldChange)

lfcA_BC_OVERSHOOT = lfcA[lfcA$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "OVERSHOOT",]),]
lfcE_BC_OVERSHOOT = lfcE[lfcE$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "OVERSHOOT",]),]
lfcD_BC_OVERSHOOT = lfcD[lfcD$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "OVERSHOOT",]),]
lfcH_BC_OVERSHOOT = lfcH[lfcH$Names %in% rownames(norm_counts_resBresC[norm_counts_resBresC$Type == "OVERSHOOT",]),]
resids_countsA_BC_OVERSHOOT = abs(lfcA_BC_OVERSHOOT$log2FoldChange)
resids_countsE_BC_OVERSHOOT = abs(lfcE_BC_OVERSHOOT$log2FoldChange)
resids_countsD_BC_OVERSHOOT = abs(lfcD_BC_OVERSHOOT$log2FoldChange)
resids_countsH_BC_OVERSHOOT = abs(lfcH_BC_OVERSHOOT$log2FoldChange)

norm_counts_noPC_resBresC = norm_counts_noPC[rownames(norm_counts_noPC) %in% resBresC_sig_overlap_samedir$Names,]
dim(norm_counts_noPC_resBresC)
lfcA_BC_NOANCPLAST = lfcA[lfcA$Names %in% rownames(norm_counts_noPC_resBresC),]
lfcE_BC_NOANCPLAST = lfcE[lfcE$Names %in% rownames(norm_counts_noPC_resBresC),]
lfcD_BC_NOANCPLAST = lfcD[lfcD$Names %in% rownames(norm_counts_noPC_resBresC),]
lfcH_BC_NOANCPLAST = lfcH[lfcH$Names %in% rownames(norm_counts_noPC_resBresC),]

norm_counts_noGC_resBresC = norm_counts_noGC[rownames(norm_counts_noGC) %in% resBresC_sig_overlap_samedir$Names,]
dim(norm_counts_noGC_resBresC)

resids_countsA_BC_NOANCPLAST = abs(lfcA_BC_NOANCPLAST$log2FoldChange)
resids_countsE_BC_NOANCPLAST = abs(lfcE_BC_NOANCPLAST$log2FoldChange)
resids_countsD_BC_NOANCPLAST = abs(lfcD_BC_NOANCPLAST$log2FoldChange)
resids_countsH_BC_NOANCPLAST = abs(lfcH_BC_NOANCPLAST$log2FoldChange)
resids_countsH_BC_NOANCPLAST = as.data.frame(resids_countsH_BC_NOANCPLAST)
colnames(resids_countsH_BC_NOANCPLAST) = "VALUE"
resids_countsH_BC_NOANCPLAST$Names = "NOANCPLAST"
dim(resids_countsH_BC_NOANCPLAST)
#SO that's more like it (woah 371 that seems like a lot)
dim(resBresC_sig_overlap_samedir)
#...that can't be right...
median(resids_countsH_BC_NOANCPLAST$VALUE)
#0.12. So the same as all the rest of them, then. 

resids_countsH_BC_RI = as.data.frame(resids_countsH_BC_RI)
colnames(resids_countsH_BC_RI) = "VALUE"
resids_countsH_BC_RI$Names = "RI"
resids_countsH_BC_REV = as.data.frame(resids_countsH_BC_REV)
colnames(resids_countsH_BC_REV) = "VALUE"
resids_countsH_BC_REV$Names = "REV"
resids_countsH_BC_OVERSHOOT = as.data.frame(resids_countsH_BC_OVERSHOOT)
colnames(resids_countsH_BC_OVERSHOOT) = "VALUE"
resids_countsH_BC_OVERSHOOT$Names = "OVERSHOOT"

Plast_FC_BC = rbind(resids_countsH_BC_OVERSHOOT, resids_countsH_BC_REV, resids_countsH_BC_RI)
Plast_FC_BC_wilcox = pairwise.wilcox.test(Plast_FC_BC$VALUE, Plast_FC_BC$Names, p.adjust.method = "BH")
Plast_FC_BC_wilcox
temp_Plast_FC_BC = Plast_FC_BC
temp_Plast_FC_BC$Names = "PLAST"
PlastNoPlast_FC_BC = rbind(temp_Plast_FC_BC, resids_countsH_BC_NOANCPLAST)
PlastNoPlast_FC_BC_wilcox = pairwise.wilcox.test(PlastNoPlast_FC_BC$VALUE, PlastNoPlast_FC_BC$Names, p.adjust.method = "BH")
PlastNoPlast_FC_BC_wilcox

PlastNoPlast_FC_BC$Set = "CEC"
PlastNoPlast_FC_KL$Set = "AP"
PlastNoPlast_FC_BCKL = rbind(PlastNoPlast_FC_BC, PlastNoPlast_FC_KL)

T_vec$title = "Transcriptome-wide"
T_barplot = ggplot(data = T_vec, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity')+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(text = element_text(size=15), axis.title.y = element_text(size = 12))+
theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
      strip.text = element_text(size = 10))+facet_grid(.~title)

T_barplot

BC_vec$title = "CEC"
BC_barplot = ggplot(data = BC_vec, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity')+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(text = element_text(size=15), axis.title.y = element_text(size = 12))+
theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
      strip.text = element_text(size = 10))+facet_grid(.~title)


BC_barplot

KL_vec$title = "AP"
KL_barplot = ggplot(data = KL_vec, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity')+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(text = element_text(size=15), axis.title.y = element_text(size = 12))+
theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
      strip.text = element_text(size = 10))+facet_grid(.~title)

KL_barplot

BCKL_PlastNoPlast_Boxplot = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_BCKL,fill = another_plot_purple, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 10))+
  scale_x_discrete(labels = c("P", "NP"))+
  facet_grid(.~Set)+
  xlab("")+ylab("|FC|")+ylim(0,0.625)+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = c("solid", "solid")), 
        text = element_text(size = 5), strip.text = element_text(size = 10))
BCKL_PlastNoPlast_Boxplot

plot_grid(T_barplot, KL_barplot, BC_barplot, BCKL_PlastNoPlast_Boxplot)

median(resids_countsH_BC_RI$VALUE)
median(resids_countsH_BC_REV$VALUE)
median(resids_countsH_BC_OVERSHOOT$VALUE)
median(resids_countsE_BC_RI)
median(resids_countsE_BC_REV)
median(resids_countsE_BC_OVERSHOOT)

##########################################################################################
##Parametric Bootstrapping (a la Ho and Zhang 2019)
##########################################################################################
plasticity_parametric_bootstrap = function(Lo_row, Lp_ow, La_row){
  Lo_mean = rowMeans(Lo_row)
  Lo_sd = rowSds(as.matrix(Lo_row))/sqrt(length(Lo_row))
  Lp_mean = rowMeans(Lp_row)
  Lp_sd = rowSds(as.matrix(Lp_row))/sqrt(length(Lp_row))
  La_mean = rowMeans(La_row)
  La_sd = rowSds(as.matrix(La_row))/sqrt(length(La_row))
  #So we have our stats for that gene row...
  #Then we need to 1,000 times calculate the actual numbers
  result_vec = rep("NA", 100)
  for (i in seq(1,100,1)){
    Lo = rnorm(1, mean = Lo_mean, sd = Lo_sd)
    Lp = rnorm(1, mean = Lp_mean, sd = Lp_sd)
    La = rnorm(1, mean = La_mean, sd = La_sd)
    
    #Assume we've already done the cutoff on these I guess...
    PC = Lp-Lo
    GC = La-Lp
    if (PC*GC > 1){
      result = "RI"
    }else{
      if (abs(GC) < 0.5*abs(PC)){
        result = "OVERSHOOT"
      }else{
        result = "REV"
      }
    }
    result_vec[i] = result 
    result_table = sort(table(result_vec), decreasing = T)
    if (result_table[1] > 95){
      final_result = names(result_table[1])
    }else{
      final_result = "UNCLASSIFIABLE"
    }
  }
  return(final_result)
}

#classified_vec = rep(NA, nrow(norm_counts_cut))

###############################################
#Running the simulations. Note - this will take about an hour to run. 
#for (i in 1:nrow(norm_counts_cut)){
#  Lo_row = norm_counts[i,c("BD11C", "BD12C", "BD1C", "SA6C", "SA7C", "SA8C")]
#  Lp_row = norm_counts[i,c("BD11Z", "BD12Z", "BD1Z", "SA6Z", "SA7Z", "SA8Z")]
#  La_row = norm_counts[i,c("GR10Z", "GR2Z", "GR8Z", "PP1Z", "PP4Z", "PP8Z")]
#  classified_vec[i] = plasticity_parametric_bootstrap(Lo_row, Lp_row, La_row)
#}

norm_counts_cut$temp_perm_calc = classified_vec
perm_table = table(norm_counts_cut$temp_perm_calc)
classifiable_perm = perm_table[1]+perm_table[2]+perm_table[3]
perm_table[1]/classifiable_perm
perm_table[2]/classifiable_perm
perm_table[3]/classifiable_perm

norm_counts_cut_param = norm_counts_cut[norm_counts_cut$temp_perm_calc != "UNCLASSIFIABLE",]
dim(norm_counts_cut_param)

norm_counts_resKresL_param = norm_counts_cut_param[rownames(norm_counts_cut_param) %in% resKresL_sig_overlap_samedir$Names,]
dim(norm_counts_resKresL_param)

table(norm_counts_resKresL_param$Type)[1]/nrow(norm_counts_resKresL_param) #% overshooting
table(norm_counts_resKresL_param$Type)[2]/nrow(norm_counts_resKresL_param) #% reversion
table(norm_counts_resKresL_param$Type)[3]/nrow(norm_counts_resKresL_param) #% reinforcement

norm_counts_resBresC_param = norm_counts_cut_param[rownames(norm_counts_cut_param) %in% resBresC_sig_overlap_samedir$Names,]

table(norm_counts_resBresC_param$Type)[1]/nrow(norm_counts_resBresC_param) #% overshooting
table(norm_counts_resBresC_param$Type)[2]/nrow(norm_counts_resBresC_param) #% reversion
table(norm_counts_resBresC_param$Type)[3]/nrow(norm_counts_resBresC_param) #% reinforcement

#############


########################################################################
################ Analysis from genotypes (popgen stats, phylogeny)
########################################################################

#SNP phylogeny libraries

library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(phytools)

#Figure 1A - SNP Phylogeny
snp_tree = read.tree(file = "C:/Users/Danie/DEA_SNPhylo.bs.tree")
ggtree(snp_tree)+geom_tiplab()+geom_text(aes(label=node))+geom_treescale()
snp_tree$edge.length
#To find which branch length to halve...
plotTree(snp_tree)
edgelabels(snp_tree$eedge.length)
snp_tree = reroot(snp_tree, 18, (0.03134)/2)
ggtree(snp_tree)+geom_tiplab()+geom_text(aes(label=node))
snp_tree2 = groupClade(snp_tree, .node=c(22, 20, 17, 15))

woof = ggtree(snp_tree2, aes(color=group), size = 1.5)+
  scale_color_manual(values=c("black", col_mine, col_coast, col_coast, col_mine))+
  theme(legend.position = "none")+
  geom_cladelabel(node=22, label = "T2", offset = 0.02, align = T, col = col_PPBD, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_cladelabel(node=20, label = "S2", offset = 0.02, align = T, col = col_PPBD, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_cladelabel(node=17, label = "S1", offset = 0.02, align = T, col = col_GRSA, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_cladelabel(node=15, label = "T1", offset = 0.02, align = T, col = col_GRSA, barsize = 2, offset.text = 0.005, fontsize = 5)+
  geom_treescale()+
  
  xlim(0,0.5)
#  geom_text(aes(label=node))+
#  geom_tiplab()
woof
gridExtra::grid.arrange(flip(woof, 19, 14) %>% flip(22, 20), ncol=1)

pdf(file = woof, "C:/Users/Danie/Documents/RNA_Seq/Figure1A_270321.pdf")
dev.off()


#Genotypes
library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)

whatever = read.csv(file = "C:/Users/Danie/0000.vcf.shared.BD.vcf.Tajima.D", sep = "\t", header = T)
par(mfrow = c(1,1))
hist(whatever$TajimaD)
plot(1,1)

par(mfrow = c(2,2))
PP_TD = read.csv(file = "C:/Users/Danie/0000.vcf.gz.DEA.DEA.biallelic.tajd.Tajima.D", sep = "\t", header = T)
GR_TD = read.csv(file = "C:/Users/Danie/0001.vcf.gz.DEA.DEA.biallelic.tajd.Tajima.D", sep = "\t", header = T)
BD_TD = read.csv(file = "C:/Users/Danie/0002.vcf.gz.DEA.DEA.biallelic.tajd.Tajima.D", sep = "\t", header = T)
SA_TD = read.csv(file = "C:/Users/Danie/0003.vcf.gz.DEA.DEA.biallelic.tajd.Tajima.D", sep = "\t", header = T)

mean(PP_TD$TajimaD)
var(PP_TD$TajimaD)
hist(BD_TD$TajimaD)
hist(GR_TD$TajimaD)
hist(PP_TD$TajimaD)
hist(SA_TD$TajimaD)
PP_TD$Group = "PP"
BD_TD$Group = "BD"
GR_TD$Group = "GR"
SA_TD$Group = "SA"
PP_TD = PP_TD[,c("Group", "TajimaD")]
BD_TD = BD_TD[,c("Group", "TajimaD")]
GR_TD = GR_TD[,c("Group", "TajimaD")]
SA_TD = SA_TD[,c("Group", "TajimaD")]
all = rbind(PP_TD, BD_TD, GR_TD, SA_TD)
all = na.omit(all)
pairwise.t.test(all$TajimaD, all$Group)

#Need to read in those pi values...I guess? I'm not sure how to test significant differneces here given how many are 0. 

t.test(value ~ Group, data = Data, var.equal=FALSE, conf.level=0.95) #So this is a Welch's t-test.
#Pi_stuff

GR_Pi = read.csv("C:/Users/Danie/0000.vcf.gz.DEA.gz.site.pi.sites.pi.total_adjust", sep = "\t")
SA_Pi = read.csv("C:/Users/Danie/0001.vcf.gz.DEA.gz.site.pi.sites.pi.total_adjust", sep = "\t")
PP_Pi = read.csv("C:/Users/Danie/0002.vcf.gz.DEA.gz.site.pi.sites.pi.total_adjust", sep = "\t")
BD_Pi = read.csv("C:/Users/Danie/0003.vcf.gz.DEA.gz.site.pi.sites.pi.total_adjust", sep = "\t")
GR_Pi$Pop = "GR"
SA_Pi$Pop = "SA"
PP_Pi$Pop = "PP"
BD_Pi$Pop = "BD"
t.test(GR_Pi$PI, SA_Pi$PI, var.equal=FALSE)
t.test(GR_Pi$PI, PP_Pi$PI, var.equal=FALSE)


All_Pi = rbind(GR_Pi, SA_Pi, PP_Pi, BD_Pi)
remove(GR_Pi)
remove(SA_Pi)
remove(PP_Pi)
remove(BD_Pi)
#So...is this the way to do it? 
pairwise.t.test(All_Pi$PI, All_Pi$Pop, var.equal=FALSE, conf.level=0.95) #So this is a Welch's t-test.
#They are all significantly different from each other...

par(mfrow=c(2,1))
hist(GR_Pi$PI)
hist(SA_Pi$PI)

##################################################
###################### Even more Miscellaneous
#################################################
#Excessive DGE in Coastal Zinc problems...
#A) What % of genes are DE - 
length(resI_sig$Names)/length(ddTxi2)
length(resJ_sig$Names)/length(ddTxi2)
#Approximately 48% Seems like a lot?

#Is the GC content of the DE genes different? 
par(mfrow = c(3,1))
GC_resI = GC_content_Filter2[GC_content_Filter2$V1 %in% resI_sig$Names,]
GC_resJ = GC_content_Filter2[GC_content_Filter2$V1 %in% resJ_sig$Names,]
GC_resK = GC_content_Filter2[GC_content_Filter2$V1 %in% resK_sig$Names,]
GC_resL = GC_content_Filter2[GC_content_Filter2$V1 %in% resL_sig$Names,]
hist(GC_content_Filter2$V2, breaks = 100)
hist(GC_resI$V2, breaks = 50)
hist(GC_resJ$V2, breaks = 50)
hist(GC_resK$V2, breaks = 50)
hist(GC_resL$V2, breaks = 50)
t.test(GC_content_Filter2$V2, GC_resI$V2)
t.test(GC_content_Filter2$V2, GC_resJ$V2)
t.test(GC_content_Filter2$V2, GC_resK$V2)
t.test(GC_content_Filter2$V2, GC_resL$V2)
t.test(GC_resI$V2, GC_resK$V2) #There are a variety of things where they're DE. 
t.test(GC_resJ$V2, GC_resL$V2) # 

#Are they all going one way? 
dim(resI_sig[resI_sig$log2FoldChange > 0,]) #So 16,235 upregulated (59%)
dim(resI_sig[resI_sig$log2FoldChange < 0,]) #11,372 are downregulated (41%)
dim(resJ_sig[resJ_sig$log2FoldChange > 0,]) #16,311 are upregulated  (59%)
dim(resJ_sig[resJ_sig$log2FoldChange < 0,]) #11,372 are downregulated (41%)

#Seems fine...
