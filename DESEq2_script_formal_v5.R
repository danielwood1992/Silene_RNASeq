#Happy tidy code as of 20:25 13/08/20

#0 LOADING IN THE DATA, FILTERING
remove(list=ls())
`%notin%` = Negate(`%in%`)

setwd(dir = "/Users/Danie/Documents/RNA_Seq")
#load(".Rdata")

#Tables etc.:
Fisher_Table = data.frame("Comparison", "All_pvalue (0 = <2.2e-16)", "All_oddsratio", "Samedir_pvalue", "Samedir_oddsratio", stringsAsFactors = FALSE)

#Colours
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

tmp = col2rgb(col_mine)
mine_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.2)
tmp = col2rgb(col_coast)
coast_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.2)
tmp = col2rgb(col_GRSA)
GRSA_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.6)
tmp = col2rgb(col_PPBD)
PPBD_plot = rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha = 0.6)



#Load relevant libraries
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
#library(rcompanion)

#Get sample name 
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
#Current approach: Either a >50bp match to latifolia genome, and/or a >50bp match to an Embryophyte

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

###########################################################################################################

#1 OVERALL LEVELS OF VARIATION

###########################################################################################################

#Examine the overall patterns of variation
#DECISION: Can use vst() or rlog() to transform data for visualisation etc. - rlog seems to be more useful for when
#library sizes/sequencing depths vary significantly. Examine these.
par(mfrow  = c(1,2))
boxplot(log2(counts(dds, normalized = FALSE)))
boxplot(log2(counts(dds, normalized = TRUE)))
#Still a bit weird some have such a narrow IQR

#20/05/20 - Library sizes appear mostly similar - using vst() transformation is appropriate. 
vst_all = vst(dds)

#DECISION: Given the very high levels of variance due to Coastal C->Z effects, the top 500 most variable genes
  #are likely to be dominated by the set of genes that vary between these populations, making them less informative
  #for the other populations. Therefore decided to visualise PCAs using the entire set of genes (200,000+)
#Function to plot first 2 principal components, from correct colours. 

object = vst_all
intgroup1 = "Cond"
intgroup2 = "Geog"
intgroup3 = "Ecotype"
ntop = 1000000

#PCA figure functions...

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
      geom_point(size = 10, stroke=4, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                          100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  + 
      scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) + 
      coord_fixed() + theme(legend.position = "none")+
      theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1) #, panel.grid.major=element_line(colour="lightgrey"))
    
  }

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
  
  
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[3,4]
    return(d)
  }
  return(pca_thing)
  ggplot(data = d, aes_string(x = "PC1", y = "PC2")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 10, stroke=4, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                                                                                                                             100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) + 
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1) #, panel.grid.major=element_line(colour="lightgrey"))
  
}



plotPCA34 <- function (object, intgroup1 = "condition", ntop = 200000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "Geog") 
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
    attr(d, "percentVar") <- percentVar[3,4]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4")) + 
    geom_point(size = 6, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC3: ", round(percentVar[3] * 
                                                  100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC4: ", round(percentVar[4]*100), "% variance")) + 
    coord_fixed() + theme(legend.position = "none")
}

#Test plots...


#Some separation by geography, some separation by ecotype...
plotPCA34(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000) #Control and zinc explains a lot of variance...

##############################################################################################################################################

#2 COUNTS IN CONTROL/ZINC

##############################################################################################################################################

#a) PCAs

#PCA just for control plants

dds_C = dds[,dds$Cond == "C"]
vst_C = vst(dds_C)
dds_Z = dds[,dds$Cond == "Z"]
vst_Z = vst(dds_Z)

pdf()

#Seems fine
plotPCA12(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-300, 200), my_ylims = c(-150, 150), my_background_col = "lightgray")
woof = plotPCA12_centroids(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-300, 200), my_ylims = c(-150, 150), my_background_col = "lightgray")
blork.x = as.data.frame(woof$x)
blork.x$groups = samples$Pop_Cond
pca.centroids = aggregate(blork.x[,1:2], list(Type = blork.x$groups), mean)
#Ok great, so we have these now? So need to a) plot? b) Calculate? 
sqrt((pca.centroids[3,2]-pca.centroids[5,2])^2 + (pca.centroids[3,3] - pca.centroids[5,3])^2)
pca.centroids[3,2]-pca.centroids[5,2] #PC1
pca.centroids[3,3]-pca.centroids[5,3] #PC2
sqrt((pca.centroids[4,2]-pca.centroids[6,2])^2 + (pca.centroids[4,3] - pca.centroids[6,3])^2)
pca.centroids[4,2]-pca.centroids[6,2] #PC1
pca.centroids[4,3]-pca.centroids[6,3] #PC2

#If we actually wanted to see if they got more similar...probably needed to do PCA of just the mine pops? 

#Not 100% sure how best to report this? 

#https://stackoverflow.com/questions/36208909/how-to-calculate-centroids-in-pca - assuming this is right?
#This seems to weight each PC equally though...but who knows what's right. 

#Great...so a) how do we add these to PCAs, and b) how do we calculate? 

Fig2E = plotPCA12(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-280, 120), my_ylims = c(-95, 110), my_background_col = "lightgray")
Fig2E
pdf(file = "Fig2E_1908_1700.pdf")
grid.draw(Fig2E)
dev.off()

plotPCA12(vst_C, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_xlims = c(-100,100), my_ylims = c(-100,75), my_background_col = ctrl_col) #Control and zinc explains a lot of variance...
Fig2A = plotPCA12(vst_C, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_xlims = c(-75,90), my_ylims = c(-75,50), my_background_col = ctrl_col) #Control and zinc explains a lot of variance...
Fig2A


pdf(file = "Fig2A_1908_1700.pdf")
grid.draw(Fig2A)
dev.off()
plotPCA34(vst_C, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000) #Control and zinc explains a lot of variance...

plotPCA12(vst_Z, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_xlims = c(-250, 250), my_ylims = c(-150, 150), my_background_col = col_zinc)
Fig2C = plotPCA12(vst_Z, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_xlims = c(-250, 220), my_ylims = c(-120, 110), my_background_col = col_zinc) #Control and zinc explains a lot of variance...
Fig2C
pdf(file = "Fig2C_1908_1700.pdf")
grid.draw(Fig2C)
dev.off()

plotPCA34(vst_Z, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000) #Control and zinc explains a lot of variance...


##################################################################################

#b) Within zinc/control comparisons: identify DGE genes. 

ddTxi3 = ddTxi2
design(ddTxi3) <- ~Pop_Cond
dds2 = DESeq(ddTxi3)
resultsNames(dds2)

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

#Adjust p-value < 0.05. 
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

write.csv(resA_sig$Names, "resA_sig_150820")
write.csv(resB_sig$Names, "resB_sig_150820")
write.csv(resC_sig$Names, "resC_sig_150820")
write.csv(resD_sig$Names, "resD_sig_150820")
write.csv(resE_sig$Names, "resE_sig_150820")
write.csv(resF_sig$Names, "resF_sig_150820")
write.csv(resG_sig$Names, "resG_sig_150820")
write.csv(resH_sig$Names, "resH_sig_150820")

dim(resA_sig)
dim(resB_sig)
dim(resC_sig)
dim(resD_sig)
dim(resE_sig)
dim(resF_sig)
dim(resG_sig)
dim(resH_sig)

#Significance tests for overlapping genes in sets of interest...
resBresC_sig_overlap = merge(as.data.frame(resB_sig), as.data.frame(resC_sig), by = "Names")
resFresG_sig_overlap = merge(as.data.frame(resF_sig), as.data.frame(resG_sig), by = "Names")
resAresE_sig_overlap = merge(as.data.frame(resA_sig), as.data.frame(resE_sig), by = "Names")
resDresH_sig_overlap = merge(as.data.frame(resD_sig), as.data.frame(resH_sig), by = "Names")
resAresD_sig_overlap = merge(as.data.frame(resA_sig), as.data.frame(resD_sig), by = "Names")

dim(resBresC_sig_overlap)
dim(resFresG_sig_overlap)
dim(resAresE_sig_overlap)
dim(resDresH_sig_overlap)

#Signs in same direction in each comparison

#Within-comparison Venn diagrams
#Fig 1B
venn.diagram(x= list(resB_sig$Names, resC_sig$Names), category.names = c("", ""), filename = "VennDiagram_BvC_140820.png", output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_control), col = c(col_GRSA, col_PPBD), rotation.degree = 180, height = 2000, width = 2000, imagetype = "png") 
Fig2B = venn.diagram(x= list(resB_sig$Names, resC_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_control), col = c(col_GRSA, col_PPBD), rotation.degree = 180, height = 2000, width = 2000, imagetype = "png") 
library(grDevices)
pdf(file = "Fig2B_1908.pdf")
grid.draw(Fig2B)
dev.off()

#Fig 1D
venn.diagram(x= list(resF_sig$Names, resG_sig$Names), category.names = c("", ""), filename = "VennDiagram_FvG_140820.png", output = TRUE, cat.pos = c(0,0), fill = c(col_zinc, col_zinc), col = c(col_GRSA, col_PPBD), rotation.degree = 180, width = 2000, height = 2000) 
Fig2D = venn.diagram(x= list(resF_sig$Names, resG_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), fill = c(col_zinc, col_zinc), col = c(col_GRSA, col_PPBD), rotation.degree = 180, width = 2000, height = 1000) 
pdf(file = "Fig2D_1908.pdf")
grid.draw(Fig2D)
dev.off()

#Additional venn diagrams for supplement
venn.diagram(x= list(resA_sig$Names, resE_sig$Names), category.names = c("A", "E"), filename = "VennDiagram_AvE.png", output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_zinc), col = c(col_coast, col_coast), rotation.degree = 180, height = 1000, width = 2000) 
venn.diagram(x= list(resD_sig$Names, resH_sig$Names), category.names = c("D", "H"), filename = "VennDiagram_DvH.png", output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_zinc), col = c(col_mine, col_mine), height = 1000, width = 2000, rotation.degree=180) 
venn.diagram(x= list(resA_sig$Names, resB_sig$Names, resC_sig$Names, resD_sig$Names), category.names = c("A", "B", "C", "D"), filename = "VennDiagram_ABCD.png", output = TRUE, cat.pos = c(0,0,0,0), height = 2000, width = 2000, rotation.degree=180) 
venn.diagram(x= list(resE_sig$Names, resF_sig$Names, resG_sig$Names, resH_sig$Names), category.names = c("E", "F", "G", "H"), filename = "VennDiagram_EFGH.png", output = TRUE, cat.pos = c(0,0,0,0), height = 2000, width = 2000, rotation.degree=180) 

#To find resKresL/resIresJ/resFresG Venn Diagram, search for #I2P9

#Genes upregulated/downregulated - print to file...

#Note: UP = higher in coast
#DOWN = higher in mine. 

#resBresC
resBresC_sig_overlap$sign = resBresC_sig_overlap$log2FoldChange.x*resBresC_sig_overlap$log2FoldChange.y
resBresC_sig_overlap_samedir = resBresC_sig_overlap[resBresC_sig_overlap$sign > 0,]
dim(resBresC_sig_overlap_samedir)

write.csv(resBresC_sig_overlap_samedir$Names, "resBresC_sig_overlap_samedir_150820.csv")
resBresC_sig_overlap_samedir_UP = resBresC_sig_overlap_samedir[resBresC_sig_overlap_samedir$log2FoldChange.x > 0,]
resBresC_sig_overlap_samedir_DOWN = resBresC_sig_overlap_samedir[resBresC_sig_overlap_samedir$log2FoldChange.x < 0,]
dim(resBresC_sig_overlap_samedir_UP)
dim(resBresC_sig_overlap_samedir_DOWN)
write.csv(resBresC_sig_overlap_samedir_UP$Names, "resBresC_sig_overlap_samedir_up_150820.csv")
write.csv(resBresC_sig_overlap_samedir_DOWN$Names, "resBresC_sig_overlap_samedir_down_150820.csv")

#resFresG
resFresG_sig_overlap$sign = resFresG_sig_overlap$log2FoldChange.x*resFresG_sig_overlap$log2FoldChange.y
resFresG_sig_overlap_samedir = resFresG_sig_overlap[resFresG_sig_overlap$sign > 0,]
dim(resFresG_sig_overlap_samedir)
write.csv(resFresG_sig_overlap_samedir$Names, "resFresG_sig_overlap_samedir_150820.csv")

resFresG_sig_overlap_samedir_UP = resFresG_sig_overlap_samedir[resFresG_sig_overlap_samedir$log2FoldChange.x > 0,]
resFresG_sig_overlap_samedir_DOWN = resFresG_sig_overlap_samedir[resFresG_sig_overlap_samedir$log2FoldChange.x < 0,]
write.csv(resFresG_sig_overlap_samedir_UP$Names, "resFresG_sig_overlap_samedir_up_150820.csv")
write.csv(resFresG_sig_overlap_samedir_DOWN$Names, "resFresG_sig_overlap_samedir_down_150820.csv")


#Others for supplement etc. 
resAresD_sig_overlap$sign = resAresD_sig_overlap$log2FoldChange.x*resAresD_sig_overlap$log2FoldChange.y
resAresD_sig_overlap_samedir = resAresD_sig_overlap[resAresD_sig_overlap$sign > 0,]
dim(resAresD_sig_overlap_samedir)
#write.csv(resAresD_sig_overlap_samedir$Names, "resAresD_sig_overlap_samedir_1607")

resDresH_sig_overlap$Sign = resDresH_sig_overlap$log2FoldChange.x*resDresH_sig_overlap$log2FoldChange.y
resDresH_sig_overlap_samedir = resDresH_sig_overlap[resDresH_sig_overlap$Sign > 0,]
dim(resDresH_sig_overlap_samedir)

resAresE_sig_overlap$Sign = resAresE_sig_overlap$log2FoldChange.x*resAresE_sig_overlap$log2FoldChange.y
resAresE_sig_overlap_samedir = resAresE_sig_overlap[resAresE_sig_overlap$Sign > 0,]
dim(resAresE_sig_overlap_samedir)

#Function - fisher.test where you can specify the overlap (i.e. when you want the "overlap" to be those in the same direction only)
#Do both...
my_fisher_test = function(overlap, set1, set2, background) {
  test_m = matrix(c(overlap, set1-overlap, set2-overlap, background-(set1+set2-overlap)), nrow = 2)
  x = fisher.test(test_m, alternative = "greater")
  return(x)
}

####19:15 13/08/20

#fish_AE_samedir = my_fisher_test(length(resAresE_sig_overlap_samedir$Names), length(resA_sig$Names), length(resE_sig$Names), length(rownames(ddTxi2)))
#fish_AE = my_fisher_test(length(resAresE_sig_overlap$Names), length(resA_sig$Names), length(resE_sig$Names), length(rownames(ddTxi2)))
#table_tmp = data.frame("AE", fish_AE$p.value, fish_AE$estimate, fish_AE_samedir$p.value, fish_AE_samedir$estimate, stringsAsFactors = FALSE)
#names(table_tmp) = names(Fisher_Table)
#Fisher_Table = rbind(Fisher_Table, table_tmp)
#Fisher_Table

#fish_AD_samedir = my_fisher_test(length(resAresD_sig_overlap_samedir$Names), length(resA_sig$Names), length(resD_sig$Names), length(rownames(ddTxi2)))
#fish_AD = my_fisher_test(length(resAresD_sig_overlap$Names), length(resA_sig$Names), length(resD_sig$Names), length(rownames(ddTxi2)))
#table_tmp = data.frame("AD", fish_AD$p.value, fish_AD$estimate, fish_AD_samedir$p.value, fish_AD_samedir$estimate, stringsAsFactors = FALSE)
#names(table_tmp) = names(Fisher_Table)
#Fisher_Table = rbind(Fisher_Table, table_tmp)

fish_BC_samedir = my_fisher_test(length(resBresC_sig_overlap_samedir$Names), length(resB_sig$Names), length(resC_sig$Names), length(rownames(ddTxi2)))
fish_BC = my_fisher_test(length(resBresC_sig_overlap$Names), length(resB_sig$Names), length(resC_sig$Names), length(rownames(ddTxi2)))
table_tmp = data.frame("BC", fish_BC$p.value, fish_BC$estimate, fish_BC_samedir$p.value, fish_BC_samedir$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)

fish_FG_samedir = my_fisher_test(length(resFresG_sig_overlap_samedir$Names), length(resF_sig$Names), length(resG_sig$Names), length(rownames(ddTxi2)))
fish_FG = my_fisher_test(length(resFresG_sig_overlap$Names), length(resF_sig$Names), length(resG_sig$Names), length(rownames(ddTxi2)))
table_tmp = data.frame("FG", fish_FG$p.value, fish_FG$estimate, fish_FG_samedir$p.value, fish_FG_samedir$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)

#fish_DH_samedir = my_fisher_test(length(resDresH_sig_overlap_samedir$Names), length(resD_sig$Names), length(resH_sig$Names), length(rownames(ddTxi2)))
#fish_DH = my_fisher_test(length(resDresH_sig_overlap$Names), length(resD_sig$Names), length(resH_sig$Names), length(rownames(ddTxi2)))
#table_tmp = data.frame("DH", fish_DH$p.value, fish_DH$estimate, fish_DH_samedir$p.value, fish_DH_samedir$estimate, stringsAsFactors = FALSE)
#names(table_tmp) = names(Fisher_Table)
#Fisher_Table = rbind(Fisher_Table, table_tmp)

Fisher_Table

######################################################################################

#c) Pairwise correlation of estimated counts in the control/zinc

#Note - these are called "resids" for historical reasons...

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

2^median(resids_countsA)
2^median(resids_countsB)
2^median(resids_countsC)
2^median(resids_countsD)
2^median(resids_countsE)
2^median(resids_countsF)
2^median(resids_countsG)
2^median(resids_countsH)

both_CZ = melt(data.frame(resids_countsA, resids_countsE, resids_countsB, resids_countsF, resids_countsC, resids_countsG, resids_countsD, resids_countsH), na.rm = TRUE)
dim(both_CZ)

both_C = melt(data.frame(resids_countsA, resids_countsB, resids_countsC, resids_countsD), na.rm = TRUE)
both_Z = melt(data.frame(resids_countsE, resids_countsF, resids_countsG, resids_countsH), na.rm = TRUE)

both_CZ_wilcox = pairwise.wilcox.test(both_CZ$value, both_CZ$variable, p.adjust.method = "BH")
pairwiseMedianTest(value ~ variable, data=both_CZ) #Could also include this? 
write.csv(as.data.frame(both_CZ_wilcox$p.value), "both_CZ_wilcox.csv")

#Fig 5A

Fig5A_old = ggplot()+
  #  geom_half_violin(data = both_CZ, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  geom_rect(aes(xmin = 0.43, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA)+
  geom_rect(aes(xmin = 4.5, xmax = 6.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD)+
  geom_rect(aes(xmin = 6.5, xmax = 8.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
    scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
#  geom_half_violin(data = both_CZ, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, col_zinc2, ctrl_col, col_zinc2, ctrl_col, col_zinc2, ctrl_col, col_zinc2))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 vs S2 (Z)", "T1 vs S1 (C)", "T1 vs S1 (Z)", "T2 vs S2 (C)", "T2 vs S2 (Z)", "T1 vs T2 (C)", "T1 vs T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,5.5)
Fig5A_old
pdf(file = "Fig5A_old.pdf")
grid.draw(Fig5A_old)
dev.off()

Fig5A_all = ggplot()+
  #  geom_half_violin(data = both_CZ, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  geom_rect(aes(xmin = 0.47, xmax = 2.4, ymin = 0, ymax = 5.4),alpha = 0.6, fill = NA, col = col_coast)+
  geom_rect(aes(xmin = 2.6, xmax = 4.4, ymin = 0, ymax = 5.4),alpha = 0.6, fill = NA, col = col_GRSA)+
  geom_rect(aes(xmin = 4.6, xmax = 6.4, ymin = 0, ymax = 5.4),alpha = 0.6, fill = NA, col = col_PPBD)+
  geom_rect(aes(xmin = 6.6, xmax = 8.4, ymin = 0, ymax = 5.4),alpha = 0.6, fill = NA, col = col_mine)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, col_zinc2, ctrl_col, col_zinc2, ctrl_col, col_zinc2, ctrl_col, col_zinc2))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S1 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("abs(log2(Shrunken FC))")+ylim(0,5.5)

Fig5A_all
#No this just looks rubbish...
#I think plot separately is the only way...

Fig5A_V0 = ggplot()+
  #  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
  geom_rect(aes(xmin = 0.43, xmax = 1.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast, col = col_coast)+
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA, col = col_GRSA)+
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD, col = col_PPBD)+
  geom_rect(aes(xmin = 3.5, xmax = 4.575, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine, col = col_mine)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_V0
pdf(file = "Fig5A_V0.pdf")
grid.draw(Fig5A_V0)
dev.off()

Fig5A_nocol = ggplot()+
  #  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
#  geom_rect(aes(xmin = 0.43, xmax = 1.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast, col = col_coast)+
#  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA, col = col_GRSA)+
#  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD, col = col_PPBD)+
#  geom_rect(aes(xmin = 3.5, xmax = 4.575, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine, col = col_mine)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_nocol
pdf(file = "Fig5A_nocol.pdf")
grid.draw(Fig5A_nocol)
dev.off()

Fig5B_nocol = ggplot()+
  geom_boxplot(data = both_Z, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (Z)", "T1 v S1 (Z)", "T2 v S2 (Z)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,5.5)
Fig5B_nocol
pdf(file = "Fig5B_nocol.pdf")
grid.draw(Fig5B_nocol)
dev.off()

Fig5A_V1 = ggplot()+
#  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
  geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_coast)+
  geom_rect(aes(xmin = 1.6, xmax = 2.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_GRSA)+
  geom_rect(aes(xmin = 2.6, xmax = 3.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_PPBD)+
  geom_rect(aes(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_mine)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
 # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_V1
pdf(file = "Fig5A_V1.pdf")
grid.draw(Fig5A_V1)
dev.off()

Fig5A_V2 = ggplot()+
  #  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  geom_rect(aes(xmin = 0.4, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_control, col = "black")+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  #scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_fill_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_V2
pdf(file = "Fig5A_V2.pdf")
grid.draw(Fig5A_V2)
dev.off()

Fig5A_V3 = ggplot()+
  #  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  #geom_rect(aes(xmin = 0.4, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_control, col = "black")+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  #scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_fill_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_V3
pdf(file = "Fig5A_V3.pdf")
grid.draw(Fig5A_V3)
dev.off()

Fig5A_V4 = ggplot()+
  #  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+

  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_V4
pdf(file = "Fig5A_V4.pdf")
grid.draw(Fig5A_V4)
dev.off()

Fig5A_V0 = ggplot()+
  #  geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
  geom_rect(aes(xmin = 0.43, xmax = 1.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast, col = col_coast)+
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA, col = col_GRSA)+
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD, col = col_PPBD)+
  geom_rect(aes(xmin = 3.5, xmax = 4.575, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine, col = col_mine)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_C, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_C, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (C)", "S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5A_V0
pdf(file = "Fig5A_V0.pdf")
grid.draw(Fig5A_V0)
dev.off()





Fig5B_V1 = ggplot()+
  #  geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
  geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_coast)+
  geom_rect(aes(xmin = 1.6, xmax = 2.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_GRSA)+
  geom_rect(aes(xmin = 2.6, xmax = 3.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_PPBD)+
  geom_rect(aes(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0.5),alpha = 0.6, fill = NA, col = col_mine)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_Z, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_ZZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(zinc_col, zinc_col, zinc_col, zinc_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (Z)", "S1 v S2 (Z)", "T1 v S1 (Z)", "T2 v S2 (Z)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5B_V1
pdf(file = "Fig5B_V1.pdf")
grid.draw(Fig5B_V1)
dev.off()

Fig5B_V2 = ggplot()+
  #  geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  geom_rect(aes(xmin = 0.4, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_zinc, col = "black")+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_Z, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_ZZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  #scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_fill_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (Z)", "S1 v S2 (Z)", "T1 v S1 (Z)", "T2 v S2 (Z)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5B_V2
pdf(file = "Fig5B_V2.pdf")
grid.draw(Fig5B_V2)
dev.off()

Fig5B_V3 = ggplot()+
  #  geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  #geom_rect(aes(xmin = 0.4, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_zinc, col = "black")+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_Z, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_ZZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  #scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  scale_fill_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (Z)", "S1 v S2 (Z)", "T1 v S1 (Z)", "T2 v S2 (Z)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5B_V3
pdf(file = "Fig5B_V3.pdf")
grid.draw(Fig5B_V3)
dev.off()

Fig5B_V4 = ggplot()+
  #  geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  # geom_half_violin(data = both_Z, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_Z, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_ZZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(zinc_col, zinc_col, zinc_col, zinc_col))+
  scale_colour_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (Z)", "S1 v S2 (Z)", "T1 v S1 (Z)", "T2 v S2 (Z)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,0.5)
Fig5B_V4
pdf(file = "Fig5B_V4.pdf")
grid.draw(Fig5B_V4)
dev.off()


Fig5B = ggplot()+
  #  geom_half_violin(data = both_CZ, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_Z, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  #geom_boxplot(data = both_CZ, aes(x = variable, y=value, fill=variable), width = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, GRSA_plot, PPBD_plot, col_mine))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15), panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  scale_x_discrete(labels = c("S1 vs S2 (Z)", "S1 v S2 (Z)", "T1 v S1 (Z)", "T2 v S2 (Z)", "T1 v T2 (Z)"))+scale_y_continuous(minor_breaks = seq(0,5,0.1), breaks = seq(0,5,0.1))+
  xlab("")+ylab("|FC|")+ylim(0,5.5)
Fig5B


#resBresC_sig_overlap_samedir; "resids" of these
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

#Fig 5B

#Shared BC genes
p = ggplot()+
  geom_rect(aes(xmin = 0.42, xmax = 1.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA)+
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD)+
  geom_rect(aes(xmin = 3.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
#  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
#  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_resBresC, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("abs (shrunken log2FC)")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
p

Fig5C_nocol = ggplot()+
#  geom_rect(aes(xmin = 0.42, xmax = 1.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
 # geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA)+
#  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD)+
#  geom_rect(aes(xmin = 3.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_resBresC, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
 # scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
Fig5C_nocol
pdf(file = "Fig5C_nocol.pdf")
grid.draw(Fig5C_nocol)
dev.off()
#Supplementary Fig 1A

Supplementary_Figure_3A = ggplot()+
#  geom_rect(aes(xmin = 0.45, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
#  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA)+
#  geom_rect(aes(xmin = 4.5, xmax = 6.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD)+
#  geom_rect(aes(xmin = 6.5, xmax = 8.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_resBresC_all, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
#  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S2 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
Supplementary_Figure_3A
pdf(file = "Supplementary_Figure_3A.pdf")
grid.draw(Supplementary_Figure_3A)
dev.off()

vst_resBresC = vst_all[rownames(vst_all) %in% resBresC_sig_overlap_samedir$Names,]
plotPCA12(vst_resBresC, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-50, 50),my_ylims = c(-25, 25), my_background_col = gold_plot_light) 

#Genes DE in B or C but not both. 

#NOT Including genes that do "overlap" but in opposite directions...
resB_notC = resB_sig[resB_sig$Names %notin% resBresC_sig_overlap$Names,]
dim(resB_notC)
resC_notB = resC_sig[resC_sig$Names %notin% resBresC_sig_overlap$Names,]
dim(resC_notB)
notshared_resBresC = c(resB_notC$Names, resC_notB$Names)
length(notshared_resBresC)
notshared_resBresC = unique(notshared_resBresC)
length(notshared_resBresC) #Cool. So then what are these guys doing? 
#write.csv(notshared_resBresC, "notshared_resBresC_200720")

resids_countsA_notshared_resBresC = abs(lfcA[lfcA$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsB_notshared_resBresC = abs(lfcB[lfcB$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsC_notshared_resBresC = abs(lfcC[lfcC$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsD_notshared_resBresC = abs(lfcD[lfcD$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsE_notshared_resBresC = abs(lfcE[lfcE$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsF_notshared_resBresC = abs(lfcF[lfcF$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsG_notshared_resBresC = abs(lfcG[lfcG$Names %in% notshared_resBresC,]$log2FoldChange)
resids_countsH_notshared_resBresC = abs(lfcH[lfcH$Names %in% notshared_resBresC,]$log2FoldChange)

median(resids_countsA_notshared_resBresC)
median(resids_countsB_notshared_resBresC)
median(resids_countsC_notshared_resBresC)
median(resids_countsD_notshared_resBresC)
median(resids_countsE_notshared_resBresC)
median(resids_countsF_notshared_resBresC)
median(resids_countsG_notshared_resBresC)
median(resids_countsH_notshared_resBresC)

2^median(resids_countsA_notshared_resBresC)
2^median(resids_countsB_notshared_resBresC)
2^median(resids_countsC_notshared_resBresC)
2^median(resids_countsD_notshared_resBresC)
2^median(resids_countsE_notshared_resBresC)
2^median(resids_countsF_notshared_resBresC)
2^median(resids_countsG_notshared_resBresC)
2^median(resids_countsH_notshared_resBresC)

both_CZ_notshared_resBresC = melt(data.frame(resids_countsA_notshared_resBresC, resids_countsE_notshared_resBresC, resids_countsB_notshared_resBresC, resids_countsF_notshared_resBresC, resids_countsC_notshared_resBresC, resids_countsG_notshared_resBresC, resids_countsD_notshared_resBresC, resids_countsH_notshared_resBresC), na.rm = TRUE)
both_CZ_notshared_resBresC_wilcox = pairwise.wilcox.test(both_CZ_notshared_resBresC$value, both_CZ_notshared_resBresC$variable, p.adjust.method = "BH")
both_CZ_notshared_resBresC_wilcox
write.csv(as.data.frame(both_CZ_notshared_resBresC_wilcox$p.value), "both_CZ_notshared_resBresC_wilcox.csv")
pairwiseMedianTest(value ~ variable, data=both_CZ_notshared_resBresC) #Could also include this? 

both_CZ_notshared_resBresC_withineco = melt(data.frame(resids_countsA_notshared_resBresC, resids_countsE_notshared_resBresC, resids_countsD_notshared_resBresC, resids_countsH_notshared_resBresC), na.rm = TRUE)

#Shared BC genes
p = ggplot()+
 # geom_half_violin(data = both_CZ_notshared_resBresC_withineco, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_coast, col_mine, col_mine))+
#  geom_half_violin(data = both_CZ_notshared_resBresC_withineco, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resBresC_withineco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, ))+
  scale_x_discrete(labels = c("SACvBDC (A)", "SAZvBDZ (E)", "GRCvSAC (B)", "GRZvSAZ (F)", "PPCvBDC (C)", "PPZvBDZ (G)", "GRCvPPC (D)", "GRZvPPZ (H)"))+
  xlab("")+ylab("Relative difference (y-x)/(y+x)")+theme(panel.background = element_rect(fill = "thistle2"))+ylim(0,1.5)
p

p = ggplot()+
  geom_rect(aes(xmin = 0.43, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA)+
  geom_rect(aes(xmin = 4.5, xmax = 6.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD)+
  geom_rect(aes(xmin = 6.5, xmax = 8.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resBresC, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S2 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("abs (shrunken log2FC)")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
p

Supplementary_Figure_3B = ggplot()+
#  geom_rect(aes(xmin = 0.43, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
#  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_GRSA)+
#  geom_rect(aes(xmin = 4.5, xmax = 6.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_PPBD)+
#  geom_rect(aes(xmin = 6.5, xmax = 8.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resBresC, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
#  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,6.5)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S2 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
Supplementary_Figure_3B
pdf(file = "Supplementary_Figure_3B.pdf")
grid.draw(Supplementary_Figure_3B)
dev.off()

vst_notshared_resBresC = vst_all[rownames(vst_all) %in% notshared_resBresC,]
plotPCA12(vst_notshared_resBresC, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-100, 100),my_ylims = c(-100, 100), my_background_col = purple_plot) 


##############################################################
#I/J/K/L

#3B) LFC between control and zinc treatments

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

#Graphs with KL overlapping FCs...

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

write.csv(resI_sig$Names, "resI_sig_150820.csv")
write.csv(resJ_sig$Names, "resJ_sig_150820.csv")
write.csv(resK_sig$Names, "resK_sig_150820.csv")
write.csv(resL_sig$Names, "resL_sig_150820.csv")

#Overlapping sets

resIresJ_sig_overlap = merge(as.data.frame(resI_sig), as.data.frame(resJ_sig), by = "Names")
resKresL_sig_overlap = merge(as.data.frame(resK_sig), as.data.frame(resL_sig), by = "Names")

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

write.csv(resIresJ_samedir_up$Names, "resIresJ_samedir_up_150820.csv")
write.csv(resIresJ_samedir_down$Names, "resIresJ_samedir_down_150820.csv")
write.csv(resKresL_samedir_up$Names, "resKresL_samedir_up_150820.csv")
write.csv(resKresL_samedir_down$Names, "resKresL_samedir_down_150820.csv")

write.csv(resIresJ_sig_overlap_samedir, "resIresJ_sig_overlap_samedir_150820.csv")
write.csv(resKresL_sig_overlap_samedir, "resKresL_sig_overlap_samedir_150820.csv")

resIresJ_resKresL = merge(resKresL_sig_overlap_samedir, resIresJ_sig_overlap_samedir, by = "Names")
write.csv(resIresJ_resKresL$Names, "resIresJ_resKresL_150820.csv" )
resIresJ_resKresL_UP = resIresJ_resKresL[resIresJ_resKresL$log2FoldChange.x.x > 0,]
dim(resIresJ_resKresL_UP)
resIresJ_resKresL_DOWN = resIresJ_resKresL[resIresJ_resKresL$log2FoldChange.x.x < 0,]
dim(resIresJ_resKresL_DOWN)
write.csv(resIresJ_resKresL_UP$Names, "resIresJ_resKresL_UP_150820.csv")
write.csv(resIresJ_resKresL_DOWN$Names, "resIresJ_resKresL_DOWN_150820.csv")

#Code - I2P9
#Genes overlapping in resIresJ, resKresL, resFresG
#Fig 1E

resIresJ_resFresG = merge(resIresJ_sig_overlap_samedir, resFresG_sig_overlap_samedir, by = "Names")
resIresJ_resFresG_notKnotL = resIresJ_resFresG[resIresJ_resFresG$Names %notin% resKresL_sig_overlap_samedir$Names,]
dim(resIresJ_resFresG_notKnotL)
venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = "VennDiagram_FG_KL_IJ_180820.png", output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, fill = c(col_mine, col_coast, col_zinc), width = 2000, height = 2000) 
Fig2F = venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = NULL, output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, fill = c(col_mine, col_coast, col_zinc), width = 2000, height = 2000) 
pdf(file = "Fig2F_1908.pdf")
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

#Fig 2A
venn.diagram(x= list(resI_sig$Names, resJ_sig$Names,resK_sig$Names, resL_sig$Names), category.names = c("S1", "S2", "T1", "T2"), filename = "VennDiagram_5_140820.png", output = TRUE, cat.pos = c(0,0,0,0), col = c(col_coast, col_coast, col_mine, col_mine)) 
Fig3A = venn.diagram(x= list(resI_sig$Names, resJ_sig$Names,resK_sig$Names, resL_sig$Names), category.names = c("S1", "S2", "T1", "T2"), filename = NULL, output = TRUE, cat.pos = c(0,0,0,0), col = c(col_coast, col_coast, col_mine, col_mine)) 
Fig3A
pdf(file = "Fig3A_1908.pdf")
grid.draw(Fig3A)
dev.off()

##################################################################################################################################################

# 3) Fold Changes

##########################################################################################

#3A) LFC within treatment 

lfcA_resBresC_bothdir = lfcA[lfcA$Names %in% resBresC_sig_overlap$Names,]
lfcB_resBresC_bothdir = lfcB[lfcB$Names %in% resBresC_sig_overlap$Names,]
lfcC_resBresC_bothdir = lfcC[lfcC$Names %in% resBresC_sig_overlap$Names,]
lfcD_resBresC_bothdir = lfcD[lfcD$Names %in% resBresC_sig_overlap$Names,]
lfcE_resBresC_bothdir = lfcE[lfcE$Names %in% resBresC_sig_overlap$Names,]
lfcF_resBresC_bothdir = lfcF[lfcF$Names %in% resBresC_sig_overlap$Names,]
lfcG_resBresC_bothdir = lfcG[lfcG$Names %in% resBresC_sig_overlap$Names,]
lfcH_resBresC_bothdir = lfcH[lfcH$Names %in% resBresC_sig_overlap$Names,]

lfcA_resBresC_samedir = lfcA[lfcA$Names %in% resBresC_sig_overlap_samedir$Names,]
lfcB_resBresC_samedir = lfcB[lfcB$Names %in% resBresC_sig_overlap_samedir$Names,]
lfcC_resBresC_samedir = lfcC[lfcC$Names %in% resBresC_sig_overlap_samediar$Names,]
lfcD_resBresC_samedir = lfcD[lfcD$Names %in% resBresC_sig_overlap_samedir$Names,]
lfcE_resBresC_samedir = lfcE[lfcE$Names %in% resBresC_sig_overlap_samedir$Names,]
lfcF_resBresC_samedir = lfcF[lfcF$Names %in% resBresC_sig_overlap_samedir$Names,]
lfcG_resBresC_samedir = lfcG[lfcG$Names %in% resBresC_sig_overlap_samedir$Names,]
lfcH_resBresC_samedir = lfcH[lfcH$Names %in% resBresC_sig_overlap_samedir$Names,]

lfcB_resB = lfcB[lfcB$Names %in% resB_sig$Names,]
lfcB_resC = lfcB[lfcB$Names %in% resC_sig$Names,]
lfcC_resB = lfcC[lfcC$Names %in% resB_sig$Names,]
lfcC_resC = lfcC[lfcC$Names %in% resC_sig$Names,]
lfcF_resB = lfcF[lfcF$Names %in% resB_sig$Names,]
lfcF_resC = lfcF[lfcF$Names %in% resC_sig$Names,]
lfcG_resB = lfcG[lfcG$Names %in% resB_sig$Names,]
lfcG_resC = lfcG[lfcG$Names %in% resC_sig$Names,]

lfcA_resA = lfcA[lfcB$Names %in% resA_sig$Names,]
lfcA_resD = lfcA[lfcB$Names %in% resD_sig$Names,]
lfcD_resA = lfcD[lfcC$Names %in% resA_sig$Names,]
lfcD_resD = lfcD[lfcC$Names %in% resD_sig$Names,]
lfcE_resA = lfcE[lfcE$Names %in% resA_sig$Names,]
lfcE_resD = lfcE[lfcE$Names %in% resD_sig$Names,]
lfcH_resA = lfcH[lfcH$Names %in% resA_sig$Names,]
lfcH_resD = lfcH[lfcH$Names %in% resD_sig$Names,]

#Graphs etc?

#3B) LFC between treatments
lfcI = lfcShrink(dds3, contrast=list("PopSA.CondZ"), type = "ashr")
lfcI$Names = rownames(lfcI)
lfcJ = lfcShrink(dds3, contrast=list("PopBD.CondZ"), type = "ashr")
lfcJ$Names = rownames(lfcJ)
lfcK = lfcShrink(dds3, contrast=list("PopGR.CondZ"), type = "ashr")
lfcK$Names = rownames(lfcK)
lfcL = lfcShrink(dds3, contrast=list("PopPP.CondZ"), type = "ashr")
lfcL$Names = rownames(lfcL)

lfcI_sig = lfcI[lfcI$Names %in% resI_sig$Names,]
lfcJ_sig = lfcJ[lfcJ$Names %in% resJ_sig$Names,]
lfcK_sig = lfcK[lfcK$Names %in% resK_sig$Names,]
lfcL_sig = lfcL[lfcL$Names %in% resL_sig$Names,]

lfcIJKL = rbind(data.frame(value = lfcI_sig$log2FoldChange, variable = "lfcI"), 
                data.frame(value = lfcJ_sig$log2FoldChange, variable = "lfcJ"),
                data.frame(value = lfcK_sig$log2FoldChange, variable = "lfcK"),
                data.frame(value = lfcL_sig$log2FoldChange, variable = "lfcL"))

ggplot()+
  geom_boxplot(data = lfcIJKL, aes(x = variable, y = value))+xlab("Comparison")+ylab("log2FC")

mean(lfcI_sig$log2FoldChange)
mean(lfcJ_sig$log2FoldChange)
mean(lfcK_sig$log2FoldChange)
mean(lfcL_sig$log2FoldChange)

lfcKvI = merge(as.data.frame(lfcK), as.data.frame(lfcI), by = "Names")
lfcIvJ = merge(as.data.frame(lfcI), as.data.frame(lfcJ), by = "Names")
lfcLvJ = merge(as.data.frame(lfcL), as.data.frame(lfcJ), by = "Names")
lfcKvL = merge(as.data.frame(lfcK), as.data.frame(lfcL), by = "Names")

#For each get lists for K, L and KvL genes. 
lfcIvJ_resK = lfcIvJ[lfcIvJ$Names %in% resK_sig$Names,]
lfcKvI_resK = lfcKvI[lfcKvI$Names %in% resK_sig$Names,]
lfcLvJ_resK = lfcLvJ[lfcLvJ$Names %in% resK_sig$Names,]
lfcKvL_resK = lfcKvL[lfcKvL$Names %in% resK_sig$Names,]

lfcIvJ_resL = lfcIvJ[lfcIvJ$Names %in% resL_sig$Names,]
lfcKvI_resL = lfcKvI[lfcKvI$Names %in% resL_sig$Names,]
lfcLvJ_resL = lfcLvJ[lfcLvJ$Names %in% resL_sig$Names,]
lfcKvL_resL = lfcKvL[lfcKvL$Names %in% resL_sig$Names,]

lfcIvJ_resKresL = lfcIvJ[lfcIvJ$Names %in% resKresL_sig_overlap_samedir$Names,]
lfcKvI_resKresL = lfcKvI[lfcKvI$Names %in% resKresL_sig_overlap_samedir$Names,]
lfcLvJ_resKresL = lfcLvJ[lfcLvJ$Names %in% resKresL_sig_overlap_samedir$Names,]
lfcKvL_resKresL = lfcKvL[lfcKvL$Names %in% resKresL_sig_overlap_samedir$Names,]

lfcKvI$diff = lfcKvI$log2FoldChange.x - lfcKvI$log2FoldChange.y
lfcLvJ$diff = lfcLvJ$log2FoldChange.x - lfcLvJ$log2FoldChange.y
lfcKvI_IJKL = lfcKvI[lfcKvI$Names %in% resIresJ_resKresL$Names,]
lfcLvJ_IJKL = lfcLvJ[lfcLvJ$Names %in% resIresJ_resKresL$Names,]
mean(lfcKvI_IJKL$diff)
mean(lfcLvJ_IJKL$diff)

FCdiffs_IJKL = rbind(data.frame(value = lfcKvI_IJKL$diff, variable = "T1 - S1"), data.frame(value = lfcLvJ_IJKL$diff, variable = "T2 - S2"))

ggplot(data = FCdiffs_IJKL, aes(x = variable, y = value, fill = variable))+
 geom_boxplot(outlier.shape = NA)+xlab("")+
  scale_fill_manual(values = c(GRSA_plot, PPBD_plot))+
  ylim(-7, 4) + ylab("Difference in log2(Fold Change)")+
  theme(panel.border = element_rect(fill = NA, color = "lightgray", size = 2), panel.background = element_blank(), legend.position = "none")+
geom_hline(yintercept = 0, linetype = "dashed", color = "black")

Fig3E = ggplot(data = FCdiffs_IJKL, aes(x = variable, y = value, fill = variable))+
  geom_boxplot(outlier.shape = NA)+xlab("")+
  scale_fill_manual(values = c(GRSA_plot, PPBD_plot))+
  ylim(-7, 4) + ylab("Difference in log2(Fold Change)")+
  theme(panel.border = element_rect(fill = NA, color = "lightgray", size = 2), panel.background = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

pdf(file = "Fig3E_1908.pdf")
grid.draw(Fig3E)
dev.off()
geom_boxplot(data = resids_CZ_resKresL_eco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  


lfcIvJ_resKresL_bothdir = lfcIvJ[lfcIvJ$Names %in% resKresL_sig_overlap$Names,]
lfcKvI_resKresL_bothdir = lfcKvI[lfcKvI$Names %in% resKresL_sig_overlap$Names,]
lfcLvJ_resKresL_bothdir = lfcLvJ[lfcLvJ$Names %in% resKresL_sig_overlap$Names,]
lfcKvL_resKresL_bothdir = lfcKvL[lfcKvL$Names %in% resKresL_sig_overlap$Names,]

#3C) LFCs of shared gene sets

lfc_BvC_all = merge(as.data.frame(lfcB), as.data.frame(lfcC), by = "Names")
lfc_BvC_resB = merge(as.data.frame(lfcB_resB), as.data.frame(lfcC_resB), by = "Names")
lfc_BvC_resC = merge(as.data.frame(lfcB_resC), as.data.frame(lfcC_resC), by = "Names")
lfc_BvC_resBresC = merge(as.data.frame(lfcB_resBresC_bothdir), as.data.frame(lfcC_resBresC_bothdir), by = "Names")

lfc_AvD_all = merge(as.data.frame(lfcA), as.data.frame(lfcD), by = "Names")
lfc_AvD_resA = merge(as.data.frame(lfcA_resA), as.data.frame(lfcD_resA), by = "Names")
lfc_AvD_resD = merge(as.data.frame(lfcA_resD), as.data.frame(lfcD_resD), by = "Names")
lfc_AvD_resAresD = merge(as.data.frame(lfcA_resAresD_bothdir), as.data.frame(lfcD_resAresD_bothdir), by = "Names")

lfc_FvG_all = merge(as.data.frame(lfcF), as.data.frame(lfcG), by = "Names")
lfc_FvG_resB = merge(as.data.frame(lfcF_resB), as.data.frame(lfcG_resB), by = "Names")
lfc_FvG_resC = merge(as.data.frame(lfcF_resC), as.data.frame(lfcG_resC), by = "Names")
lfc_FvG_resBresC = merge(as.data.frame(lfcF_resBresC_bothdir), as.data.frame(lfcG_resBresC_bothdir), by = "Names")

lfc_EvH_all = merge(as.data.frame(lfcE), as.data.frame(lfcH), by = "Names")
lfc_EvH_resA = merge(as.data.frame(lfcE_resA), as.data.frame(lfcH_resA), by = "Names")
lfc_EvH_resD = merge(as.data.frame(lfcE_resD), as.data.frame(lfcH_resD), by = "Names")
lfc_EvH_resAresD = merge(as.data.frame(lfcE_resAresD_bothdir), as.data.frame(lfcH_resAresD_bothdir), by = "Names")

#####

lfc_BvC_resB = lfc_BvC_all[lfc_BvC_all$Names %in% resB_sig$Names & lfc_BvC_all$Names %notin% resBresC_sig_overlap$Names,]
dim(lfc_BvC_resB)
lfc_BvC_resB$Sign = lfc_BvC_resB$log2FoldChange.x*lfc_BvC_resB$log2FoldChange.y
lfc_BvC_resB_PC = 100*length(lfc_BvC_resB[lfc_BvC_resB$Sign > 0,]$Sign)/length(lfc_BvC_resB$Sign)
lfc_BvC_resC = lfc_BvC_all[lfcB$Names %in% resC_sig$Names & lfcB$Names %notin% resBresC_sig_overlap$Names,]
lfc_BvC_resC$Sign = lfc_BvC_resC$log2FoldChange.x*lfc_BvC_resC$log2FoldChange.y
lfc_BvC_resC_PC = 100*length(lfc_BvC_resC[lfc_BvC_resC$Sign > 0,]$Sign)/length(lfc_BvC_resC$Sign)
lfc_BvC_resBresC = lfc_BvC_all[lfc_BvC_all$Names %in% resBresC_sig_overlap$Names,]
lfc_BvC_resBresC$Sign = lfc_BvC_resBresC$log2FoldChange.x*lfc_BvC_resBresC$log2FoldChange.y
lfc_BvC_resBresC_PC = 100*length(lfc_BvC_resBresC[lfc_BvC_resBresC$Sign > 0,]$Sign)/length(lfc_BvC_resBresC$Sign)
lfc_BvC_resB_PC
lfc_BvC_resC_PC
lfc_BvC_resBresC_PC

lfc_BvC_resKresL = lfc_BvC_all[lfc_BvC_all$Names %in% resKresL_sig_overlap$Names,]

lfc_AvD_resA = lfc_AvD_all[lfc_AvD_all$Names %in% resA_sig$Names & lfc_AvD_all$Names %notin% resAresD_sig_overlap$Names,]
dim(lfc_AvD_resA)
lfc_AvD_resA$Sign = lfc_AvD_resA$log2FoldChange.x*lfc_AvD_resA$log2FoldChange.y
lfc_AvD_resA_PC = 100*length(lfc_AvD_resA[lfc_AvD_resA$Sign > 0,]$Sign)/length(lfc_AvD_resA$Sign)
lfc_AvD_resD = lfc_AvD_all[lfc_AvD_all$Names %in% resD_sig$Names & lfc_AvD_all$Names %notin% resAresD_sig_overlap$Names,]
lfc_AvD_resD$Sign = lfc_AvD_resD$log2FoldChange.x*lfc_AvD_resD$log2FoldChange.y
lfc_AvD_resD_PC = 100*length(lfc_AvD_resD[lfc_AvD_resD$Sign > 0,]$Sign)/length(lfc_AvD_resD$Sign)
lfc_AvD_resAresD = lfc_AvD_all[lfc_AvD_all$Names %in% resAresD_sig_overlap$Names,]
lfc_AvD_resAresD$Sign = lfc_AvD_resAresD$log2FoldChange.x*lfc_AvD_resAresD$log2FoldChange.y
lfc_AvD_resAresD_PC = 100*length(lfc_AvD_resAresD[lfc_AvD_resAresD$Sign > 0,]$Sign)/length(lfc_AvD_resAresD$Sign)
lfc_AvD_resA_PC
lfc_AvD_resD_PC
lfc_AvD_resAresD_PC

lfc_AvD_resKresL = lfc_AvD_all[lfc_AvD_all$Names %in% resKresL_sig_overlap$Names,]

#IJKL LFCs

lfcIvJ_resI = lfcIvJ[lfcIvJ$Names %in% resI_sig$Names & lfcIvJ$Names %notin% resIresJ_sig_overlap$Names,]
dim(lfcIvJ_resI)
lfcIvJ_resI$Sign = lfcIvJ_resI$log2FoldChange.x*lfcIvJ_resI$log2FoldChange.y
lfcIvJ_resI_PC = 100*length(lfcIvJ_resI[lfcIvJ_resI$Sign > 0,]$Sign)/length(lfcIvJ_resI$Sign)
lfcIvJ_resJ = lfcIvJ[lfcIvJ$Names %in% resJ_sig$Names & lfcIvJ$Names %notin% resIresJ_sig_overlap$Names,]
lfcIvJ_resJ$Sign = lfcIvJ_resJ$log2FoldChange.x*lfcIvJ_resJ$log2FoldChange.y
lfcIvJ_resJ_PC = 100*length(lfcIvJ_resJ[lfcIvJ_resJ$Sign > 0,]$Sign)/length(lfcIvJ_resJ$Sign)
lfcIvJ_resIresJ_bothdir = lfcIvJ[lfcIvJ$Names %in% resIresJ_sig_overlap$Names,]
lfcIvJ_resIresJ_bothdir$Sign = lfcIvJ_resIresJ_bothdir$log2FoldChange.x*lfcIvJ_resIresJ_bothdir$log2FoldChange.y
lfcIvJ_resIresJ_PC = 100*length(lfcIvJ_resIresJ_bothdir[lfcIvJ_resIresJ_bothdir$Sign > 0,]$Sign)/length(lfcIvJ_resIresJ_bothdir$Sign)
lfcIvJ_resI_PC #NaN
lfcIvJ_resJ_PC #100
lfcIvJ_resIresJ_PC #97.23

lfcKvL_resK = lfcKvL[lfcKvL$Names %in% resK_sig$Names & lfcKvL$Names %notin% resKresL_sig_overlap$Names,]
dim(lfcKvL_resK)
lfcKvL_resK$Sign = lfcKvL_resK$log2FoldChange.x*lfcKvL_resK$log2FoldChange.y
lfcKvL_resK_PC = 100*length(lfcKvL_resK[lfcKvL_resK$Sign > 0,]$Sign)/length(lfcKvL_resK$Sign)
lfcKvL_resL = lfcKvL[lfcKvL$Names %in% resL_sig$Names & lfcKvL$Names %notin% resKresL_sig_overlap$Names,]
lfcKvL_resL$Sign = lfcKvL_resL$log2FoldChange.x*lfcKvL_resL$log2FoldChange.y
lfcKvL_resL_PC = 100*length(lfcKvL_resL[lfcKvL_resL$Sign > 0,]$Sign)/length(lfcKvL_resL$Sign)
lfcKvL_resKresL_bothdir = lfcKvL[lfcKvL$Names %in% resKresL_sig_overlap$Names,]
lfcKvL_resKresL_bothdir$Sign = lfcKvL_resKresL_bothdir$log2FoldChange.x*lfcKvL_resKresL_bothdir$log2FoldChange.y
lfcKvL_resKresL_PC = 100*length(lfcKvL_resKresL_bothdir[lfcKvL_resKresL_bothdir$Sign > 0,]$Sign)/length(lfcKvL_resKresL_bothdir$Sign)
lfcKvL_resK_PC
lfcKvL_resL_PC
lfcKvL_resKresL_PC

lfcKvL_resBresC = lfcKvL[lfcKvL$Names %in% resBresC_sig_overlap$Names,]
lfcKvL_resBresC$Sign = lfcKvL_resBresC$log2FoldChange.x*lfcKvL_resBresC$log2FoldChange.y
lfcKvL_resBresC_PC = 100*length(lfcKvL_resBresC[lfcKvL_resBresC$Sign > 0,]$Sign)/length(lfcKvL_resBresC$Sign)
lfcKvL_resBresC_PC

PC_shared = melt(c(lfc_BvC_resBresC_PC, lfc_AvD_resAresD_PC, lfcIvJ_resIresJ_PC, lfcKvL_resKresL_PC))
PC_shared$Names = c("lfc_BvC_resBresC_PC", "lfc_AvD_resAresD_PC", "lfcIvJ_resIresJ_PC", "lfcKvL_resKresL_PC")

ggplot(data=PC_shared, aes(x=Names, y=value, fill=Names))+geom_bar(stat="identity")+
  scale_x_discrete(limits = PC_shared$Names)+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))+labs(y = "% shared DE genes changing in same direction")+
  ylim(0,100)+scale_fill_manual(values=c("lightgreen", "yellow1", col_coast, col_mine))+guides(fill=FALSE)

#Plotting FCs

#lfc plots for resBresC in resB and resC, and resAresD in resA and resD

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
#So obviously there's not a relationship here. 

ggplot()+
  xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.5, data = lfc_AvD_resAresD)+
  ylab("SAC vs BDC (A) log2(Fold Change)")+
  xlab("GRC vs PPC (D) log2(Fold Change)")+
  theme(panel.background = element_rect(fill = purple_plot))


#resKresL overlaps

ggplot()+
  #xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcIvJ, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.3, data = lfcIvJ_resIresJ_bothdir)+
  scale_fill_viridis_c()+

    ylab("S2 - Control to Zinc log2(Fold Change)")+
  xlab("S1 - Control to Zinc log2(Fold Change)")+
#  theme(panel.background = element_rect(fill = coast_plot))
  theme(panel.border = element_rect(colour = col_coast, fill = NA, size = 2), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))

Fig3B = ggplot()+
  #xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcIvJ, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.3, data = lfcIvJ_resIresJ_bothdir)+
  scale_fill_viridis_c()+
  
  ylab("S2 - Control to Zinc log2(Fold Change)")+
  xlab("S1 - Control to Zinc log2(Fold Change)")+
  #  theme(panel.background = element_rect(fill = coast_plot))
  theme(panel.border = element_rect(colour = col_coast, fill = NA, size = 2), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))


#7.5 v 7 = 581 v 622

581/622
7.0*581/622

pdf(file = "Fig3B_1908.pdf", width = 7.5, height = 6.538)
grid.draw(Fig3B)
dev.off()

  

woof = summary(lm(lfcIvJ_resIresJ_bothdir$log2FoldChange.y ~ lfcIvJ_resIresJ_bothdir$log2FoldChange.x))
woof
woof$r.squared

ggplot()+
  xlim(-10.2,10.2)+ ylim(-10.2,10.2)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcKvL, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.3, data = lfcKvL_resKresL_bothdir)+
    scale_fill_viridis_c()+
  ylab("T2 - Control to Zinc log2(Fold Change)")+
  xlab("T1 - Control to Zinc log2(Fold Change)")+
  #theme(panel.background = element_rect(fill = mine_plot))
  theme(panel.border = element_rect(colour = col_mine, fill = NA, size = 2), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))

Fig3C = ggplot()+
  xlim(-10.2,10.2)+ ylim(-10.2,10.2)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcKvL, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.3, data = lfcKvL_resKresL_bothdir)+
  scale_fill_viridis_c()+
  ylab("T2 - Control to Zinc log2(Fold Change)")+
  xlab("T1 - Control to Zinc log2(Fold Change)")+
  #theme(panel.background = element_rect(fill = mine_plot))
  theme(panel.border = element_rect(colour = col_mine, fill = NA, size = 2), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))

pdf(file = "Fig3C_1908.pdf", width = 7.5, height = 6.538)
grid.draw(Fig3C)
dev.off()



summary(lm(lfcKvL_resKresL_bothdir$log2FoldChange.y ~ lfcKvL_resKresL_bothdir$log2FoldChange.x))


#############################################
#resKresL/resBresC overlaps


#Fig 4A
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

#Fig 4B
Supplementary_Figure_2B = ggplot()+
  xlim(-7,7)+ ylim(-7,7)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_fill_viridis_c()+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth =0.3, data = lfc_BvC_resKresL)+
  ylab("T2 to S2 - log2(Fold Change)")+
  xlab("T1 to S1 - log2(Fold Change)")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
Supplementary_Figure_2B
pdf(file = "Supplementary_Figure_2B.pdf")
grid.draw(Supplementary_Figure_2B)
dev.off()
#Fig 4C
Supplementary_Figure_2C = ggplot()+
  xlim(-6,6)+ ylim(-6,6)+  
  # geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcKvL, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_fill_viridis_c()+
    # geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcKvL_resB, col = col_GRSA)+
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcKvL_resC, col = col_PPBD)+
  geom_bin2d(aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcKvL_resBresC, binwidth=0.3)+
  ylab("T2 - Control to Zinc log2(Fold Change)")+
  xlab("T1 - Control to Zinc log2(Fold Change)")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
Supplementary_Figure_2C
pdf(file = "Supplementary_Figure_2C.pdf")
grid.draw(Supplementary_Figure_2C)
dev.off()  
####################################################################################

#resKresL overlaps, Venn diagrams, etc. 

#a)resKresL overlap Venn diagram...

venn.diagram(x= list(resB_sig$Names, resC_sig$Names), category.names = c("", ""), filename = "VennDiagram_BvC_140820.png", output = TRUE, cat.pos = c(0,0), fill = c(col_control, col_control), col = c(col_GRSA, col_PPBD), rotation.degree = 180, height = 1000, width = 2000, imagetype = "png") 


#b) resKresL overlapping genes

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
  geom_boxplot(data = resids_CZ_resKresL_eco, aes(x = variable, y=value), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
#  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete(labels = c("S1C v S2 (C)", "S1 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+ylim(0,1)+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2), text = element_text(size = 15))
Supplementary_Figure_5A
pdf(file = "Supplementary_Figure_5A.pdf")
grid.draw(Supplementary_Figure_5A)
dev.off()


vst_resKresL = vst_all[rownames(vst_all) %in% resKresL_sig_overlap_samedir$Names,]
plotPCA12(vst_resKresL, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-100, 100),my_ylims = c(-40, 40), my_background_col = col_mine, thicko = 2)

Fig3D_orange = plotPCA12(vst_resKresL, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-80, 80),my_ylims = c(-40, 40), my_background_col = col_mine)
Fig3D_gray = plotPCA12(vst_resKresL, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-80, 80),my_ylims = c(-40, 40), my_background_col = "lightgray")
Fig3D_gray
Fig3D_orange

pdf(file = "Fig3D_1908_1700_orange.pdf")
grid.draw(Fig3D_orange)
dev.off()
pdf(file = "Fig3D_19_08_1700_gray.pdf")
grid.draw(Fig3D_gray)
dev.off()

#################################
#notshared resKresL
################################

resK_notL = resK_sig[resK_sig$Names %notin% resKresL_sig_overlap$Names,]
dim(resK_notL)
resL_notK = resL_sig[resL_sig$Names %notin% resKresL_sig_overlap$Names,]
dim(resL_notK)
notshared_resKresL = c(resK_notL$Names, resL_notK$Names)
length(notshared_resKresL)
notshared_resKresL = unique(notshared_resKresL)
length(notshared_resKresL) #Cool. So then what are these guys doing? 
#write.csv(notshared_resKresL, "notshared_resKresL_2007")

resids_countsA_notshared_resKresL = abs(lfcA[lfcA$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsB_notshared_resKresL = abs(lfcB[lfcB$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsC_notshared_resKresL = abs(lfcC[lfcC$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsD_notshared_resKresL = abs(lfcD[lfcD$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsE_notshared_resKresL = abs(lfcE[lfcE$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsF_notshared_resKresL = abs(lfcF[lfcF$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsG_notshared_resKresL = abs(lfcG[lfcG$Names %in% notshared_resKresL,]$log2FoldChange)
resids_countsH_notshared_resKresL = abs(lfcH[lfcH$Names %in% notshared_resKresL,]$log2FoldChange)

median(resids_countsA_notshared_resKresL)
median(resids_countsE_notshared_resKresL)
median(resids_countsB_notshared_resKresL)
median(resids_countsF_notshared_resKresL)
median(resids_countsC_notshared_resKresL)
median(resids_countsG_notshared_resKresL)
median(resids_countsD_notshared_resKresL)
median(resids_countsH_notshared_resKresL)

#So there are actually a bunch of things that are NA (so not plotted/counted) in some, but not other, comparisons...

both_CZ_notshared_resKresL = melt(data.frame(resids_countsA_notshared_resKresL, resids_countsE_notshared_resKresL, resids_countsB_notshared_resKresL, resids_countsF_notshared_resKresL, resids_countsC_notshared_resKresL, resids_countsG_notshared_resKresL, resids_countsD_notshared_resKresL, resids_countsH_notshared_resKresL), na.rm = TRUE)
both_CZ_notshared_resKresL_wilcox = pairwise.wilcox.test(both_CZ_notshared_resKresL$value, both_CZ_notshared_resKresL$variable, p.adjust.method = "BH")
both_CZ_notshared_resKresL_wilcox
write.csv(as.data.frame(both_CZ_notshared_resKresL_wilcox$p.value), "both_CZ_notshared_resKresL_wilcox.csv")
both_CZ_notshared_resKresL_eco = melt(data.frame(resids_countsA_notshared_resKresL, resids_countsE_notshared_resKresL, resids_countsD_notshared_resKresL, resids_countsH_notshared_resKresL), na.rm = TRUE)

p = ggplot()+
#  geom_half_violin(data = both_CZ_notshared_resKresL, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_coast, col_coast, col_mine, col_mine))+
 # geom_half_violin(data = both_CZ_notshared_resKresL, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resKresL_eco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete(labels = c("S1C v S2C", "S1Z v S2Z", "T1C v T2C", "T1Z v T2Z"))+
  xlab("")+ylab("sqrt(% Residual)")+
  theme(panel.background = element_rect(fill = rgb(1,0,1,0.05)))+ylim(0,1)
p



p = ggplot()+
  geom_rect(aes(xmin = 0.42, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resKresL_eco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S2 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("abs (shrunken log2FC)")+ylim(0,1)+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
p

Supplementary_Figure_5B = ggplot()+
#  geom_rect(aes(xmin = 0.42, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
 # geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resKresL_eco, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
#  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+ylim(0,1)+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
Supplementary_Figure_5B
pdf(file = "Supplementary_Figure_5B.pdf")
grid.draw(Supplementary_Figure_5B)
dev.off()

vst_notshared_resKresL = vst_all[rownames(vst_all) %in% notshared_resKresL,]
plotPCA12(vst_notshared_resKresL, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-100, 100),my_ylims = c(-100, 100), my_background_col = purple_col2) 

#13/08/20 20:12

##########################################################################################################################
######resMresN filtering
##########################################################################################################################

#Filtering...
#First, those that are DE in T1 but not S1...
resM_base_names = resK_sig[resK_sig$Names %notin% resI_sig$Names,]$Names
resN_base_names = resL_sig[resL_sig$Names %notin% resJ_sig$Names,]$Names

length(resM_base_names)
length(resN_base_names)

#Second, those that are DE in both, but significantly different...

resM_KIsig = resK_sig[resK_sig$Names %in% resI_sig$Names,]
resN_LJsig = resL_sig[resL_sig$Names %in% resJ_sig$Names,]

resM = results(dds3, contrast = list("PopGR.CondZ", "PopSA.CondZ"))
resN = results(dds3, contrast = list("PopPP.CondZ", "PopBD.CondZ"))

resM_sig = resM[is.na(resM$padj) == FALSE,]
resM_sig = resM_sig[resM_sig$padj < 0.05,]
resN_sig = resN[is.na(resN$padj) == FALSE,]
resN_sig = resN_sig[resN_sig$padj < 0.05,]
resM_sig$Names = rownames(resM_sig)
resN_sig$Names = rownames(resN_sig)

resM_KIsig_sigdiff = resM_KIsig[resM_KIsig$Names %in% resM_sig$Names,]
resN_LJsig_sigdiff = resN_LJsig[resN_LJsig$Names %in% resN_sig$Names,]
dim(resM_KIsig_sigdiff)
dim(resN_LJsig_sigdiff)

#So then of these, we need those in the same direction where the magnitude of K/L is greater than I/J...

resM_lfcs = merge(as.data.frame(lfcK[lfcK$Names %in% resM_KIsig_sigdiff$Names,]), as.data.frame(lfcI[lfcI$Names %in% resM_KIsig_sigdiff$Names, ]), by = "Names")
resN_lfcs = merge(as.data.frame(lfcL[lfcL$Names %in% resN_LJsig_sigdiff$Names,]), as.data.frame(lfcJ[lfcJ$Names %in% resN_LJsig_sigdiff$Names, ]), by = "Names")
resM_lfcs$Abs =  abs(resM_lfcs$log2FoldChange.x) - abs(resM_lfcs$log2FoldChange.y) #log2x here is lfcK, then...
resN_lfcs$Abs =  abs(resN_lfcs$log2FoldChange.x) - abs(resN_lfcs$log2FoldChange.y) #log2x here is lfcI
dim(resM_lfcs)
dim(resN_lfcs)

#Then we need to find ones that are in the opposite direction...
resM_lfcs$dir = resM_lfcs$log2FoldChange.x*resM_lfcs$log2FoldChange.y
resN_lfcs$dir = resN_lfcs$log2FoldChange.x*resN_lfcs$log2FoldChange.y
resM_lfcs_oppo = resM_lfcs[resM_lfcs$dir < 0,]
resN_lfcs_oppo = resN_lfcs[resN_lfcs$dir < 0,]
dim(resM_lfcs_oppo)
dim(resN_lfcs_oppo)

#Of those in the same direction, we need to find those where K/L has a bigger change...
resM_KIsame = resM_lfcs[resM_lfcs$dir > 0,]
resN_LJsame = resN_lfcs[resN_lfcs$dir > 0,]
dim(resM_KIsame)
dim(resN_LJsame)

resM_KIsame$Abs = abs(resM_KIsame$log2FoldChange.x) - abs(resM_KIsame$log2FoldChange.y)
resN_LJsame$Abs = abs(resN_LJsame$log2FoldChange.x) - abs(resN_LJsame$log2FoldChange.y)

resM_magnitude = resM_KIsame[resM_KIsame$Abs > 0,]
resN_magnitude = resN_LJsame[resN_LJsame$Abs > 0,]
dim(resM_magnitude) #98 and 239, then...
dim(resN_magnitude)

resM_final_isig = unique(c(c(resM_magnitude$Names), c(resM_lfcs_oppo$Names)))
resN_final_isig = unique(c(c(resN_magnitude$Names), c(resN_lfcs_oppo$Names)))
length(resM_final_isig)
length(resN_final_isig)


resM_final = unique(c(resM_final_isig, c(resM_base_names)))
resN_final = unique(c(resN_final_isig, c(resN_base_names)))
#Genes with significantly different expression responses to zinc...

length(resM_final)
length(resN_final)

#5B) From these, get K/I and L/J FCs...keeping only those DGE in K/L respectively...
resM_filter = lfcKvI[lfcKvI$Names %in% resM_final,]
resN_filter = lfcLvJ[lfcLvJ$Names %in% resN_final,]
resM_filter_isig = lfcKvI[lfcKvI$Names %in% resM_final_isig,]
resN_filter_isig = lfcLvJ[lfcLvJ$Names %in% resN_final_isig,]

write.csv(resM_filter$Names, "resM_filter_150820.csv")
write.csv(resN_filter$Names, "resN_filter_150820.csv")

resMresN = merge(as.data.frame(resM_filter), as.data.frame(resN_filter), by = "Names")
dim(resMresN)

resMresN_resBresC = intersect(resMresN$Names,resBresC_sig_overlap$Names)
length(resMresN_resBresC)

lfcKvL_resMresN = lfcKvL[lfcKvL$Names %in% resMresN$Names,]
lfcKvI_resMresN = lfcKvI[lfcKvI$Names %in% resMresN$Names,]
lfcLvJ_resMresN = lfcLvJ[lfcLvJ$Names %in% resMresN$Names,]
lfcIvJ_resMresN = lfcIvJ[lfcIvJ$Names %in% resMresN$Names,]


#See https://stackoverflow.com/questions/43963293/how-to-define-color-of-intersection-in-a-venn-diagram

library(polyclip)
par(mfrow = c(1,1))
vp = venn.diagram(x= list(resM_filter$Names, resN_filter$Names), alpha = 1, fill = c("white", "red"), category.names = c("", ""), filename = NULL, rotation.degree = 180, cat.pos = c(0,0)) 
A <- list(list(x = as.vector(vp[[3]][[1]]), y = as.vector(vp[[3]][[2]])))
B <- list(list(x = as.vector(vp[[4]][[1]]), y = as.vector(vp[[4]][[2]])))

AintB = polyclip(A,B)
ix <- sapply(vp, function(x) grepl("text", x$name, fixed = TRUE))
labs <- do.call(rbind.data.frame, lapply(vp[ix], `[`, c("x", "y", "label")))
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
polygon(B[[1]], col = NA, border = col_GRSA)
polygon(A[[1]], col = NA, border = col_PPBD)
#polygon(AintB[[1]], col ="white")
polygon(AintB[[1]], col =purple_plot, border = NA)
text(x = labs$x, y = labs$y, labels = labs$label)



#K vs. L : resMresN
#Fig3D
Fig4D = ggplot()+
  xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  scale_fill_viridis_c()+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.5, data = lfcKvL_resMresN)+
  ylab("T2 - Control to Zinc log2(Fold Change)")+
  xlab("T1- Control to Zinc log2(Fold Change)")+
  theme(panel.background = element_rect(fill = NA, colour = col_mine, size = 2))
Fig4D  

pdf(file = "Fig4D_1908.pdf", width = 7.5, height = 6.538)
grid.draw(Fig4D)
dev.off()

  
summary(lm(lfcKvL_resMresN$log2FoldChange.y ~ lfcKvL_resMresN$log2FoldChange.x))

#K v I: resM significant resI
Fig4A = ggplot()+
  xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_fill_viridis_c()+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.5, data = resM_filter_isig)+
  ylab("S2 - Control to Zinc log2(Fold Change)")+
  xlab("T2 - Control to Zinc log2(Fold Change)")+
  theme(panel.background = element_rect(fill = NA, colour = col_GRSA, size = 2))+
  geom_abline(slope = 1, intercept = 0)
Fig4A

pdf(file = "Fig4A_1908.pdf", width = 7.5, height = 6.538)
grid.draw(Fig4A)
dev.off()


#L vs. J : resN significant resJ

Fig4B = ggplot()+
  xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_fill_viridis_c()+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.5, data = resN_filter_isig)+
  ylab("S2 - Control to Zinc log2(Fold Change)")+
  xlab("T2 - Control to Zinc log2(Fold Change)")+
  theme(panel.background = element_rect(fill = NA, color = col_PPBD, size = 2))+
  geom_abline(slope = 1, intercept = 0)

Fig4B

pdf(file = "Fig4B_1908.pdf", width = 7.5, height = 6.538)
grid.draw(Fig4B)
dev.off()



#resMresN: samedir
resMresN$SignBoth = resMresN$log2FoldChange.x.x*resMresN$log2FoldChange.x.y
resMresN_samedir = resMresN[resMresN$SignBoth > 0,]
dim(resMresN_samedir)

#resMresN: overlap with resBresC
write.csv(resMresN_samedir$Names, "resMresN_samedir_2007")
resMresN_samedir_resBresC = merge(resMresN_samedir, resBresC_sig_overlap_samedir, by = "Names")
dim(resMresN_samedir_resBresC)
resMresN_samedir_resBresC$Names
write.csv(resMresN_samedir_resBresC$Names, "resMresN_samedir_resBresC_150820.csv")

#Includes
#DN17476_c0_g1: GSTUM_ARATH
#DN6068_c0_g1: GSTXA_ARATH

resMresN_samedir_up = resMresN_samedir[resMresN_samedir$log2FoldChange.x.x > 0,] #so this is presumably...wait...is this right?
resMresN_samedir_down = resMresN_samedir[resMresN_samedir$log2FoldChange.x.x < 0,]
dim(resMresN_samedir_up)
dim(resMresN_samedir_down)
write.csv(resMresN_samedir_up$Names, "resMresN_samedir_up_150820.csv")
write.csv(resMresN_samedir_down$Names, "resMresN_samedir_down_150820.csv")

fish_MN_samedir = my_fisher_test(length(resMresN_samedir$Names), length(resM_filter$Names), length(resN_filter$Names), length(rownames(ddTxi2)))
fish_MN = my_fisher_test(length(resMresN$Names), length(resM_filter$Names), length(resN_filter$Names), length(rownames(ddTxi2)))
table_tmp = data.frame("MN", fish_MN$p.value, fish_MN$estimate, fish_MN_samedir$p.value, fish_MN_samedir$estimate, stringsAsFactors = FALSE)
names(table_tmp) = names(Fisher_Table)
Fisher_Table = rbind(Fisher_Table, table_tmp)
Fisher_Table


# Venn diagram

#resMresN : "resids" stuff

resids_countsA_resMresN = abs(lfcA[lfcA$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsB_resMresN = abs(lfcB[lfcB$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsC_resMresN = abs(lfcC[lfcC$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsD_resMresN = abs(lfcD[lfcD$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsE_resMresN = abs(lfcE[lfcE$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsF_resMresN = abs(lfcF[lfcF$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsG_resMresN = abs(lfcG[lfcG$Names %in% resMresN_samedir$Names,]$log2FoldChange)
resids_countsH_resMresN = abs(lfcH[lfcH$Names %in% resMresN_samedir$Names,]$log2FoldChange)

both_resMresN = melt(data.frame(resids_countsA_resMresN, resids_countsE_resMresN, resids_countsB_resMresN, resids_countsF_resMresN, resids_countsC_resMresN, resids_countsG_resMresN, resids_countsD_resMresN, resids_countsH_resMresN), na.rm = TRUE)
both_resMresN_wilcox = pairwise.wilcox.test(both_resMresN$value, both_resMresN$variable, p.adjust.method = "BH")
write.csv(as.data.frame(both_resMresN_wilcox$p.value), "both_resMresN_wilcox.csv")
both_resMresN_eco = melt(data.frame(resids_countsA_resMresN, resids_countsE_resMresN, resids_countsD_resMresN, resids_countsH_resMresN), na.rm = TRUE)

ggplot()+
  #geom_half_violin(data = both_resMresN_eco, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_coast, col_coast, col_mine, col_mine))+
  #geom_half_violin(data = both_resMresN_eco, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_resMresN_eco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete(labels = c("S1C v S2C", "S1Z v S2Z", "T1C v T1Z", "T2C v T2Z", "PPCvBDC (C)", "PPZvBDZ (G)", "GRCvPPC (D)", "GRZvPPZ (H)"))+
  xlab("")+ylab("sqrt(% Residual)")+ylim(0,1)
  theme(panel.background = element_rect(fill = "lemonchiffon"))

  
p = ggplot()+
    geom_rect(aes(xmin = 0.42, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
    geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
    
    #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
    scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
    new_scale_fill()+
    scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
    #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
    new_scale_fill()+
    geom_boxplot(data = both_resMresN_eco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
    scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
    scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S2 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
    xlab("")+ylab("abs (shrunken log2FC)")+ylim(0,1)+
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
  p

Fig5D_nocol = ggplot()+
#    geom_rect(aes(xmin = 0.42, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
 #   geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
    
    #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
    scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
    new_scale_fill()+
    scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
    #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
    new_scale_fill()+
    geom_boxplot(data = both_resMresN_eco, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
    #scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
    scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v T2 (C)", "T1 v S2 (Z)"))+
    xlab("")+ylab("|FC|")+ylim(0,0.75)+
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))

Fig5D_nocol
pdf(file = "Fig5D_nocol.pdf")  
grid.draw(Fig5D_nocol)    
dev.off()

median(resids_countsD_resMresN)
median(resids_countsH_resMresN)
median(resids_countsA_resMresN)
median(resids_countsE_resMresN)
  
2^median(resids_countsD_resMresN)
2^median(resids_countsH_resMresN)
2^median(resids_countsA_resMresN)
2^median(resids_countsE_resMresN)

vst_resMresN = vst_all[rownames(vst_all) %in% resMresN_samedir$Names,]
Fig4C = plotPCA12(vst_resMresN, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-30, 30),my_ylims = c(-25, 25), my_background_col = "lightgray") 
Fig4C
pdf(file = "Fig4C_1908_1700.pdf")
grid.draw(Fig4C)
dev.off()

############################
#resMresN - not shared
###########################

resM_notN = resM_filter[resM_filter$Names %notin% resN_filter$Names,]
dim(resM_notN)
resN_notM = resN_filter[resN_filter$Names %notin% resM_filter$Names,]
dim(resN_notM)
notshared_resMresN = c(resM_notN$Names, resN_notM$Names)
length(notshared_resMresN)
notshared_resMresN = unique(notshared_resMresN)
length(notshared_resMresN) #Cool. So then what are these guys doing? 
write.csv(notshared_resMresN, "notshared_resMresN_2007")

resids_countsA_notshared_resMresN = abs(lfcA[lfcA$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsB_notshared_resMresN = abs(lfcB[lfcB$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsC_notshared_resMresN = abs(lfcC[lfcC$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsD_notshared_resMresN = abs(lfcD[lfcD$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsE_notshared_resMresN = abs(lfcE[lfcE$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsF_notshared_resMresN = abs(lfcF[lfcF$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsG_notshared_resMresN = abs(lfcG[lfcG$Names %in% notshared_resMresN,]$log2FoldChange)
resids_countsH_notshared_resMresN = abs(lfcH[lfcH$Names %in% notshared_resMresN,]$log2FoldChange)

median(resids_countsA_resMresN)
median(resids_countsB_resMresN)
median(resids_countsC_resMresN)
median(resids_countsD_resMresN)
median(resids_countsE_resMresN)
median(resids_countsF_resMresN)
median(resids_countsG_resMresN)
median(resids_countsH_resMresN)

2^median(resids_countsA_resMresN)
2^median(resids_countsB_resMresN)
2^median(resids_countsC_resMresN)
2^median(resids_countsD_resMresN)
2^median(resids_countsE_resMresN)
2^median(resids_countsF_resMresN)
2^median(resids_countsG_resMresN)
2^median(resids_countsH_resMresN)

both_CZ_notshared_resMresN = melt(data.frame(resids_countsA_notshared_resMresN, resids_countsE_notshared_resMresN, resids_countsB_notshared_resMresN, resids_countsF_notshared_resMresN, resids_countsC_notshared_resMresN, resids_countsG_notshared_resMresN, resids_countsD_notshared_resMresN, resids_countsH_notshared_resMresN), na.rm = TRUE)
both_CZ_notshared_resMresN_wilcox = pairwise.wilcox.test(both_CZ_notshared_resMresN$value, both_CZ_notshared_resMresN$variable, p.adjust.method = "BH")
both_CZ_notshared_resMresN_wilcox
write.csv(as.data.frame(both_CZ_notshared_resMresN_wilcox$p.value), "both_CZ_notshared_resMresN_wilcox.csv")
both_CZ_notshared_resMresN_eco = melt(data.frame(resids_countsA_notshared_resMresN, resids_countsE_notshared_resMresN, resids_countsD_notshared_resMresN, resids_countsH_notshared_resMresN), na.rm = TRUE)

Supplementary_Figure_4 = ggplot()+
#  geom_half_violin(data = both_CZ_notshared_resMresN_eco, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_coast, col_coast, col_mine, col_mine))+
 # geom_half_violin(data = both_CZ_notshared_resMresN_eco, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
#  geom_boxplot(data = both_CZ_notshared_resMresN_eco, aes(x = variable, y=value, fill=variable), width = 0.1, outlier.shape = NA, show.legend = FALSE)+
  geom_boxplot(data = both_CZ_notshared_resMresN, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
#  scale_fill_manual(values = c(zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 V S1 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2))+ylim(0,6)
Supplementary_Figure_4
pdf(file = "Supplementary_Figure_4.pdf")
grid.draw(Supplementary_Figure_4)
dev.off()

p = ggplot()+
  geom_rect(aes(xmin = 0.42, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_coast, col_coast, col_coast, col_coast, col_mine, col_mine))+
  #  geom_half_violin(data = both_CZ_resBresC, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(data = both_CZ_notshared_resMresN_eco, aes(x = variable, y=value, fill=variable), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))+ylim(0,10)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "S1 v S2 (Z)", "T1 v S1 (C)", "T1 v S2 (Z)", "T2 v S2 (C)", "T2 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("abs (shrunken log2FC)")+ylim(0,1)+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2))
p


vst_notshared_resMresN = vst_all[rownames(vst_all) %in% notshared_resMresN,]
plotPCA12(vst_notshared_resMresN, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-50, 50),my_ylims = c(-30, 30), my_background_col = purple_plot) 

#Export the fisher table
write.csv(Fisher_Table, "Fisher_Table.csv")



###################################

#Potentially for supplement
#Si) resD not H ("convergent")

#Excessive DGE in Coastal Zinc problems...
#A) What % of genes are DE - 
length(resI_sig$Names)/length(ddTxi2)
length(resJ_sig$Names)/length(ddTxi2)
#Approximately 48%

#Is the GC contnet of the DE genes different? 
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

save.image()

