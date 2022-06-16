#0 LOADNG IN THE DATA, FILTERING
remove(list=ls())
`%notin%` = Negate(`%in%`)
setwd(dir = "~/RNA_Seq")

###########################################
#1 FUNCTION
###########################################-----
#install.packages(c("wesanderson", "RColorBrewer", "ggsci"))

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
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

plotPCA12_revcol <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
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
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) + scale_fill_manual(values = c(col_mine, col_coast))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

#Plot first 2 PCs (centroids)

plotPCA12_centroids <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 4) 
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

plotPCA34 <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
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
  ggplot(data = d, aes_string(x = "PC3", y = "PC4")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 8, stroke=2, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC3: ", round(percentVar[3] * 
                                                                                                                                                            100), "% variance")) + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC4: ", round(percentVar[4]*100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

plotPCA34_revcol <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
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
  ggplot(data = d, aes_string(x = "PC3", y = "PC4")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 8, stroke=2, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC3: ", round(percentVar[3] * 
                                                                                                                                                            100), "% variance")) + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC4: ", round(percentVar[4]*100), "% variance")) + scale_fill_manual(values = c(col_mine, col_coast))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

plotPCA14_revcol <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
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
    attr(d, "percentVar") <- percentVar[1,4]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC4")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 8, stroke=2, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC3: ", round(percentVar[1] * 
                                                                                                                                                            100), "% variance")) + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC4: ", round(percentVar[4]*100), "% variance")) + scale_fill_manual(values = c(col_mine, col_coast))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

plotPCA13 <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
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
    attr(d, "percentVar") <- percentVar[1,3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 8, stroke=2, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                                                                                                                            100), "% variance")) + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC3: ", round(percentVar[3]*100), "% variance")) + scale_fill_manual(values = c(col_coast, col_mine))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

plotPCA13_revcol <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
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
    attr(d, "percentVar") <- percentVar[1,3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3")) + xlim(my_xlims[1], my_xlims[2]) + ylim(my_ylims[1], my_ylims[2])+ 
    #    ggplot(data = d, aes_string(x = "PC1", y = "PC2"))+ 
    geom_point(size = 8, stroke=2, aes(shape = group1, col = group3, fill = group2)) + scale_shape_manual(values = c(21,24)) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                                                                                                                            100), "% variance")) + 
    scale_colour_manual(values = c(col_PPBD, col_GRSA)) + ylab(paste0("PC3: ", round(percentVar[3]*100), "% variance")) + scale_fill_manual(values = c(col_mine, col_coast))  +
    
    coord_fixed() + theme(legend.position = "none")+
    theme(panel.border = element_rect(colour = my_background_col, fill = NA, size = thicko), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15)) #, panel.grid.major=element_line(colour="lightgrey"))
  
}

plotPCA_3D <- function (object, intgroup1 = "Geog", ntop = 2000000, returnData = FALSE, intgroup2 = "Ecotype", intgroup3 = "condition", my_xlims = c(-200, 200), my_ylims = c(-200, 200), my_background_col = col_control, thicko = 8) 
{
  #TEMP
  object=vst_all
  #
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
  
  get_colors <- function(groups, group.col = palette()){
    groups <- as.factor(groups)
    ngrps <- length(levels(groups))
    if(ngrps > length(group.col)) 
      group.col <- rep(group.col, ngrps)
    color2 <- group.col[as.numeric(groups)]
    names(color2) <- as.vector(groups)
    return(color2)
  }
    
  
  intgroup3.df$Ecotype = as.factor(intgroup3.df$Ecotype)
  intgroup3.df
  get_colors(intgroup3.df$Ecotype)
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], 
                  group1 = group1, group2 = group2, group3 = group3, name = colnames(object))
  intgrop3d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group1 = group1, 
                          intgroup1.df, group2 = group2, intgroup2.df, name = colnames(object))
  
  
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1,3]
    return(d)
  }
  plot3d(d[,1:3], size=10, type='p'
      , xlim = c(-50,50), ylim=c(-50,50), zlim=c(-50,50), )
  color_vector = get_colors(intgroup3.df$Ecotype)
  color_vector
  rgl.spheres(d[,1], d[,2], d[,3], color = color_vector)  
  
  #plot3d(d[,1:3], color = color_vector)
  text3d(d[,1]+2, d[,2]+10, d[,3]+2,
         texts=c(rownames(d)), cex= 0.7, pos=3)
#  rgl.postscript("3DPCA-Merge.pdf", "pdf")
  color_vector
  vector_thing = color_vector
  with(d, plot3d(d[,1],d[,2],d[,3], col = color_vector, type ='p'))
  rglwidget()
  
  
  
}

heatmap_samples <-function(vst_all) {
  sampleDists <- dist(t(assay(vst_all)))
  library("RColorBrewer")
  #install.packages("pheatmap")
  library(pheatmap)
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vst_all$Pop_Cond)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
}



#permutation stuff: does plotting a/c vs. b/c give a spurious correlation?
a = rnorm(100, mean = 10, sd = 1)
b = rnorm(100, mean = 10, sd = 1)
c = rnorm(100, mean = 10, sd = 1)
woof1 = a/c
woof2 = b/c
plot(woof1, woof2, xlim = c(-3,3), ylim = c(-3,3))
summary(lm(woof2 ~ woof1)) #So there is a correlation...


###################################
#2 Fisher Table setup
###################################
Fisher_Table = data.frame("Comparison", "All_pvalue (0 = <2.2e-16)", "All_oddsratio", "Samedir_pvalue", "Samedir_oddsratio", stringsAsFactors = FALSE)

##################################
#3 Colours
##################################

library(wesanderson)
library(RColorBrewer)
#install.packages("ggsci")
library(ggsci)
library(RColorBrewer)

color_list = c("#75EAA9","#B6EA50", "#F5FF74")

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

tmp = col2rgb(col_mine)
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

#install.packages(c("DESeq2", "tximport", "ggplot2", "reshape2", "VennDiagram", "gghalves", "ggnewscale", "cowplot"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("DESeq2")
#install.packages("tximport")
#BiocManager::install("tximport")
#BiocManager::install("rhdf5")
#install.packages("DESeq2")

library(DESeq2) #1.26.0
library(tximport) #1.14.2
library(ggplot2)
library(reshape2)
library(VennDiagram)
library(gghalves)
library(ggnewscale)
#library(rcompanion)
library(cowplot)

##########################################
#5 Setup and filtering
##########################################

samples = read.table("~/RNA_Seq/samples_file.txt_forslueth.tsv", sep = "\t", header = TRUE)                     

files = file.path("~/RNA_Seq/Formal_RNASeq/", samples$Sample, "abundance.h5")
samples
#Ok well backup this PC too

#Adjust samples table to add "Pop_Cond" variable
samples$Pop_Cond = paste(samples$Pop, samples$Cond, sep="")
samples$Ecotype = c("C", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M","M", "C", "C", "C", "C", "C", "C")
samples$Geog = c("E", "E", "E", "E", "E", "E", "W", "W", "W", "W", "W", "W", "E", "E", "E", "E", "E", "E", "W", "W", "W", "W", "W", "W")
samples$Ecotype = as.factor(samples$Ecotype)
samples$Geog = as.factor(samples$Geog)
samples


#Read in Trinity gene-transcript map
tx2gene = read.csv("~/RNA_Seq/Trinity-clusters.txt", sep = "\t", header = FALSE)
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

#blast transcripts against SwissProt/UniProt10 database [version 290]
blast_info = read.csv("~/RNA_Seq/blastp.outfmt6.emb.longesti", header = FALSE, sep = "\t")
#Use blat to search against Silene uniflora genome (GCA_018983105.1_ASM1898310v1_genomic.fna)
blat_uniflora_info = read.csv("~/SR29_sorted_merged_outfile_best.psl.longesti", sep = "\t", header = FALSE)
#GC content caluclated for each contigs longest isoform
GC_content = read.csv("~/RNA_Seq/SU_RedSamp_k25mincov2.fasta.trans_filtered.GClongesti", header=FALSE, sep="\t")
emb_yes = blast_info[blast_info$V13 == "emb_yes",]
emb_no = blast_info[blast_info$V13 == "emb_no",]
dim(emb_yes)
dim(emb_no)

dim(blast_info) #73k out of...121k. So 50k don't have a match to anything.
dim(ddTxi_backup) 

blat_uniflora_info$match_length = blat_uniflora_info$V1+blat_uniflora_info$V2
blat_uniflora_info$pc = blat_uniflora_info$V1/(blat_uniflora_info$V1+blat_uniflora_info$V2)
#i) Match length to uniflora genome must be >200 and %ID > 90%
blat_uniflora_good = blat_uniflora_info[blat_uniflora_info$match_length > 200 & blat_uniflora_info$pc > 0.9,]
blat_200plus_list = blat_uniflora_good$V10
length(blat_200plus_list)
#ii) Match length must be >200 and %ID > 70%
emb_200plus70pc = emb_yes[emb_yes$V4 > 200 & emb_yes$V3 > 70,]

#Filtering 1: only genes either matching i) or ii) are included, provisionally. 
GC_content_Filter2 = GC_content[GC_content$V1 %in% blat_200plus_list | GC_content$V1 %in% emb_200plus70pc$V1,]
dim(GC_content_Filter2)
#Filtering 2: any of these genes whose top hit is to a non-embryophyte (regardless of length etc.) excluded
GC_content_Filter2 = GC_content_Filter2[GC_content_Filter2$V1 %notin% emb_no$V1,]
dim(GC_content_Filter2)
hist(GC_content_Filter2$V2, breaks = 100)

ddTxi = subset(ddTxi, rownames(ddTxi) %in% GC_content_Filter2$V1)
dim(ddTxi)

#Filtering 3: removes genes with <10 counts across all samples
keep = rowSums(counts(ddTxi)) >= 10
hist(rowSums(counts(ddTxi)), breaks = 5000, xlim = c(0,10000))
summary(rowSums(counts(ddTxi)))
head(rowSums(counts(ddTxi)))

ddTxi2 = ddTxi[keep,]
GC_filter3 = GC_content_Filter2[GC_content_Filter2$V1 %in% rownames(ddTxi2),]
dim(GC_filter3)
dim(ddTxi2) #So this is the final set of genes to be included.

par(mfrow = c(2,1))
hist(GC_content$V2, breaks = 100, xlim = c(20,70), main = "Unfiltered")
hist(GC_filter3$V2, breaks = 100, xlim = c(20,70), main = "New filtering") #So this still has a load of crap

#############
# 5B Transcriptome quality metrics
#############

#Length/ExN50
trinity_quality = read.table("~/SU_RedSamp_k25mincov2.fasta.stats", header = T)
head(trinity_quality)
plot(trinity_quality$Ex, trinity_quality$ExN50)
trinity_quality2 = read.table("~/SU_RedSamp_k25mincov2.fasta.trans_filtered.sub.stats", header = T)
head(trinity_quality2)

plot(trinity_quality2$Ex, trinity_quality2$ExN50)

#Multimapping transcripts?
sr30 = read.csv("~/SR29B.sorted_merged_outfile_best.psl.filt.re.sorted.not.SR30", header = FALSE, sep = "\t")
sr30 = read.csv("~/sorted_merged_outfile_best.psl.filt.re.sorted.not.SR30", header = FALSE, sep = "\t")
dim(sr30)
hist(sr30$V2, xlab = "Percentage match length on different scaffold", main = "nice genome")
head(sr30)
sr30 = na.omit(sr30)
dim(sr30)

plot(sr30$V2, sr30$V3)
sum(sr30$V2 > 0.95)

write.csv(rownames(ddTxi3), "background_ddTxi3_100622")

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
#6 Counts, PCAs/heatmaps in control/zinc
####################################################

vst_all = vst(dds)

dds_C = dds[,dds$Cond == "C"]
vst_C = vst(dds_C)
dds_Z = dds[,dds$Cond == "Z"]
vst_Z = vst(dds_Z)

dds_S = dds[,dds$Ecotype %in% "C"]
vst_S = vst(dds_S)
dds_T = dds[,dds$Ecotype %in% "M"]
vst_T = vst(dds_T)

#Heatmaps...
heatmap_samples(vst_all)
heatmap_samples(vst_C)
heatmap_samples(vst_Z)
heatmap_samples(vst_S)
heatmap_samples(vst_T)

#PCA all individuals, all conditions (with centoirds and arrows)

#https://stackoverflow.com/questions/36208909/how-to-calculate-centroids-in-pca - based on this
PCA_all = plotPCA12_centroids(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-300, 200), my_ylims = c(-150, 150), my_background_col = "black")
PCA_all

plotPCA_3D(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-300, 200), my_ylims = c(-150, 150), my_background_col = "black")
dev.off()

PCA_34 = plotPCA34(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-150, 150), my_ylims = c(-150, 150), my_background_col = "black")
PCA_34

PCA_13 = plotPCA13(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop= 1000000, my_xlims = c(-350, 350), my_ylims = c(-150, 150), my_background_col = "black")
PCA_13


blork.x = as.data.frame(PCA_all$x)
blork.x$groups = samples$Pop_Cond
pca.centroids = aggregate(blork.x[,1:2], list(Type = blork.x$groups), mean)
pca_centroids_c = pca.centroids[c(1,3,5,7),]

PCA_all = plotPCA12(vst_all, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-200, 100),my_ylims = c(-100, 100), my_background_col = "black", thicko = 4)
PCA_all


#PCA of whole transcriptome, across all populations and treatments
#KR12
PCA_all_plot = PCA_all+
  geom_segment(aes(x = pca.centroids[7,2], y = pca.centroids[7,3], xend = pca.centroids[8,2], yend = pca.centroids[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[7,2], y = pca.centroids[7,3], xend = pca.centroids[8,2], yend = pca.centroids[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[1,2], y = pca.centroids[1,3], xend = pca.centroids[2,2], yend = pca.centroids[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5,col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[1,2], y = pca.centroids[1,3], xend = pca.centroids[2,2], yend = pca.centroids[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot,  lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[5,2], y = pca.centroids[5,3], xend = pca.centroids[6,2], yend = pca.centroids[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[5,2], y = pca.centroids[5,3], xend = pca.centroids[6,2], yend = pca.centroids[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[3,2], y = pca.centroids[3,3], xend = pca.centroids[4,2], yend = pca.centroids[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids[3,2], y = pca.centroids[3,3], xend = pca.centroids[4,2], yend = pca.centroids[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_point(data =  pca_centroids_c, size = 4, shape = 21, fill = c(PPBD_plot, GRSA_plot, PPBD_plot, GRSA_plot), alpha = 1)
PCA_all_plot
#Try a 3D plot?
#install.packages("rgl")
library(rgl)

PCA_S = plotPCA12(vst_S, my_xlims = c(-200,150), my_ylims = c(-100,100), intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_background_col = "black", thicko = 4)
PCA_S

PCA_T = plotPCA12_revcol(vst_T, my_xlims = c(-100,100), my_ylims = c(-75,75), intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_background_col = "black", thicko = 4)
PCA_T
PCA_T = plotPCA13_revcol(vst_T, my_xlims = c(-200,150), my_ylims = c(-100,100), intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_background_col = "black", thicko = 4)
PCA_T
PCA_T = plotPCA14_revcol(vst_T, my_xlims = c(-200,150), my_ylims = c(-100,100), intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000, my_background_col = "black", thicko = 4)
PCA_T

#KR7
PCA_control_plot = plotPCA12(vst_C, intgroup1 = "Cond", my_xlims = c(-75, 75), my_ylims = c(-75,75), intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 2000000,  my_background_col = "black", thicko = 4)
PCA_control_plot

#pdf(file = "Fig2A_200721_ControlPCA.pdf")
grid.draw(PCA_control_plot)
#dev.off()

#################################################
#7 Within condition comparisons
#################################################

#7.1 : Look at each populations separately 
ddTxi3 = ddTxi2
design(ddTxi3) <- ~Pop_Cond
dds2 = DESeq(ddTxi3)
resultsNames(dds2)
colnames(dds2) = paste(dds2$IND, dds2$Cond, sep="")
blork = counts(dds2, normalized=T)
head(blork)
#write.csv(blork, "~/Normalised_Counts,tsv", sep = "\t")

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

#write.csv(resA_sig$Names, "resA_sig_210721")
#write.csv(resB_sig$Names, "resB_sig_210721")
#write.csv(resC_sig$Names, "resC_sig_210721")
#write.csv(resD_sig$Names, "resD_sig_210721")
#write.csv(resE_sig$Names, "resE_sig_210721")
#write.csv(resF_sig$Names, "resF_sig_210721")
#write.csv(resG_sig$Names, "resG_sig_210721")
#write.csv(resH_sig$Names, "resH_sig_210721")

dim(resA_sig)
dim(resB_sig)
dim(resC_sig)
dim(resD_sig)
dim(resE_sig)
dim(resF_sig)
dim(resG_sig)
dim(resH_sig)


############
# 7.2 - Genes that have a significant difference when combining both populations. 
############

#MA6.2 1/N - calculating within ecotype values, so that we can look at combined genes. 
ddTxi_ecowithin = ddTxi2
paste(ddTxi_ecowithin$Ecotype, ddTxi_ecowithin$Cond, sep = "")
ddTxi_ecowithin$Eco_Cond = as.factor(paste(ddTxi_ecowithin$Ecotype, ddTxi_ecowithin$Cond, sep = ""))
#design(ddTxi_ecowithin) <- ~0+Eco_Cond
design(ddTxi_ecowithin) <- ~Eco_Cond
ddTxi_ecowithin$Eco_Cond
dds5 = DESeq(ddTxi_ecowithin)
resultsNames(dds5)

resY = results(dds5, contrast=c("Eco_Cond", "CZ", "MZ"))
resY$Names = rownames(resY)
dim(resY)
resY = resY[is.na(resY$padj) == FALSE,]
resY_sig = resY[resY$padj < 0.05,]
dim(resY_sig) #14,929...ok...
resY_sig$Names = rownames(resY_sig)

###############################################
#####8 between condition comparisons
###############################################

#I/J/K/L
#Includes individual as factor (woo clones)

ddTxi_ind = ddTxi2
ddTxi_ind$Ind_plant = as.factor(c(1, 2, 2, 1, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3))
#See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
design(ddTxi_ind) <- ~Pop + Pop:Ind_plant + Pop:Cond
#Double check you're happy with this?
dds3 = DESeq(ddTxi_ind) #This must have gotten deleted?

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

#write.csv(resI_sig$Names, "resI_sig_210721.csv")
#write.csv(resJ_sig$Names, "resJ_sig_210721.csv")
#write.csv(resK_sig$Names, "resK_sig_210721.csv")
#write.csv(resL_sig$Names, "resL_sig_210721.csv")

######
# 8.2 - Combining across ecotypes (to get resX)
#MA6.2 (2/N) - combining the plants across ecotype, i.e. so genes where there's a general ecotype effect.
ddTxi_acrosspopind = ddTxi2
ddTxi2$IND
ddTxi2$Ecotype
ddTxi_acrosspopind$Ind_plant = as.factor(c(1, 2, 2, 1, 3, 3, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 4, 4, 5, 5, 6, 6))
ddTxi_acrosspopind$Ind_plant
design(ddTxi_acrosspopind) <- ~Ecotype + Ecotype:Ind_plant + Ecotype:Cond
dds4 = DESeq(ddTxi_acrosspopind)

resultsNames(dds4)

resX = results(dds4, contrast = list("EcotypeC.CondZ"))
resX$Names = rownames(resX)
resX_sig = resX[resX$padj < 0.05,]
resX_sig$Names = rownames(resX_sig)
dim(resX_sig) #Ok so 17k genes in this category.
#Ok cool. 

############################################
#9 Overlapping sets (includes venn diagrams)
############################################

#Overlapping gene sets of interest

#Within-condition

resBresC_sig_overlap = merge(as.data.frame(resB_sig), as.data.frame(resC_sig), by = "Names")
#KR8(1/2)
dim(resBresC_sig_overlap)
resFresG_sig_overlap = merge(as.data.frame(resF_sig), as.data.frame(resG_sig), by = "Names")
write.csv(resFresG_sig_overlap$Names, "resFresG_sig_overlap_samedir_100622.csv")

resAresE_sig_overlap = merge(as.data.frame(resA_sig), as.data.frame(resE_sig), by = "Names")
resDresH_sig_overlap = merge(as.data.frame(resD_sig), as.data.frame(resH_sig), by = "Names")
resAresD_sig_overlap = merge(as.data.frame(resA_sig), as.data.frame(resD_sig), by = "Names")

dim(resBresC_sig_overlap)
#KR2: Number of genes DE between S1/T1 and S2/T2 in the zinc. KR_2
dim(resFresG_sig_overlap)
dim(resAresE_sig_overlap)
dim(resDresH_sig_overlap)

#resBresC overlaps in same direction - CEC genes
resBresC_sig_overlap$sign = resBresC_sig_overlap$log2FoldChange.x*resBresC_sig_overlap$log2FoldChange.y
resBresC_sig_overlap_samedir = resBresC_sig_overlap[resBresC_sig_overlap$sign > 0,]
#KR9
dim(resBresC_sig_overlap_samedir)
write.csv(resBresC_sig_overlap_samedir$Names, "resBresC_sig_overlap_samedir_100622.csv")

dim(resB_sig)
dim(resC_sig)
dim(resBresC_sig_overlap)
dim(resBresC_sig_overlap_samedir) 

#KR10 - greater overlap resBresC than expected by chance.
fish_BC_samedir = my_fisher_test(length(resBresC_sig_overlap_samedir$Names), length(resB_sig$Names), length(resC_sig$Names), length(rownames(ddTxi2)))
fish_BC_samedir

#Upregulated/downregulated CEC genes.
#write.csv(resBresC_sig_overlap_samedir$Names, "resBresC_sig_overlap_samedir_210721.csv")
resBresC_sig_overlap_samedir_UP = resBresC_sig_overlap_samedir[resBresC_sig_overlap_samedir$log2FoldChange.x > 0,]
resBresC_sig_overlap_samedir_DOWN = resBresC_sig_overlap_samedir[resBresC_sig_overlap_samedir$log2FoldChange.x < 0,]
dim(resBresC_sig_overlap_samedir_UP)
dim(resBresC_sig_overlap_samedir_DOWN)
#write.csv(resBresC_sig_overlap_samedir_UP$Names, "resBresC_sig_overlap_samedir_up_210721.csv")
#write.csv(resBresC_sig_overlap_samedir_DOWN$Names, "resBresC_sig_overlap_samedir_down_210721.csv")

#Between conditions

#KR1: % of transcriptome DE in both sensitive populations
resIresJ_sig_overlap = merge(as.data.frame(resI_sig), as.data.frame(resJ_sig), by = "Names")
dim(resIresJ_sig_overlap)
dim(resIresJ_sig_overlap)[1]/dim(resI)[1] #% of transcriptome DE in sensitive populations between treatments

resIresJ_samedir_up = resIresJ_sig_overlap_samedir[resIresJ_sig_overlap_samedir$log2FoldChange.x > 0,]
write.csv(resIresJ_samedir_up$Names, "resIresJ_samedir_up_100622.csv")


#T
resKresL_sig_overlap = merge(as.data.frame(resK_sig), as.data.frame(resL_sig), by = "Names")
dim(resK_sig)
dim(resL_sig)
dim(resKresL_sig_overlap)
fish_KL = my_fisher_test(length(resKresL_sig_overlap$Names), length(resK_sig$Names), length(resL_sig$Names), length(rownames(ddTxi_ind)))
fish_KL

resKresL_notInotJ = resKresL_sig_overlap[resKresL_sig_overlap$Names %notin% resIresJ_sig_overlap$Names,]
1-(dim(resKresL_notInotJ)[1]/dim(resKresL_sig_overlap)[1])
#So 80% of the resKresL genes are actually also in resI and resJs

#Multiply FCs; those with positive value are in same direction. 
resIresJ_sig_overlap$Sign = resIresJ_sig_overlap$log2FoldChange.x*resIresJ_sig_overlap$log2FoldChange.y
resIresJ_sig_overlap_samedir = resIresJ_sig_overlap[resIresJ_sig_overlap$Sign > 0,]

resKresL_sig_overlap$Sign = resKresL_sig_overlap$log2FoldChange.x*resKresL_sig_overlap$log2FoldChange.y
resKresL_sig_overlap_samedir = resKresL_sig_overlap[resKresL_sig_overlap$Sign > 0,]
dim(resKresL_sig_overlap_samedir)

fish_KL_samedir = my_fisher_test(length(resKresL_sig_overlap_samedir$Names), length(resK_sig$Names), length(resL_sig$Names), length(rownames(ddTxi_ind)))
fish_KL = my_fisher_test(length(resKresL_sig_overlap$Names), length(resK_sig$Names), length(resL_sig$Names), length(rownames(ddTxi_ind)))

#Defining DP(N)

#ALTERNATE: So here we have resKresL_sig_overlap_samedir. However, we will want to exclude things that are just identical to the S ones...
#S1/T1 - a list of genes that aren't different in either condition.
temp_list = rownames(dds3)
#Pop1
#Note: previous way of defining these I think doesn't make sense. I think they need to be DE in Z, or in C, in both. 
#Otherwise they're not reasonably doing the same thing I think. 

BorF = c(resB_sig$Names, resF_sig$Names)
length(BorF)
#Pop2
CorG = c(resC_sig$Names, resG_sig$Names)
length(CorG)
#So we will want things that are in both of these sets: between ecotypes in at least 1 condition, in both sets.
BorF_CorG = intersect(CorG, BorF)
length(BorF_CorG) #So there are 10k genes fitting this description
dim(resKresL_sig_overlap_samedir[resKresL_sig_overlap_samedir$Names %in% BorF_CorG,])
#So that's now 142. Why are these numbers changing slightly?

#So it needs to be resBresC and resFresG I guess. 
BandC_FandG = unique(c(resBresC_sig_overlap$Names, resFresG_sig_overlap$Names))
length(BandC_FandG)
#So they are different in both. But it doesn't necessarily have to be in same direction...
#worth discussing with Alex explicitly. 

#This is the defintion of the DP(N) genes. 

#Defining DPN genes. 

#So we need the DP(N) genes for each population to compare this to...
BorF = unique(BorF)
resK_DPN = resK_sig[resK_sig$Names %in% BorF,]

#KR14
dim(resK_sig)
dim(resK_DPN)

CorG = unique(CorG)
resL_DPN = resL_sig[resL_sig$Names %in% CorG,]

dim(resL_sig)
dim(resL_DPN)
length(unique(intersect(resK_DPN$Names, resL_DPN$Names)))
resKresL_DPN_sig_overlap_samedir = resKresL_sig_overlap_samedir[resKresL_sig_overlap_samedir$Names %in% BandC_FandG,]
dim(resKresL_DPN_sig_overlap_samedir) #So this is 138. So it's the same, essentially. Fine, good. 
write.csv(resKresL_DPN_sig_overlap_samedir$Names, "resKresL_DPN_sig_overlap_samedir_100622.csv")
fish_DPN = my_fisher_test(length(resKresL_DPN_sig_overlap_samedir$Names), length(resK_DPN$Names), length(resL_DPN$Names), length(rownames(ddTxi_ind)))
fish_DPN

resKresL_DPN_sig_overlap_samedir_ancplast = resKresL_DPN_sig_overlap_samedir[resKresL_DPN_sig_overlap_samedir$Names %in% resIresJ_sig_overlap$Names,]
resKresL_DPN_sig_overlap_samedir_noancplast = resKresL_DPN_sig_overlap_samedir[resKresL_DPN_sig_overlap_samedir$Names %notin% resIresJ_sig_overlap$Names,]
dim(resKresL_DPN_sig_overlap_samedir_ancplast)
dim(resKresL_DPN_sig_overlap_samedir_noancplast)
write.csv(resKresL_DPN_sig_overlap_samedir_ancplast$Names, "resKresL_DPN_sig_overlap_samedir_ancplast_100622.csv")
write.csv(resKresL_DPN_sig_overlap_samedir_noancplast$Names, "resKresL_DPN_sig_overlap_samedir_noancplast_100622.csv")

resIresJ_resKresLDPN = merge(resKresL_DPN_sig_overlap_samedir, resIresJ_sig_overlap_samedir, by = "Names")

#KR8 - % DP also DE in S1/S2
dim(resIresJ_resKresLDPN)[1]/dim(resKresL_DPN_sig_overlap_samedir)[1]
#83.3%

Supplementary_Figure_2A = venn.diagram(x= list(resKresL_sig_overlap_samedir$Names, resBresC_sig_overlap_samedir$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), rotation.degree=180, cex= 2) 
Supplementary_Figure_2A
#pdf(file = "Supplementary_Figure_2A.pdf")
grid.draw(Supplementary_Figure_2A)
#dev.off()

#KR17: Overlap between resBresC and resKresL
resBresC_resKresLDPN = intersect(resBresC_sig_overlap_samedir$Names, resKresL_DPN_sig_overlap_samedir$Names)
length(resBresC_resKresLDPN) #21 genes...

#Total number of DP(N)+CEC genes KR23
resBresC_AND_resKresLDPN = unique(c(resBresC_sig_overlap_samedir$Names, resKresL_DPN_sig_overlap_samedir$Names))
length(resBresC_AND_resKresLDPN)/dim(ddTxi3)[1]

#"The adaptive components of evolutuonary change (AP and CEC genes) make up only
resKresLDPN_ANDOR_resBresC = unique(c(resKresL_DPN_sig_overlap_samedir$Names, resBresC_sig_overlap$Names))
length(resKresLDPN_ANDOR_resBresC)  
length(resKresLDPN_ANDOR_resBresC)/dim(ddTxi)[1]
#% of the transcriptome"

#Venn diagrams
#Figure S2D- BvC
#KR8(2/2)
Fig2B = venn.diagram(x= list(resB_sig$Names, resC_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), rotation.degree = 180, height = 2000, width = 2000, imagetype = "png", cex = 3) 
library(grDevices)
pdf(file = "S2D_VennDiagram_BC.pdf")
grid.draw(Fig2B)
dev.off()

#FigS2B
#KR4
venn.diagram(x= list(resF_sig$Names, resG_sig$Names), category.names = c("", ""), filename = "VennDiagram_FvG_140820.png", output = TRUE, cat.pos = c(0,0), rotation.degree = 180, width = 2000, height = 2000) 
Fig2D = venn.diagram(x= list(resF_sig$Names, resG_sig$Names), category.names = c("", ""), filename = NULL, output = TRUE, cat.pos = c(0,0), rotation.degree = 180, width = 2000, height = 1000, cex = 3) 
pdf(file = "S2B_VennDiagram_FG.pdf")
grid.draw(Fig2D)
dev.off()

#KR2
#FigS2A
venn.diagram(x= list(resI_sig$Names, resJ_sig$Names,resK_sig$Names, resL_sig$Names), category.names = c("S1", "S2", "T1", "T2"), filename = "VennDiagram_5_140820.png", output = TRUE, cat.pos = c(0,0,0,0), col = c(col_coast, col_coast, col_mine, col_mine)) 
Fig3A = venn.diagram(x= list(resI_sig$Names, resJ_sig$Names, resK_sig$Names, resL_sig$Names), category.names = c("", "", "", ""), filename = NULL, output = TRUE, cat.pos = c(0,0,0,0), cex = 2) 
pdf(file = "S2A_VennDiagram_IJKL.pdf")
grid.draw(Fig3A)
dev.off()

#FigS2C
venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = "VennDiagram_FG_KL_IJ_180820.png", output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, fill = c(col_mine, col_coast, col_zinc), width = 2000, height = 2000) 
Fig2F = venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = NULL, output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, fill = c(col_mine, col_coast, col_zinc), width = 2000, height = 2000) 
Fig2F = venn.diagram(x= list(resKresL_sig_overlap$Names, resIresJ_sig_overlap$Names, resFresG_sig_overlap$Names), category.names = c("", "", ""), filename = NULL, output = TRUE, cat.pos = c(0,0,0), rotation.degree = 180, width = 2000, height = 2000, cex = 3) 
pdf(file = "S2C_Venn_IJKLFG.pdf")
grid.draw(Fig2F)
dev.off()

#KR6
resFresG_resIresJ = resFresG_sig_overlap[resFresG_sig_overlap$Names %in% resIresJ_sig_overlap$Names,]
dim(resFresG_resIresJ)[1]/dim(resFresG_sig_overlap)[1]
resFresG_resKresL = resFresG_sig_overlap[resFresG_sig_overlap$Names %in% resKresL_sig_overlap$Names,]
dim(resFresG_resKresL)[1]/dim(resFresG_sig_overlap)[1]


###############################################
#10 PCAs of resKresL [AP] / resBresC [CEC]
###############################################

vst_resKresLDPN = vst_all[rownames(vst_all) %in% resKresL_DPN_sig_overlap_samedir$Names,]
PCA_resKresLDPN = plotPCA12_centroids(vst_resKresLDPN, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-45, 45),my_ylims = c(-20, 20), my_background_col = col_mine)
blork_resKresLDPN.x = as.data.frame(PCA_resKresLDPN$x)
blork_resKresLDPN.x$groups = samples$Pop_Cond
pca.centroids_resKresLDPN = aggregate(blork_resKresLDPN.x[,1:2], list(Type = blork_resKresLDPN.x$groups), mean)
PCA_resKresLDPN = plotPCA12(vst_resKresLDPN, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-45, 45),my_ylims = c(-20, 20), my_background_col = "lightgray")
pca.centroids_resKresLDPN

#KR16 (1/2)
PCA_resKresLDPN_centroids_SupplementaryFigure = PCA_resKresLDPN+
  geom_segment(aes(x = pca.centroids_resKresLDPN[7,2], y = pca.centroids_resKresLDPN[7,3], xend = pca.centroids_resKresLDPN[8,2], yend = pca.centroids_resKresLDPN[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[7,2], y = pca.centroids_resKresLDPN[7,3], xend = pca.centroids_resKresLDPN[8,2], yend = pca.centroids_resKresLDPN[8,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[1,2], y = pca.centroids_resKresLDPN[1,3], xend = pca.centroids_resKresLDPN[2,2], yend = pca.centroids_resKresLDPN[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5,col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[1,2], y = pca.centroids_resKresLDPN[1,3], xend = pca.centroids_resKresLDPN[2,2], yend = pca.centroids_resKresLDPN[2,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = coast_plot,  lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[5,2], y = pca.centroids_resKresLDPN[5,3], xend = pca.centroids_resKresLDPN[6,2], yend = pca.centroids_resKresLDPN[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[5,2], y = pca.centroids_resKresLDPN[5,3], xend = pca.centroids_resKresLDPN[6,2], yend = pca.centroids_resKresLDPN[6,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = 'round', linejoin = 'round')+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[3,2], y = pca.centroids_resKresLDPN[3,3], xend = pca.centroids_resKresLDPN[4,2], yend = pca.centroids_resKresLDPN[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 2.5, col = "black", lineend = "round", linejoin = "round")+ #Not sure exactly how to add these? 
  geom_segment(aes(x = pca.centroids_resKresLDPN[3,2], y = pca.centroids_resKresLDPN[3,3], xend = pca.centroids_resKresLDPN[4,2], yend = pca.centroids_resKresLDPN[4,3]), arrow = arrow(length = unit(0.03, "npc")), size = 1.5, col = mine_plot, lineend = "round", linejoin = "round") #Not sure exactly how to add these? 
#  geom_point(data = pca.centroids_resKresLDPN, size = 4, shape = 21, fill = c(PPBD_plot, GRSA_plot, PPBD_plot, GRSA_plot), alpha = 1)
PCA_resKresLDPN_centroids_SupplementaryFigure


vst_resBresC = vst_all[rownames(vst_all) %in% resBresC_sig_overlap_samedir$Names,]
plotPCA12(vst_resBresC, intgroup1 = "Cond", intgroup2 = "Ecotype", intgroup3 = "Geog", ntop = 1000000, my_xlims = c(-50, 50),my_ylims = c(-25, 25), my_background_col = gold_plot_light) 


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

#Figure S3A
Supplementary_Figure_1 = ggplot()+
  xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfc_BvC_all, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_fill_viridis_c()+
  geom_bin2d(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), binwidth=0.5, data = lfc_BvC_resBresC)+
  ylab("T2 v S2 log2(Fold Change)")+
  xlab("T1 v S1 log2(Fold Change)")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), panel.background = element_blank(), axis.text=element_text(size=20), axis.title = element_text(size=15))
Supplementary_Figure_1
#pdf(file = "~/RNA_Seq/Supplementary_Figure_1_1.pdf")
grid.draw(Supplementary_Figure_1)
#dev.off()

summary(lm(lfc_BvC_resBresC$log2FoldChange.y ~ lfc_BvC_resBresC$log2FoldChange.x))
#So not a very positive correlaiton...

#3B) LFC between treatments

lfcKvI = merge(as.data.frame(lfcK), as.data.frame(lfcI), by = "Names")
lfcIvJ = merge(as.data.frame(lfcI), as.data.frame(lfcJ), by = "Names")
lfcLvJ = merge(as.data.frame(lfcL), as.data.frame(lfcJ), by = "Names")
lfcKvL = merge(as.data.frame(lfcK), as.data.frame(lfcL), by = "Names")

#For each get lists for K, L and KvL genes. 
lfcIvJ_resKresLDPN = lfcIvJ[lfcIvJ$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]
lfcKvI_resKresLDPN = lfcKvI[lfcKvI$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]
lfcLvJ_resKresLDPN = lfcLvJ[lfcLvJ$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]
lfcKvL_resKresLDPN = lfcKvL[lfcKvL$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]

lfcIvJ_resIresJ = lfcIvJ[lfcIvJ$Names %in% resIresJ_sig_overlap$Names,]

#KR15
resKresLDPN_LFC_Plot = ggplot()+
  #xlim(-10,10)+ ylim(-10,10)+  
  #  geom_point(mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y), data = lfcIvJ, col = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_bin2d(mapping = aes(x = lfcKvL_resKresLDPN$log2FoldChange.x, y = lfcKvL_resKresLDPN$log2FoldChange.y), binwidth=0.3)+
  scale_fill_viridis_c()+
  xlim(-10, 10)+
  ylim(-10, 10)+
  geom_abline()+
  
  ylab("T2 Control to Zinc - log2(Fold Change)")+
  xlab("T1 Control to Zinc - log2(Fold Change)")+
  #  theme(panel.background = element_rect(fill = coast_plot))
  theme(text= element_text(size = 20), panel.border = element_rect(colour = "black", fill = NA, size = 4), panel.background = element_blank()) #, panel.grid.major=element_line(colour="lightgrey"))
resKresLDPN_LFC_Plot

summary(lm(lfcKvL_resKresLDPN$log2FoldChange.y ~ lfcKvL_resKresLDPN$log2FoldChange.x))

#1.077, adjusted R^2 = 0.94
summary(lm(lfcKvI_resKresLDPN$log2FoldChange.x ~ lfcKvI_resKresLDPN$log2FoldChange.y))
#Slope = 0.886, adj R^2 = 0.857
summary(lm(lfcLvJ_resKresLDPN$log2FoldChange.x ~ lfcLvJ_resKresLDPN$log2FoldChange.y))
#Slope = 0.911, adj. R^2 = 0.828
summary(lm(lfcIvJ_resKresLDPN$log2FoldChange.y ~ lfcIvJ_resKresLDPN$log2FoldChange.x))
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

###12A - transcriptome wide resids

resids_countsA = abs(lfcA$log2FoldChange)
resids_countsB = abs(lfcB$log2FoldChange)
resids_countsC = abs(lfcC$log2FoldChange)
resids_countsD = abs(lfcD$log2FoldChange)
resids_countsE = abs(lfcE$log2FoldChange)
resids_countsF = abs(lfcF$log2FoldChange)
resids_countsG = abs(lfcG$log2FoldChange)
resids_countsH = abs(lfcH$log2FoldChange)

#KR11 (1/4)
median(resids_countsA)
quantile(resids_countsA) #Well that's obviously very big. 
median(resids_countsB)
median(resids_countsC)
#KR11 (2/4)
median(resids_countsD)
quantile(resids_countsD)
median(resids_countsE)
median(resids_countsF)
median(resids_countsG)
median(resids_countsH)
quantile(resids_countsH)

both_CZ = melt(data.frame(resids_countsA, resids_countsE, resids_countsB, resids_countsF, resids_countsC, resids_countsG, resids_countsD, resids_countsH), na.rm = TRUE)
both_CZ_reduced = melt(data.frame(resids_countsA, resids_countsE, resids_countsD, resids_countsH), na.rm = TRUE)
both_C = melt(data.frame(resids_countsA, resids_countsB, resids_countsC, resids_countsD), na.rm = TRUE)
both_Z = melt(data.frame(resids_countsE, resids_countsF, resids_countsG, resids_countsH), na.rm = TRUE)

#KR11 (3/4)
C_wilcox = wilcox.test(resids_countsA, resids_countsD, paired = TRUE, alternative = "two.sided")
C_wilcox

C_Mood = mood.medtest(resids_countsA, resids_countsB)
both_CZ_wilcox = pairwise.wilcox.test(both_CZ$value, both_CZ$variable, p.adjust.method = "BH", paired = TRUE)
both_CZ_wilcox
both_C_wilcox = pairwise.wilcox.test(both_C$value, both_C$variable, p.adjust.method = "BH", paired = TRUE)
both_C_wilcox

pairwiseMedianTest(value ~ variable, data=both_CZ) #Alternative stats test
#write.csv(as.data.frame(both_CZ_wilcox$p.value), "both_CZ_wilcox.csv")

#KR11 (4/4)
#Supplementary Figure - Transcriptome-wide |FC| in control
Fig5A_nocol = ggplot()+
  geom_boxplot(fill = ctrl_green, data = both_C, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col))+
  theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15),  panel.background = element_rect(fill = "white", colour = "black", size = 2), )+
  xlab("")+ylab("|FC|")+coord_cartesian(ylim=c(0,1))
Fig5A_nocol
#pdf(file = "TranscriptomeWide_FC_Control.pdf")
#grid.draw(Fig5A_nocol)
#dev.off()

#resBresC_sig_overlap_samedir; "resids" of these

resids_countsA_resBresC = abs(lfcA[lfcA$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsB_resBresC = abs(lfcB[lfcB$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsC_resBresC = abs(lfcC[lfcC$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsD_resBresC = abs(lfcD[lfcD$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsE_resBresC = abs(lfcE[lfcE$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsF_resBresC = abs(lfcF[lfcF$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsG_resBresC = abs(lfcG[lfcG$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsH_resBresC = abs(lfcH[lfcH$Names %in% resBresC_sig_overlap_samedir$Names,]$log2FoldChange)

#KR7 (1/3)- CEC |FC|S1-S2 = 0.055
median(resids_countsA_resBresC)
quantile(resids_countsA_resBresC)
median(resids_countsB_resBresC)
median(resids_countsC_resBresC)
#KR7 (2/3) = CEC |FC|T1-T2 = 0.12
median(resids_countsD_resBresC)
quantile(resids_countsD_resBresC)
median(resids_countsE_resBresC)
median(resids_countsF_resBresC)
median(resids_countsG_resBresC)
median(resids_countsH_resBresC)
quantile(resids_countsH_resBresC)

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
#write.csv(as.data.frame(both_CZ_resBresC_wilcox$p.value), "both_CZ_resBresC_wilcox.csv")
pairwiseMedianTest(value ~ variable, data=both_CZ_resBresC_all) #Could also include this? 
both_CZ_resBresC = melt(data.frame(resids_countsA_resBresC, resids_countsB_resBresC, resids_countsC_resBresC, resids_countsD_resBresC), na.rm = TRUE)

#KR12 - Graph showing CEC |FC| values
boxplot_resBresC_FC_220721 = ggplot()+
  geom_boxplot(fill = ctrl_green, data = both_CZ_resBresC, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  # scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 23), axis.text.y = element_text(size = 15))+coord_cartesian(ylim=c(0,10))+
  #scale_x_discrete(labels = c("S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+0
  xlab("")+ylab("Absolute log2 fold change (|FC|)")+
  #  theme(panel.background = element_rect(fill = NA, colour = "black", size = 4))
  theme(panel.background = element_rect(fill = NA), panel.border=element_rect(color="black", fill = NA, size = 4))
boxplot_resBresC_FC_220721
#pdf(file = "boxplot_resBresC_FC_220721.pdf")
grid.draw(boxplot_resBresC_FC_220721)
#dev.off()

#I guess this is just the resA and resD comparison. Not used but useful to have. 
boxplot_resBresC_FC_220721 = ggplot()+
  geom_boxplot(fill = ctrl_green, data = both_CZ_resBresC_two, aes(x = variable, y=value), width = 0.6, outlier.shape = NA, show.legend = FALSE)+
  # scale_fill_manual(values = c(ctrl_col, ctrl_col, ctrl_col, ctrl_col, ctrl_col, zinc_col, ctrl_col, zinc_col))+
  theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15))+ylim(0,0.4)+
  scale_x_discrete(labels = c("S1 v S2 (C)", "T1 v S1 (C)", "T2 v S2 (C)", "T1 v T2 (C)"))+
  xlab("")+ylab("|FC|")+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1))
boxplot_resBresC_FC_220721

########## resKresL resids

#########Compare across conditions

resids_countsA_resKresLDPN = abs(lfcA[lfcA$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsB_resKresLDPN = abs(lfcB[lfcB$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsC_resKresLDPN = abs(lfcC[lfcC$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsD_resKresLDPN = abs(lfcD[lfcD$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsE_resKresLDPN = abs(lfcE[lfcE$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsF_resKresLDPN = abs(lfcF[lfcF$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsG_resKresLDPN = abs(lfcG[lfcG$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)
resids_countsH_resKresLDPN = abs(lfcH[lfcH$Names %in% resKresL_DPN_sig_overlap_samedir$Names,]$log2FoldChange)

median(resids_countsA_resKresLDPN)
median(resids_countsB_resKresLDPN)
median(resids_countsC_resKresLDPN)
median(resids_countsD_resKresLDPN)
median(resids_countsE_resKresLDPN)
median(resids_countsF_resKresLDPN)
median(resids_countsG_resKresLDPN)
median(resids_countsH_resKresLDPN)
quantile(resids_countsH_resKresLDPN)

2^median(resids_countsA_resKresLDPN)
2^median(resids_countsB_resKresLDPN)
2^median(resids_countsC_resKresLDPN)
2^median(resids_countsD_resKresLDPN)
2^median(resids_countsE_resKresLDPN)
2^median(resids_countsF_resKresLDPN)
2^median(resids_countsG_resKresLDPN)
2^median(resids_countsH_resKresLDPN)

resids_CZ_resKresLDPN = melt(data.frame(resids_countsA_resKresLDPN, resids_countsE_resKresLDPN, resids_countsB_resKresLDPN, resids_countsF_resKresLDPN, resids_countsC_resKresLDPN, resids_countsG_resKresLDPN, resids_countsD_resKresLDPN, resids_countsH_resKresLDPN), na.rm = TRUE)
resids_CZ_resKresLDPN_wilcox = pairwise.wilcox.test(resids_CZ_resKresLDPN$value, resids_CZ_resKresLDPN$variable, p.adjust.method = "BH")
resids_CZ_resKresLDPN_wilcox
#write.csv(as.data.frame(resids_CZ_resKresLDPN_wilcox$p.value), "resids_CZ_resKresLDPN_wilcox.csv")
pairwiseMedianTest(value ~ variable, data=resids_CZ_resKresLDPN) #Could also include this? 

resids_CZ_resKresLDPN_eco = melt(data.frame(resids_countsA_resKresLDPN, resids_countsE_resKresLDPN, resids_countsD_resKresLDPN, resids_countsH_resKresLDPN), na.rm  = TRUE)

#KR16: DP(N) |FC|s between ecotypes/conditions: T1/T2 get more similar in zinc
Supplementary_Figure_5A = ggplot()+
  #  geom_rect(aes(xmin = 0.43, xmax = 2.5, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_coast)+
  # geom_rect(aes(xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf),alpha = 0.6, fill = col_mine)+
  
  #geom_half_violin(data = resids_CZ_resKresLDPN_eco, aes(x = variable, y=value, fill = variable), side = "l", show.legend = FALSE)+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_mine, col_mine, col_mine, col_mine))+
  new_scale_fill()+
  scale_fill_manual(values = c(col_coast, col_coast, col_mine, col_mine, col_coast, col_coast, col_mine, col_mine))+
  #geom_half_violin(data = resids_CZ_resKresLDPN_eco, aes(x = variable, y=value, fill = variable), side = "r", show.legend = FALSE)+
  new_scale_fill()+
  geom_boxplot(fill = c(ctrl_green, another_plot_purple, ctrl_green, another_plot_purple), data = resids_CZ_resKresLDPN_eco, aes(x = variable, y=value), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("S1C v S2 (C)", "S1 v S2 (Z)", "T1 v T2 (C)", "T1 v T2 (Z)"))+
  xlab("")+ylab("|FC|")+coord_cartesian(ylim=c(0,0.9))+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2), text = element_text(size = 15))

Supplementary_Figure_5A #So a pretty similar look result, they become more similar in zinc.
#pdf(file = "Supplementary_Figure_5A.pdf")
grid.draw(Supplementary_Figure_5A)
#dev.off()

#Some permutations?

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

#############################################################################################

##################################################################################################
#13 Plasticity - REV/RI classification, |FC| of plastic vs. nonplastic, parametric bootstrapping
##################################################################################################

#Need to calculate the Lo, Lp and La for various sets of genes. 
#Then include those where Lp and La both are above a certain cutoff...

#Previous authiors used TMM normalised counts measured in edgeR. 

######################################
#### 13A - Simulating (a/b) ~ (c/b)
#######################################
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

#######################################
#### 13B.1 - calculating Lo/Lp/La, EC/GC, REV/RI/OVER; entire dataset.
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

#New
#So here we have resX, significant in C between conditions, and resY, significant in Z between ecotypes.
resXY = resX[resX_sig$Names %in% resY_sig$Names,]
norm_counts_cut_XY = norm_counts[rownames(norm_counts) %in% resXY$Names,]
#KR18
dim(norm_counts_cut_XY)[1]/dim(norm_counts)[1]
#MA6.2 Significance cutoff.

norm_counts_cut_XY$DirStat = norm_counts_cut_XY$GC*norm_counts_cut_XY$PC
norm_counts_cut_XY$SameDir = ifelse(norm_counts_cut_XY$DirStat > 0, "Y", "N")
norm_counts_cut_XY$logPC = norm_counts_cut_XY$PC*log2(abs(norm_counts_cut_XY$PC))/abs(norm_counts_cut_XY$PC)
head(norm_counts_cut_XY$logPC)
head(norm_counts_cut_XY$PC)
norm_counts_cut_XY$logGC = norm_counts_cut_XY$GC*log2(abs(norm_counts_cut_XY$GC))/abs(norm_counts_cut_XY$GC)
#KR12
dim(norm_counts_cut_XY)[1]/dim(norm_counts)[1]
norm_counts_noPC_XY = norm_counts[rownames(norm_counts) %in% resY_sig$Names & rownames(norm_counts) %notin% resX_sig$Names,]

dim(norm_counts_noPC_XY)
dim(resX_sig) #~17k
dim(resY_sig) #~15k
length(intersect(resX_sig$Names, resY_sig$Names)) #

#KR12: 76.7% of genes have substnatial PC and EC
dim(norm_counts_cut_XY)[1]/dim(norm_counts)[1] #...well...this is now 73% of all genes. OK. 


norm_counts_cut_XY$Type = ifelse(norm_counts_cut_XY$DirStat > 0, "RI", ifelse(abs(norm_counts_cut_XY$GC) > 0.5*abs(norm_counts_cut_XY$PC), "REV", "OVERSHOOT"))



#MA6.1: Significance cutoff only. (1/N)
#So for these I think you'd want things that display this behaviour in both replicates
#aka are in I, J, F, G
#resIJFG = Reduce(intersect, list(resI_sig$Names, resJ_sig$Names, resF_sig$Names, resG_sig$Names))
#length(resIJFG)
#resIJ = Reduce(intersect, list(resI_sig$Names, resJ_sig$Names))
#length(resIJ)
#resFG = Reduce(intersect, list(resF_sig$Names, resG_sig$Names))
#length(resFG)
#length(resIJFG) #Still only 9k genes. That's not that many.

#norm_counts_cut = norm_counts[rownames(norm_counts) %in% resIJFG,]

#So...get those that are DE in these conditions regardless of ecotype.

#MA6.2: Signficance cutoff just with ecoytpe differences...
#Should increase number of genes identified.

####################################################################################


########################################
### 13B.2 - considering each individual geographic pair individually. 
########################################

#Individual geographic pairs...

#Original way...
norm_counts$Lo_G1 = rowMeans(norm_counts[,c("BD11C", "BD12C", "BD1C")])
norm_counts$Lp_G1 = rowMeans(norm_counts[,c("BD11Z", "BD12Z", "BD1Z")])
norm_counts$La_G1 = rowMeans(norm_counts[,c("PP1Z", "PP4Z", "PP8Z")])

norm_counts$Lo_G2 = rowMeans(norm_counts[,c("SA6C", "SA7C", "SA8C")])
norm_counts$Lp_G2 = rowMeans(norm_counts[,c("SA6Z", "SA7Z", "SA8Z")])
norm_counts$La_G2 = rowMeans(norm_counts[,c("GR10Z", "GR2Z", "GR8Z")])

norm_counts$PC_G1 = norm_counts$Lp_G1 - norm_counts$Lo_G1
norm_counts$GC_G1 = norm_counts$La_G1 - norm_counts$Lp_G1
norm_counts$PC_G2 = norm_counts$Lp_G2 - norm_counts$Lo_G2
norm_counts$GC_G2 = norm_counts$La_G2 - norm_counts$Lp_G2

norm_counts$Cutoff_GC_G1 = ifelse(abs(norm_counts$GC_G1) > cutoff*norm_counts$Lo_G1, "Y", "N")
norm_counts$Cutoff_PC_G1 = ifelse(abs(norm_counts$PC_G1) > cutoff*norm_counts$Lo_G1, "Y", "N")
norm_counts$Cutoff_GC_G2 = ifelse(abs(norm_counts$GC_G2) > cutoff*norm_counts$Lo_G2, "Y", "N")
norm_counts$Cutoff_PC_G2 = ifelse(abs(norm_counts$PC_G2) > cutoff*norm_counts$Lo_G2, "Y", "N")

#MA6 New (XY)
rownames(norm_counts)
#So here we will just need the ones that are...a) in resI and resF and b) in resJ and resG
norm_counts_cut_XY_G1 = norm_counts[rownames(norm_counts) %in% resI_sig$Names & rownames(norm_counts) %in% resF_sig$Names,]
norm_counts_cut_XY_G2 = norm_counts[rownames(norm_counts) %in% resJ_sig$Names & rownames(norm_counts) %in% resG_sig$Names,]

norm_counts_cut_XY_G1$DirStat = norm_counts_cut_XY_G1$GC_G1*norm_counts_cut_XY_G1$PC_G1
norm_counts_cut_XY_G2$DirStat = norm_counts_cut_XY_G2$GC_G2*norm_counts_cut_XY_G2$PC_G2

norm_counts_cut_XY_G1$Type = ifelse(norm_counts_cut_XY_G1$DirStat > 0, "RI", ifelse(abs(norm_counts_cut_XY_G1$GC_G1) > 0.5*abs(norm_counts_cut_XY_G1$PC_G1), "REV", "OVERSHOOT"))
norm_counts_cut_XY_G2$Type = ifelse(norm_counts_cut_XY_G2$DirStat > 0, "RI", ifelse(abs(norm_counts_cut_XY_G2$GC_G2) > 0.5*abs(norm_counts_cut_XY_G2$PC_G2), "REV", "OVERSHOOT"))

100*table(norm_counts_cut_XY_G1$Type)[1]/nrow(norm_counts_cut_XY_G1) #% overshooting
100*table(norm_counts_cut_XY_G1$Type)[2]/nrow(norm_counts_cut_XY_G1) #% reversion
100*table(norm_counts_cut_XY_G1$Type)[3]/nrow(norm_counts_cut_XY_G1) #% reinforcement

100*table(norm_counts_cut_XY_G2$Type)[1]/nrow(norm_counts_cut_XY_G2) #% overshooting
100*table(norm_counts_cut_XY_G2$Type)[2]/nrow(norm_counts_cut_XY_G2) #% reversion
100*table(norm_counts_cut_XY_G2$Type)[3]/nrow(norm_counts_cut_XY_G2) #% reinforcement

100*dim(norm_counts_cut_XY_G1)[1]/dim(norm_counts)[1]
100*dim(norm_counts_cut_XY_G2)[1]/dim(norm_counts)[1]
#...so this is a much lower % of genes. But that is what we expect? A higher degree of evidence required here. 

####################################################
#13C - Whole transcriptome - combined data
####################################################

#New way (XY)
#
100*dim(norm_counts_cut_XY)[1]/dim(norm_counts)[1] 
#KR19
100*table(norm_counts_cut_XY$Type)[1]/nrow(norm_counts_cut_XY) #% overshooting
100*table(norm_counts_cut_XY$Type)[2]/nrow(norm_counts_cut_XY) #% reversion
100*table(norm_counts_cut_XY$Type)[3]/nrow(norm_counts_cut_XY) #% reinforcement

#Yup ok so....

table(norm_counts_cut_XY$Type)[1]
table(norm_counts_cut_XY$Type)[2]
table(norm_counts_cut_XY$Type)[3]

head(norm_counts_cut_XY)
T_table = table(norm_counts_cut_XY$Type)
T_vec = c(0,0,0)
T_vec[2] = table(norm_counts_cut_XY$Type)[1] #% overshooting
T_vec[1] = table(norm_counts_cut_XY$Type)[2] #% reversion
T_vec[3] = table(norm_counts_cut_XY$Type)[3] #% reinforcement
T_vec = data.frame(val = T_vec)
T_vec = cbind(T_vec, type = c("REV", "OVER", "RI"))
T_vec$type= factor(T_vec$type, levels = T_vec$type)

#|FC| values for whole transcriptome
lfcA_T_REV = lfcA[lfcA$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "REV",]),]
lfcE_T_REV = lfcE[lfcE$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "REV",]),]
lfcD_T_REV = lfcD[lfcD$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "REV",]),]
lfcH_T_REV = lfcH[lfcH$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "REV",]),]
resids_countsA_T_REV = abs(lfcA_T_REV$log2FoldChange)
resids_countsE_T_REV = abs(lfcE_T_REV$log2FoldChange)
resids_countsD_T_REV = abs(lfcD_T_REV$log2FoldChange)
resids_countsH_T_REV = abs(lfcH_T_REV$log2FoldChange)

lfcA_T_RI = lfcA[lfcA$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "RI",]),]
lfcE_T_RI = lfcE[lfcE$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "RI",]),]
lfcD_T_RI = lfcD[lfcD$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "RI",]),]
lfcH_T_RI = lfcH[lfcH$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "RI",]),]
resids_countsA_T_RI = abs(lfcA_T_RI$log2FoldChange)
resids_countsE_T_RI = abs(lfcE_T_RI$log2FoldChange)
resids_countsD_T_RI = abs(lfcD_T_RI$log2FoldChange)
resids_countsH_T_RI = abs(lfcH_T_RI$log2FoldChange)

lfcA_T_OVERSHOOT = lfcA[lfcA$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "OVERSHOOT",]),]
lfcE_T_OVERSHOOT = lfcE[lfcE$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "OVERSHOOT",]),]
lfcD_T_OVERSHOOT = lfcD[lfcD$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "OVERSHOOT",]),]
lfcH_T_OVERSHOOT = lfcH[lfcH$Names %in% rownames(norm_counts_cut_XY[norm_counts_cut_XY$Type == "OVERSHOOT",]),]
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
#So this is using the XY_method. I guess we will also need a version of this with the older method?

##########
#MA6 old method (for supplement)
###
#########


#########

#####################################################################

######################################################
#13D.1 -reKresL_samedir genes (AP genes)
######################################################

#...right. 
#a) Change to resKresL_DPN
#b) change to the XY method. 

dim(resKresL_DPN_sig_overlap_samedir)
norm_counts_resKresLDPN_XY = norm_counts_cut_XY[rownames(norm_counts_cut_XY) %in% resKresL_DPN_sig_overlap_samedir$Names,]
norm_counts_resKresLDPN_XY = norm_counts_cut_XY[rownames(norm_counts_cut_XY) %in% resKresL_DPN_sig_overlap_samedir$Names,]

#KR24: DP(N) genes that have "substantial" |PC| and |EC|
dim(norm_counts_resKresLDPN_XY)[1]/dim(resKresL_DPN_sig_overlap_samedir)[1]

#KR25: DP(N) gene REV, RI, OVER
100*table(norm_counts_resKresLDPN_XY$Type)[1]/nrow(norm_counts_resKresLDPN_XY) #% overshooting- 14%
100*table(norm_counts_resKresLDPN_XY$Type)[2]/nrow(norm_counts_resKresLDPN_XY) #% reversion - 66%
100*table(norm_counts_resKresLDPN_XY$Type)[3]/nrow(norm_counts_resKresLDPN_XY) #% reinforcement - 18.7% 

table(norm_counts_resKresLDPN_XY$Type)[2]
table(norm_counts_resKresLDPN_XY$Type)[1]
table(norm_counts_resKresLDPN_XY$Type)[3]

library(cowplot)

KL_Table = table(norm_counts_resKresLDPN_XY$Type)
KL_vec = c(0,0,0)
KL_vec[2] = table(norm_counts_resKresLDPN_XY$Type)[1] #% overshooting
KL_vec[1] = table(norm_counts_resKresLDPN_XY$Type)[2] #% reversion
KL_vec[3] = table(norm_counts_resKresLDPN_XY$Type)[3] #% reinforcement
KL_vec = data.frame(val = KL_vec)
KL_vec = cbind(KL_vec, type = c("REV", "OVER", "RI"))
KL_vec$type= factor(KL_vec$type, levels = KL_vec$type)


woof = binom.test(x = c((KL_Table[1]+KL_Table[3]), KL_Table[2]), p = (T_table[1]+T_table[3])/(T_table[1]+T_table[2]+T_table[3]), alternative = "two.sided")
#Note - no test statistic here.

#Are "adaptive" plastic changes overrepresented in this gene set?M 

#|FC| values for whole transcriptome
lfcA_KL_REV = lfcA[lfcA$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "REV",]),]
lfcE_KL_REV = lfcE[lfcE$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "REV",]),]
lfcD_KL_REV = lfcD[lfcD$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "REV",]),]
lfcH_KL_REV = lfcH[lfcH$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "REV",]),]
resids_countsA_KL_REV = abs(lfcA_KL_REV$log2FoldChange)
resids_countsE_KL_REV = abs(lfcE_KL_REV$log2FoldChange)
resids_countsD_KL_REV = abs(lfcD_KL_REV$log2FoldChange)
resids_countsH_KL_REV = abs(lfcH_KL_REV$log2FoldChange)

lfcA_KL_RI = lfcA[lfcA$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "RI",]),]
lfcE_KL_RI = lfcE[lfcE$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "RI",]),]
lfcD_KL_RI = lfcD[lfcD$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "RI",]),]
lfcH_KL_RI = lfcH[lfcH$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "RI",]),]
resids_countsA_KL_RI = abs(lfcA_KL_RI$log2FoldChange)
resids_countsE_KL_RI = abs(lfcE_KL_RI$log2FoldChange)
resids_countsD_KL_RI = abs(lfcD_KL_RI$log2FoldChange)
resids_countsH_KL_RI = abs(lfcH_KL_RI$log2FoldChange)

lfcA_KL_OVERSHOOT = lfcA[lfcA$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "OVERSHOOT",]),]
lfcE_KL_OVERSHOOT = lfcE[lfcE$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "OVERSHOOT",]),]
lfcD_KL_OVERSHOOT = lfcD[lfcD$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "OVERSHOOT",]),]
lfcH_KL_OVERSHOOT = lfcH[lfcH$Names %in% rownames(norm_counts_resKresLDPN_XY[norm_counts_resKresLDPN_XY$Type == "OVERSHOOT",]),]
resids_countsA_KL_OVERSHOOT = abs(lfcA_KL_OVERSHOOT$log2FoldChange)
resids_countsE_KL_OVERSHOOT = abs(lfcE_KL_OVERSHOOT$log2FoldChange)
resids_countsD_KL_OVERSHOOT = abs(lfcD_KL_OVERSHOOT$log2FoldChange)
resids_countsH_KL_OVERSHOOT = abs(lfcH_KL_OVERSHOOT$log2FoldChange)

norm_counts_resKresLDPN_XY$PC_pc = norm_counts_resKresLDPN_XY$PC/norm_counts_resKresLDPN_XY$Lp
norm_counts_resKresLDPN_XY$GC_pc = norm_counts_resKresLDPN_XY$GC/norm_counts_resKresLDPN_XY$Lp

norm_counts_noPC_resKresLDPN_XY = norm_counts_noPC_XY[rownames(norm_counts_noPC_XY) %in% resKresL_DPN_sig_overlap_samedir$Names,]
dim(norm_counts_noPC_resKresLDPN_XY) #Note: only 18 genes.
lfcA_KL_NOANCPLAST = lfcA[lfcA$Names %in% rownames(norm_counts_noPC_resKresLDPN_XY),]
lfcE_KL_NOANCPLAST = lfcE[lfcE$Names %in% rownames(norm_counts_noPC_resKresLDPN_XY),]
lfcD_KL_NOANCPLAST = lfcD[lfcD$Names %in% rownames(norm_counts_noPC_resKresLDPN_XY),]
lfcH_KL_NOANCPLAST = lfcH[lfcH$Names %in% rownames(norm_counts_noPC_resKresLDPN_XY),]

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
Plast_FC_KL_wilcox
median(resids_countsH_KL_RI$VALUE)
median(resids_countsH_KL_REV$VALUE)
median(resids_countsH_KL_OVERSHOOT$VALUE)

temp_Plast_FC_KL = Plast_FC_KL
temp_Plast_FC_KL$Names = "PLAST"
median(temp_Plast_FC_KL$VALUE)
quantile(temp_Plast_FC_KL$VALUE)
median(resids_countsH_KL_NOANCPLAST$VALUE)
quantile(resids_countsH_KL_NOANCPLAST$VALUE)


PlastNoPlast_FC_KL = rbind(temp_Plast_FC_KL, resids_countsH_KL_NOANCPLAST)
median(temp_Plast_FC_KL$VALUE) #Plast
median(resids_countsH_KL_NOANCPLAST$VALUE) #No Plast... Ahhhh dear, this is not the right way round...

#KR36 (1/2)
PlastNoPlast_FC_KL_wilcox = wilcox.test(temp_Plast_FC_KL$VALUE, resids_countsH_KL_NOANCPLAST$VALUE, alternative = "two.sided")
PlastNoPlast_FC_KL_wilcox #no significant difference (insufficient numbers I guess...)

KL_PlastNoPlast_Boxplot = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_KL, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("NP", "P"))+
  xlab("")+ylab("|FC|")+ylim(0,0.625)+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2), text = element_text(size = 15))
KL_PlastNoPlast_Boxplot
#Should note there aren't enough here.

median(resids_countsE_KL_RI)
median(resids_countsE_KL_REV)
median(resids_countsE_KL_OVERSHOOT)

#Ok...I think this probably needs some more thought. 

pairwiseMedianTest(value ~ variable, data=both_CZ) #Might be better to do this, really.

###################################################################
#13D.2 AP genes only in one of the two replicates...
###################################################################

#So will need to tweak this based on...well DP(N) really. Which is what. MA5
#Old
#resK_notL = resK_sig[resK_sig$Names %notin% resKresL_sig_overlap_samedir$Names,] 
#dim(resK_notL) #485 in G1 not G2
#resL_notK = resL_sig[resL_sig$Names %notin% resKresL_sig_overlap_samedir$Names,]
#dim(resL_notK) #2365 in G2 not G1

length(BorF)
length(CorG)
dim(resK_sig)
dim(resL_sig)
#MA5
resKDPN_notLDPN = resK_DPN[resK_DPN$Names %notin% resL_DPN$Names,]
dim(resKDPN_notLDPN)
resLDPN_notKDPN = resL_DPN[resL_DPN$Names %notin% resK_DPN$Names,]
dim(resLDPN_notKDPN) #2365 in G2 not G1

#Ok so for the individaul things...
norm_counts_cut_XY_G1_resKDPN_notLDPN = norm_counts_cut_XY_G1[rownames(norm_counts_cut_XY_G1) %in% resKDPN_notLDPN$Names,]
norm_counts_cut_XY_G2_resLDPN_notKDPN = norm_counts_cut_XY_G2[rownames(norm_counts_cut_XY_G2) %in% resLDPN_notKDPN$Names,]

#i) what proportion of these genes don't meet the EC/PC cutoffs?
dim(norm_counts_cut_XY_G1_resKDPN_notLDPN) #405 meet cutoff
dim(norm_counts_cut_XY_G1_resKDPN_notLDPN)[1]/nrow(resKDPN_notLDPN) #So 83.5% pass both thresholds

dim(norm_counts_cut_XY_G2_resLDPN_notKDPN) #1564
dim(norm_counts_cut_XY_G2_resLDPN_notKDPN)/nrow(resLDPN_notKDPN) #So only 66% pass both thresholds.

#ii) What proportion show substantial ancestral plasticity?
#For this we will need to go back and look at the pre-cutoff files...
norm_counts_G1_resKDPN_notLDPN = norm_counts[rownames(norm_counts) %in% resKDPN_notLDPN$Names,]
norm_counts_G2_resLDPN_notKDPN = norm_counts[rownames(norm_counts) %in% resLDPN_notKDPN$Names,]
dim(norm_counts_G1_resKDPN_notLDPN) #Same as before, good
dim(norm_counts_G2_resLDPN_notKDPN) #Same as before, good

#So for G1...
table(norm_counts_G1_resKDPN_notLDPN$Cutoff_GC_G1) #So 437/485 pass the GC cutoff
table(norm_counts_G1_resKDPN_notLDPN$Cutoff_GC_G1)[2]/nrow(resKDPN_notLDPN) #so that's 90%
table(norm_counts_G1_resKDPN_notLDPN$Cutoff_PC_G1) #So 448/485.
table(norm_counts_G1_resKDPN_notLDPN$Cutoff_PC_G1)[2]/nrow(resKDPN_notLDPN) #so 92% pass PC as well. 

table(norm_counts_G2_resLDPN_notKDPN$Cutoff_GC_G2) #So 683/2365
table(norm_counts_G2_resLDPN_notKDPN$Cutoff_GC_G2)[2]/nrow(resLDPN_notKDPN) #...71.2% show substantial GC (i.e. derived genetic change)
table(norm_counts_G2_resLDPN_notKDPN$Cutoff_PC_G2) #So 3461/4035, or...
table(norm_counts_G2_resLDPN_notKDPN$Cutoff_PC_G2)[2]/nrow(resLDPN_notKDPN) #but still, 92.6% show substantail PC

#Ok so I think the real question is...for genes in resKDPN_notLDPN [G1 only}...do they show substantially lower anc. plast. in G2? I.e. is that why they're not DGE?
table(norm_counts_G1_resKDPN_notLDPN$Cutoff_PC_G2)[2]/nrow(resKDPN_notLDPN) #Nope...
table(norm_counts_G2_resLDPN_notKDPN$Cutoff_PC_G1)[2]/nrow(resLDPN_notKDPN) #Nope...

#So there isn't really strong evidence here that a lack of ancestral plasticity is driving recruitment.

norm_counts_cut_XY_G1_resKresLDPN = norm_counts_cut_XY_G1[rownames(norm_counts_cut_XY_G1) %in% resKresL_DPN_sig_overlap_samedir$Names,]
norm_counts_cut_XY_G2_resKresLDPN = norm_counts_cut_XY_G2[rownames(norm_counts_cut_XY_G2) %in% resKresL_DPN_sig_overlap_samedir$Names,]
dim(norm_counts_cut_XY_G1) #very similar
dim(norm_counts_cut_XY_G2) #to each other

#So 112 in both. I guess..conceivable?
dim(norm_counts_cut_XY_G1_resKresLDPN)[1]
dim(norm_counts_cut_XY_G2_resKresLDPN)[1]

dim(norm_counts_cut_XY_G1_resKresLDPN)[1]/nrow(resKresL_DPN_sig_overlap_samedir) #291:
dim(norm_counts_cut_XY_G2_resKresLDPN)[1]/nrow(resKresL_DPN_sig_overlap_samedir) #300:
#So similar across both pairs...

#SO probably some pertinent questions are:
#Are the shared genes showing similar patterns in the individual ones.

#KR34: More reversion in single population "DP(N)" genes (89.9% and 77.6%) compared to shared ones (70.8%). 
100*table(norm_counts_cut_XY_G1_resKDPN_notLDPN$Type)[1]/nrow(norm_counts_cut_XY_G1_resKDPN_notLDPN) #% overshooting: 4.4
100*table(norm_counts_cut_XY_G1_resKDPN_notLDPN$Type)[2]/nrow(norm_counts_cut_XY_G1_resKDPN_notLDPN) #% reversion: 89.9%
100*table(norm_counts_cut_XY_G1_resKDPN_notLDPN$Type)[3]/nrow(norm_counts_cut_XY_G1_resKDPN_notLDPN) #% reinforcement; 5.7%
table(norm_counts_cut_XY_G1_resKDPN_notLDPN$Type)[3] #So there are 23 reinforced genes in here...

100*table(norm_counts_cut_XY_G2_resLDPN_notKDPN$Type)[1]/nrow(norm_counts_cut_XY_G2_resLDPN_notKDPN) #% overshooting 13.1%
100*table(norm_counts_cut_XY_G2_resLDPN_notKDPN$Type)[2]/nrow(norm_counts_cut_XY_G2_resLDPN_notKDPN) #% reversion 77.6%
100*table(norm_counts_cut_XY_G2_resLDPN_notKDPN$Type)[3]/nrow(norm_counts_cut_XY_G2_resLDPN_notKDPN) #% reinforcement; 9.27%

table(norm_counts_cut_XY_G2_resLDPN_notKDPN$Type)[3] #And 290 reinforced genes in here.

100*table(norm_counts_cut_XY_G1_resKresLDPN$Type)[1]/nrow(norm_counts_cut_XY_G1_resKresLDPN) #% overshooting: 11.7%
100*table(norm_counts_cut_XY_G1_resKresLDPN$Type)[2]/nrow(norm_counts_cut_XY_G1_resKresLDPN) #% reversion: 69.8%
100*table(norm_counts_cut_XY_G1_resKresLDPN$Type)[3]/nrow(norm_counts_cut_XY_G1_resKresLDPN) #% reinforcement: 18.5%: so way more???

100*table(norm_counts_cut_XY_G2_resKresLDPN$Type)[1]/nrow(norm_counts_cut_XY_G2_resKresLDPN) #% overshooting: 11.8%
100*table(norm_counts_cut_XY_G2_resKresLDPN$Type)[2]/nrow(norm_counts_cut_XY_G2_resKresLDPN) #% reversion: 71.6%
100*table(norm_counts_cut_XY_G2_resKresLDPN$Type)[3]/nrow(norm_counts_cut_XY_G2_resKresLDPN) #% reinforcement: 16.6%...

#So this is done. I guess it needs to be redone for the 0.2x thing, right?

##################################################################
#13E.1 resBresC_samedir genes (CEC genes)
##################################################################
norm_counts_XY_resBresC = norm_counts_cut_XY[rownames(norm_counts_cut_XY) %in% resBresC_sig_overlap_samedir$Names,]
norm_counts_XY_resBresC = norm_counts_cut_XY[rownames(norm_counts_cut_XY) %in% resBresC_sig_overlap_samedir$Names,]

#KR28 - CEC genes % showing "substantial" EC and PC. 

dim(norm_counts_XY_resBresC)[1]/nrow(resBresC_sig_overlap_samedir)
dim(resBresC_sig_overlap_samedir)
#KR29 - CEC genes, OVER=3.3%, REV=62.7%, RI=34%. See graphs. 
100*table(norm_counts_XY_resBresC$Type)[1]/(table(norm_counts_XY_resBresC$Type)[1]+table(norm_counts_XY_resBresC$Type)[2]+table(norm_counts_XY_resBresC$Type)[3]) #% overshooting
100*table(norm_counts_XY_resBresC$Type)[2]/(table(norm_counts_XY_resBresC$Type)[1]+table(norm_counts_XY_resBresC$Type)[2]+table(norm_counts_XY_resBresC$Type)[3]) #% reversion
100*table(norm_counts_XY_resBresC$Type)[3]/(table(norm_counts_XY_resBresC$Type)[1]+table(norm_counts_XY_resBresC$Type)[2]+table(norm_counts_XY_resBresC$Type)[3]) #% reinforcement

BC_Table = table(norm_counts_XY_resBresC$Type)
BC_Table

BC_vec = c(0,0,0)
BC_vec[2] = table(norm_counts_XY_resBresC$Type)[1] #% overshooting
BC_vec[1] = table(norm_counts_XY_resBresC$Type)[2] #% reversion
BC_vec[3] = table(norm_counts_XY_resBresC$Type)[3] #% reinforcement
BC_vec = data.frame(val = BC_vec)
BC_vec = cbind(BC_vec, type = c("REV", "OVER", "RI"))
BC_vec$type= factor(BC_vec$type, levels = BC_vec$type)

binom.test(x = c((BC_Table[1]+BC_Table[3]), BC_Table[2]), p = (KL_Table[1]+KL_Table[3])/(KL_Table[1]+KL_Table[2]+KL_Table[3]))
binom.test(x = c((BC_Table[3]), BC_Table[2]+BC_Table[1]), p = (KL_Table[3])/(KL_Table[1]+KL_Table[2]+KL_Table[3]))
binom.test(x = c((BC_Table[3]), BC_Table[2]+BC_Table[1]), p = (T_table[3])/(T_table[1]+T_table[2]+T_table[3]))

###################################################################
### 13 E.2 resBnotC, resCnotB: are there differences?
###################################################################

#13.E.2.1 : WHole transcriptome vs. BnotC, CnotB

resB_notC = resB_sig[resB_sig$Names %notin% resBresC_sig_overlap_samedir$Names,]
dim(resB_notC) #1785
resC_notB = resC_sig[resC_sig$Names %notin% resBresC_sig_overlap_samedir$Names,]
dim(resC_notB) #4035

norm_counts_cut_G1_resBnotC = norm_counts_cut_XY_G1[rownames(norm_counts_cut_XY_G1) %in% resB_notC$Names,]
norm_counts_cut_G2_resCnotB = norm_counts_cut_XY_G2[rownames(norm_counts_cut_XY_G2) %in% resC_notB$Names,]

#i) what proportion of these genes don't meet the EC/PC cutoffs?
dim(norm_counts_cut_G1_resBnotC) #1231
dim(norm_counts_cut_G1_resBnotC)[1]/nrow(resB_notC) #So 68.9% pass both thresholds...
#Update: only 37% of these genes pass the threshold, why?
#I guess it's again because the evidence is more stringent, this is what we were seeing before.

dim(norm_counts_cut_G2_resCnotB) #2630 
dim(norm_counts_cut_G2_resCnotB)/nrow(resC_notB) #So 65.2% pass both thresholds...
#So again now a much smaller proportion of these pass the threshold.

#ii) What proportion show substantial ancestral plasticity?
#For this we will need to go back and look at the pre-cutoff files...
norm_counts_G1_resBnotC = norm_counts[rownames(norm_counts) %in% resB_notC$Names,]
norm_counts_G2_resCnotB = norm_counts[rownames(norm_counts) %in% resC_notB$Names,]
dim(norm_counts_G1_resBnotC) #1785
dim(norm_counts_G2_resCnotB) #4035

#So for G1...
table(norm_counts_G1_resBnotC$Cutoff_GC_G1) #So 1456/1785, or...
table(norm_counts_G1_resBnotC$Cutoff_GC_G1)[2]/nrow(resB_notC) #...81.5% show substantial GC, and...
table(norm_counts_G1_resBnotC$Cutoff_PC_G1) #So 1456/1785, or...
table(norm_counts_G1_resBnotC$Cutoff_PC_G1)[2]/nrow(resB_notC) #...81.5% show substantial PC? But they're not the same genes I think. Weird.
#So in terms of the raw numbers (PC vs. GC still calculated in the same way) they're the same...
#but a much smaller subset are now being included because of the way we do the cutoffs.

table(norm_counts_G2_resCnotB$Cutoff_GC_G2) #So 3034/4035, or...
table(norm_counts_G2_resCnotB$Cutoff_GC_G2)[2]/nrow(resC_notB) #...75.2% show subantial GC (i.e. derived genetic change)
table(norm_counts_G2_resCnotB$Cutoff_PC_G2) #So 3461/4035, or...
table(norm_counts_G2_resCnotB$Cutoff_PC_G2)[2]/nrow(resC_notB) #...85.8% show substantial PC (i.e. ancestral plastic change)
#So in terms of the raw numbers (PC vs. GC still calculated in the same way) they're the same...
#but a much smaller subset are now being included because of the way we do the cutoffs.

#Ok so I think the real question is...for genes in resBnotC [G1 only}...do they show substantially lower anc. plast. in G2? I.e. is that why they're not DGE?
table(norm_counts_G1_resBnotC$Cutoff_PC_G2)[2]/nrow(resB_notC) #SO 81.4% of these are ancestrally plastic.
table(norm_counts_G2_resCnotB$Cutoff_PC_G1)[2]/nrow(resC_notB) #SO 85.2% of these are ancestrally plastic.

#I guess it could be the magnitude of plasticity is different? But maybe we should stop here. Does any of this add up to much?

#So there isn't really strong evidence here that a lack of ancestral plasticity is driving recruitment.

norm_counts_cut_G1_resBresC = norm_counts_cut_XY_G1[rownames(norm_counts_cut_XY_G1) %in% resBresC_sig_overlap_samedir$Names,]
norm_counts_cut_G2_resBresC = norm_counts_cut_XY_G2[rownames(norm_counts_cut_XY_G2) %in% resBresC_sig_overlap_samedir$Names,]
dim(norm_counts_cut_G1_resBresC)[1]/nrow(resBresC_sig_overlap_samedir)#291:
dim(norm_counts_cut_G2_resBresC)[1]/nrow(resBresC_sig_overlap_samedir) #300:
#So similar across both pairs...
#So you're still losing most genes here...apparently most of the reinforced ones.

dim(norm_counts_cut_G1_resBnotC)[1]/nrow(resB_notC)
dim(norm_counts_cut_G2_resCnotB)[1]/nrow(resC_notB)

#SO probably some pertinent questions are:
#Are the shared genes showing similar patterns in the individual ones.

#KR33: More reversion in one sample "CEC" genes 
100*table(norm_counts_cut_G1_resBnotC$Type)[1]/nrow(norm_counts_cut_G1_resBnotC) #% overshooting: 6.4% 
100*table(norm_counts_cut_G1_resBnotC$Type)[2]/nrow(norm_counts_cut_G1_resBnotC) #% reversion: 77.5%
100*table(norm_counts_cut_G1_resBnotC$Type)[3]/nrow(norm_counts_cut_G1_resBnotC) #% reinforcement; 16%
table(norm_counts_cut_G1_resBnotC$Type)[3] #So there are 197 reinforced genes in here...

100*table(norm_counts_cut_G2_resCnotB$Type)[1]/nrow(norm_counts_cut_G2_resCnotB) #% overshooting 7.14%
100*table(norm_counts_cut_G2_resCnotB$Type)[2]/nrow(norm_counts_cut_G2_resCnotB) #% reversion 81.8%
100*table(norm_counts_cut_G2_resCnotB$Type)[3]/nrow(norm_counts_cut_G2_resCnotB) #% reinforcement; 11.0%
table(norm_counts_cut_G2_resCnotB$Type)[3] #And 290 reinforced genes in here.

100*table(norm_counts_cut_G1_resBresC$Type)[1]/nrow(norm_counts_cut_G1_resBresC) #% overshooting: 4.46% 
100*table(norm_counts_cut_G1_resBresC$Type)[2]/nrow(norm_counts_cut_G1_resBresC) #% reversion: 60.1%
100*table(norm_counts_cut_G1_resBresC$Type)[3]/nrow(norm_counts_cut_G1_resBresC) #% reinforcement: 35.4%

100*table(norm_counts_cut_G2_resBresC$Type)[1]/nrow(norm_counts_cut_G2_resBresC) #% overshooting: 1.7%
100*table(norm_counts_cut_G2_resBresC$Type)[2]/nrow(norm_counts_cut_G2_resBresC) #% reversion: 67%
100*table(norm_counts_cut_G2_resBresC$Type)[3]/nrow(norm_counts_cut_G2_resBresC) #% reinforcement: 31.3%

#So a bit of a problem here somewhere...
#So to spell it out...
#For the combined set of genes (XY), you get most passing that threshold as before.
#Considering each population separately, there's a smaller number of genes passing both thresholds.
#That's because the power you have decreases. Only ~37-40% of genes in BnotC, CnotB and BandC, are able to be included.
#I guess an issue is that this might be a slightly skewed set of genes.

lfcA_BC_REV = lfcA[lfcA$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "REV",]),]
lfcE_BC_REV = lfcE[lfcE$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "REV",]),]
lfcD_BC_REV = lfcD[lfcD$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "REV",]),]
lfcH_BC_REV = lfcH[lfcH$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "REV",]),]
resids_countsA_BC_REV = abs(lfcA_BC_REV$log2FoldChange)
resids_countsE_BC_REV = abs(lfcE_BC_REV$log2FoldChange)
resids_countsD_BC_REV = abs(lfcD_BC_REV$log2FoldChange)
resids_countsH_BC_REV = abs(lfcH_BC_REV$log2FoldChange)

lfcA_BC_RI = lfcA[lfcA$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "RI",]),]
lfcE_BC_RI = lfcE[lfcE$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "RI",]),]
lfcD_BC_RI = lfcD[lfcD$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "RI",]),]
lfcH_BC_RI = lfcH[lfcH$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "RI",]),]
resids_countsA_BC_RI = abs(lfcA_BC_RI$log2FoldChange)
resids_countsE_BC_RI = abs(lfcE_BC_RI$log2FoldChange)
resids_countsD_BC_RI = abs(lfcD_BC_RI$log2FoldChange)
resids_countsH_BC_RI = abs(lfcH_BC_RI$log2FoldChange)

lfcA_BC_OVERSHOOT = lfcA[lfcA$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "OVERSHOOT",]),]
lfcE_BC_OVERSHOOT = lfcE[lfcE$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "OVERSHOOT",]),]
lfcD_BC_OVERSHOOT = lfcD[lfcD$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "OVERSHOOT",]),]
lfcH_BC_OVERSHOOT = lfcH[lfcH$Names %in% rownames(norm_counts_XY_resBresC[norm_counts_XY_resBresC$Type == "OVERSHOOT",]),]
resids_countsA_BC_OVERSHOOT = abs(lfcA_BC_OVERSHOOT$log2FoldChange)
resids_countsE_BC_OVERSHOOT = abs(lfcE_BC_OVERSHOOT$log2FoldChange)
resids_countsD_BC_OVERSHOOT = abs(lfcD_BC_OVERSHOOT$log2FoldChange)
resids_countsH_BC_OVERSHOOT = abs(lfcH_BC_OVERSHOOT$log2FoldChange)

norm_counts_noPC_XY_resBresC = norm_counts_noPC_XY[rownames(norm_counts_noPC_XY) %in% resBresC_sig_overlap_samedir$Names,]
dim(norm_counts_noPC_XY_resBresC)# So Out of 371, 

lfcA_BC_NOANCPLAST = lfcA[lfcA$Names %in% rownames(norm_counts_noPC_XY_resBresC),]
lfcE_BC_NOANCPLAST = lfcE[lfcE$Names %in% rownames(norm_counts_noPC_XY_resBresC),]
lfcD_BC_NOANCPLAST = lfcD[lfcD$Names %in% rownames(norm_counts_noPC_XY_resBresC),]
lfcH_BC_NOANCPLAST = lfcH[lfcH$Names %in% rownames(norm_counts_noPC_XY_resBresC),]

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
quantile(resids_countsH_BC_NOANCPLAST$VALUE)
#0.12. So the same as all the rest of them, then. 

resids_countsH_BC_RI = as.data.frame(resids_countsH_BC_RI)
colnames(resids_countsH_BC_RI) = "VALUE"
resids_countsH_BC_RI$Names = "Reinforcement"
resids_countsH_BC_REV = as.data.frame(resids_countsH_BC_REV)
colnames(resids_countsH_BC_REV) = "VALUE"
resids_countsH_BC_REV$Names = "REV"
resids_countsH_BC_OVERSHOOT = as.data.frame(resids_countsH_BC_OVERSHOOT)
colnames(resids_countsH_BC_OVERSHOOT) = "VALUE"
resids_countsH_BC_OVERSHOOT$Names = "OVER"

Plast_FC_BC = rbind(resids_countsH_BC_OVERSHOOT, resids_countsH_BC_REV, resids_countsH_BC_RI)

#KR37 (1/2)
Plast_FC_BC_wilcox = wilcox.test(resids_countsH_BC_NOANCPLAST$VALUE, Plast_FC_BC$VALUE, alternative = "two.sided")
Plast_FC_BC_wilcox

temp_Plast_FC_BC = Plast_FC_BC
temp_Plast_FC_BC$Names = "PLAST"

PlastNoPlast_FC_BC = rbind(resids_countsH_BC_NOANCPLAST, temp_Plast_FC_BC)
median(temp_Plast_FC_BC$VALUE)
median(resids_countsH_BC_NOANCPLAST$VALUE)


PlastNoPlast_FC_KL_wilcox = wilcox.test(temp_Plast_FC_KL$VALUE, resids_countsH_KL_NOANCPLAST$VALUE, alternative = "two.sided")

#KR36 (1/2)
PlastNoPlast_FC_BC_wilcox = wilcox.test(temp_Plast_FC_BC$VALUE, resids_countsH_BC_NOANCPLAST$VALUE, p.adjust.method = "BH")
PlastNoPlast_FC_BC_wilcox #So there isn't any difference. Fine.

################################################
#### |FC| values for different classes - AP, CEC, PLAST, NOPLAST, etc.
###############################################

T_vec$title = "Transcriptome-wide"
T_vec$type = c("Reversion", "Overshooting", "Reinforcement")
T_vec$type = factor(T_vec$type, levels = c("Reinforcement", "Overshooting", "Reversion"), ordered = T)

#KR20 (1/2)
T_barplot = ggplot(data = T_vec, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity', color = "black")+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(axis.text.y = element_text(size = 25),
                               axis.text.x = element_text(size = 25), axis.title.x = element_text(size=25))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
        strip.text = element_text(size = 30))+facet_grid(.~title)+coord_flip()+scale_fill_manual(values = c(color_list[1], color_list[2], color_list[3]))

T_barplot

#KR30 (1/2) - Barpolot for CEC genes in OVER, REV, RI
BC_vec$title = "CEC genes"
BC_vec$type = c("Reversion", "Overshooting", "Reinforcement")
BC_vec$type = factor(BC_vec$type, levels = c("Reinforcement", "Overshooting", "Reversion"), ordered = T)
BC_barplot = ggplot(data = BC_vec, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity', color = "black")+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(axis.text.y = element_text(size = 25),
                               axis.text.x = element_text(size = 25), axis.title.x = element_text(size=25))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
        strip.text = element_text(size = 30))+facet_grid(.~title)+coord_flip()+scale_fill_manual(values = color_list)

BC_barplot



KL_vec$title = "AP genes"
KL_vec$type = c("Reversion", "Overshooting", "Reinforcement")
KL_vec$type = factor(KL_vec$type, levels = c("Reinforcement", "Overshooting", "Reversion"), ordered = T)

#KR26 (1/2) - Fig3B
KL_barplot = ggplot(data = KL_vec, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity', color = "black")+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(axis.text.y = element_text(size = 25),
                               axis.text.x = element_text(size = 25), axis.title.x = element_text(size=25))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
        strip.text = element_text(size = 30))+facet_grid(.~title)+coord_flip()+
  scale_fill_manual(values = color_list)

KL_barplot

#: Comparable levels of |FC| in genes lacking vs. having ancestral plasticity in CEC genes

#KR37 (2/2)
BC_PlastNoPlast_Boxplot = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_BC,fill = another_plot_purple, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size =20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("NP", "P"))+coord_cartesian(ylim = c(0,1.1))+
  #facet_grid(.~Set)+
  xlab("")+ylab("Absolute log2 fold change (|FC|)")+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = c("solid", "solid")), 
        text = element_text(size = 30),
        strip.text = element_text(size = 30))
BC_PlastNoPlast_Boxplot

#KR38 (2/2): Comparable levels of |FC| in zinc in genes lacking vs having ancestral plasticity in DP(0) genes. 

KL_PlastNoPlast_Boxplot = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_KL,fill = another_plot_purple, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size =20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("NP", "P"))+
  #facet_grid(.~Set)+
  xlab("")+ylab("Absolute log2 fold change (|FC|)")+coord_cartesian(ylim=c(0,0.7))+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = c("solid", "solid")), 
        text = element_text(size = 30),
        strip.text = element_text(size = 30))
KL_PlastNoPlast_Boxplot

median(resids_countsH_BC_RI$VALUE)
median(resids_countsH_BC_REV$VALUE)
median(resids_countsH_BC_OVERSHOOT$VALUE)
median(resids_countsE_BC_RI)
median(resids_countsE_BC_REV)
median(resids_countsE_BC_OVERSHOOT)


###########################################
### Plotting the green/yellow diagrams
############################################

#EC/GC plots

ncc_lim = max(c(quantile(norm_counts_cut_XY$PC)[2], quantile(norm_counts_cut_XY$PC)[4], quantile(norm_counts_cut_XY$GC)[2], quantile(norm_counts_cut_XY$GC)[4]))
ap_lim = max(c(quantile(norm_counts_resKresLDPN_XY$PC)[2], quantile(norm_counts_resKresLDPN_XY$PC)[4], quantile(norm_counts_resKresLDPN_XY$GC)[2], quantile(norm_counts_resKresLDPN_XY$GC)[4]))
cec_lim = max(c(quantile(norm_counts_XY_resBresC$PC)[2], quantile(norm_counts_XY_resBresC$PC)[4], quantile(norm_counts_XY_resBresC$GC)[2], quantile(norm_counts_XY_resBresC$GC)[4]))

biggest_lim = max(c(ncc_lim, ap_lim, cec_lim))
biggest_lim


color_list = c("#75EAA9","#B6EA50", "#F5FF74")
color_list
over1 = data.frame(x = c(0, -Inf, -Inf), y = c(0,0,0.5*biggest_lim*1.1))
over2 = data.frame(x = c(0, Inf, Inf), y = c(0,0,-0.5*biggest_lim*1.1))
over3 = data.frame(x = c(0, -Inf, -Inf, 0), y = c(0, 0.5*biggest_lim*1.1, Inf, Inf)) 
over4 = data.frame(x = c(0, Inf, Inf, 0), y = c(0, -0.5*biggest_lim*1.1, -Inf, -Inf)) 

over1 = data.frame(x = c(0, -Inf, -Inf), y = c(0,0,0.5*biggest_lim*1.1))
over2 = data.frame(x = c(0, Inf, Inf), y = c(0,0,-0.5*biggest_lim*1.1))
over3 = data.frame(x = c(0, -Inf, -Inf, 0), y = c(0, 0.5*biggest_lim*1.1, Inf, Inf)) 
over4 = data.frame(x = c(0, Inf, Inf, 0), y = c(0, -0.5*biggest_lim*1.1, -Inf, -Inf)) 

#KR21
Whole_ECGC = ggplot(data = data.frame())+geom_point()+
  #geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  #geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  #geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  #geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_cut_XY, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
Whole_ECGC

#KR27
WholeAP_ECGC = ggplot(data = data.frame())+geom_point()+
  #geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  #geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  #geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  #geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_resKresLDPN_XY, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
WholeAP_ECGC

#KR31: Fig 5C
WholeCEC_ECGC = ggplot(data = data.frame())+geom_point()+
  #geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  #geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  #geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  #geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_XY_resBresC, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
WholeCEC_ECGC

#KR20 (2/2)
Middle_ECGC = ggplot(data = data.frame())+geom_point()+xlim(-biggest_lim, biggest_lim)+ylim(-biggest_lim, biggest_lim)+
  geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_cut_XY, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "white")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "white")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
Middle_ECGC

#KR26 (2/2) - Figure 3B
AP_ECGC = ggplot(data = data.frame())+geom_point()+xlim(-biggest_lim, biggest_lim)+ylim(-biggest_lim, biggest_lim)+
  geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+geom_abline(slope = -0.5, intercept = 0)+
  geom_abline(slope = -1, intercept = 0, lty = 2)+geom_bin2d(data = norm_counts_resKresLDPN_XY, aes(PC, GC), bins=50)+
  geom_hline(yintercept = -biggest_lim*1.0)+geom_hline(yintercept = biggest_lim)+
  geom_vline(xintercept = -biggest_lim)+geom_vline(xintercept = biggest_lim)+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
AP_ECGC

# KR30 (2/2)- CEC genes green/yellow plot showing EC vs. P
CEC_ECGC = ggplot(data = data.frame())+geom_point()+xlim(-biggest_lim, biggest_lim)+ylim(-biggest_lim, biggest_lim)+
  geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_XY_resBresC, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 22))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
CEC_ECGC

##########################################################################################
##Parametric Bootstrapping (a la Ho and Zhang 2019)
##########################################################################################
#KR

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
#write.csv(classified_vec, file = "final_param_classifications.csv")


perm_table = table(norm_counts_cut$temp_perm_calc)
classifiable_perm = perm_table[1]+perm_table[2]+perm_table[3]
perm_table[1]/classifiable_perm
perm_table[2]/classifiable_perm
perm_table[3]/classifiable_perm

perm_table[1]
perm_table[2]
perm_table[3]

norm_counts_cut_param = norm_counts_cut[norm_counts_cut$temp_perm_calc != "UNCLASSIFIABLE",]
dim(norm_counts_cut_param)

norm_counts_resKresL_param = norm_counts_cut_param[rownames(norm_counts_cut_param) %in% resKresL_sig_overlap_samedir$Names,]
dim(norm_counts_resKresL_param)

table(norm_counts_resKresL_param$Type)[1]/nrow(norm_counts_resKresL_param) #% overshooting
table(norm_counts_resKresL_param$Type)[2]/nrow(norm_counts_resKresL_param) #% reversion
table(norm_counts_resKresL_param$Type)[3]/nrow(norm_counts_resKresL_param) #% reinforcement

table(norm_counts_resKresL_param$Type)[1]
table(norm_counts_resKresL_param$Type)[2]
table(norm_counts_resKresL_param$Type)[3]

norm_counts_resBresC_param = norm_counts_cut_param[rownames(norm_counts_cut_param) %in% resBresC_sig_overlap_samedir$Names,]

table(norm_counts_resBresC_param$Type)[1]/nrow(norm_counts_resBresC_param) #% overshooting
table(norm_counts_resBresC_param$Type)[2]/nrow(norm_counts_resBresC_param) #% reversion
table(norm_counts_resBresC_param$Type)[3]/nrow(norm_counts_resBresC_param) #% reinforcement

#############


###############
# Old method of considering if genes are substantially DE

#1) Defining
norm_counts_cut_0.2 = norm_counts[norm_counts$Cutoff_GC == "Y" & norm_counts$Cutoff_PC == "Y",]
#KR38
dim(norm_counts_cut_0.2)[1]/dim(norm_counts)[1] #21,150
dim(norm_counts_cut_XY[rownames(norm_counts_cut_XY) %in% rownames(norm_counts_cut_0.2),])[1]/dim(norm_counts_cut_XY)[1]
length(unique(c(rownames(norm_counts_cut_XY), c(rownames(norm_counts_cut_0.2)))))/dim(norm_counts)[1]

norm_counts_noPC_0.2 = norm_counts[norm_counts$Cutoff_GC == "Y" & norm_counts$Cutoff_PC == "N",]
norm_counts_noGC_0.2 = norm_counts[norm_counts$Cutoff_GC == "N" & norm_counts$Cutoff_PC == "Y",]

norm_counts_cut_0.2$DirStat = norm_counts_cut_0.2$GC*norm_counts_cut_0.2$PC
norm_counts_cut_0.2$SameDir = ifelse(norm_counts_cut_0.2$DirStat > 0, "Y", "N")
norm_counts_cut_0.2$logPC = norm_counts_cut_0.2$PC*log2(abs(norm_counts_cut_0.2$PC))/abs(norm_counts_cut_0.2$PC)
head(norm_counts_cut_0.2$logPC)
head(norm_counts_cut_0.2$PC)
norm_counts_cut_0.2$logGC = norm_counts_cut_0.2$GC*log2(abs(norm_counts_cut_0.2$GC))/abs(norm_counts_cut_0.2$GC)

dim(norm_counts_cut_0.2)[1]/dim(norm_counts)[1] #...well...this is now 73% of all genes. OK. 
norm_counts_cut_0.2$Type = ifelse(norm_counts_cut_0.2$DirStat > 0, "RI", ifelse(abs(norm_counts_cut_0.2$GC) > 0.5*abs(norm_counts_cut_0.2$PC), "REV", "OVERSHOOT"))

#Looking at each population separately

#2) Looking at each population individually
norm_counts_cut_0.2_G1 = norm_counts[norm_counts$Cutoff_GC_G1 == "Y" & norm_counts$Cutoff_PC_G1 == "Y",]
norm_counts_cut_0.2_G2 = norm_counts[norm_counts$Cutoff_GC_G2 == "Y" & norm_counts$Cutoff_PC_G2 == "Y",]

norm_counts_cut_0.2_G1$DirStat = norm_counts_cut_0.2_G1$GC_G1*norm_counts_cut_0.2_G1$PC_G1
norm_counts_cut_0.2_G2$DirStat = norm_counts_cut_0.2_G2$GC_G2*norm_counts_cut_0.2_G2$PC_G2

norm_counts_cut_0.2_G1$Type = ifelse(norm_counts_cut_0.2_G1$DirStat > 0, "RI", ifelse(abs(norm_counts_cut_0.2_G1$GC_G1) > 0.5*abs(norm_counts_cut_0.2_G1$PC_G1), "REV", "OVERSHOOT"))
norm_counts_cut_0.2_G2$Type = ifelse(norm_counts_cut_0.2_G2$DirStat > 0, "RI", ifelse(abs(norm_counts_cut_0.2_G2$GC_G2) > 0.5*abs(norm_counts_cut_0.2_G2$PC_G2), "REV", "OVERSHOOT"))

100*table(norm_counts_cut_0.2_G1$Type)[1]/nrow(norm_counts_cut_0.2_G1) #% overshooting
100*table(norm_counts_cut_0.2_G1$Type)[2]/nrow(norm_counts_cut_0.2_G1) #% reversion
100*table(norm_counts_cut_0.2_G1$Type)[3]/nrow(norm_counts_cut_0.2_G1) #% reinforcement

100*table(norm_counts_cut_0.2_G2$Type)[1]/nrow(norm_counts_cut_0.2_G2) #% overshooting
100*table(norm_counts_cut_0.2_G2$Type)[2]/nrow(norm_counts_cut_0.2_G2) #% reversion
100*table(norm_counts_cut_0.2_G2$Type)[3]/nrow(norm_counts_cut_0.2_G2) #% reinforcement

100*dim(norm_counts_cut_0.2_G1)[1]/dim(norm_counts)[1]
100*dim(norm_counts_cut_0.2_G2)[1]/dim(norm_counts)[1]

#3) Combined transcriptome
100*dim(norm_counts_cut_0.2)[1]/dim(norm_counts)[1] #So actually...super low? 

#KR40
100*table(norm_counts_cut_0.2$Type)[1]/nrow(norm_counts_cut_0.2) #% overshooting
100*table(norm_counts_cut_0.2$Type)[2]/nrow(norm_counts_cut_0.2) #% reversion
100*table(norm_counts_cut_0.2$Type)[3]/nrow(norm_counts_cut_0.2) #% reinforcement

#Yup ok so....

table(norm_counts_cut_0.2$Type)[1]
table(norm_counts_cut_0.2$Type)[2]
table(norm_counts_cut_0.2$Type)[3]

head(norm_counts_cut_0.2)
T_table_0.2 = table(norm_counts_cut_0.2$Type)
T_vec_0.2 = c(0,0,0)
T_vec_0.2[2] = table(norm_counts_cut_0.2$Type)[1] #% overshooting
T_vec_0.2[1] = table(norm_counts_cut_0.2$Type)[2] #% reversion
T_vec_0.2[3] = table(norm_counts_cut_0.2$Type)[3] #% reinforcement
T_vec_0.2 = data.frame(val = T_vec_0.2)
T_vec_0.2 = cbind(T_vec_0.2, type = c("REV", "OVER", "RI"))
T_vec_0.2$type= factor(T_vec_0.2$type, levels = T_vec_0.2$type)

#|FC| values for whole transcriptome
lfcA_T_REV_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "REV",]),]
lfcE_T_REV_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "REV",]),]
lfcD_T_REV_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "REV",]),]
lfcH_T_REV_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "REV",]),]
resids_countsA_T_REV_0.2 = abs(lfcA_T_REV_0.2$log2FoldChange)
resids_countsE_T_REV_0.2 = abs(lfcE_T_REV_0.2$log2FoldChange)
resids_countsD_T_REV_0.2 = abs(lfcD_T_REV_0.2$log2FoldChange)
resids_countsH_T_REV_0.2 = abs(lfcH_T_REV_0.2$log2FoldChange)

lfcA_T_RI_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "RI",]),]
lfcE_T_RI_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "RI",]),]
lfcD_T_RI_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "RI",]),]
lfcH_T_RI_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "RI",]),]
resids_countsA_T_RI_0.2 = abs(lfcA_T_RI_0.2$log2FoldChange)
resids_countsE_T_RI_0.2 = abs(lfcE_T_RI_0.2$log2FoldChange)
resids_countsD_T_RI_0.2 = abs(lfcD_T_RI_0.2$log2FoldChange)
resids_countsH_T_RI_0.2 = abs(lfcH_T_RI_0.2$log2FoldChange)

lfcA_T_OVERSHOOT_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "OVERSHOOT",]),]
lfcE_T_OVERSHOOT_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "OVERSHOOT",]),]
lfcD_T_OVERSHOOT_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "OVERSHOOT",]),]
lfcH_T_OVERSHOOT_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_cut_0.2[norm_counts_cut_0.2$Type == "OVERSHOOT",]),]
resids_countsA_T_OVERSHOOT_0.2 = abs(lfcA_T_OVERSHOOT_0.2$log2FoldChange)
resids_countsE_T_OVERSHOOT_0.2 = abs(lfcE_T_OVERSHOOT_0.2$log2FoldChange)
resids_countsD_T_OVERSHOOT_0.2 = abs(lfcD_T_OVERSHOOT_0.2$log2FoldChange)
resids_countsH_T_OVERSHOOT_0.2 = abs(lfcH_T_OVERSHOOT_0.2$log2FoldChange)

resids_countsH_T_RI_0.2 = as.data.frame(resids_countsH_T_RI_0.2)
colnames(resids_countsH_T_RI_0.2) = "VALUE"
resids_countsH_T_RI_0.2$Names = "RI"
resids_countsH_T_REV_0.2 = as.data.frame(resids_countsH_T_REV_0.2)
colnames(resids_countsH_T_REV_0.2) = "VALUE"
resids_countsH_T_REV_0.2$Names = "REV"
resids_countsH_T_OVERSHOOT_0.2 = as.data.frame(resids_countsH_T_OVERSHOOT_0.2)
colnames(resids_countsH_T_OVERSHOOT_0.2) = "VALUE"
resids_countsH_T_OVERSHOOT_0.2$Names = "OVERSHOOT"

median(resids_countsH_T_RI_0.2$VALUE) #Reinforcement |FC|
median(resids_countsH_T_REV_0.2$VALUE) #Reversion |FC|
median(resids_countsH_T_OVERSHOOT_0.2$VALUE) #Overshooting |FC|


Plast_FC_T_0.2 = rbind(resids_countsH_T_OVERSHOOT_0.2, resids_countsH_T_REV_0.2, resids_countsH_T_RI_0.2)

median(resids_countsH_T_RI_0.2$VALUE)
median(resids_countsH_T_REV_0.2$VALUE)
median(resids_countsH_T_OVERSHOOT_0.2$VALUE)

median(resids_countsE_T_RI_0.2)
median(resids_countsE_T_REV_0.2)
median(resids_countsE_T_OVERSHOOT_0.2)

Plast_FC_T_wilcox_0.2 = pairwise.wilcox.test(Plast_FC_T_0.2$VALUE, Plast_FC_T_0.2$Names, p.adjust.method = "BH")
Plast_FC_T_wilcox_0.2
pairwiseMedianTest(value ~ variable, data=both_CZ_0.2) #Might be better to do this, really. 

#4) reskresLDPN genes: combined

dim(resKresL_DPN_sig_overlap_samedir)
norm_counts_resKresLDPN_0.2 = norm_counts_cut_0.2[rownames(norm_counts_cut_0.2) %in% resKresL_DPN_sig_overlap_samedir$Names,]
norm_counts_resKresLDPN_0.2 = norm_counts_cut_0.2[rownames(norm_counts_cut_0.2) %in% resKresL_DPN_sig_overlap_samedir$Names,]

#KR43: 39% of DP(0) genes have "substantial" |PC| and |EC|
dim(norm_counts_resKresLDPN_0.2)[1]/dim(resKresL_DPN_sig_overlap_samedir)[1]
#KR44:
dim(norm_counts_resKresLDPN_0.2[rownames(norm_counts_resKresLDPN_0.2) %in% rownames(norm_counts_resKresLDPN_XY),])[1]/dim(norm_counts_resKresLDPN_XY)[1]
length(unique(rownames(norm_counts_resKresLDPN_0.2), rownames(norm_counts_resKresLDPN_XY)))/dim(resKresL_DPN_sig_overlap_samedir)[1]

#Whereas with the old method we do say 87%. 

#KR45 (1/3) - DP(0) genes OVER=12.8%, REV=70.8%, RI=16.4%. 
100*table(norm_counts_resKresLDPN_0.2$Type)[1]/nrow(norm_counts_resKresLDPN_0.2) #% overshooting- 14%
100*table(norm_counts_resKresLDPN_0.2$Type)[2]/nrow(norm_counts_resKresLDPN_0.2) #% reversion - 66%
100*table(norm_counts_resKresLDPN_0.2$Type)[3]/nrow(norm_counts_resKresLDPN_0.2) #% reinforcement - 18.7% 

table(norm_counts_resKresLDPN_0.2$Type)[2]
table(norm_counts_resKresLDPN_0.2$Type)[1]
table(norm_counts_resKresLDPN_0.2$Type)[3]

KL_Table_0.2 = table(norm_counts_resKresLDPN_0.2$Type)
KL_vec_0.2 = c(0,0,0)
KL_vec_0.2[2] = table(norm_counts_resKresLDPN_0.2$Type)[1] #% overshooting
KL_vec_0.2[1] = table(norm_counts_resKresLDPN_0.2$Type)[2] #% reversion
KL_vec_0.2[3] = table(norm_counts_resKresLDPN_0.2$Type)[3] #% reinforcement
KL_vec_0.2 = data.frame(val = KL_vec_0.2)
KL_vec_0.2 = cbind(KL_vec_0.2, type = c("REV", "OVER", "RI"))
KL_vec_0.2$type= factor(KL_vec_0.2$type, levels = KL_vec_0.2$type)

woof_0.2 = binom.test(x = c((KL_Table_0.2[1]+KL_Table_0.2[3]), KL_Table_0.2[2]), p = (T_table_0.2[1]+T_table_0.2[3])/(T_table_0.2[1]+T_table_0.2[2]+T_table_0.2[3]), alternative = "two.sided")
woof_0.2
#Note - no test statistic here.

#Are "adaptive" plastic changes overrepresented in this gene set?M 

#|FC| values for whole transcriptome
lfcA_KL_REV_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "REV",]),]
lfcE_KL_REV_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "REV",]),]
lfcD_KL_REV_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "REV",]),]
lfcH_KL_REV_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "REV",]),]
resids_countsA_KL_REV_0.2 = abs(lfcA_KL_REV_0.2$log2FoldChange)
resids_countsE_KL_REV_0.2 = abs(lfcE_KL_REV_0.2$log2FoldChange)
resids_countsD_KL_REV_0.2 = abs(lfcD_KL_REV_0.2$log2FoldChange)
resids_countsH_KL_REV_0.2 = abs(lfcH_KL_REV_0.2$log2FoldChange)

lfcA_KL_RI_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "RI",]),]
lfcE_KL_RI_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "RI",]),]
lfcD_KL_RI_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "RI",]),]
lfcH_KL_RI_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "RI",]),]
resids_countsA_KL_RI_0.2 = abs(lfcA_KL_RI_0.2$log2FoldChange)
resids_countsE_KL_RI_0.2 = abs(lfcE_KL_RI_0.2$log2FoldChange)
resids_countsD_KL_RI_0.2 = abs(lfcD_KL_RI_0.2$log2FoldChange)
resids_countsH_KL_RI_0.2 = abs(lfcH_KL_RI_0.2$log2FoldChange)

lfcA_KL_OVERSHOOT_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "OVERSHOOT",]),]
lfcE_KL_OVERSHOOT_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "OVERSHOOT",]),]
lfcD_KL_OVERSHOOT_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "OVERSHOOT",]),]
lfcH_KL_OVERSHOOT_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_resKresLDPN_0.2[norm_counts_resKresLDPN_0.2$Type == "OVERSHOOT",]),]
resids_countsA_KL_OVERSHOOT_0.2 = abs(lfcA_KL_OVERSHOOT_0.2$log2FoldChange)
resids_countsE_KL_OVERSHOOT_0.2 = abs(lfcE_KL_OVERSHOOT_0.2$log2FoldChange)
resids_countsD_KL_OVERSHOOT_0.2 = abs(lfcD_KL_OVERSHOOT_0.2$log2FoldChange)
resids_countsH_KL_OVERSHOOT_0.2 = abs(lfcH_KL_OVERSHOOT_0.2$log2FoldChange)

norm_counts_resKresLDPN_0.2$PC_pc = norm_counts_resKresLDPN_0.2$PC/norm_counts_resKresLDPN_0.2$Lp
norm_counts_resKresLDPN_0.2$GC_pc = norm_counts_resKresLDPN_0.2$GC/norm_counts_resKresLDPN_0.2$Lp

norm_counts_noPC_resKresLDPN_0.2 = norm_counts_noPC_0.2[rownames(norm_counts_noPC_0.2) %in% resKresL_DPN_sig_overlap_samedir$Names,]
dim(norm_counts_noPC_resKresLDPN_0.2) #Note: only 18 genes.
lfcA_KL_NOANCPLAST_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_noPC_resKresLDPN_0.2),]
lfcE_KL_NOANCPLAST_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_noPC_resKresLDPN_0.2),]
lfcD_KL_NOANCPLAST_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_noPC_resKresLDPN_0.2),]
lfcH_KL_NOANCPLAST_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_noPC_resKresLDPN_0.2),]

resids_countsA_KL_NOANCPLAST_0.2 = abs(lfcA_KL_NOANCPLAST_0.2$log2FoldChange)
resids_countsE_KL_NOANCPLAST_0.2 = abs(lfcE_KL_NOANCPLAST_0.2$log2FoldChange)
resids_countsD_KL_NOANCPLAST_0.2 = abs(lfcD_KL_NOANCPLAST_0.2$log2FoldChange)
resids_countsH_KL_NOANCPLAST_0.2 = abs(lfcH_KL_NOANCPLAST_0.2$log2FoldChange)
resids_countsH_KL_NOANCPLAST_0.2 = as.data.frame(resids_countsH_KL_NOANCPLAST_0.2)
colnames(resids_countsH_KL_NOANCPLAST_0.2) = "VALUE"
resids_countsH_KL_NOANCPLAST_0.2$Names = "NOANCPLAST"
dim(resids_countsH_KL_NOANCPLAST_0.2)
median(resids_countsH_KL_NOANCPLAST_0.2$VALUE)
#0.12. So the same as all the rest of them, then. 

resids_countsH_KL_RI_0.2 = as.data.frame(resids_countsH_KL_RI_0.2)
colnames(resids_countsH_KL_RI_0.2) = "VALUE"
resids_countsH_KL_RI_0.2$Names = "RI"
resids_countsH_KL_REV_0.2 = as.data.frame(resids_countsH_KL_REV_0.2)
colnames(resids_countsH_KL_REV_0.2) = "VALUE"
resids_countsH_KL_REV_0.2$Names = "REV"
resids_countsH_KL_OVERSHOOT_0.2 = as.data.frame(resids_countsH_KL_OVERSHOOT_0.2)
colnames(resids_countsH_KL_OVERSHOOT_0.2) = "VALUE"
resids_countsH_KL_OVERSHOOT_0.2$Names = "OVERSHOOT"

Plast_FC_KL_0.2 = rbind(resids_countsH_KL_OVERSHOOT_0.2, resids_countsH_KL_REV_0.2, resids_countsH_KL_RI_0.2)

median(resids_countsH_KL_RI_0.2$VALUE)
median(resids_countsH_KL_REV_0.2$VALUE)
median(resids_countsH_KL_OVERSHOOT_0.2$VALUE)

temp_Plast_FC_KL_0.2 = Plast_FC_KL_0.2
temp_Plast_FC_KL_0.2$Names = "PLAST"
median(temp_Plast_FC_KL_0.2$VALUE)
quantile(temp_Plast_FC_KL_0.2$VALUE)
median(resids_countsH_KL_NOANCPLAST_0.2$VALUE)
quantile(resids_countsH_KL_NOANCPLAST_0.2$VALUE)

PlastNoPlast_FC_KL_0.2 = rbind(temp_Plast_FC_KL_0.2, resids_countsH_KL_NOANCPLAST_0.2)
median(temp_Plast_FC_KL_0.2$VALUE) #Plast
median(resids_countsH_KL_NOANCPLAST_0.2$VALUE) #No Plast... Ahhhh dear, this is not the right way round...

PlastNoPlast_FC_KL_wilcox_0.2 = wilcox.test(temp_Plast_FC_KL_0.2$VALUE, resids_countsH_KL_NOANCPLAST_0.2$VALUE, alternative = "two.sided")
PlastNoPlast_FC_KL_wilcox_0.2 #no significant difference (insufficient numbers I guess...)

KL_PlastNoPlast_Boxplot_0.2 = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_KL_0.2, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("NP", "P"))+
  xlab("")+ylab("|FC|")+ylim(0,0.625)+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 2), text = element_text(size = 15))
KL_PlastNoPlast_Boxplot_0.2
#Should note there aren't enough here.

###################################################################
#0.2_13D.2 AP genes only in one of the two replicates...
###################################################################

#We already have KnotL and LnotK

#Ok so for the individaul things...
norm_counts_cut_0.2_G1_resKDPN_notLDPN = norm_counts_cut_0.2_G1[rownames(norm_counts_cut_0.2_G1) %in% resKDPN_notLDPN$Names,]
norm_counts_cut_0.2_G2_resLDPN_notKDPN = norm_counts_cut_0.2_G2[rownames(norm_counts_cut_0.2_G2) %in% resLDPN_notKDPN$Names,]

#i) what proportion of these genes don't meet the EC/PC cutoffs?
dim(norm_counts_cut_0.2_G1_resKDPN_notLDPN) #405 meet cutoff
dim(norm_counts_cut_0.2_G1_resKDPN_notLDPN)[1]/nrow(resKDPN_notLDPN) #So 83.5% pass both thresholds

dim(norm_counts_cut_0.2_G2_resLDPN_notKDPN) #1564
dim(norm_counts_cut_0.2_G2_resLDPN_notKDPN)/nrow(resLDPN_notKDPN) #So only 66% pass both thresholds.

#ii) What proportion show substantial ancestral plasticity?
#For this we will need to go back and look at the pre-cutoff files...
norm_counts_G1_resKDPN_notLDPN_0.2 = norm_counts[rownames(norm_counts) %in% resKDPN_notLDPN$Names,]
norm_counts_G2_resLDPN_notKDPN_0.2 = norm_counts[rownames(norm_counts) %in% resLDPN_notKDPN$Names,]
dim(norm_counts_G1_resKDPN_notLDPN) #Same as before, good
dim(norm_counts_G2_resLDPN_notKDPN) #Same as before, good

#So for G1...
table(norm_counts_G1_resKDPN_notLDPN_0.2$Cutoff_GC_G1) #So 437/485 pass the GC cutoff
table(norm_counts_G1_resKDPN_notLDPN_0.2$Cutoff_GC_G1)[2]/nrow(resKDPN_notLDPN) #so that's 90%
table(norm_counts_G1_resKDPN_notLDPN_0.2$Cutoff_PC_G1) #So 448/485.
table(norm_counts_G1_resKDPN_notLDPN_0.2$Cutoff_PC_G1)[2]/nrow(resKDPN_notLDPN) #so 92% pass PC as well. 

table(norm_counts_G2_resLDPN_notKDPN_0.2$Cutoff_GC_G2) #So 683/2365
table(norm_counts_G2_resLDPN_notKDPN_0.2$Cutoff_GC_G2)[2]/nrow(resLDPN_notKDPN) #...71.2% show substantial GC (i.e. derived genetic change)
table(norm_counts_G2_resLDPN_notKDPN_0.2$Cutoff_PC_G2) #So 3461/4035, or...
table(norm_counts_G2_resLDPN_notKDPN_0.2$Cutoff_PC_G2)[2]/nrow(resLDPN_notKDPN) #but still, 92.6% show substantail PC

#Ok so I think the real question is...for genes in resKDPN_notLDPN [G1 only}...do they show substantially lower anc. plast. in G2? I.e. is that why they're not DGE?
table(norm_counts_G1_resKDPN_notLDPN_0.2$Cutoff_PC_G2)[2]/nrow(resKDPN_notLDPN) #Nope...
table(norm_counts_G2_resLDPN_notKDPN_0.2$Cutoff_PC_G1)[2]/nrow(resLDPN_notKDPN) #Nope...

#So there isn't really strong evidence here that a lack of ancestral plasticity is driving recruitment.

norm_counts_cut_0.2_G1_resKresLDPN = norm_counts_cut_0.2_G1[rownames(norm_counts_cut_0.2_G1) %in% resKresL_DPN_sig_overlap_samedir$Names,]
norm_counts_cut_0.2_G2_resKresLDPN = norm_counts_cut_0.2_G2[rownames(norm_counts_cut_0.2_G2) %in% resKresL_DPN_sig_overlap_samedir$Names,]
dim(norm_counts_cut_0.2_G1) #very similar
dim(norm_counts_cut_0.2_G2) #to each other

#So 112 in both. I guess..conceivable?
dim(norm_counts_cut_0.2_G1_resKresLDPN)[1]
dim(norm_counts_cut_0.2_G2_resKresLDPN)[1]

dim(norm_counts_cut_0.2_G1_resKresLDPN)[1]/nrow(resKresL_DPN_sig_overlap_samedir) #291:
dim(norm_counts_cut_0.2_G2_resKresLDPN)[1]/nrow(resKresL_DPN_sig_overlap_samedir) #300:
#So similar across both pairs...

#SO probably some pertinent questions are:
#Are the shared genes showing similar patterns in the individual ones.

#KR55: More reversion in single population "DP(N)" genes (89.9% and 77.6%) compared to shared ones (70.8%). 
100*table(norm_counts_cut_0.2_G1_resKDPN_notLDPN$Type)[1]/nrow(norm_counts_cut_0.2_G1_resKDPN_notLDPN) #% overshooting: 4.4
100*table(norm_counts_cut_0.2_G1_resKDPN_notLDPN$Type)[2]/nrow(norm_counts_cut_0.2_G1_resKDPN_notLDPN) #% reversion: 89.9%
100*table(norm_counts_cut_0.2_G1_resKDPN_notLDPN$Type)[3]/nrow(norm_counts_cut_0.2_G1_resKDPN_notLDPN) #% reinforcement; 5.7%
table(norm_counts_cut_0.2_G1_resKDPN_notLDPN$Type)[3] #So there are 23 reinforced genes in here...

100*table(norm_counts_cut_0.2_G2_resLDPN_notKDPN$Type)[1]/nrow(norm_counts_cut_0.2_G2_resLDPN_notKDPN) #% overshooting 13.1%
100*table(norm_counts_cut_0.2_G2_resLDPN_notKDPN$Type)[2]/nrow(norm_counts_cut_0.2_G2_resLDPN_notKDPN) #% reversion 77.6%
100*table(norm_counts_cut_0.2_G2_resLDPN_notKDPN$Type)[3]/nrow(norm_counts_cut_0.2_G2_resLDPN_notKDPN) #% reinforcement; 9.27%

table(norm_counts_cut_0.2_G2_resLDPN_notKDPN$Type)[3] #And 290 reinforced genes in here.

100*table(norm_counts_cut_0.2_G1_resKresLDPN$Type)[1]/nrow(norm_counts_cut_0.2_G1_resKresLDPN) #% overshooting: 11.7%
100*table(norm_counts_cut_0.2_G1_resKresLDPN$Type)[2]/nrow(norm_counts_cut_0.2_G1_resKresLDPN) #% reversion: 69.8%
100*table(norm_counts_cut_0.2_G1_resKresLDPN$Type)[3]/nrow(norm_counts_cut_0.2_G1_resKresLDPN) #% reinforcement: 18.5%: so way more???

100*table(norm_counts_cut_0.2_G2_resKresLDPN$Type)[1]/nrow(norm_counts_cut_0.2_G2_resKresLDPN) #% overshooting: 11.8%
100*table(norm_counts_cut_0.2_G2_resKresLDPN$Type)[2]/nrow(norm_counts_cut_0.2_G2_resKresLDPN) #% reversion: 71.6%
100*table(norm_counts_cut_0.2_G2_resKresLDPN$Type)[3]/nrow(norm_counts_cut_0.2_G2_resKresLDPN) #% reinforcement: 16.6%...

#5) resBresC

##################################################################
#13E.1 resBresC_samedir genes (CEC genes)
##################################################################
norm_counts_0.2_resBresC = norm_counts_cut_0.2[rownames(norm_counts_cut_0.2) %in% resBresC_sig_overlap_samedir$Names,]
norm_counts_0.2_resBresC = norm_counts_cut_0.2[rownames(norm_counts_cut_0.2) %in% resBresC_sig_overlap_samedir$Names,]

#KR48 - CEC genes, 73% show "substantial" EC and PC. 
dim(norm_counts_0.2_resBresC)[1]/nrow(resBresC_sig_overlap_samedir)
dim(resBresC_sig_overlap_samedir)
dim(norm_counts_XY_resBresC[rownames(norm_counts_XY_resBresC) %in% rownames(norm_counts_0.2_resBresC),])[1]/dim(norm_counts_XY_resBresC)[1]
length(unique(c(rownames(norm_counts_0.2_resBresC), rownames(norm_counts_XY_resBresC))))/length(resBresC_sig_overlap_samedir$Names)


#KR50 (1/3) - CEC genes, OVER=3.3%, REV=62.7%, RI=34%. See graphs. 
100*table(norm_counts_0.2_resBresC$Type)[1]/(table(norm_counts_0.2_resBresC$Type)[1]+table(norm_counts_0.2_resBresC$Type)[2]+table(norm_counts_0.2_resBresC$Type)[3]) #% overshooting
100*table(norm_counts_0.2_resBresC$Type)[2]/(table(norm_counts_0.2_resBresC$Type)[1]+table(norm_counts_0.2_resBresC$Type)[2]+table(norm_counts_0.2_resBresC$Type)[3]) #% reversion
100*table(norm_counts_0.2_resBresC$Type)[3]/(table(norm_counts_0.2_resBresC$Type)[1]+table(norm_counts_0.2_resBresC$Type)[2]+table(norm_counts_0.2_resBresC$Type)[3]) #% reinforcement

BC_Table_0.2 = table(norm_counts_0.2_resBresC$Type)
BC_Table_0.2

BC_vec_0.2 = c(0,0,0)
BC_vec_0.2[2] = table(norm_counts_0.2_resBresC$Type)[1] #% overshooting
BC_vec_0.2[1] = table(norm_counts_0.2_resBresC$Type)[2] #% reversion
BC_vec_0.2[3] = table(norm_counts_0.2_resBresC$Type)[3] #% reinforcement
BC_vec_0.2 = data.frame(val = BC_vec_0.2)
BC_vec_0.2 = cbind(BC_vec_0.2, type = c("REV", "OVER", "RI"))
BC_vec_0.2$type= factor(BC_vec_0.2$type, levels = BC_vec_0.2$type)

#binom.test(x = c((BC_Table_0.2[1]+BC_Table_0.2[3]), BC_Table_0.2[2]), p = (KL_Table_0.2[1]+KL_Table_0.2[3])/(KL_Table_0.2[1]+KL_Table_0.2[2]+KL_Table_0.2[3]))
#KR53: test for proportions of CEC vs. DP genes. 
binom.test(x = c((BC_Table_0.2[3]), BC_Table_0.2[2]+BC_Table_0.2[1]), p = (KL_Table_0.2[3])/(KL_Table_0.2[1]+KL_Table_0.2[2]+KL_Table_0.2[3]))
binom.test(x = c((BC_Table_0.2[3]), BC_Table_0.2[2]+BC_Table_0.2[1]), p = (T_table_0.2[3])/(T_table_0.2[1]+T_table_0.2[2]+T_table_0.2[3]))

###################################################################
### 13 E.2 resBnotC, resCnotB: are there differences?
###################################################################

#13.E.2.1 : WHole transcriptome vs. BnotC, CnotB

norm_counts_cut_G1_resBnotC_0.2 = norm_counts_cut_0.2_G1[rownames(norm_counts_cut_0.2_G1) %in% resB_notC$Names,]
norm_counts_cut_G2_resCnotB_0.2 = norm_counts_cut_0.2_G2[rownames(norm_counts_cut_0.2_G2) %in% resC_notB$Names,]

#i) what proportion of these genes don't meet the EC/PC cutoffs?
dim(norm_counts_cut_G1_resBnotC_0.2) #1231
dim(norm_counts_cut_G1_resBnotC_0.2)[1]/nrow(resB_notC) #So 68.9% pass both thresholds...
#Update: only 37% of these genes pass the threshold, why?
#I guess it's again because the evidence is more stringent, this is what we were seeing before.

dim(norm_counts_cut_G2_resCnotB_0.2) #2630 
dim(norm_counts_cut_G2_resCnotB_0.2)/nrow(resC_notB) #So 65.2% pass both thresholds...
#So again now a much smaller proportion of these pass the threshold.

#ii) What proportion show substantial ancestral plasticity?
#For this we will need to go back and look at the pre-cutoff files...
norm_counts_G1_resBnotC = norm_counts[rownames(norm_counts) %in% resB_notC$Names,]
norm_counts_G2_resCnotB = norm_counts[rownames(norm_counts) %in% resC_notB$Names,]
dim(norm_counts_G1_resBnotC) #1785
dim(norm_counts_G2_resCnotB) #4035

#So for G1...

table(norm_counts_G1_resBnotC$Cutoff_GC_G1) #So 1456/1785, or...
table(norm_counts_G1_resBnotC$Cutoff_GC_G1)[2]/nrow(resB_notC) #...81.5% show substantial GC, and...
table(norm_counts_G1_resBnotC$Cutoff_PC_G1) #So 1456/1785, or...
table(norm_counts_G1_resBnotC$Cutoff_PC_G1)[2]/nrow(resB_notC) #...81.5% show substantial PC? But they're not the same genes I think. Weird.
#So in terms of the raw numbers (PC vs. GC still calculated in the same way) they're the same...
#but a much smaller subset are now being included because of the way we do the cutoffs.

table(norm_counts_G2_resCnotB$Cutoff_GC_G2) #So 3034/4035, or...
table(norm_counts_G2_resCnotB$Cutoff_GC_G2)[2]/nrow(resC_notB) #...75.2% show subantial GC (i.e. derived genetic change)
table(norm_counts_G2_resCnotB$Cutoff_PC_G2) #So 3461/4035, or...
table(norm_counts_G2_resCnotB$Cutoff_PC_G2)[2]/nrow(resC_notB) #...85.8% show substantial PC (i.e. ancestral plastic change)
#So in terms of the raw numbers (PC vs. GC still calculated in the same way) they're the same...
#but a much smaller subset are now being included because of the way we do the cutoffs.

#Ok so I think the real question is...for genes in resBnotC [G1 only}...do they show substantially lower anc. plast. in G2? I.e. is that why they're not DGE?
table(norm_counts_G1_resBnotC$Cutoff_PC_G2)[2]/nrow(resB_notC) #SO 81.4% of these are ancestrally plastic.
table(norm_counts_G2_resCnotB$Cutoff_PC_G1)[2]/nrow(resC_notB) #SO 85.2% of these are ancestrally plastic.

#I guess it could be the magnitude of plasticity is different? But maybe we should stop here. Does any of this add up to much?

#So there isn't really strong evidence here that a lack of ancestral plasticity is driving recruitment.

norm_counts_cut_G1_resBresC_0.2 = norm_counts_cut_0.2_G1[rownames(norm_counts_cut_0.2_G1) %in% resBresC_sig_overlap_samedir$Names,]
norm_counts_cut_G2_resBresC_0.2 = norm_counts_cut_0.2_G2[rownames(norm_counts_cut_0.2_G2) %in% resBresC_sig_overlap_samedir$Names,]
dim(norm_counts_cut_G1_resBresC_0.2)[1]/nrow(resBresC_sig_overlap_samedir)#291:
dim(norm_counts_cut_G2_resBresC_0.2)[1]/nrow(resBresC_sig_overlap_samedir) #300:
#So similar across both pairs...
#So you're still losing most genes here...apparently most of the reinforced ones.

dim(norm_counts_cut_G1_resBnotC_0.2)[1]/nrow(resB_notC)
dim(norm_counts_cut_G2_resCnotB_0.2)[1]/nrow(resC_notB)

#SO probably some pertinent questions are:
#Are the shared genes showing similar patterns in the individual ones.

#KR54: More reversion in one sample "CEC" genes 
100*table(norm_counts_cut_G1_resBnotC_0.2$Type)[1]/nrow(norm_counts_cut_G1_resBnotC_0.2) #% overshooting: 6.4% 
100*table(norm_counts_cut_G1_resBnotC_0.2$Type)[2]/nrow(norm_counts_cut_G1_resBnotC_0.2) #% reversion: 77.5%
100*table(norm_counts_cut_G1_resBnotC_0.2$Type)[3]/nrow(norm_counts_cut_G1_resBnotC_0.2) #% reinforcement; 16%
table(norm_counts_cut_G1_resBnotC_0.2$Type)[3] #So there are 197 reinforced genes in here...

100*table(norm_counts_cut_G2_resCnotB_0.2$Type)[1]/nrow(norm_counts_cut_G2_resCnotB_0.2) #% overshooting 7.14%
100*table(norm_counts_cut_G2_resCnotB_0.2$Type)[2]/nrow(norm_counts_cut_G2_resCnotB_0.2) #% reversion 81.8%
100*table(norm_counts_cut_G2_resCnotB_0.2$Type)[3]/nrow(norm_counts_cut_G2_resCnotB_0.2) #% reinforcement; 11.0%
table(norm_counts_cut_G2_resCnotB_0.2$Type)[3] #And 290 reinforced genes in here.

100*table(norm_counts_cut_G1_resBresC$Type)[1]/nrow(norm_counts_cut_G1_resBresC) #% overshooting: 4.46% 
100*table(norm_counts_cut_G1_resBresC$Type)[2]/nrow(norm_counts_cut_G1_resBresC) #% reversion: 60.1%
100*table(norm_counts_cut_G1_resBresC$Type)[3]/nrow(norm_counts_cut_G1_resBresC) #% reinforcement: 35.4%

100*table(norm_counts_cut_G2_resBresC$Type)[1]/nrow(norm_counts_cut_G2_resBresC) #% overshooting: 1.7%
100*table(norm_counts_cut_G2_resBresC$Type)[2]/nrow(norm_counts_cut_G2_resBresC) #% reversion: 67%
100*table(norm_counts_cut_G2_resBresC$Type)[3]/nrow(norm_counts_cut_G2_resBresC) #% reinforcement: 31.3%

#So a bit of a problem here somewhere...
#So to spell it out...
#For the combined set of genes (0.2), you get most passing that threshold as before.
#Considering each population separately, there's a smaller number of genes passing both thresholds.
#That's because the power you have decreases. Only ~37-40% of genes in BnotC, CnotB and BandC, are able to be included.
#I guess an issue is that this might be a slightly skewed set of genes.

lfcA_BC_REV_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "REV",]),]
lfcE_BC_REV_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "REV",]),]
lfcD_BC_REV_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "REV",]),]
lfcH_BC_REV_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "REV",]),]
resids_countsA_BC_REV_0.2 = abs(lfcA_BC_REV_0.2$log2FoldChange)
resids_countsE_BC_REV_0.2 = abs(lfcE_BC_REV_0.2$log2FoldChange)
resids_countsD_BC_REV_0.2 = abs(lfcD_BC_REV_0.2$log2FoldChange)
resids_countsH_BC_REV_0.2 = abs(lfcH_BC_REV_0.2$log2FoldChange)

lfcA_BC_RI_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "RI",]),]
lfcE_BC_RI_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "RI",]),]
lfcD_BC_RI_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "RI",]),]
lfcH_BC_RI_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "RI",]),]
resids_countsA_BC_RI_0.2 = abs(lfcA_BC_RI_0.2$log2FoldChange)
resids_countsE_BC_RI_0.2 = abs(lfcE_BC_RI_0.2$log2FoldChange)
resids_countsD_BC_RI_0.2 = abs(lfcD_BC_RI_0.2$log2FoldChange)
resids_countsH_BC_RI_0.2 = abs(lfcH_BC_RI_0.2$log2FoldChange)

lfcA_BC_OVERSHOOT_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "OVERSHOOT",]),]
lfcE_BC_OVERSHOOT_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "OVERSHOOT",]),]
lfcD_BC_OVERSHOOT_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "OVERSHOOT",]),]
lfcH_BC_OVERSHOOT_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_0.2_resBresC[norm_counts_0.2_resBresC$Type == "OVERSHOOT",]),]
resids_countsA_BC_OVERSHOOT_0.2 = abs(lfcA_BC_OVERSHOOT_0.2$log2FoldChange)
resids_countsE_BC_OVERSHOOT_0.2 = abs(lfcE_BC_OVERSHOOT_0.2$log2FoldChange)
resids_countsD_BC_OVERSHOOT_0.2 = abs(lfcD_BC_OVERSHOOT_0.2$log2FoldChange)
resids_countsH_BC_OVERSHOOT_0.2 = abs(lfcH_BC_OVERSHOOT_0.2$log2FoldChange)

norm_counts_noPC_0.2_resBresC = norm_counts_noPC_0.2[rownames(norm_counts_noPC_0.2) %in% resBresC_sig_overlap_samedir$Names,]
dim(norm_counts_noPC_0.2_resBresC)# So Out of 371, 

lfcA_BC_NOANCPLAST_0.2 = lfcA[lfcA$Names %in% rownames(norm_counts_noPC_0.2_resBresC),]
lfcE_BC_NOANCPLAST_0.2 = lfcE[lfcE$Names %in% rownames(norm_counts_noPC_0.2_resBresC),]
lfcD_BC_NOANCPLAST_0.2 = lfcD[lfcD$Names %in% rownames(norm_counts_noPC_0.2_resBresC),]
lfcH_BC_NOANCPLAST_0.2 = lfcH[lfcH$Names %in% rownames(norm_counts_noPC_0.2_resBresC),]

resids_countsA_BC_NOANCPLAST_0.2 = abs(lfcA_BC_NOANCPLAST_0.2$log2FoldChange)
resids_countsE_BC_NOANCPLAST_0.2 = abs(lfcE_BC_NOANCPLAST_0.2$log2FoldChange)
resids_countsD_BC_NOANCPLAST_0.2 = abs(lfcD_BC_NOANCPLAST_0.2$log2FoldChange)
resids_countsH_BC_NOANCPLAST_0.2 = abs(lfcH_BC_NOANCPLAST_0.2$log2FoldChange)
resids_countsH_BC_NOANCPLAST_0.2 = as.data.frame(resids_countsH_BC_NOANCPLAST_0.2)
colnames(resids_countsH_BC_NOANCPLAST_0.2) = "VALUE"

resids_countsH_BC_NOANCPLAST_0.2$Names = "NOANCPLAST"
dim(resids_countsH_BC_NOANCPLAST_0.2)
#SO that's more like it (woah 371 that seems like a lot)
dim(resBresC_sig_overlap_samedir)
#...that can't be right...
median(resids_countsH_BC_NOANCPLAST_0.2$VALUE)
quantile(resids_countsH_BC_NOANCPLAST_0.2$VALUE)
#0.12. So the same as all the rest of them, then. 

resids_countsH_BC_RI_0.2 = as.data.frame(resids_countsH_BC_RI_0.2)
colnames(resids_countsH_BC_RI_0.2) = "VALUE"
resids_countsH_BC_RI_0.2$Names = "Reinforcement"
resids_countsH_BC_REV_0.2 = as.data.frame(resids_countsH_BC_REV_0.2)
colnames(resids_countsH_BC_REV_0.2) = "VALUE"
resids_countsH_BC_REV_0.2$Names = "REV"
resids_countsH_BC_OVERSHOOT_0.2 = as.data.frame(resids_countsH_BC_OVERSHOOT_0.2)
colnames(resids_countsH_BC_OVERSHOOT_0.2) = "VALUE"
resids_countsH_BC_OVERSHOOT_0.2$Names = "OVER"

Plast_FC_BC_0.2 = rbind(resids_countsH_BC_OVERSHOOT_0.2, resids_countsH_BC_REV_0.2, resids_countsH_BC_RI_0.2)

Plast_FC_BC_wilcox_0.2 = wilcox.test(resids_countsH_BC_NOANCPLAST_0.2$VALUE, Plast_FC_BC_0.2$VALUE, alternative = "two.sided")
Plast_FC_BC_wilcox_0.2

temp_Plast_FC_BC_0.2 = Plast_FC_BC_0.2
temp_Plast_FC_BC_0.2$Names = "PLAST"

PlastNoPlast_FC_BC_0.2 = rbind(resids_countsH_BC_NOANCPLAST_0.2, temp_Plast_FC_BC_0.2)
median(temp_Plast_FC_BC_0.2$VALUE)
median(resids_countsH_BC_NOANCPLAST_0.2$VALUE)

PlastNoPlast_FC_KL_wilcox_0.2 = wilcox.test(temp_Plast_FC_KL_0.2$VALUE, resids_countsH_KL_NOANCPLAST_0.2$VALUE, alternative = "two.sided")
#KR56 (1/2)
PlastNoPlast_FC_KL_wilcox_0.2
PlastNoPlast_FC_BC_wilcox_0.2 = wilcox.test(temp_Plast_FC_BC_0.2$VALUE, resids_countsH_BC_NOANCPLAST_0.2$VALUE, p.adjust.method = "BH")
#KR57 (1/2)
PlastNoPlast_FC_BC_wilcox_0.2 #So there isn't any difference. Fine.

#Ok so that's the BC stuff done.

#6 Plotting the graphs

T_vec_0.2$title = "Transcriptome-wide"
T_vec_0.2$type = c("Reversion", "Overshooting", "Reinforcement")
T_vec_0.2$type = factor(T_vec$type, levels = c("Reinforcement", "Overshooting", "Reversion"), ordered = T)

#KR41 (1/2)
T_barplot_0.2 = ggplot(data = T_vec_0.2, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity', color = "black")+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(axis.text.y = element_text(size = 25),
                               axis.text.x = element_text(size = 25), axis.title.x = element_text(size=25))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
        strip.text = element_text(size = 30))+facet_grid(.~title)+coord_flip()+scale_fill_manual(values = c(color_list[1], color_list[2], color_list[3]))

T_barplot_0.2

#KR51 - Barpolot for CEC genes in OVER, REV, RI
BC_vec_0.2$title = "CEC genes"
BC_vec_0.2$type = c("Reversion", "Overshooting", "Reinforcement")
BC_vec_0.2$type = factor(BC_vec_0.2$type, levels = c("Reinforcement", "Overshooting", "Reversion"), ordered = T)
BC_barplot_0.2 = ggplot(data = BC_vec_0.2, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity', color = "black")+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(axis.text.y = element_text(size = 25),
                               axis.text.x = element_text(size = 25), axis.title.x = element_text(size=25))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
        strip.text = element_text(size = 30))+facet_grid(.~title)+coord_flip()+scale_fill_manual(values = color_list)

BC_barplot_0.2

KL_vec_0.2$title = "AP genes"
KL_vec_0.2$type = c("Reversion", "Overshooting", "Reinforcement")
KL_vec_0.2$type = factor(KL_vec_0.2$type, levels = c("Reinforcement", "Overshooting", "Reversion"), ordered = T)

#KR46 (1/2) - Barplot showing DP(0) genes classified into OVER, REV, RI
KL_barplot_0.2 = ggplot(data = KL_vec_0.2, aes(x = type, y = val, fill = type))+
  geom_bar(width = 0.6, show.legend = FALSE, stat = 'identity', color = "black")+xlab("")+ylab("Number of genes\n")+
  scale_fill_viridis_d()+theme(axis.text.y = element_text(size = 25),
                               axis.text.x = element_text(size = 25), axis.title.x = element_text(size=25))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"), 
        strip.text = element_text(size = 30))+facet_grid(.~title)+coord_flip()+
  scale_fill_manual(values = color_list)

KL_barplot_0.2

#KR21: Comparable levels of |FC| in genes lacking vs. having ancestral plasticity in CEC genes
#KR57 (2/2)
BC_PlastNoPlast_Boxplot_0.2 = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_BC_0.2,fill = another_plot_purple, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size =20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("NP", "P"))+coord_cartesian(ylim = c(0,1.1))+
  #facet_grid(.~Set)+
  xlab("")+ylab("Absolute log2 fold change (|FC|)")+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = c("solid", "solid")), 
        text = element_text(size = 30),
        strip.text = element_text(size = 30))
BC_PlastNoPlast_Boxplot_0.2

#KR56 (2/2)
KL_PlastNoPlast_Boxplot_0.2 = ggplot()+
  geom_boxplot(data = PlastNoPlast_FC_KL_0.2,fill = another_plot_purple, aes(x = Names, y=VALUE), width = 0.7, outlier.shape = NA, show.legend = FALSE)+
  #  scale_fill_manual(values = c(ctrl_col, zinc_col, ctrl_col, zinc_col, zinc_col, ctrl_col, zinc_col, ctrl_col))+
  theme(axis.text.x = element_text(size =20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20))+
  scale_x_discrete(labels = c("NP", "P"))+
  #facet_grid(.~Set)+
  xlab("")+ylab("Absolute log2 fold change (|FC|)")+coord_cartesian(ylim=c(0,0.7))+
  #theme(panel.background = element_rect(fill = "lemonchiffon1"))
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = c("solid", "solid")), 
        text = element_text(size = 30),
        strip.text = element_text(size = 30))
KL_PlastNoPlast_Boxplot_0.2


###########################################
### Plotting the green/yellow diagrams
############################################

#EC/GC plots

ncc_lim_0.2 = max(c(quantile(norm_counts_cut_0.2$PC)[2], quantile(norm_counts_cut_0.2$PC)[4], quantile(norm_counts_cut_0.2$GC)[2], quantile(norm_counts_cut_0.2$GC)[4]))
ap_lim_0.2 = max(c(quantile(norm_counts_resKresLDPN_0.2$PC)[2], quantile(norm_counts_resKresLDPN_0.2$PC)[4], quantile(norm_counts_resKresLDPN_0.2$GC)[2], quantile(norm_counts_resKresLDPN_0.2$GC)[4]))
cec_lim_0.2 = max(c(quantile(norm_counts_0.2_resBresC$PC)[2], quantile(norm_counts_0.2_resBresC$PC)[4], quantile(norm_counts_0.2_resBresC$GC)[2], quantile(norm_counts_0.2_resBresC$GC)[4]))

biggest_lim_0.2 = max(c(ncc_lim_0.2, ap_lim_0.2, cec_lim_0.2))
biggest_lim_0.2


color_list = c("#75EAA9","#B6EA50", "#F5FF74")
color_list
over1_0.2 = data.frame(x = c(0, -Inf, -Inf), y = c(0,0,0.5*biggest_lim_0.2*1.1))
over2_0.2 = data.frame(x = c(0, Inf, Inf), y = c(0,0,-0.5*biggest_lim_0.2*1.1))
over3_0.2 = data.frame(x = c(0, -Inf, -Inf, 0), y = c(0, 0.5*biggest_lim_0.2*1.1, Inf, Inf)) 
over4_0.2 = data.frame(x = c(0, Inf, Inf, 0), y = c(0, -0.5*biggest_lim_0.2*1.1, -Inf, -Inf)) 

over1_0.2 = data.frame(x = c(0, -Inf, -Inf), y = c(0,0,0.5*biggest_lim_0.2*1.1))
over2_0.2 = data.frame(x = c(0, Inf, Inf), y = c(0,0,-0.5*biggest_lim_0.2*1.1))
over3_0.2 = data.frame(x = c(0, -Inf, -Inf, 0), y = c(0, 0.5*biggest_lim_0.2*1.1, Inf, Inf)) 
over4_0.2 = data.frame(x = c(0, Inf, Inf, 0), y = c(0, -0.5*biggest_lim_0.2*1.1, -Inf, -Inf)) 

#KR42
Whole_ECGC_0.2 = ggplot(data = data.frame())+geom_point()+
  #geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  #geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  #geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  #geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_cut_0.2, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
Whole_ECGC

#KR47 
WholeAP_ECGC_0.2 = ggplot(data = data.frame())+geom_point()+
  #geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  #geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  #geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  #geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_resKresLDPN_0.2, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
WholeAP_ECGC_0.2

#KR52
WholeCEC_ECGC_0.2 = ggplot(data = data.frame())+geom_point()+
  #geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  #geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  #geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  #geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  #geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_0.2_resBresC, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
WholeCEC_ECGC_0.2

#KR41 (2/2) - yellow/green diagram for DP(0) genes, EC vs PC
Middle_ECGC_0.2 = ggplot(data = data.frame())+geom_point()+xlim(-biggest_lim_0.2, biggest_lim_0.2)+ylim(-biggest_lim_0.2, biggest_lim_0.2)+
  geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_cut_0.2, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "white")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "white")+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
Middle_ECGC_0.2

#KR46 (2/2)
AP_ECGC_0.2 = ggplot(data = data.frame())+geom_point()+xlim(-biggest_lim_0.2, biggest_lim_0.2)+ylim(-biggest_lim_0.2, biggest_lim_0.2)+
  geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+geom_abline(slope = -0.5, intercept = 0)+
  geom_abline(slope = -1, intercept = 0, lty = 2)+geom_bin2d(data = norm_counts_resKresLDPN_0.2, aes(PC, GC), bins=50)+
  geom_hline(yintercept = -biggest_lim_0.2*1.0)+geom_hline(yintercept = biggest_lim_0.2)+
  geom_vline(xintercept = -biggest_lim_0.2)+geom_vline(xintercept = biggest_lim_0.2)+
  theme(text = element_text(size = 20))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
AP_ECGC_0.2

#KR51 (2/2) - CEC genes green/yellow plot showing EC vs. P
CEC_ECGC_0.2 = ggplot(data = data.frame())+geom_point()+xlim(-biggest_lim_0.2, biggest_lim_0.2)+ylim(-biggest_lim_0.2, biggest_lim_0.2)+
  geom_rect(aes(xmin = -Inf, xmax  = 0, ymin = -Inf, ymax = 0), fill=color_list[1])+
  geom_rect(aes(xmin = Inf, xmax  = 0, ymin = Inf, ymax = 0), fill=color_list[1])+
  geom_polygon(data = over1, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over2, aes(x=x, y=y), fill = color_list[2])+
  geom_polygon(data = over3, aes(x=x, y=y), fill = color_list[3])+
  geom_polygon(data = over4, aes(x=x, y=y), fill = color_list[3])+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_bin2d(data = norm_counts_0.2_resBresC, aes(PC, GC), bins=50)+
  geom_abline(slope = -0.5, intercept = 0, col = "black")+
  geom_abline(slope = -1, intercept = 0, lty = 2, col = "black")+
  theme(text = element_text(size = 22))+xlab("Plastic Change (PC)")+ylab("Evolutionary Change (EC)")
CEC_ECGC_0.2




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
#snp_tree = read.tree(file = "C:/Users/Danie/DEA_SNPhylo.bs.tree")
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

pdf(file = woof, "C:~/RNA_Seq/Figure1A_270321.pdf")
dev.off()

##################################################
###################### Even more Miscellaneous
