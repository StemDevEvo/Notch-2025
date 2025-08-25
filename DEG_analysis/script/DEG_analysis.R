## Libraries
library(edgeR)
library(rgl)
library(ggfortify)
library("FactoMineR")
library("factoextra")
library("smacof")
library("dplyr")
library(ggplot2)
library(gplots)
library(statmod)
library(magrittr)
library(EnhancedVolcano)
library(Seurat)

#### Functions

plotPCA <- function(x, groups) {
  plot3d(x, col=c("red","blue","green","yellow"), type="s", size=1, axes=F)
  axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
  grid3d("x")
  grid3d("y")
  grid3d("z")
}



# Plot a sample correlation matrix using the log CPM count and a list
# giving column indexes.
plot_individual_cor_heatmap <- function(matrix, list_pos){
  
  nbr_of_rep = length(list_pos)
  head_names = colnames(matrix[,list_pos])
  submatrix = matrix[,list_pos]
  cor_matrix = matrix(nrow = nbr_of_rep, ncol = nbr_of_rep)
  rownames(cor_matrix) = head_names
  colnames(cor_matrix) = head_names
  
  for (i in 1:nbr_of_rep) {
    for (y in 1:nbr_of_rep) {
      cor_matrix[i,y] = cor(x = submatrix[,i], y = submatrix[,y])
    }
  }
  colors = c(seq(0.5,0.7,length=100),seq(0.71,0.79,length=100), seq(0.80,1,length=100))
  mypalette = colorRampPalette(c("purple","black","yellow"))(n = 299)
  
  par(mar = c(5.1, 6.1, 6.1, 10))
  heatmap.2(cor_matrix, col=mypalette,breaks = colors, key=T, keysize=1.5,
            density.info="none", trace="none",cexCol=0.9, cexRow = 0.9,
            main = "Sample correlation matrix", margins = c(6.5,7),
            dendrogram = 'both', symkey = F, symm = F, symbreaks = F,
            scale = "none", lwid = c(1,5), lhei = c(1.5,5)
  )
  
}


shorter <- function(str2short){
  if(grepl("^", str2short, fixed = T)){
    return(unlist(strsplit(str2short,split = "^", fixed = T))[[1]])
  }else{
    return(unlist(strsplit(str2short,split = "`", fixed = T))[[1]])
  }
}


# Concatenate a list into a single string.
agglom <- function(str2aggl){
  aggl = ""
  for (y in 1:length(unlist(str2aggl))) {
    aggl = paste(aggl, shorter(unlist(str2aggl)[y]))
  }
  return(aggl)
}

# Add annotation to a Edger DEG output matrix with the trinotate report.
add_annotation <- function(table, annotation){
  
  blastx = c()
  blastp = c()
  GOterms = c()
  kegg = c()
  
  for (i in 1:length(table[,1])) {
    line = annotation[which(annotation$X.gene_id == rownames(table)[i]),]
    blastx = c(blastx, agglom(line[3]))
    blastp = c(blastp, agglom(line[7]))
    GOterms = c(GOterms, agglom(line[13]))
    kegg = c(kegg, agglom(line[12])) 
  }
  table$blastx = blastx
  table$blastp = blastp
  table$GOterms = GOterms
  table$kegg = kegg
  return(table)
}

##


#### Main

# Data Loading:

data = read.table("kallisto/kallisto.gene.counts.matrix",
                  header = T, row.names = 1, com='', sep = "\t")

trinotate = read.table("trinotate_annotation_report.xls",
                       header=T, stringsAsFactors = F,sep = "\t", quote="", com='')

cluster_correspondance = read.table("cluster_to_geneID_correspondance.csv",
                                    sep = ",", header = T)


DMSO_1dpa_1tmm = read.table("Tables/DMSO_1dpa/DMSO_1dpa_1TMM_threshold_matrix.csv",
                       header = T, row.names = 1, com='', sep = ",")

DMSO_2dpa_1tmm = read.table("Tables/DMSO_2dpa/DMSO_2dpa_1TMM_threshold_matrix.csv",
                            header = T, row.names = 1, com='', sep = ",")

LY_1dpa_1tmm = read.table("Tables/LY_1dpa/LY_1dpa_1TMM_threshold_matrix.csv",
                            header = T, row.names = 1, com='', sep = ",")

LY_2dpa_1tmm = read.table("Tables/LY_2dpa/LY_2dpa_1TMM_threshold_matrix.csv",
                            header = T, row.names = 1, com='', sep = ",")


# Contrast preparation:
rep = factor(c("rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3"))
cond = factor(c("DMSO_1dpa","DMSO_1dpa","DMSO_1dpa","DMSO_2dpa","DMSO_2dpa","DMSO_2dpa","LY_1dpa","LY_1dpa","LY_1dpa","LY_2dpa","LY_2dpa","LY_2dpa"))
cond = relevel(cond, ref="DMSO_1dpa")

# Create EdgeR object and filter out lowly or non expressed genes.


keep = unique(c(rownames(DMSO_1dpa_1tmm),rownames(DMSO_2dpa_1tmm),
                rownames(LY_1dpa_1tmm),rownames(LY_2dpa_1tmm)))
all_gene = rownames(data)
keep2 = rep(FALSE,length(all_gene))
names(keep2) = all_gene
for (i in 1:length(all_gene)) {
  if(all_gene[i] %in% keep){
    keep2[i] = TRUE
  }
}


data = DGEList(counts = data, group = cond)
keep <- filterByExpr(data)
table(keep2)
data <- data[keep2, , keep.lib.sizes = F]


data = calcNormFactors(data)

### Analysis

## Data quality

# Plot MDS

data_mds = cpm(data, log=TRUE, prior.count = 1)
plotMDS(data_mds, col=rep(1:4, each = 3), main = "MDS plot before batch correction")
plotMD(data_mds, column = 1)

data_bc = removeBatchEffect(cpm(data, log=TRUE, prior.count = 1), batch = rep, design = model.matrix(~0+cond))
plotMDS(data_bc, col=rep(1:4, each = 3), main = "MDS plot after batch correction")

## Apply batch effect, run:
data_mds = data_bc
##


design = model.matrix(~rep+rep:cond)
logFC <- predFC(data, design, prior.count = 1, dispersion = 0.05)
cor(logFC[,6:6])


# compute top 1000 variable gene

data_sd = apply(data_mds, MARGIN = 1, FUN = sd)

data_sd = data.frame(data_mds, data_sd)

data_sd = data_sd[order(-data_sd$data_sd),]

top1000 = head(data_sd, n = 1000)

# ACP generation.
# 2D PCA
res.pca <- PCA(t(data_mds), graph = F)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0,100))
var <- get_pca_var(res.pca)
fviz_pca_ind(res.pca,axes = c(1,2), geom = c("point","text"),
             col.var = cond,
             title = "PCA from the initial matrix", col.ind = cond, addEllipses = T) + 
  theme(panel.background = element_rect(fill= "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_shape_manual(values=c(15,15,15,15,15,15,15,15,15,15,15,15))

# 3D PCA
pca = prcomp(t(data_mds), scale. = T)
plotPCA(pca$x[,1:3], cond)

#

# Dendrogram
dd <- dist(t(top1000[1:length(top1000[1,])-1]), method = "euclidean")
hc <- hclust(dd, method = "complete")

plot(hc, main = "Dendrogram clustering on top 1000 variable genes", xlab = "Condition")

dd <- dist(t(data_mds), method = "euclidean")
hc <- hclust(dd, method = "complete")

plot(hc, main = "Dendrogram clustering on all genes", xlab = "Condition")


# Heatmap generation

# Basic heatmap
heatmap(x = as.matrix(top1000[1:12]), main = "Heatmap",
        ylab = "Genes", xlab = "Condition")

# Heatmap with heatmap2
# matrix arrangement for heatmap.
mat_df =as.matrix(top1000[,1:12])
rownames(mat_df) = rownames(top1000)

colors = c(seq(-10,-2,length=100),seq(-1.9,1.9,length=100), seq(2,10,length=100))
mypalette = colorRampPalette(c("white","orange","red"))(n = 299)


heatmap.2(mat_df, col=mypalette,breaks = colors, key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9,
          main = "Heatmap of the top 1000 variable transcripts", labRow = F, margins = c(6.5,1.5),
          dendrogram = 'column', symkey = F, symm = F, symbreaks = F,
          scale = "none", xlab = "Samples", ylab = "Transcripts")

## Correlation matrix

cor_matrix = matrix(nrow = length(rep), ncol = length(rep))
rownames(cor_matrix) = paste(cond,rep,sep="_")
colnames(cor_matrix) = paste(cond,rep,sep="_")
cor_matrix
for (i in 1:length(cor_matrix[,1])) {
  for (y in 1:length(cor_matrix[,1])) {
    cor_matrix[i,y] = cor(x = data_bc[,i], y = data_bc[,y])
  }
}

colors = c(seq(0.4,0.7,length=100),seq(0.71,0.79,length=100), seq(0.80,1,length=100))
mypalette = colorRampPalette(c("purple","black","yellow"))(n = 299)

hclust2 = function(x){
  return(hclust(x, method = "mcquitty"))
}

par(mar = c(5.1, 6.1, 6.1, 3.1))
heatmap.2(cor_matrix,Rowv = T, Colv = T, col=mypalette,breaks = colors, key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9,
          main = "Sample correlation matrix", margins = c(8,10),
          dendrogram = 'both', symkey = F, symm = F, symbreaks = F,
          scale = "none", lwid = c(1,5), lhei = c(1.5,5),
          RowSideColors = c(
            rep("#81061A",3),
            rep("#A3B0E4",3),
            rep("#A3D0E4",3),
            rep("#B0E4A3",3)
          ),
          ColSideColors = c(
            rep("#81061A",3),
            rep("#A3B0E4",3),
            rep("#A3D0E4",3),
            rep("#B0E4A3",3)
          )
)



##
# Differential expression

design = model.matrix(~0+cond+rep)
rownames(design) <- colnames(data)

data <- estimateDisp(data, design, robust = T)

# tests
fit <- glmFit(data, design)

cond_fixed = unique(paste("cond", cond, sep =""))
nbr_of_cond = length(cond_fixed)
for (i in 1:(nbr_of_cond-1)){
  for(y in 1:(nbr_of_cond-i)){
    cond1 = cond_fixed[i]
    cond2 = cond_fixed[i+y]
    actual_contrast = makeContrasts(paste(cond2, "-", cond1,sep=""), levels = colnames(design))
    
    print(actual_contrast)
    lrt <- glmLRT(fit, contrast = actual_contrast)
    DET_matrix = topTags(lrt, n = length(lrt$table[,1]))
    cond1 = substr(cond1, 5,100)
    cond2 = substr(cond2, 5,100)
    volcano = EnhancedVolcano(DET_matrix$table,
                    lab = NA,
                    x = 'logFC',
                    y = 'FDR',
                    ylim = c(0,20),
                    title = paste("Volcano Plot of ",cond1, " vs ", cond2, " conditions",sep=""),
                    subtitle = bquote(italic("Log2 FC threshold = 1 and p_value (FDR) threshold = 0.05")),
                    FCcutoff = 1,
                    pCutoff = 0.05,
                    caption = paste0("total = ", nrow(DET_matrix$table), " genes"))
    ggsave(filename = paste("DEG/",cond1,"vs",cond2,"/volcano_plot.pdf", sep = ""), plot = volcano, device = "pdf")
    
    thresholded_matrix = DET_matrix$table[which(abs(DET_matrix$table$logFC) > 1 & DET_matrix$table$FDR < 0.05),]
    pdf(file = paste("DEG/",cond1,"vs",cond2,"/MA_plot.pdf", sep = ""))
    plotSmear(lrt, de.tags = rownames(thresholded_matrix))
    dev.off()
    write.table(x = thresholded_matrix, file = paste("DEG/",cond1,"vs",cond2,"/DEG_result_",cond1,"vs",cond2,"_FC_and_FDR_filter.tsv", sep = ""),
                sep = "\t", row.names = T, col.names = NA)
    
    write.table(x = DET_matrix$table, file = paste("DEG/",cond1,"vs",cond2,"/DEG_result_",cond1,"vs",cond2,".tsv", sep = ""),
                sep = "\t", row.names = T, col.names = NA)
    
    annotated_matrix = add_annotation(thresholded_matrix, trinotate)
    
    write.table(x = annotated_matrix,
                file = paste0("DEG/",cond1,"vs",cond2,"/DEG_result_",cond1,"vs",cond2,"_FC_and_FDR_filter_trinotate.tsv"),
                sep = "\t", row.names = T, col.names = NA)
  }
}
##