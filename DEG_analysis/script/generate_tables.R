# Libraries
library("VennDiagram") 
library("UpSetR")
library(ComplexUpset)
library(ggvenn)
library(ggVennDiagram)
library(ggplot2)
library("RColorBrewer")


### Functions

save_matrix <- function(matrix, filename){
  return()
  write.table(x = matrix, file = filename, sep = ",", row.names = T, col.names = NA)
}


####

### Main

# Loading data
TPM_table = read.table("kallisto/kallisto.gene.TMM.EXPR.matrix",
                       header=T, row.names=1, com='', check.names=F)

trinotate = read.table("trinotate_annotation_report.xls",
                       header=T, stringsAsFactors = F,sep = "\t", quote="", com='')
#

# Add median information
TPM_table$median_DMSO_1dpa = apply(TPM_table[,c(1:3)],1,median)
TPM_table$median_DMSO_2dpa = apply(TPM_table[,c(4:6)],1,median)
TPM_table$median_LY_1dpa = apply(TPM_table[,c(7:9)],1,median)
TPM_table$median_LY_2dpa = apply(TPM_table[,c(10:12)],1,median)

#

## DMSO_1dpa matrix generation
# All genes
matrix_DMSO_1dpa = TPM_table[,c(1,2,3,13)]
save_matrix(matrix_DMSO_1dpa, "Tables/DMSO_1dpa/DMSO_1dpa_TMM_matrix.csv")

# 1TPM threshold
matrix_DMSO_1dpa_1TPM = matrix_DMSO_1dpa[which(matrix_DMSO_1dpa$median_DMSO_1dpa>1),]
save_matrix(matrix_DMSO_1dpa_1TPM, "Tables/DMSO_1dpa/DMSO_1dpa_1TMM_threshold_matrix.csv")

# 5TPM threshold
matrix_DMSO_1dpa_5TPM = matrix_DMSO_1dpa[which(matrix_DMSO_1dpa$median_DMSO_1dpa>5),]
save_matrix(matrix_DMSO_1dpa_5TPM, "Tables/DMSO_1dpa/DMSO_1dpa_5TMM_threshold_matrix.csv")

# 10TPM threshold
matrix_DMSO_1dpa_10TPM = matrix_DMSO_1dpa[which(matrix_DMSO_1dpa$median_DMSO_1dpa>10),]
save_matrix(matrix_DMSO_1dpa_10TPM, "Tables/DMSO_1dpa/DMSO_1dpa_10TMM_threshold_matrix.csv")



## DMSO_2dpa matrix generation
# All genes
matrix_DMSO_2dpa = TPM_table[,c(4,5,6,14)]
save_matrix(matrix_DMSO_2dpa, "Tables/DMSO_2dpa/DMSO_2dpa_TMM_matrix.csv")

# 1TPM threshold
matrix_DMSO_2dpa_1TPM = matrix_DMSO_2dpa[which(matrix_DMSO_2dpa$median_DMSO_2dpa>1),]
save_matrix(matrix_DMSO_2dpa_1TPM, "Tables/DMSO_2dpa/DMSO_2dpa_1TMM_threshold_matrix.csv")

# 5TPM threshold
matrix_DMSO_2dpa_5TPM = matrix_DMSO_2dpa[which(matrix_DMSO_2dpa$median_DMSO_2dpa>5),]
save_matrix(matrix_DMSO_2dpa_5TPM, "Tables/DMSO_2dpa/DMSO_2dpa_5TMM_threshold_matrix.csv")

# 10TPM threshold
matrix_DMSO_2dpa_10TPM = matrix_DMSO_2dpa[which(matrix_DMSO_2dpa$median_DMSO_2dpa>10),]
save_matrix(matrix_DMSO_2dpa_10TPM, "Tables/DMSO_2dpa/DMSO_2dpa_10TMM_threshold_matrix.csv")



## LY_1dpa matrix generation
# All genes
matrix_LY_1dpa = TPM_table[,c(7,8,9,15)]
save_matrix(matrix_LY_1dpa, "Tables/LY_1dpa/LY_1dpa_TMM_matrix.csv")

# 1TPM threshold
matrix_LY_1dpa_1TPM = matrix_LY_1dpa[which(matrix_LY_1dpa$median_LY_1dpa>1),]
save_matrix(matrix_LY_1dpa_1TPM, "Tables/LY_1dpa/LY_1dpa_1TMM_threshold_matrix.csv")

# 5TPM threshold
matrix_LY_1dpa_5TPM = matrix_LY_1dpa[which(matrix_LY_1dpa$median_LY_1dpa>5),]
save_matrix(matrix_LY_1dpa_5TPM, "Tables/LY_1dpa/LY_1dpa_5TMM_threshold_matrix.csv")

# 10TPM threshold
matrix_LY_1dpa_10TPM = matrix_LY_1dpa[which(matrix_LY_1dpa$median_LY_1dpa>10),]
save_matrix(matrix_LY_1dpa_10TPM, "Tables/LY_1dpa/LY_1dpa_10TMM_threshold_matrix.csv")





## LY_2dpa matrix generation
# All genes
matrix_LY_2dpa = TPM_table[,c(10,11,12,16)]
save_matrix(matrix_LY_2dpa, "Tables/LY_2dpa/LY_2dpa_TMM_matrix.csv")

# 1TPM threshold
matrix_LY_2dpa_1TPM = matrix_LY_2dpa[which(matrix_LY_2dpa$median_LY_2dpa>1),]
save_matrix(matrix_LY_2dpa_1TPM, "Tables/LY_2dpa/LY_2dpa_1TMM_threshold_matrix.csv")

# 5TPM threshold
matrix_LY_2dpa_5TPM = matrix_LY_2dpa[which(matrix_LY_2dpa$median_LY_2dpa>5),]
save_matrix(matrix_LY_2dpa_5TPM, "Tables/LY_2dpa/LY_2dpa_5TMM_threshold_matrix.csv")

# 10TPM threshold
matrix_LY_2dpa_10TPM = matrix_LY_2dpa[which(matrix_LY_2dpa$median_LY_2dpa>10),]
save_matrix(matrix_LY_2dpa_10TPM, "Tables/LY_2dpa/LY_2dpa_10TMM_threshold_matrix.csv")


#### Venn diagram

y <- list(
  DMSO_1dpa = rownames(matrix_DMSO_1dpa_1TPM),
  DMSO_2dpa = rownames(matrix_DMSO_2dpa_1TPM),
  LY_1dpa = rownames(matrix_LY_1dpa_1TPM),
  LY_2dpa = rownames(matrix_LY_2dpa_1TPM)
)

UpSetR::upset(fromList(y), order.by = "freq",number.angles = 0,
      line.size = 1.5,point.size = 3.5,
      mainbar.y.label = "Gene intersections", sets.x.label = "Genes Per Samples",
      sets = c("DMSO_1dpa","DMSO_2dpa","LY_1dpa","LY_2dpa"),
      text.scale = c(1.3,1.3,1,1,2,1.1), keep.order = T, nsets = 4,
      matrix.color = "black", main.bar.color = "black", show.numbers = "yes")

brewer_palette<-brewer.pal(4,"Set2")

venn.diagram(y,"Figures/Venn_1TMM.tiff",
             category = c("DMSO_1dpa", "DMSO_2dpa", "LY_1dpa", "LY_2dpa"),
             scaled = F,
             ext.text = T,
             lwd = 2,
             fill = brewer_palette,
             cex = 1,
             #cat.pos = c(-20, 20),
             cat.cex = 1,
             cat.dist = 0.1)

####