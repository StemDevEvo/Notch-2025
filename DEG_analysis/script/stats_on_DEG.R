#Librairies
library(stringr)
library(plyr)
library("VennDiagram") 
library("RColorBrewer")
library(ggplot2)
library(ggvenn)

### Function

removequotes = function(list_to_remove){
  return(str_replace_all(list_to_remove, '\"', ""))
}

get_percent_gene_with_mouse_annotation= function(DEG_output){
  nbr_without = length(which(DEG_output$Mouse == "."))
  percentage = (nrow(DEG_output)-nbr_without) / nrow(DEG_output) *100
  return(percentage)
}

get_percentage_cluster = function(DEG_output, cluster){
  nbr_of_gene_with_cluster = length(which(DEG_output$Cluster == cluster))
  percentage = (nbr_of_gene_with_cluster) / nrow(DEG_output)*100
  return(percentage)
}

##


### Parameters

setwd("DEG")

liste_directories = list.dirs()

##



### Main

# Compute stats on all conditions

for (i in liste_directories) {
  if (i != ".") {
    current_table = read.table(file = paste0(i,"/DEG_result_",substr(i,3,50),"_FC_and_FDR_filter_annotated.tsv"),header = T, sep = "\t", quote = "", stringsAsFactors = F)
    current_table$X.. = removequotes(current_table$X..)
    rownames(current_table) = current_table$X..
    list_of_cluster = 0:11
    percent_with_mouse_annotation = get_percent_gene_with_mouse_annotation(current_table)
    percentage_of_clusters = list()
    for (y in list_of_cluster) {
      percentage_of_clusters[[paste0("Cluster_",as.character(y))]] = get_percentage_cluster(current_table, y)
    }
    up_regulated_genes = length(which(current_table$X.logFC. > 0))
    down_regulated_genes = length(which(current_table$X.logFC. < 0))
    output_dataframe = data.frame(up_regulated_genes, down_regulated_genes,
               Gene_with_mouse_annotation = percent_with_mouse_annotation,
               as.data.frame(percentage_of_clusters))
    write.table(output_dataframe,
                file =  paste0(i,"/statistics.tsv"),
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
}


# unique gene from the DMSO vs LY conditions

DMSO_vs_LY_table = read.table(file = paste0("DMSOvsLY/DEG_result_DMSOvsLY_FC_and_FDR_filter_annotated.tsv"),header = T, sep = "\t", quote = "", stringsAsFactors = F)

DMSO_1dpavsLY_1dpa = read.table(file = paste0("DMSO_1dpavsLY_1dpa/DEG_result_DMSO_1dpavsLY_1dpa_FC_and_FDR_filter_annotated.tsv"),header = T, sep = "\t", quote = "", stringsAsFactors = F)

DMSO_2dpavsLY_2dpa = read.table(file = paste0("DMSO_2dpavsLY_2dpa/DEG_result_DMSO_2dpavsLY_2dpa_FC_and_FDR_filter_annotated.tsv"),header = T, sep = "\t", quote = "", stringsAsFactors = F)


list_gene_unique_separate = unique(c(DMSO_1dpavsLY_1dpa$X.., DMSO_2dpavsLY_2dpa$X..))

list_gene_unique_regroup = DMSO_vs_LY_table[which(!(DMSO_vs_LY_table$X.. %in% list_gene_unique_separate)),]


write.table(x = list_gene_unique_regroup, file = "DMSOvsLY/DEG_result_gene_unique.tsv",
            col.names = T,row.names = list_gene_unique_regroup$X.., sep = "\t")

y = list(
  genes_DMSOvsLY = DMSO_vs_LY_table$X..,
  genes_addition = list_gene_unique_separate
)

brewer_palette<-brewer.pal(3,"Set2")

ggvenn::ggvenn(y, show_percentage = F,
               fill_color = brewer_palette[1:2], fill_alpha = 0.7, stroke_size = 0.5,
               text_size = 5)

##