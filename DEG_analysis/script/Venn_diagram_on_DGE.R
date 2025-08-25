# Libraries
library("VennDiagram") 
library("RColorBrewer")
library(ggplot2)
library(ggvenn)


######## Parameters

DMSO_vs_LY_1dpa = read.table(file = "DEG/DMSO_1dpavsLY_1dpa/DEG_result_DMSO_1dpavsLY_1dpa_FC_and_FDR_filter.tsv",
                             sep = "\t", header = T, row.names = 1)

DMSO_vs_LY_2dpa = read.table(file = "DEG/DMSO_2dpavsLY_2dpa/DEG_result_DMSO_2dpavsLY_2dpa_FC_and_FDR_filter.tsv",
                             sep = "\t", header = T, row.names = 1)


##



######## Main

y <- list(
  dpa1 = rownames(DMSO_vs_LY_1dpa),
  dpa2 = rownames(DMSO_vs_LY_2dpa)
)

brewer_palette<-brewer.pal(3,"Set2")

# Using VennDiagram

VennDiagram::venn.diagram(y,"figures/Venn_1dpa_vs_2dpa.tiff",
             category = c("1 dpa","2 dpa"),
             scaled = F,
             ext.text = T,
             lwd = 2,
             fill = brewer_palette[1:2],
             cex = 1,
             cat.cex = 1,
             cat.dist = 0.05)

# Using ggvenn

dpa1_up = length(which(DMSO_vs_LY_1dpa[which(!(rownames(DMSO_vs_LY_1dpa) %in% rownames(DMSO_vs_LY_2dpa))),1] > 0))
dpa1_down = length(which(DMSO_vs_LY_1dpa[which(!(rownames(DMSO_vs_LY_1dpa) %in% rownames(DMSO_vs_LY_2dpa))),1] < 0))

dpa2_up = length(which(DMSO_vs_LY_2dpa[which(!(rownames(DMSO_vs_LY_2dpa) %in% rownames(DMSO_vs_LY_1dpa))),1] > 0))
dpa2_down = length(which(DMSO_vs_LY_2dpa[which(!(rownames(DMSO_vs_LY_2dpa) %in% rownames(DMSO_vs_LY_1dpa))),1] < 0))

ggvenn::ggvenn(y, show_percentage = F,
               fill_color = brewer_palette[1:2], fill_alpha = 0.7, stroke_size = 0.5,
               text_size = 5) +
  annotate("text", x=0.75,y=0.5,label = as.character(dpa2_up), colour = "darkgreen")+ # Top right
  annotate("text", x=0.75,y=-0.5,label = as.character(dpa2_down), colour = "red")+ # Bottom right
  annotate("text", x=-0.75,y=0.5,label = as.character(dpa1_up), colour = "darkgreen")+ # Top left
  annotate("text", x=-0.75,y=-0.5,label = as.character(dpa1_down), colour = "red") # Bottom left



##