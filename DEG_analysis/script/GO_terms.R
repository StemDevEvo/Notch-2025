# Script to generate GOterms annotation / enrichment from a list of gene.

library(clusterProfiler)
library(AnnotationHub)
library(enrichplot)

### Function


prepare_DEG_gene = function(blast_res){
  if(blast_res == "."){
    return()
  }else{
    return(unlist(strsplit(blast_res, " "))[1])
  }
}

return_uniprot_id = function(full_name){
  splitted_id = unlist(strsplit(x=full_name, split = "|", fixed = T))
  return(splitted_id[2])
}

match_id = function(id, blast_table){
  return(blast_table$BlastResult[which(blast_table$TranscriptID == paste(id,"t1",sep=""))])
}


####



#########################

# OrgDb : Mouse.

hub <- AnnotationHub()

q <- query(hub, "org.Mm.eg.db")
id <- q$ah_id[length(q)]
OrgDb = hub[[id]]

##########################

# Loading TMM matrices to determine universe (all genes initially compared)

DMSO_1dpa_1tmm = read.table("Tables/DMSO_1dpa/DMSO_1dpa_1TMM_threshold_matrix.csv",
                            header = T, row.names = 1, com='', sep = ",")

DMSO_2dpa_1tmm = read.table("Tables/DMSO_2dpa/DMSO_2dpa_1TMM_threshold_matrix.csv",
                            header = T, row.names = 1, com='', sep = ",")

LY_1dpa_1tmm = read.table("Tables/LY_1dpa/LY_1dpa_1TMM_threshold_matrix.csv",
                          header = T, row.names = 1, com='', sep = ",")

LY_2dpa_1tmm = read.table("Tables/LY_2dpa/LY_2dpa_1TMM_threshold_matrix.csv",
                          header = T, row.names = 1, com='', sep = ",")


universe = unique(c(rownames(DMSO_1dpa_1tmm),rownames(DMSO_2dpa_1tmm),
                    rownames(LY_1dpa_1tmm),rownames(LY_2dpa_1tmm)))

###########################

# Load mouse annotation :

ann_table = read.table(file = "Mus_musculus_table.tsv", 
                       sep = "\t", header = T)


##########################

# List to compare

liste_directories = list.dirs(path = "DEG/")

universe_specie = lapply(X = universe, FUN = match_id, blast_table = ann_table)

df = data.frame(universe_transcript = universe,
                universe_specie = unlist(universe_specie))

universe_gene = df$universe_specie[df$universe_specie != "."]
universe_gene = unlist(lapply(X = universe_gene, FUN = return_uniprot_id))

for (i in liste_directories){
  if (i != "DEG/") {
    i = substr(i,6,50)
    DEG_result = read.table(file = paste0("DEG/",i,"/twisted.tsv"),
                            sep = "\t", header = T, row.names = 1)
    
    
    DEG_gene = unlist(lapply(DEG_result$Mouse, prepare_DEG_gene))
    
    DEG_gene = unlist(lapply(X = DEG_gene, FUN = return_uniprot_id))
    
    
    ## GO terms enrichment (one condition only)
    
    file_path = paste("GOterms/",i,"/",sep ="")
    ont = c("BP","MF", "CC")
    for (y in ont) {
      ego <- enrichGO(gene = DEG_gene,
                      universe = universe_gene,
                      OrgDb = OrgDb,
                      keyType = "UNIPROT",
                      ont = y,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = T)
      
      write.table(x = ego, file = paste(file_path,"GO_terms_enrichment_",y,".tsv", sep = ""), sep = "\t",quote = F)
      
      
      if(length(which(ego@result$p.adjust < 0.05)) > 5){
        
        pdf(file = paste(file_path,"GO_terms_plot_",y,".pdf",sep=""))
        print(goplot(ego, showCategory = 10))
        dev.off()
        
        pdf(file = paste(file_path,"barplot_top20_",y,".pdf",sep=""))
        print(barplot(ego, showCategory = 20))
        dev.off()
        
        pdf(file = paste(file_path,"dotplot_top20_",y,".pdf", sep =""))
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        
        egox = setReadable(ego, OrgDb, "UNIPROT")
        
        egox2 = pairwise_termsim(egox)
        pdf(file = paste(file_path,"Tree_plot_clustering_",y,".pdf",sep=""), width = 17, height = 8)
        print(treeplot(egox2, hclust_method = "average"))
        dev.off()
        
        ego_enr = pairwise_termsim(ego)
        
        pdf(file = paste(file_path,"Network_plot",y,".pdf",sep=""), width = 12, height = 12)
        print(emapplot(ego_enr))
        dev.off()
      }
    }
  }
}


########### Compare cluster

list_DEG_genes = list()
for (i in liste_directories){
  if (i != "DEG/") {
    i = substr(i,6,50)
    gene_list = read.table(
      file = paste("DEG/",i,"/twisted.tsv",sep = ""),
      sep = "\t", header = T, row.names = 1)
    
    unordered_list_mouse = gene_list$Mouse
    
    DEG_gene = unlist(lapply(X = unlist(lapply(unordered_list_mouse, prepare_DEG_gene)), FUN = return_uniprot_id))
    
    list_DEG_genes[[i]] = DEG_gene
  }
}

## Compare 2 cluster

ck <- compareCluster(geneCluster = list_DEG_genes, fun = enrichGO, keyType = "UNIPROT",
                     OrgDb = OrgDb, universe = universe_gene, ont = "MF",pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)


ck <- setReadable(ck, OrgDb = OrgDb, keyType = "UNIPROT")

dotplot(ck, showCategory = 10, label_format = 20)

################# Kegg pathways
library(pathview)

get_logFC = function(DEG_result){
  if(DEG_result$Mouse == "."){
    return()
  }else{
    return(DEG_result$X.logFC.)
  }
}

for (i in liste_directories){
  i = substr(i,6,50)
  if (i != "") {
    if(!dir.exists(paste0("Kegg/",i))){
      dir.create(paste0("Kegg/",i))
    }
    file_path = paste("Kegg/",i,"/",sep ="")
    DEG_result = read.table(file = paste0("DEG/",i,"/twisted.tsv"),
                            sep = "\t", header = T, row.names = 1)
    DEG_gene = unlist(lapply(DEG_result$Mouse, prepare_DEG_gene))
    DEG_gene = unlist(lapply(X = DEG_gene, FUN = return_uniprot_id))
    
    logFC = ifelse(DEG_result$Mouse ==".",NA,DEG_result$X.logFC.)
    logFC = logFC[!is.na(logFC)]
    names(logFC) = DEG_gene
    kk <- enrichKEGG(gene = DEG_gene,
                     universe = universe_gene,
                     keyType = "uniprot",
                     organism = "mmu",
                     pvalueCutoff = 0.05)
    write.table(x = kk, file = paste0(file_path,"Kegg_pathways_enrichment.tsv"), sep = "\t",quote = F, col.names = NA)
  
    
    ENTREZ_name_logfc = id2eg(ids = names(logFC), category = "UNIPROT", org = "Mm")
    
    names(logFC) = ENTREZ_name_logfc[,2]
    
    significant_ids = kk@result$ID[1:length(which(kk@result$p.adjust < 0.05))]
    
    for (y in 1:5) {
      current_directory = getwd()
      setwd(file_path)
      current_pathway_id = significant_ids[y]
      if(!is.na(current_pathway_id)){
        pathview(gene.data = logFC,
               pathway.id = current_pathway_id,
               gene.idtype = "ENTREZ",
               species = "mmu",
               limit = list(gene=max(abs(logFC)), cpd=1))
        setwd(current_directory)
      }
    }
  }
}
