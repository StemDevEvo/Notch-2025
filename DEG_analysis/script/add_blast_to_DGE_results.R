# Libraries
library(stringr)
library(plyr)


### Functions

add_info = function(current_table, blast_result, specie_name){
  info = c()
  for (i in current_table$X..) {
    i = paste(i, "t1", sep = "")
    info = c(info,as.character(blast_result$Name[which(blast_result$TranscriptID == i)]))
  }
  current_table[specie_name] = info
  return(current_table)
}

add_cluster_info = function(current_table, cluster_correspondance){
  clusters = c()
  for (i in current_table$X..) {
    if(i %in% cluster_correspondance$GeneID){
      clusters = c(clusters,cluster_correspondance$cluster[which(cluster_correspondance$GeneID == i)])
    }else{
      clusters = c(clusters,0)
    }
  }
  current_table["Cluster"] = clusters
  return(current_table)
}

removequotes = function(list_to_remove){
  return(str_replace_all(list_to_remove, '\"', ""))
}

##

### Parameters

cluster_correspondance = read.table("cluster_to_geneid.tsv",
                                    sep = ",", header = T)

# Blast result table with Name and TranscriptID equivalent.
mouse.path = "Mus_musculus_table.tsv"
human.path = "Hsap_table.tsv"
aplysia.path = "Aplysia_clifornica_table.tsv"
drosophila.path = "Drosophila_table.tsv"
saccoglossus.path = "Saccoglossus_kowalevskii_table.tsv"

liste_directories = list.dirs()

mouse = read.table(file = mouse.path, sep = "\t", header = T, quote = "")
human = read.table(file = human.path, sep = "\t", header = T, quote = "")
aplysia = read.table(file = aplysia.path, sep = "\t", header = T, quote = "")
drosophila = read.table(file = drosophila.path, sep = "\t", header = T, quote = "")
saccoglossus = read.table(file = saccoglossus.path, sep = "\t", header = T, quote = "")

##




### Main



for (i in liste_directories) {
  if (i != ".") {
    current_table = read.table(file = paste0(i,"/DEG_result_",substr(i,3,50),"_FC_and_FDR_filter.tsv"),header = T, sep = "\t", quote = "", stringsAsFactors = F)
    current_table$X.. = removequotes(current_table$X..)
    rownames(current_table) = current_table$X..
    current_table = add_info(current_table, mouse, "Mouse")
    current_table = add_info(current_table, human, "Human")
    current_table = add_info(current_table, drosophila, "Drosophila")
    current_table = add_info(current_table, aplysia, "Aplysia clifornica")
    current_table = add_info(current_table, saccoglossus, "Saccoglossus kowalevskii")
    current_table = add_cluster_info(current_table, cluster_correspondance)
    write.table(x = current_table,
                file = paste0(i,"/DEG_result_",substr(i,3,50),"_FC_and_FDR_filter_annotated.tsv"),
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
}





##