rm(list = ls())
options(warn=-1)

list_of_packages <- c("tidyverse", "data.table")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) { suppressPackageStartupMessages(install.packages(s))}
}

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE)) }

KEEG_list<-function(x ,sep){
  my_list <- str_split(x, ",")%>% unlist()
  element<-my_list[grep("[K]\\d",  my_list)]
  text<-str_split(element," ")%>% unlist()
  if (length(text) > 0) {
    result <- c(text[grep("[K]\\d", text)])[1]
  }
  else {
    result <- NA
  }
  return(result)
}

path_to_db_p = "DATASET/AA/proteins.faa"
path_to_db = "DATASET/DNA/genes.fna"
path_to_taxonomy_db = "DATASET/taxonomy.tsv"
path_to_annot_db ="DATASET/annotations.tsv"
path_to_metadata = "DATASET/BS_MAGv2_sample_metadata.tsv"

if ( !file.exists("DATASET/big_tbl.rds") ) {

  tax <- read_tsv(path_to_taxonomy_db, show_col_types = FALSE, col_names=F)
  tax <- tax[-1,]
  colnames(tax) <- c("GeneID", "assigned_taxonomy")

  ann <- read_tsv(path_to_annot_db, show_col_types = FALSE, col_names=F)
  colnames(ann) <-ann[1,]
  colnames(ann)[1] <- "GeneID"
  ann <-ann[-1,]

  ann <- ann %>% rowwise() %>% mutate(KEGG = KEEG_list(EggNOG_annotation, ","))
  big_tbl <- tax %>%  full_join(ann, by = "GeneID")

  saveRDS(big_tbl, file="DATASET/big_tbl.rds")

}

if ( !file.exists(paste(path_to_db, ".ndb", sep=""))) {
  makeblastdb(path_to_db, dbtype = "nucl")
}

if ( !file.exists(paste(path_to_db_p, ".pto", sep=""))) {
  makeblastdb(path_to_db_p, dbtype = "prot")
}


#get the coordinates
if (!file.exists("DATASET/ref_coords.rds")) {
ref_coords<-as.data.frame(t(read.table(path_to_metadata)))%>%
  rownames_to_column(.,var="ref_id")%>%
  select(ref_id, Lon, Lat,Sal,Depth,Temp ) %>%
  mutate(Lon=as.numeric(Lon), Lat=as.numeric(Lat))

saveRDS(ref_coords, file="DATASET/ref_coords.rds")
}
