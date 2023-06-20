rm(list = ls())
options(warn=-1)

list_of_packages <- c("tidyverse", "data.table", "rBLAST", "leaflet", "sp", "magrittr", "KEGGREST")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) { suppressPackageStartupMessages(install.packages(s))}
}

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE)) }

#parameters defined by user (from Rshiny)
Nhits=10
minPerIden = 90
minPeralg= 90
E_value=1e-3
Btype="blastn" #blastp
C_pus=8
usingtext = F
#TEMPORAL
text=">User1\nGGAGAGAGGGAGACAGGCAGGGAGACAACCAGGCAGGGAGAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGGAGACAGACAGGCAGGGAGAGAGGGAGACAGGCAGGGAGACAACCAGGCAGGGAGAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGG"

#REQUIRED DATA
big_tbl <- readRDS("DATASET/big_tbl.rds")
ref_coords<-readRDS("DATASET/ref_coords.rds")
#Functions
dash_split <- function(x){
  unlist(strsplit(x,"_k"))[1]
}
run_blast<-function(bl, query,cpus, hits, evalue, minIden,minalg, type) {
  lista=query@ranges@NAMES
  if (type == "blastn") {blast_arg=paste("-num_threads", cpus,"-max_target_seqs", hits, "-evalue", evalue, "-soft_masking false","-dust no") }
  if (type == "blastp") {blast_arg=paste("-num_threads", cpus,"-max_target_seqs", hits, "-evalue", evalue, "-soft_masking false") }
  blast_table <- vector("list", length =length(query))
  names(blast_table) <- lista
  for (i in 1:length(query)) {
    nam=lista[i]
    blast_table[[nam]] <- predict(bl, query[i,], BLAST_args = blast_arg)
    LEN=query[i,]@ranges@width
    blast_table[[nam]] <- blast_table[[nam]] %>% mutate(Perc.Alignment = Alignment.Length*100/LEN) %>% filter(Perc.Ident > minIden & Perc.Alignment > minalg  )
  }
  return(blast_table)
}
Hit_blast <- function(B_type,using_text = FALSE,text="", Cpus, N_hits, Evalue, minPer_Iden, minPer_alg) {
  #REQUIRED INFORMATION - Defaults (no shiny input required)
  path_to_db_p = "DATASET/AA/proteins.faa"
  path_to_db = "DATASET/DNA/genes.fna"
  path_to_query_p = "DATASET/USER_queryp.faa"
  path_to_query_n = "DATASET/USER_query.fna"
  
if (B_type == "blastn") { path_to_query=path_to_query_n }
if (B_type == "blastp") { path_to_query=path_to_query_p}

if (using_text){
  fwrite(as.list(text), file=path_to_query, quote = F)
} 


##add code that check query input: 1) it is fasta 2) is AA when selecting blastp 3) no more than ?10 sequences
if (B_type == "blastp") {
  blp <- blast(db=path_to_db_p, type = "blastp")
query_p<-readAAStringSet(path_to_query)
RESULTS=run_blast(blp, query_p, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
}


if (B_type == "blastn") {
  bln <- blast(db=path_to_db, type = "blastn")
  query_n<-readDNAStringSet(path_to_query)
  RESULTS=run_blast(bln, query_n, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
}
return(RESULTS)  
}
get_annotation <- function(selected_hit) {
Annotation_and_taxonomy_result <- big_tbl %>% filter(GeneID %in% selected_hit) 
}
get_localisation <- function(selected_hit) {
  #selected_hit=Selected_hit
  selectes=sapply(selected_hit, dash_split) 
 hit_seqs <- ref_coords %>% filter(ref_id %in% selectes) 
 coordinates(hit_seqs) <- ~Lon + Lat
  
  #For time serie data points, only one is shown since lon and lat are the same
  map<-leaflet(hit_seqs) %>% addTiles() %>% addCircles() %>% 
    addMarkers(data=hit_seqs, 
               label=paste("Ref.: ",hit_seqs$ref_id,
                           "Sal.: ",hit_seqs$Sal,
                           "Temp.: ",hit_seqs$Temp,
                           "Depth: ",hit_seqs$Depth, sep="\n")) 

  return(map)
}
expand_kegg <- function(kegg_id){
  kegg_info <- keggGet(kegg_id)
  try(kegg_BRITE <- kegg_info[[1]]$BRITE %>% as_tibble(), silent=T)
  try(kegg_PATHWAY <-kegg_info[[1]]$PATHWAY %>% as_tibble(), silent = T)
  df=list()
  if (exists("kegg_BRITE")) { df[["BRITE"]] = kegg_BRITE }
  if (exists("kegg_PATHWAY")) { df[["PATHWAY"]] = kegg_PATHWAY }
  return(df)
}


#Code to be integrated in Shiny with corresponding outputs
Query_Hits <- Hit_blast(Btype,usingtext,text, C_pus, Nhits, E_value, minPerIden, minPeralg)
#example Selected_hit=Query_Hits[[1]]$SubjectID[1] #one or several hits per Query
Annotation_per_hit <- get_annotation(Selected_hit) 
Localisation_in_the_Baltic_Sea<-get_localisation(Selected_hit)

KEGG_id="K19584" #"K09490" #"K19584"
Output_KEG<-expand_kegg(KEGG_id)
