rm(list = ls())
options(warn=-1)

library(tidyverse)
library(rBLAST)
library(ggplot2)
library(leaflet)
library(sp)
library(data.table)

#parameters defined by user
N_hits=10  #Rshiny input
minPer_Iden = 90 #Rshiny input
minPer_alg= 90 #Rshiny input
Evalue=1e-3 #Rshiny input
B_type="blastn" #blastp, #in the future blastx, tblastn #Rshiny input in a if (DNA, Protein)
Cpus=4

path_to_query_p="DATASET/USER_queryp.faa" # Path TO save when user provide the file 
path_to_query="DATASET/USER_query.fna" #"DATASET/queryn.fna" # Path TO save when user provide the file 

using_text = TRUE # dependend on user option
if(B_type == "blastn" & using_text){
  #here we read the input from Shiny input text
  #Function checking text is in fasta format, type, N sequences 
  text=">User1\nATCTATGATTTGGCGCAAGACGACTGGTATGCAAGCTCCCAATTAGAAGAAATATACAACCTTATATCAACATCAGAAAAGATACCTCAATTGGACGGTTGGGCAATTAACAAGACTGACAGCAGATGTGCCAAGTCTTGA\n>User2\nCGCGCGCGGTGCTACGGCGGGGCGGGGCCGGTGGGGACCGGTCCCATACAACGCCCTACTTACACAAACATAGTGTCCGTCATGAAGTGTGCTGTTCCAGGTACAGCTGGCGTGGAGCAGCGGCACATGCCATTTAAAGTACATGCGTCCATAAGTAGCATGAGCAAATTGTATGTACAGCTAAAGCCGCGTATCATAGGTAACCCGTATGAGGTATGTACATGA"
  fwrite(as.list(text), file=path_to_query, quote = F)
} else { #save input file uploaded by the user to path_to_query 
  } 
 
if(B_type == "blastp" & using_text){
  #here we read the input from Shiny input text
  text=">User1\nMPDTGDIRTSSGTFLSRHQDPDGVLEKIERRIADATHVPYTHGEPFNVLRYTPGQKYDSHYDTFDPVSYGVQTSQRVASFLLYLTDVEEGGETHFPLEGRNGLERLKNIDYKSCDGGLLVRPRAGDALLFWNVFPNATFDKHALHGGCPVTKGEKWVATKWIRDKSFGKA*\n>User2\nMRKRLAQFRKLGRTPQHKWAMLRNMVTSLIKHERIKTTLPKAKELRHLADKVIGYAKKPNQTYGKQLALQVVTEKPMVTKLMEVLGPRYADRDGGYTRVLKLSRPRRGDNAPMAIIEYVDRPGEVRAARVPE"
  fwrite(as.list(text), file=path_to_query_p, quote = F)
} else { #save input file uploaded by the user to path_to_query_p 
} 

# Remove this part when data is alredy available
## WE should provide these paths. TO be DONE only once
path_to_db_p = "BAGS/AA/bags_v1_proteins.faa" #"DATASET/AA/proteins.faa"
path_to_db = "BAGS/DNA/bags_v1_genes.fna" # "DATASET/DNA/genes.fna"
path_to_taxonomy_db = "BAGS/bags_v1_taxonomy.tsv" #DATASET/taxonomy.tsv
path_to_annot_db ="BAGS/bags_v1_annotations.tsv" # DATASET/annotations.tsv
path_to_metadata = "DATASET/BS_MAGv2_sample_metadata.tsv"

#Functions
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
#dash_split
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


#if we go for "DATASET/big_tbl.rds", then remove !file.exists("DATASET/ann.rds") & !file.exists("DATASET/tax.rds" & and corresponding save/readRDS
#if ( !file.exists("DATASET/ann.rds") & !file.exists("DATASET/tax.rds" & !file.exists("DATASET/big_tbl.rds") ) ) {
if ( !file.exists("DATASET/big_tbl.rds") ) {
    
  tax <- read_tsv(path_to_taxonomy_db, show_col_types = FALSE, col_names=F) 
  tax <- tax[-1,]
  colnames(tax) <- c("GeneID", "assigned_taxonomy")
  
  ann <- read_tsv(path_to_annot_db, show_col_types = FALSE, col_names=F)
  colnames(ann) <-ann[1,]
  colnames(ann)[1] <- "GeneID"
  ann <-ann[-1,]
  # KEEG_list(ann$EggNOG_annotation, ",")
  ann <- ann %>% rowwise() %>% mutate(KEGG = KEEG_list(EggNOG_annotation, ","))
  big_tbl <- tax %>%  full_join(ann, by = "GeneID")
  
#  saveRDS(ann, file="DATASET/ann.rds")
#  saveRDS(tax, file="DATASET/tax.rds")
  saveRDS(big_tbl, file="DATASET/big_tbl.rds")

} else {
#  ann <- readRDS("DATASET/ann.rds")
#  tax <- readRDS("DATASET/tax.rds")
  big_tbl <- readRDS("DATASET/big_tbl.rds")
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
} else {
  ref_coords<-readRDS("DATASET/ref_coords.rds")
}
####


##add code that check query input: 1) it is fasta 2) is AA when selecting blastp 3) no more than ?10 sequences
if (B_type == "blastp") {
  blp <- blast(db=path_to_db_p, type = "blastp")
query_p<-readAAStringSet(path_to_query_p)
#RESULTS_p=run_blast(blp, query_p, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg)
RESULTS=run_blast(blp, query_p, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
}


if (B_type == "blastn") {
  bln <- blast(db=path_to_db, type = "blastn")
  query_n<-readDNAStringSet(path_to_query)
  RESULTS=run_blast(bln, query_n, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
}


###ADD code handling if no hit

query_name=names(RESULTS)
Selected_hit=vector(length = length(query_name))
for (q in 1:length(query_name)){
Selected_hit[q] <- RESULTS[[query_name[q]]]$SubjectID #Selected by user after blast results - one or several
}
Annotation_and_taxonomy_result <- big_tbl %>% filter(GeneID %in% Selected_hit) #OUTPUT Annotation PER HIT/sequence


###LOCALISATION

selectes=sapply(Selected_hit, dash_split) 

hit_seqs <- ref_coords %>% filter(ref_id %in% selectes) 
link=as.data.frame(matrix(nrow=length(query_name), ncol=2))
names(link)=c("ref_id", "query_id")
for (j in 1:length(query_name) ) {
  k=query_name[j]
 link[j,] = c(dash_split(RESULTS[[k]]$SubjectID), RESULTS[[k]]$QueryID)
}

hit_seqs <- hit_seqs %>% left_join(link)
coordinates(hit_seqs) <- ~Lon + Lat

#For time serie data points, only one is shown since lon and lat are the same
map<-leaflet(hit_seqs) %>% addTiles() %>% addCircles() %>% 
  addMarkers(data=hit_seqs, 
             label=paste("Query: ",hit_seqs$query_id,"\n",
                         "Ref.: ",hit_seqs$ref_id,"\n",
                         "Sal.: ",hit_seqs$Sal,"\n",
                         "Temp.: ",hit_seqs$Temp,"\n",
                         "Depth: ",hit_seqs$Depth, sep="")#,
             #labelOptions(textOnly = T) 
             )




