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
} # else { #save input file uploaded by the user to path_to_query } 
 
if(B_type == "blastp" & using_text){
  #here we read the input from Shiny input text
  text=">User1\nMPDTGDIRTSSGTFLSRHQDPDGVLEKIERRIADATHVPYTHGEPFNVLRYTPGQKYDSHYDTFDPVSYGVQTSQRVASFLLYLTDVEEGGETHFPLEGRNGLERLKNIDYKSCDGGLLVRPRAGDALLFWNVFPNATFDKHALHGGCPVTKGEKWVATKWIRDKSFGKA*\n>User2\nMRKRLAQFRKLGRTPQHKWAMLRNMVTSLIKHERIKTTLPKAKELRHLADKVIGYAKKPNQTYGKQLALQVVTEKPMVTKLMEVLGPRYADRDGGYTRVLKLSRPRRGDNAPMAIIEYVDRPGEVRAARVPE"
  fwrite(as.list(text), file=path_to_query_p, quote = F)
} #else { #save input file uploaded by the user to path_to_query_p } 

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
query_p<-readAAStringSet(path_to_query_p)  #put it earlier
#RESULTS_p=run_blast(blp, query_p, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg)
RESULTS=run_blast(blp, query_p, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
}


if (B_type == "blastn") {
  bln <- blast(db=path_to_db, type = "blastn")
  query_n<-readDNAStringSet(path_to_query)
  RESULTS=run_blast(bln, query_n, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
}


###ADD code handling if no hit
# one hit at the time
query_name=names(RESULTS)
Selected_hit=vector(length = length(query_name))
for (q in 1:length(query_name)){
Selected_hit[q] <- RESULTS[[query_name[q]]]$SubjectID #Selected by user after blast results - one or several
}
Annotation_and_taxonomy_result <- big_tbl %>% filter(GeneID %in% Selected_hit) #OUTPUT Annotation PER HIT/sequence

####KEGG pathway list
library(tidyverse)
library(magrittr)
#install.packages("KEGGREST")
library(KEGGREST)

# Get the KEGG map for K09490
kegg_id<-Annotation_and_taxonomy_result$KEGG
kegg_info <- keggGet(kegg_id)
kegg_info <- kegg_info %>%
  transpose() %>%
  as_tibble() 
kegg_map <- kegg_info$PATHWAY[[1]]%>% as_tibble()
print(kegg_map)

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

#add KEG tables (no map)
######################## plot of logo in r####################################
library(png)
library(ggplot2)
library(tidyverse)
library(data.table)
getTextImage <- function(text) {
  filename <- tempfile()
  png(filename = filename)
  plot.new()
  cex <- 1
  repeat {
    if (strwidth(text, cex = 2*cex) > 1) break
    if (strheight(text, cex = 2*cex) > 1) break 
    cex <- 2*cex
  }
  text(0.5, 0.5, text, cex = cex)
  dev.off()
  image <- readPNG(filename)
  unlink(filename)    # clean up file
  if (length(dim(image)) == 3)
    image <- image[,,1] # just keep one channel
  image
}

randomText <- function(n, text) {
  image <- getTextImage(text)
  nx <- dim(image)[1]
  ny <- dim(image)[2]
  hits <- 0
  x <- y <- numeric(n)
  while (hits < n) {
    tryx <- runif(1)
    tryy <- runif(1)
    keep <- image[round((nx-1)*tryx + 1), round((ny-1)*tryy + 1)] == 0
    if (keep) {
      hits <- hits + 1
      # Need to rotate so it looks good
      x[hits] <- tryy
      y[hits] <- 1 - tryx
    }
  }
  cbind(x, y)
}

data<-data.frame(randomText(1000, "BAGS"))


p1 <- ggplot(data) +
  geom_point(aes(x = x, y = y), colour = "green", shape=23,size = 2) +
  labs( subtitle = "BAltic Gene Set gene catalogue")

p1<-p1+theme_void()+
theme(plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="green"))
print(p1)
############################################
# Just to prove that the maths works, plot the hexagon described by unit length
path <- data.frame(Pro = c(1,1,0,0,0,1,1), Cla = c(0,1,1,1,0,0,0), Neu = c(0,0,0,1,1,1,0))
pd<-path %>%
  mutate(S1_x = Pro*cos(pi/2), S1_y = Pro*sin(pi/2), 
         S2_x = Cla*cos(pi/6), S2_y = -Cla*sin(pi/6), 
         S3_x = -Neu*cos(pi/6), S3_y = -Neu*sin(pi/6)) %>% 
  mutate(x = S1_x + S2_x + S3_x, y = S1_y + S2_y + S3_y)
pd$z<-c(1.2,0.8, 1.2,0.8,1.2, 0.8, 1.2)
pd$w<-seq(0.8,2.4,by=0.25)
#ggplot(aes(x = x, y=y))+geom_path()
dt.triangle2 <- data.table(group = c(1,1,1), polygon.x = c(1.5,2.54,2.54), polygon.y = c(0.4,0.4 ,1))
dt.triangle1 <- data.table(group = c(1,1,1), polygon.x = c(0.5,0.5,1.5), polygon.y = c(1 ,0.4 ,0.4))

data<-data.frame(randomText(1000, "RaukR"))
p2 <- ggplot(data) +
  geom_point(aes(x = x*2.9, y = y*2.9), colour = "black", shape=23,size = 1.5) 

p2<-p2+geom_bar(data=pd,aes(x=w,y=z), fill="darkcyan",stat = "identity", position = "dodge",width=.25)+
  geom_polygon(data = dt.triangle1,aes(x=polygon.x,y=polygon.y,group=group),fill="white")+
  geom_polygon(data = dt.triangle2,aes(x=polygon.x,y=polygon.y,group=group),fill="white")+
  geom_rect(data=pd,aes(xmin=min(x+1.54), xmax=max(x+1.58), ymin=0, ymax=0.4), fill="white") +
  geom_path(data=pd,aes(x=x+1.54,y=y+1.4), color="cyan3",linewidth=2)+
  annotate(geom="text", x=1.56, y=0.65, label="2023",
           color="black",size = 6)+
theme_void()
print(p2)


