rm(list = ls())
options(warn=-1)
source("build_refDB.R")
source("BAGS_blast.R")

path_to_query_p = "DATASET/USER_queryp.faa"
path_to_query_n = "DATASET/USER_query.fna"
path_to_db_p = "DATASET/AA/proteins.faa"
path_to_db = "DATASET/DNA/genes.fna"
###INPUT
#parameters defined by user (from Rshiny)
Nhits=10
minPerIden = 50
minPeralg= 50
E_value=1e-3
Btype="blastn" #blastn
C_pus=8
usingtext = T
#TEMPORAL
text=">User1\nGGAGAGAGGGAGACAGGCAGGGAGACAACCAGGCAGGGAGAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGGAGACAGACAGGCAGGGAGAGAGGGAGACAGGCAGGGAGACAACCAGGCAGGGAGAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGGAGACAGACAGGCAGGGAGACAGGGAGACAGGG"
#text=">pufM\nQVRGHTDYGVTEGTNSEPRYNLTATSHLLGRIGDAQLGPLYLGWSGVASLICGFIAIEIIGLNMWASVNWDPVQFVRQLPWLALEPPLPKYGLKVLPPLNEGGWWLI"
#text=">wzc\nMNETMTTPQNRSSVDNSSDEIDLGKLLGILLNAKWIILVTTFFFAVGGVAVALLSTPIYKADALIQIEQKSSGGISSLVGDMGELFSQESSATTEIEIIKSRMILGDTVDKFNLTTVAEPKYLPVIGKGLARIAGKVNQIEISRYTVPEYAQEMKHTLVVLDAEKKTYQLVRGDEQVVLKGVAGELAKNDGYELFVTELRSHNDQEFAIGQRSRLEAIEWLKQNLAISERGKQTGILQLSFEGENRKQIGELLNHISQTYFLQNVERNSAEAEKSLTFLKGHLPDIKTSLTTAEDTLNRFRQDNESIDLGLEAKSTLDVMVKLEAQLNELTFKESEISQRFTKDHPAYRSLLDKRETLLKERERLNQQVQKLPKTQREVLRMTRDVEVNQQIYIQLLNKVQELSIIKAGTVGNVRILDEAQSYAKPVKPKKPLIVVLATLLGGMLSVALVLVKAALHRGVENPDEIEQIGLSVYASVPKSNLQLELANKLARKKRNTDLTLLAESNPADLSIEALRGLRTSLHFAMMEAKNNVLMISGPAPGIGKSFVSTNFAAVAAKTGQKVLLIDADMRKGYLQQCFGLNWENGLSDLLSGKVTRDVAVQSAKVENLDIITRGQVPPNPSELLMHPRFKELVDWASEHYDLVIIDTPPVLAVTDPSIVGAIAGTTLMVARFGQNTVKEIDVARSRFEQAGIEVKGVILNAIEKKASSSYGYGYYNYSYGESNKA"
#REQUIRED DATA

 ###EXPECTED OUTPUT

#Code to be integrated in Shiny with corresponding outputs

#example Selected_hit=Query_Hits[[1]]$SubjectID[1] #one or several hits per Query
#Selected_hit="P4201_120_k141_443270_1"
#Selected_hit=TEST2$GeneID
#Annotation_per_hit <- get_annotation(Selected_hit) 
#Localisation_in_the_Baltic_Sea<-get_localisation(Selected_hit)

#KEGG_id=unique(Annotation_per_hit$KEGG)[1] #"K19584" #"K09490" #"K19584"
#Output_KEG<-expand_kegg(KEGG_id)


##
##add code that check query input: 1) it is fasta 2) is AA when selecting blastp 3) no more than ?10 sequences
#library(shiny)


ui=fluidPage(
  titlePanel("BAGS.v1: BAltic Gene Set gene catalogue"),
  fluidRow(
    column(12,
           includeMarkdown("abstract.md")
    )
  ),
  fluidRow(column(5,
                  # textInput("text_input",label="textInput",value=""),
                  textAreaInput("text_input",label="textAreaInput", value = ""),
                  # tableOutput("table_output"),
                  # textOutput("text_output")
                  tableOutput("table_output")#,
                  #leafletOutput("mymap")
  ))
)

server=function(input, output) {
  output$table_output <- renderTable({
    Query_Hits <- Hit_blast(Btype,usingtext,text, C_pus, Nhits, E_value, minPerIden, minPeralg)
    
    if (nrow(Query_Hits[[1]]) < 1) {print("No hit found")} else {return(Query_Hits[[1]])}
    
  })
  
  ###### map
  # output$table_output <- renderTable({
  output$mymap <- renderLeaflet({
    #LOCALISSTION TABLE
  })
  
}

shinyApp(ui=ui,server=server)