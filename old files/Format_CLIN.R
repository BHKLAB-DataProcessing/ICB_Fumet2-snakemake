library(data.table)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t")
selected_cols <-  c( "my_id_rna" , "Sex" , "Age" , "Histology" , "Tumor.stage" , "bestRECIST" , "PFS" , "evtPFS" , "OS" , "evtOS" )
clin = cbind( clin_original[ , selected_cols ] , "Lung" , "PD-1/PD-L1" , NA , NA , NA , NA )
colnames(clin) = c( "patient" , "sex" , "age"  ,"histo"  , "stage" , "recist" , "t.pfs" , "pfs" , "t.os" , "os" , "primary" , "drug_type" , "response.other.info" , "response" , "dna" , "rna")

clin$rna = "tpm"
clin$response = Get_Response( data=clin )

clin$stage = ifelse( clin$stage %in% 3 , "III" , 
				ifelse( clin$stage %in% 4 , "IV" , 
				ifelse( clin$stage %in% c(0,1,2) , "I/II" , NA ) )) 

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]
clin$patient = sapply( clin$patient , function(x){ paste( unlist( strsplit( x , "-" , fixed=TRUE)) , collapse="." ) } )
clin$patient = paste( "P" , clin$patient , sep="" )

clin_original$my_id_rna = paste( "P" , clin_original$my_id_rna , sep="" )
clin_original$my_id_rna = str_replace_all(clin_original$my_id_rna, '-', '.')
clin = format_clin_data(clin_original, "my_id_rna", selected_cols, clin)

annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Fumet.2', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, treatmentid='', .after='tissueid')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

