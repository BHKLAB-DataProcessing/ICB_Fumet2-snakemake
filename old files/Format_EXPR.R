library(data.table)
library(dplyr)
library(R.utils)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr = as.data.frame( fread( file.path(input_dir, 'EXPR.csv.gz') , stringsAsFactors=FALSE  , sep=";" , dec="," ))
colnames(expr) = expr[1,]
expr = expr[-1,]

#########################################################################################################
#########################################################################################################
## Compute Gene level TPM
transcript = expr[ , 1 ]
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

geneID = as.matrix( getBM(attributes=c("hgnc_symbol","ensembl_gene_id" , "ensembl_transcript_id"), filters="ensembl_transcript_id",values=transcript, mart=human) )
rownames(geneID) = geneID[ , "ensembl_transcript_id" ]
geneID = geneID[ !geneID[ , "hgnc_symbol" ] %in% "" , ]

intersect = intersect( transcript , rownames(geneID) )
expr = expr[ expr[ , 1 ] %in% intersect , ]
geneID = geneID[ intersect , ]

data = as.data.frame(cbind( expr[,1], geneID[ expr[,1] , "hgnc_symbol" ], expr[,-1]))
colnames(data) = c( "transcript_id" , "gene_id" , colnames(expr)[-1] )
data = data[!is.na(data$gene),]

sample = colnames(data)[-(1:2)]
tpm = NULL
for(i in 1:length(sample)){
  d = data[ , c( "transcript_id" , "gene_id" , sample[i] ) ]
  colnames(d) = c( "transcript_id" , "gene_id" ,"TPM" )
  d$TPM = as.numeric(as.character(d$TPM))
  gene_tpm <- group_by(d, gene_id) %>% summarize(TPM = sum(TPM , na.rm=TRUE))
  
  tpm = cbind( tpm , as.data.frame( gene_tpm )$TPM )
}
rownames(tpm) = unique(sort(data$gene_id))
colnames(tpm) = sample

#########################################################################################################
#########################################################################################################

colnames(tpm) = sapply( colnames(tpm) , function(x){ paste( unlist( strsplit( x , "-" , fixed=TRUE)) , collapse="." ) } )
colnames(tpm) = paste( "P" , colnames(tpm) , sep="" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
tpm = log2( tpm[ , colnames(tpm) %in% case[ case$expr %in% 1 , ]$patient ] + 0.001 )

write.table( tpm , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
