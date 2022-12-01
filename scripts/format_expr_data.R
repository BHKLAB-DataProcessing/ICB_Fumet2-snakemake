# wget 'https://zenodo.org/record/7158382/files/ICB_Fumet2_RNAseq.zip?download=1' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet2/ICB_Fumet2_RNAseq.zip
# wget 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA786565&result=read_run&fields=sample_title&format=tsv&download=true&limit=0' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet2/filereport_read_run_PRJNA786565_tsv.txt

library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr_list <- readRDS(file.path(input_dir, "expr_list.rds"))
clin <- read.table(file.path(output_dir, 'CLIN.csv'), sep=';')

for(assay_name in names(expr_list)){
  df <- expr_list[[assay_name]]
  df <- df[, colnames(df)[colnames(df) %in% rownames(clin)]]
  write.table(
    df,
    file= file.path(output_dir, paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')),
    quote=FALSE,
    sep=";",
    col.names=TRUE,
    row.names=TRUE
  )
}
