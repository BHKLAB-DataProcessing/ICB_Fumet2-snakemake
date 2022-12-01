library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
annot_dir <- args[2]

# create a sample mapping to format processed RNA-seq data.
seq_sample <- read.csv(file.path(work_dir, "filereport_read_run_PRJNA786565_tsv.txt"), sep = "\t")
seq_sample$sample_title <- str_replace_all(seq_sample$sample_title, "\\W", "_")
seq_sample$sample_title[!is.na(as.numeric(seq_sample$sample_title))] <- paste0('p', seq_sample$sample_title[!is.na(as.numeric(seq_sample$sample_title))])

# EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R")
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))
dir.create(file.path(work_dir, "rnaseq"))
unzip(file.path(work_dir, "ICB_Fumet2_RNAseq.zip"), exdir = file.path(work_dir, "rnaseq"))

unlink(file.path(work_dir, "rnaseq", "__MACOSX"), recursive = TRUE)

# save expr_list.rds
process_kallisto_output(work_dir, tx2gene)

expr_list <- readRDS(file.path(work_dir, 'expr_list.rds'))
for(name in names(expr_list)){
  df <- data.frame(expr_list[[name]])
  sampleids <- unlist(lapply(colnames(df), function(col){
    sampleid <- seq_sample[seq_sample$run_accession == col, c('sample_title')]
  }))
  colnames(df) <- sampleids
  expr_list[[name]] <- df
}

saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))