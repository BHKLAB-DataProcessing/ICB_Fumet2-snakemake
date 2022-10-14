# wget 'https://zenodo.org/record/7158382/files/ICB_Fumet2_RNAseq.zip?download=1' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet2/ICB_Fumet2_RNAseq.zip
# wget 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA786565&result=read_run&fields=sample_title&format=tsv&download=true&limit=0' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet2/filereport_read_run_PRJNA786565_tsv.txt

library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

# create a sample mapping to format processed RNA-seq data.
seq_sample <- read.csv(file.path(input_dir, "filereport_read_run_PRJNA786565_tsv.txt"), sep = "\t")
seq_sample$sample_title <- str_replace_all(seq_sample$sample_title, "\\W", "_")

# EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R")
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))
dir.create(file.path(input_dir, "rnaseq"))
unzip(file.path(input_dir, "ICB_Fumet2_RNAseq.zip"), exdir = file.path(input_dir, "rnaseq"))

unlink(file.path(input_dir, "rnaseq", "__MACOSX"), recursive = TRUE)

# save expr_list.rds
process_kallisto_output(input_dir, tx2gene)

expr_list <- readRDS(file.path(input_dir, "expr_list.rds"))
for (name in names(expr_list)) {
  expr_df <- data.frame(expr_list[[name]])
  colnames(expr_df) <- unlist(lapply(colnames(expr_df), function(run_accession) {
    sampleid <- seq_sample$sample_title[seq_sample$run_accession == run_accession]
    if (!is.na(as.numeric(sampleid))) {
      sampleid <- paste0("p", sampleid)
    }
    return(sampleid)
  }))
  write.table(expr_df, file = file.path(output_dir, paste0("ICB_Fumet2_", name, ".tsv")), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
}
