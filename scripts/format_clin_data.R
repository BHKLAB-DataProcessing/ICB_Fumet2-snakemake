# wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet2/GSE190265_series_matrix.txt.gz

library(data.table)
library(stringr)
library(GEOquery)
library(tibble)

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

gunzip(file.path(input_dir, "GSE190265_series_matrix.txt.gz"))
clin <- getGEO(filename = file.path(input_dir, "GSE190265_series_matrix.txt"), destdir = input_dir)
clin <- pData(clin)
colnames(clin)[colnames(clin) == "title"] <- "patientid"
clin$patientid <- str_replace_all(clin$patientid, "\\W", '_')
rownames(clin) <- clin$patientid

pfs_df <- read.table(file.path(input_dir, 'GSE190265_samples_info_France3.csv.gz'), header=TRUE, sep=';')
pfs_df$sample <- str_replace_all(pfs_df$sample, "\\W", '_')
rownames(pfs_df) <- pfs_df$sample
pfs_df$sample <- NULL
pfs_df <- pfs_df[rownames(pfs_df) %in% rownames(clin), ]

clin <- clin[rownames(clin) %in% rownames(pfs_df), ]
clin <- clin[order(rownames(clin)), ]
pfs_df <- pfs_df[order(rownames(pfs_df)), ]
clin <- cbind(clin, pfs_df)

# TODO: Format the clinical data with common columns etc.
selected_cols <- c('patientid', 'tissue:ch1', 'time_PFS', 'evtPFS')
remaining_cols <- colnames(clin)[!colnames(clin) %in% selected_cols]
clin <- clin[, c(selected_cols, remaining_cols)]

colnames(clin)[colnames(clin) %in% selected_cols] <- c('patientid', 'tissueid', 't.pfs', 'pfs')
clin$tissueid <- 'Lung'
clin <- add_column(clin, response=NA, response.other.info=NA, recist=NA, .after='tissueid')
clin$t.pfs <- as.numeric(clin$t.pfs)
clin$pfs <- as.numeric(clin$pfs)
clin$response <- Get_Response(clin)

clin <- add_column(
  clin, 
  sex=NA,
  age=NA,
  stage=NA,
  histo=NA,
  treatmentid="anti-PD-1/anti-PD-L1",
  drug_type=NA,
  dna=NA,
  rna=NA,
  t.os=NA,
  os=NA,
  survival_unit='month',
  .after='patientid'
)

clin <- clin[, c("patientid", "sex", "age", "treatmentid", "tissueid", "histo", "stage", "response.other.info", "recist", "response", "drug_type", "dna", "rna", "t.pfs", "pfs", "t.os", "os", 'survival_unit', remaining_cols)]
colnames(clin)[colnames(clin) %in% c('t.pfs', 'pfs', "t.os", "os")] <- c("survival_time_pfs", "event_occurred_pfs", 'survival_time_os', 'event_occurred_os')

clin$patientid[!is.na(as.numeric(clin$patientid))] <- paste0('p', clin$patientid[!is.na(as.numeric(clin$patientid))])
rownames(clin) <- clin$patientid

write.table(clin, file = file.path(output_dir, "ICB_Fumet2_metadata.tsv"), row.names = TRUE, col.names=TRUE, sep = "\t")
