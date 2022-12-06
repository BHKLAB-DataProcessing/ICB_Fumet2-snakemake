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
colnames(clin)[colnames(clin) == "title"] <- "patient"
clin$patient <- str_replace_all(clin$patient, "\\W", "_")
rownames(clin) <- clin$patient

pfs_df <- read.table(file.path(input_dir, "GSE190265_samples_info_France3.csv.gz"), header = TRUE, sep = ";")
pfs_df$sample <- str_replace_all(pfs_df$sample, "\\W", "_")
rownames(pfs_df) <- pfs_df$sample
pfs_df$sample <- NULL
pfs_df <- pfs_df[rownames(pfs_df) %in% rownames(clin), ]

clin <- clin[rownames(clin) %in% rownames(pfs_df), ]
clin <- clin[order(rownames(clin)), ]
pfs_df <- pfs_df[order(rownames(pfs_df)), ]
clin <- cbind(clin, pfs_df)

# TODO: Format the clinical data with common columns etc.
selected_cols <- c("patient", "tissue:ch1", "time_PFS", "evtPFS")
remaining_cols <- colnames(clin)[!colnames(clin) %in% selected_cols]
clin <- clin[, c(selected_cols, remaining_cols)]

colnames(clin)[colnames(clin) %in% selected_cols] <- c("patient", "tissueid", "t.pfs", "pfs")
clin$tissueid <- "Lung"
clin <- add_column(clin, response = NA, response.other.info = NA, recist = NA, .after = "tissueid")
clin$t.pfs <- as.numeric(clin$t.pfs)
clin$pfs <- as.numeric(clin$pfs)
clin$response <- Get_Response(clin)

clin <- add_column(
  clin,
  sex = NA,
  age = NA,
  primary = "Lung",
  stage = NA,
  histo = NA,
  treatmentid = "",
  drug_type = "anti-PD-1/anti-PD-L1",
  dna = NA,
  rna = NA,
  t.os = NA,
  os = NA,
  survival_unit = "month",
  survival_type = NA,
  .after = "patient"
)

clin <- clin[, c(
  "patient", "sex", "age", "primary", "treatmentid", "tissueid", "histo", "stage", "response.other.info", "recist",
  "response", "drug_type", "dna", "rna", "t.pfs", "pfs", "t.os", "os", "survival_unit", "survival_type",
  remaining_cols
)]

clin$patient[!is.na(as.numeric(clin$patient))] <- paste0("p", clin$patient[!is.na(as.numeric(clin$patient))])
rownames(clin) <- clin$patient
write.table(clin, file = file.path(output_dir, "CLIN.csv"), quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)

case <- as.data.frame(cbind(clin$patient, rep(0, length(clin$patient)), rep(0, length(clin$patient)), rep(1, length(clin$patient))))
colnames(case) <- c("patient", "snv", "cna", "expr")
rownames(case) <- clin$patient
write.table(case, file = file.path(output_dir, "cased_sequenced.csv"), quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)
