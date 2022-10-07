# wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet2/GSE190265_series_matrix.txt.gz

library(data.table)
library(stringr)
library(GEOquery)

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


write.table(clin, file = file.path(output_dir, "ICB_Fumet2_metadata.tsv"), row.names = TRUE, col.names=TRUE, sep = "\t")
