# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )
prefix = "/Users/minoru/Code/bhklab/ICBCuration/work_dir/ICB_Fumet2/"

rule all:
    input:
        prefix + "processed/ICB_Fumet2_expr_gene_tpm.tsv",
        prefix + "processed/ICB_Fumet2_expr_gene_counts.tsv",
        prefix + "processed/ICB_Fumet2_expr_isoform_tpm.tsv",
        prefix + "processed/ICB_Fumet2_expr_isoform_counts.tsv",
        prefix + "processed/ICB_Fumet2_metadata.tsv"

rule format_expr_data:
    input:
        prefix + "annotation/Gencode.v40.annotation.RData",
        prefix + "download/ICB_Fumet2_RNAseq.zip",
        prefix + "download/filereport_read_run_PRJNA786565_tsv.txt"
    output:
        prefix + "processed/ICB_Fumet2_expr_gene_tpm.tsv",
        prefix + "processed/ICB_Fumet2_expr_gene_counts.tsv",
        prefix + "processed/ICB_Fumet2_expr_isoform_tpm.tsv",
        prefix + "processed/ICB_Fumet2_expr_isoform_counts.tsv"
    shell:
        """
        Rscript scripts/format_expr_data.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation \
        """

rule format_clin_data:
    input:
        prefix + "annotation/curation_tissue.csv",
        prefix + "annotation/curation_drug.csv",
        prefix + "download/GSE190265_series_matrix.txt.gz",
        prefix + "download/GSE190265_samples_info_France3.csv.gz",
    output:
        prefix + "processed/ICB_Fumet2_metadata.tsv"
    shell:
        """
        Rscript scripts/format_clin_data.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_annotation:
    output:
        prefix + "annotation/Gencode.v40.annotation.RData",
        prefix + "annotation/curation_tissue.csv",
        prefix + "annotation/curation_drug.csv"
    shell:
        """
        wget https://github.com/BHKLAB-DataProcessing/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        wget https://github.com/BHKLAB-DataProcessing/ICB_Common/raw/main/data/curation_drug.csv -O {prefix}annotation/curation_drug.csv
        wget https://github.com/BHKLAB-DataProcessing/ICB_Common/raw/main/data/curation_tissue.csv -O {prefix}annotation/curation_tissue.csv 
        """

rule download_data:
    output:
        prefix + "download/GSE190265_series_matrix.txt.gz",
        prefix + "download/GSE190265_samples_info_France3.csv.gz",
        prefix + "download/ICB_Fumet2_RNAseq.zip",
        prefix + "download/filereport_read_run_PRJNA786565_tsv.txt"
    shell:
        """
        wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz' \
            -O {prefix}download/GSE190265_series_matrix.txt.gz

        wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190265&format=file&file=GSE190265%5Fsamples%5Finfo%5FFrance3%2Ecsv%2Egz' \
            -O {prefix}download/GSE190265_samples_info_France3.csv.gz

        wget 'https://zenodo.org/record/7158382/files/ICB_Fumet2_RNAseq.zip?download=1' \
            -O {prefix}download/ICB_Fumet2_RNAseq.zip

        wget 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA786565&result=read_run&fields=sample_accession,sample_title&format=tsv&download=true&limit=0' \
            -O {prefix}download/filereport_read_run_PRJNA786565_tsv.txt
        """
