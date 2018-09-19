library(biomaRt)
library(DT)
listMarts() 

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
all_datasets <- listDatasets(ensembl)

ensembl = useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
filters = listFilters(ensembl)
