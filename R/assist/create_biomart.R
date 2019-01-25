# This script is used to generate biomart
library(biomaRt)

listMarts()
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl)
listAttributes(ensembl)
searchAttributes(mart = ensembl, pattern = "end")[,c(1,2)]

biomart <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene', 'hgnc_symbol', 'gene_biotype', 'start_position', 'end_position', 'transcript_start','transcript_end','transcription_start_site','chromosome_name'),mart = ensembl)

saveRDS(biomart, file='./support/biomart.rds')
