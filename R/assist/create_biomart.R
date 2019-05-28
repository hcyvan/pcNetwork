# This script is used to generate biomart
library(biomaRt)

listMarts()
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl)
listAttributes(ensembl)->a
searchAttributes(mart = ensembl, pattern = "GO")[,c(1,2)]

biomart <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene', 'hgnc_symbol','external_gene_name', 'gene_biotype', 'start_position', 'end_position', 'transcript_start','transcript_end','transcription_start_site','chromosome_name','go_id','name_1006','definition_1006'),mart = ensembl)

saveRDS(biomart, file='./support/biomart.rds')
biomart<-readRDS('./support/biomart.rds')

write.csv(biomart, file='./data/biomart.csv')

########################
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
searchAttributes(mart = ensembl, pattern = "GO")[,c(1,2)]

biomart37 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene', 'hgnc_symbol','external_gene_name', 'gene_biotype', 'start_position', 'end_position', 'transcript_start','transcript_end','transcription_start_site','chromosome_name','name_1006','definition_1006'),mart = ensembl)
