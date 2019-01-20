library(biomaRt)
listMarts()
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl)
listAttributes(ensembl)
searchAttributes(mart = ensembl, pattern = "ensembl")[,c(1)]

biomart2 <- getBM(attributes = c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'gene_biotype',
                                 'start_position', 'end_position','chromosome_name'),mart = ensembl)
biomart3 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene', 'hgnc_symbol', 'gene_biotype', 'start_position', 'end_position','chromosome_name'),mart = ensembl)

sequence2 <- getSequence(id = biomart2$ensembl_gene_id,
                         type = 'ensembl_gene_id',
                         seqType = 'gene_flank',
                         upstream = 1000,
                         mart=ensembl)


 dim(biomart2)
head(biomart2)
