data.count <- read.delim('./data/data.count.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE, row.names = 'GeneID')

sample.id <- names(data.count)
cases <- substr(sample.id, 1,12)
sample.type <- substr(sample.id, 14, 15)
sample.type <- ifelse(sample.type=='11', 'N', 'T')

data.sample <- data.frame(case.id=cases, sample.type=sample.type, sample.id=sample.id)
write.csv(data.sample, file = './data/data.sample.csv', row.names = FALSE)

#--------------------------------- 
table(table(data.sample$case.id))
table(data.sample$sample.type)

case_3 <- names(sort(table(data.sample$case.id), decreasing = T)[1:2])
filter(data.sample, case.id %in% case_3)


################################################ download biomaRt data
library(biomaRt)
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
searchAttributes(mart = ensembl, pattern = "ensembl")
biomart.symbol.biotype = getBM(attributes = c('ensembl_gene_id','hgnc_symbol', 'gene_biotype'), mart = ensembl)
save(biomart.symbol.biotype, file = './cache/biomart.symbol.biotype.rda')