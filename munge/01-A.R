count.data <- count.data %>% mutate(GeneID=str_split_fixed(GeneID, '\\.', 2)[,1])
fpkm.data <- fpkm.data %>% mutate(GeneID=str_split_fixed(GeneID, '\\.', 2)[,1])

diff.all.106 <- read.delim('./reports/diff.all.106.csv', sep = ',')[,-1]
diff.lncrna.pcg.106 <- filter(diff.all.106, GeneType %in% c(config$lncRNA, config$PCGs))
diff.all.551 <- read.delim('./reports/diff.all.551.csv', sep = ',')[,-1]
diff.lncrna.pcg.551 <- filter(diff.all.551, GeneType %in% c(config$lncRNA, config$PCGs))

bioMart <- distinct(mart.export, Gene.stable.ID, .keep_all = TRUE)


load('cache/diff.pvalue.rda')

