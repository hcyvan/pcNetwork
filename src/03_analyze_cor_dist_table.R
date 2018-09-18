library(ProjectTemplate)
load.project()

load('cache/diff.cor.pairs.rda')
load('cache/diff.cor.pairs.106.rda')


pickDist <- function(data, b1, b2=0) {
  filter(data, abs(dist) < b1 & abs(dist) > b2) %>%
    filter((type1%in%config$lncRNA & type2%in%config$PCGs) | (type1%in%config$PCGs & type2%in%config$lncRNA)) %>%
    filter(significant)
}
getCisGene <- function(pairs, bound) {
  same.chromosome <- pairs %>% filter(dist != Inf)
  gene.near <- pickDist(same.chromosome, bound) %>% mutate(symbol1=helper.getGeneSymbol(Var1), symbol2=helper.getGeneSymbol(Var2))
  gene.near <- gene.near %>% mutate_if(is.factor, as.character)
  gene.near
}


cis.2k.551 <- getCisGene(diff.cor.pairs, 2000) 
write.csv(cis.2k.551, 'reports/cis.2k.551.csv')
cis.2k.106 <- getCisGene(diff.cor.pairs.106, 2000)
write.csv(cis.2k.106, 'reports/cis.2k.106.csv')

pair551 <- str_c(cis.2k.551$Var1,'-',cis.2k.551$Var2)
pair106 <- str_c(cis.2k.106$Var1,'-',cis.2k.106$Var2)
