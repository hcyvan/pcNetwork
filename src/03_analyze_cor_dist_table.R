library(ProjectTemplate)
#load.project()

load('cache/diff.cor.pairs.rda')


same.chromosome <- diff.cor.pairs %>% filter(dist != Inf)

pick_dist <- function(data, b1, b2=0) {
  filter(data, abs(dist) < b1 & abs(dist) > b2) %>%
    filter((type1%in%config$lncRNA & type2%in%config$PCGs) | (type1%in%config$PCGs & type2%in%config$lncRNA)) %>%
    filter(significant)
}
diff.1k <- pick_dist(same.chromosome, 2000)
diff.1k <- mutate(diff.1k, symbol1=helper.getGeneSymbol(Var1), symbol2=helper.getGeneSymbol(Var2))

write.csv(diff.1k, 'reports/diff.1k.csv')

