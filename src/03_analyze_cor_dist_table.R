library(ProjectTemplate)
load.project()

load('cache/diff.cor.pairs.rda')


same.chromosome <- diff.cor.pairs %>% filter(dist != Inf)

pick_dist <- function(data, b1, b2=0) {
  filter(data, abs(dist) < b1 & abs(dist) > b2) %>%
    filter((type1%in%config$lncRNA & type2%in%config$PCGs) | (type1%in%config$PCGs & type2%in%config$lncRNA)) %>%
    filter(significant)
}
pick_dist(same.chromosome, 1000)
