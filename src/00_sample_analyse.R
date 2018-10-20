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
