#07 上下调基因
rm(list = ls())
load("step4.Rdata")
library(stringr)
deg <-merge(deg,ex_deg,by=c('miRNA'),all.y = T)
k <- str_detect(deg$miRNA,'')
deg <- deg[k,]
write.table(deg,file = "ex_all.txt",
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")

up_cgs = c(deg$miRNA[deg$change.x=="up"])
down_cgs = c(deg$miRNA[deg$change.x=="down"])

save(up_cgs,down_cgs,file = "all_miRNAs.Rdata")
