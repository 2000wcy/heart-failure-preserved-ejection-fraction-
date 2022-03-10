#05差异分析
rm(list = ls()) 
load(file = "step2.Rdata")
load(file = "step3.Rdata")
#差异分析，用limma包来做
#需要表达矩阵和group_list，不需要改
library(limma)
design=model.matrix(~groups)
fit=lmFit(EXP,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
library(dplyr)
deg <- mutate(deg,miRNA=rownames(deg))
head(deg)
deg <- deg[!duplicated(deg$miRNA),]
#3.加change列,标记上下调基因
summary(deg)
logFC_t=0.3
P.Value_t= 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value< P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,
                "down",
                ifelse(k2,
                       "up",
                       "stable"))
table(change)
deg <- mutate(deg,change)
ex_deg <- deg[k1|k2,]

# write.table(ex_deg,file = "ex_deg.txt",
#             quote = F,
#             sep = "\t",
#             col.names = T,
#             row.names = T)
#  挑出rhf与对照组进行的研究

save(groups,deg,logFC_t,P.Value_t,ex_deg,file = "step4.Rdata")
