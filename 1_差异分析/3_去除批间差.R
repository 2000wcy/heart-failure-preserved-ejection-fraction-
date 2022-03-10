#03 batch effect
rm(list = ls())
library(dplyr)
load(file = "step2.Rdata")
probe_id <- row.names(exp1)
exp1 <-mutate(exp1,probe_id) 
exp1 <- merge(exp1,ids1,by="probe_id")
probe_id <- row.names(exp2)
exp2 <-mutate(exp2,probe_id) 
exp2 <- merge(exp2,ids2,by="probe_id")
probe_id <- row.names(exp3)
exp3 <-mutate(exp3,probe_id)
exp3 <- merge(exp3,ids3,by="probe_id")
exp31 <- merge(exp3,exp1,by="miRNA")
exp312 <- merge(exp31,exp2,by="miRNA")
rownames(exp312) <- exp312$miRNA
exp312 <- exp312[-1:-7,-1:-2]
exp312 <- exp312[-47:-50,]
exp312 <- exp312[,-which(colnames(exp312)=='probe_id.y')]
exp312 <- exp312[,-which(colnames(exp312)=='probe_id')]
boxplot(exp312)
group <- c(as.character(group3),as.character(group1),as.character(group2))
design<-model.matrix(~group)
dat=as.data.frame(exp312)
library(limma)
batch <- c(rep("GPL18067",41),rep("GPL11434",33),rep("GPL18066",16))
EXP <- removeBatchEffect(dat,batch=batch,design = design)
boxplot(EXP)
EXP <- as.data.frame(t(EXP))
EXP$group <- group
library(stringr)
# 选择性运行
# 挑选出hfPEF和hfREF(补充研究)
#K1 <- str_detect(EXP$group,"hfPEF");table(K1)#19
#K2 <- str_detect(EXP$group,"hfREF");table(K2)#41

# 挑选出hfREF和control(补充研究)
#K1 <- str_detect(EXP$group,"hfREF");table(K1)#41
#K2 <- str_detect(EXP$group,"control");table(K2)#30

# 挑选出hfPEF和control(文章主要研究部分)
K1 <- str_detect(EXP$group,"hfPEF");table(K2)#41
K2 <- str_detect(EXP$group,"control");table(K2)#30

EXP <- EXP[K1|K2,]
groups <- EXP$group
EXP <- EXP[,-47]
EXP <- t(EXP)
save(groups,EXP,file = "step3.Rdata")
