#02
#1.group_list(实验分组)和ids(芯片注释)，每次都需要改
rm(list = ls())  
load("../1_差异分析/step1.Rdata")
library(stringr)
gr1 = pd1$title;gr2 = pd2$title;gr3 = pd3$title;

k1 = str_detect(gr1,"hfPEF");table(k1)
k2 = str_detect(gr1,"Control");table(k2)
k3 = str_detect(gr1,"hfREF");table(k3)
exp1 = exp1[,rownames(pd1)[k1|k2|k3]]
pd1 = pd1[k1|k2|k3,]
x1 = str_detect(gr2,"hfPEF");table(x1)
x2 = str_detect(gr2,"Control");table(x2)
x3 = str_detect(gr2,"hfREF");table(x3)
exp2 = exp2[,rownames(pd2)[x1|x2|x3]]
pd2 = pd2[x1|x2|x3,]
y1 = str_detect(gr3,"hfPEF");table(y1)
y2 = str_detect(gr3,"Control");table(y2)
y3 = str_detect(gr3,"hfREF");table(y3)
exp3 = exp3[,rownames(pd3)[y1|y2|y3]]
pd3 = pd3[y1|y2|y3,]

group1 = ifelse(str_detect(pd1$title,"hfPEF"),
                "hfPEF",ifelse(str_detect(pd1$title,"hfREF"),"hfREF","control"))
table(group1)
group1 = factor(group1,levels = c("hfPEF","control","hfREF"))

group2 = ifelse(str_detect(pd2$title,"hfPEF"),
                "hfPEF",ifelse(str_detect(pd2$ title,"hfREF"),"hfREF","control"))
table(group2)
group2 = factor(group2,levels = c("hfPEF","control","hfREF"))

group3 = ifelse(str_detect(pd3$title,"hfPEF"),
                "hfPEF",ifelse(str_detect(pd3$title,"hfREF"),"hfREF","control"))
table(group3)
group3 = factor(group3,levels = c("hfPEF","control","hfREF"))
#2.探针注释----
temp1 =read.table("../1_差异分析/GPL11434.txt",
                  header = T,skip = 7,
                  sep = "\t",quote ="",fill = T)

k2 = str_detect(temp1$name,'');table(k2)
ids1 = temp1[k2,c(1,2)]
colnames(ids1) = c("probe_id","miRNA")
head(ids1)

temp2 =read.table("../1_差异分析/GPL18066.txt",
                  header = T,skip = 4,
                  sep = "\t",quote ="",fill = T)

x2 = str_detect(temp2$ID,'');table(x2)
ids2 = temp2[x2,c(1,2)]
colnames(ids2) = c("probe_id","miRNA")
head(ids2)

temp3 =read.table("../1_差异分析/GPL18067.txt",
                  header = T,skip = 4,
                  sep = "\t",quote ="",fill = T)
y2 = str_detect(temp3$ID,'');table(y2)
ids3 = temp3[y2,c(1,2)]
colnames(ids3) = c("probe_id","miRNA")
head(ids3)

save(exp1,exp2,exp3,group1,group2,group3,ids1,ids2,ids3,file = "step2.Rdata")
