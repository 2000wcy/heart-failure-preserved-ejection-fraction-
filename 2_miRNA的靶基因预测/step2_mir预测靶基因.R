rm(list = ls())
options(stringsAsFactors = F)
load("../1_差异分析/all_miRNAs.Rdata")
x = c(up_cgs,down_cgs)
write.table(x,file = "mirwalk_input.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")
#mirwalk工具预测，结果文件放在工作目录下
#http://mirwalk.umm.uni-heidelberg.de/
m_miRNA = read.csv("miRWalk_miRNA_Targets.csv")
colnames(m_miRNA)
colnames(m_miRNA)[1] = "mirnaid"
k = m_miRNA$validated!="";table(k)
m_miRNA = m_miRNA[k,]

k2 = m_miRNA$TargetScan ==1;table(k2)
k3 = m_miRNA$miRDB==1;table(k3)

m_miRNA = m_miRNA[k2|k3,]
library(dplyr)
m_miRNA = distinct(m_miRNA,mirnaid,
                 genesymbol,.keep_all = T)

#得到了多少个m/miRNA?
length(unique(m_miRNA$mirnaid))
length(unique(m_miRNA$genesymbol))

m_miRNA = m_miRNA[,c(1,3)]
save(m_miRNA,file = "mi_mRNA.Rdata")
#cytoscape可视化文件----
#网络文件
write.table(m_miRNA,file = "m_miRNA_edges.txt",
            row.names = F,
            quote = F,
            sep = "\t")
#特征表格
mRNAs = unique(m_miRNA$genesymbol)
nodes = data.frame(symbol = c(down_cgs,up_cgs,mRNAs),
                   type = rep(c("down_mi","up_mi","pc"),
                              times = c(length(down_cgs),
                                        length(up_cgs),
                                        length(mRNAs))))
write.table(nodes,"m_miRNA_nodes.txt",
            row.names = F,
            quote = F,
            sep = "\t")


