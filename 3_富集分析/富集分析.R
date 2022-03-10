library(ggplot2)
library(tidyr)
rm(list=ls())
gene <- read.table("../2_miRNA的靶基因预测/m_miRNA_edges.txt",
                    header = T)
gene <- unique(gene$genesymbol)
write.table(gene,file = "gene.txt",quote = F,
            col.names = F,row.names =F )
#https://david.ncifcrf.gov/
#KEGG
rt = read.table(file = 'KEGG.txt' ,sep = '\t',header = T,quote = '')
keggsig = rt[rt$PValue < 0.05,]
keggsig = keggsig[keggsig$Count>=3,]
keggsig = separate(keggsig,Term,sep = ":",
                   into = c("ID","Term"))
write.csv(keggsig,file = "KEGG.csv")
ggplot(keggsig,aes(x=Fold.Enrichment,y=Term))+
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="blue",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number",
    x="Fold enrichment",
    y="Pathway name",
    title="Pathway enrichment")+
  theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.3)),
        axis.title.x = element_text(size=rel(1.3)),
        axis.title.y = element_blank())
ggsave('KEGG_DAVID.png',width = 7,height = 4)

#GO
BP <- read.table(file = "GO_BP.txt",header = T,sep = "\t")
BP <- BP[BP$PValue < 0.05,]
BP <- BP[BP$Count>=3,]
BP <- BP[order(BP$Count,decreasing = T),]
BP <- BP[1:2,]

CC <- read.table(file = "GO_CC.txt",header = T,sep = "\t")
CC <- CC[CC$PValue < 0.05,]
CC <- CC[CC$Count>=3,]
CC <- CC[order(CC$Count,decreasing = T),]
CC <- CC[1:2,]

MF <- read.table(file = "GO_MF.txt",header = T,sep = "\t")
MF <- MF[MF$PValue < 0.05,]
MF <- MF[MF$Count>=3,]
MF <- MF[order(MF$Count,decreasing = T),]
MF <- MF[1:3,]
GO <- rbind(BP,CC,MF)
GO <-  separate(GO,Term,sep = "~",
                   into = c("ID","Term"))
write.csv(GO,file = "GO.csv")
GO$GO_Term <- c(rep("BP",2),rep("CC",2),rep("MF",3))
rownames(GO) <-1:7
GO_order <- factor(as.integer(rownames(GO)),labels = GO$Term)
ggplot(GO,aes(x=GO_order,y=Count,fill=GO_Term))+
  geom_bar(stat = "identity",width = 0.5)+
  coord_flip()+
  xlab(" ")+
ylab("Gene number") 
  
  

ggsave('GO_DAVID.png',width = 10,height = 4)


