rm(list= ls())
library(ggplot2)
library(pROC)
source("ROC_Plot.R")
source("Theme_Publication.R")
load(file = "../1_差异分析/deg_exp.Rdata")
load(file = "../1_差异分析/step3.Rdata")
DEG <- as.data.frame(t(DEG))
group <- model.matrix(~groups)
DEG$group <- group[,2]
DEG$group[which(DEG$group=="1")] <- 2
DEG$group[which(DEG$group=="0")] <- 1
DEG$group[which(DEG$group=="2")] <- 0
write.table(DEG,file = "DEG.txt",quote = F)

x1 <- DEG$`hsa-miR-720`
x2 <- DEG$`hsa-miR-628-3p`
x3 <- DEG$`hsa-miR-1321`
x4 <- DEG$`hsa-miR-636`
x5 <- DEG$`hsa-miR-1246`
x6 <- DEG$`hsa-miR-1827`

y <- DEG$group

R1 <- roc(y,x1,ci = T); R2 <- roc(y,x2,ci = T); R3 <- roc(y,x3,ci = T);R4 <- roc(y,x4,ci = T)
R5 <- roc(y,x5,ci = T);R6 <- roc(y,x6,ci = T)

dat1 <- data.frame(Sensitivities = R1$sensitivities, 
                   FalsePositiveRate = (1-R1$specificities), 
                   Model = rep("hsa-miR-720",length(R1$sensitivities)))
dat2 <- data.frame(Sensitivities = R2$sensitivities,
                   FalsePositiveRate = (1-R2$specificities),
                   Model = rep("hsa-miR-628-3p",length(R2$sensitivities)))
dat3 <- data.frame(Sensitivities = R3$sensitivities,
                   FalsePositiveRate = (1-R3$specificities),
                   Model = rep("hsa-miR-1321",length(R3$sensitivities)))
dat4 <- data.frame(Sensitivities = R4$sensitivities,
                   FalsePositiveRate = (1-R4$specificities),
                   Model = rep("hsa-miR-636",length(R4$sensitivities)))
dat5 <- data.frame(Sensitivities = R5$sensitivities, 
                   FalsePositiveRate = (1-R5$specificities), 
                   Model = rep("hsa-miR-1246",length(R5$sensitivities)))
dat6 <- data.frame(Sensitivities = R6$sensitivities, 
                   FalsePositiveRate = (1-R6$specificities), 
                   Model = rep("hsa-miR-1827",length(R6$sensitivities)))

dat <- rbind(dat1,dat2,dat3,dat4,dat5,dat6)

colnames(dat)[2] <- ("1-Specificities")

ggplot(dat,aes(`1-Specificities`,Sensitivities,colour = Model)) + geom_roc_plot()+
  scale_colour_Publication() + theme_Publication()+
  theme(panel.border = element_rect(colour = 'black'),legend.title = element_blank(),
        legend.position = c(0.9,0.15),legend.direction = 'vertical')

#做表
mir <- (c('hsa-miR-720','hsa-miR-628-3p','hsa-miR-1321',
          'hsa-miR-636','hsa-miR-1246','hsa-miR-1827'))

auc <- data.frame(rn= c('low','auc','up'))
roc <- list(R1,R2,R3,R4,R5,R6)

for (n in 1:6) {
  auc1 <- data.frame(roc[[n]]$ci)
  auc <- cbind(auc,auc1)
}
rownames(auc) <- auc[,1]
auc <- auc[,-1]
colnames(auc) <- mir
write.csv(auc,'auc_PEF_Con.csv')


##----------
# 如果前面挑选了hfREF、control和hfREF、hfPEF，
# 运行此部分即可
rm(list= ls())
library(ggplot2)
library(pROC)
source("ROC_Plot.R")
source("Theme_Publication.R")
load(file = "../1_差异分析/deg_exp.Rdata")
load(file = "../1_差异分析/step3.Rdata")
DEG <- as.data.frame(t(DEG))
group <- model.matrix(~groups)
DEG$group <- group[,2]
DEG$group[which(DEG$group=="1")] <- 2
DEG$group[which(DEG$group=="0")] <- 1
DEG$group[which(DEG$group=="2")] <- 0
write.table(DEG,file = "DEG.txt",quote = F)

x1 <- DEG$`hsa-miR-720`
x2 <- DEG$`hsa-miR-1246`
x3 <- DEG$`hsa-miR-1827`

y <- DEG$group


dat1 <- data.frame(Sensitivities = R1$sensitivities, 
                   FalsePositiveRate = (1-R1$specificities), 
                   Model = rep("hsa-miR-720",length(R1$sensitivities)))
dat2 <- data.frame(Sensitivities = R2$sensitivities, 
                   FalsePositiveRate = (1-R2$specificities), 
                   Model = rep("hsa-miR-1246",length(R2$sensitivities)))
dat3 <- data.frame(Sensitivities = R3$sensitivities, 
                   FalsePositiveRate = (1-R3$specificities), 
                   Model = rep("hsa-miR-1827",length(R3$sensitivities)))


dat <- rbind(dat1,dat2,dat3)
colnames(dat)[2] <- ("1-Specificities")

ggplot(dat,aes(`1-Specificities`,Sensitivities,colour = Model)) + geom_roc_plot()+
   scale_colour_Publication() + theme_Publication()+
   theme(panel.border = element_rect(colour = 'black'),legend.title = element_blank(),
         legend.position = c(0.9,0.15),legend.direction = 'vertical')

#做表

mir <- (c('hsa-miR-720','hsa-miR-1246','hsa-miR-1827'))
auc <- data.frame(rn= c('low','auc','up'))
roc <- list(R1,R2,R3)

for (n in 1:3) {
  auc1 <- data.frame(roc[[n]]$ci)
  auc <- cbind(auc,auc1)
}
rownames(auc) <- auc[,1]
auc <- auc[,-1]
colnames(auc) <- mir
# 导出对应即可
# write.csv(auc,'auc_REF_Con.csv')
# write.csv(auc,'auc_PEF_REF.csv')
