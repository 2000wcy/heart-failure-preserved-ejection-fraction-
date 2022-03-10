#04 检查数据
#PCA
rm(list = ls())
load(file = "step3.Rdata")
dat=as.data.frame(t(EXP))
library(FactoMineR)#画主成分分析图需要加载这两个包
library(factoextra) 
# pca的统一操作走起
dat.pca <- PCA(dat, graph = FALSE)
pca_plot<- fviz_pca_ind(dat.pca,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = groups, # color by groups
                        #palette = c("#00AFBB", "#E7B800"),
                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Groups"
)
pca_plot
ggsave(filename = "PCA.PNG",width = 8,height = 6)
save(pca_plot,file = "pca_plot.Rdata")

#热图
annotation_col=data.frame(groups)
rownames(annotation_col)=colnames(EXP) 
library(pheatmap)
pheatmap(EXP,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row")

dev.off()

