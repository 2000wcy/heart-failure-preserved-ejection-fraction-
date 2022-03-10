#06 差异基因的火山图和热图
rm(list = ls()) 
load(file = 'step3.Rdata')
load(file = "step4.Rdata")
DEG=rownames(ex_deg)
#aa <- c("hsa-miR-720","hsa-miR-1246","hsa-miR-1827")
DEG=EXP[DEG,]
#DEG=EXP[aa,]
#1.火山图----
library(dplyr)
library(ggplot2)
dat  = deg

#if(F){
# for_label <- dat%>% 
#  filter(symbol %in% c("TRPM3","SFRP1")) 
#}
if(F){
  for_label <- dat %>% head(10)
}
if(T) {
  x1 = dat[order(dat$logFC,decreasing = T),] %>% 
    filter(change == "up") %>% 
    head(2)
  x2 = dat[order(dat$logFC,decreasing = F),] %>% 
    filter(change == "down") %>% 
    head(4)
  for_label = rbind(x1,x2)
}
#%>%
#  head(10)
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label =miRNA ),
    data = for_label,
    color="black"
  )
volcano_plot

#2.差异基因热图----
x=deg$logFC 
names(x)=deg$miRNA
#cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))

cg = deg$miRNA[deg$change !="stable"]
exp=EXP[,order(groups)]
group = groups[order(groups)]
n=exp[cg,]
dim(n)
#作热图
library(pheatmap)
annotation_col=data.frame(group)
rownames(annotation_col)=colnames(n) 
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,cellheight = 20,cellwidth = 4,
                                   show_colnames =F,
                                   show_rownames = T,
                                   scale = "row",
                                   cluster_cols = F,
                                   cluster_rows = T, 
                                   annotation_col=annotation_col)) 
heatmap_plot

load("pca_plot.Rdata")
library(patchwork)
(heatmap_plot + volcano_plot)+ plot_annotation(tag_levels = "A")
#要根据合并的图形将绘图区拉到合适的比例，否则会报错
ggsave(filename = "HF.PNG",width = 15,height = 5)
#做ROC曲线可能要用
save(DEG,file = "deg_exp.Rdata" )
