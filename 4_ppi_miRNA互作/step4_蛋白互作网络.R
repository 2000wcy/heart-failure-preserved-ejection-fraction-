rm(list = ls())
options(stringsAsFactors = F)
load("../1_差异分析/all_miRNAs.Rdata")
load("../2_miRNA的靶基因预测/mi_mRNA.Rdata")

mRNAs = unique(m_miRNA$genesymbol)
nodes = data.frame(symbol = c(down_cgs,up_cgs,mRNAs),
                   type = rep(c("down_mi","up_mi","pc"),
                              times = c(length(down_cgs),
                                        length(up_cgs),
                                        length(mRNAs))))

write.table(nodes,"ppi_mir_nodes.txt",
            row.names = F,
            quote = F,
            sep = "\t")

#ppi输入数据
x = unique(m_miRNA$genesymbol)
write.table(x,file = "string_input.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")

#ppi输出数据和m_miRNA合并----
#https://string-db.org/cgi/input.pl?sessionId=rKtCnBvlHDFo&input_page_active_form=multiple_identifiers
pp = read.table("string_interactions.tsv")[,c(1,2)]

#只要输入的这些基因，不要其他
k1 = pp$V1 %in% mRNAs
k2 = pp$V2 %in% mRNAs
table(k1&k2)
pp = pp[k1&k2,]

#只要与其他基因有互作的基因
k3 = m_miRNA$genesymbol %in% c(pp[,1],pp[,2]);table(k3)
m_miRNA_ppi = m_miRNA[k3,]

colnames(pp) = colnames(m_miRNA_ppi)
pm = rbind(pp,m_miRNA_ppi)

#连接的类别，用于画箭头
pm$tp = rep(c("ppi","mi"),
            times = c(nrow(pp),nrow(m_miRNA_ppi)))

write.table(pm,"ppi_miRNA_edges.txt",
            row.names = F,
            quote = F,
            sep = "\t")

