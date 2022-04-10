##基因注释信息
library(tidyverse)
blastp <- read.table("D:/Graduate/空间转录组课题/数据表格/soybean.protein.1.fasta.blastp.anno",header = F,sep = '\t',quote = "")
for(i in 1:length(blastp$V1)){
  blastp$V1[i] <- unlist(strsplit(blastp$V1[i],split = '[.]'))[1]
}
colnames(blastp)[1] <- "geneID"
colnames(blastp)[2] <- "ipr_acc"
colnames(blastp)[3] <- "ipr_desc"
write.table(blastp,"D:/Graduate/空间转录组课题/数据表格/blastp.anno.txt",sep = "\t",row.names = F,col.names = T,quote = F)
##计算细胞核内外基因表达差异p-value及基因注释合并
temp$p_value <- rep(0,length(temp$geneID))
temp <- na.omit(temp)
for(i in 1:length(temp$geneID)){
  table <- matrix(c(temp[i,]$Count.x,(temp[i,]$Count.x + temp[i,]$Count.y/111478*107289)/2,temp[i,]$Count.y/111478*107289,(temp[i,]$Count.x + temp[i,]$Count.y/111478*107289)/2),nrow = 2,ncol = 2)
  temp$p_value[i] <- chisq.test(table)$p.value
}
Anno <- temp %>% left_join(blastp, by = "geneID")
write.table(Anno,"D:/Graduate/空间转录组课题/数据表格/26_02hq_anno.txt",sep = "\t",row.names = F,col.names = T,quote = F)
##找出细胞核内/外表达更显著基因
cyto_gene <- temp[which(temp$Count.x>temp$Count.y/111478*107289),]
nucl_gene <- temp[which(temp$Count.x<temp$Count.y/111478*107289),]
cyto_gene <- cyto_gene %>% left_join(blastp, by = "geneID")
nucl_gene <- nucl_gene %>% left_join(blastp, by = "geneID")
write.table(cyto_gene,"D:/Graduate/空间转录组课题/数据表格/26_02hq_cyto_gene.txt",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(nucl_gene,"D:/Graduate/空间转录组课题/数据表格/26_02hq_nucl_gene.txt",sep = "\t",row.names = F,col.names = T,quote = F)