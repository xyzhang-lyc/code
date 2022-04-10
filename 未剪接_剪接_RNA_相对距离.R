##新生RNA研究
type <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq.gem",header = T)##导入基因表达矩阵
nucl <- type[which(type$label != -1 & type$Type == "U"),]##找出表达矩阵中的未剪接RNA
cyto <- type[which(type$label != -1 & type$Type == "S"),]##找出表达矩阵中的剪接RNA
N <- type[which(type$label != -1 & type$Type == "N"),]##找出表达矩阵中的其他RNA
##合并数据
all <- rbind(nucl,cyto)
all <- rbind(all,N)
##定义每个细胞未剪接/剪接RNA距离细胞核中心点相对距离
cyto_mean <- rep(0,3589)
nucl_mean <- rep(0,3589)
##定义每个细胞未剪接RNA除以剪接RNA的过滤参数
UdS <- rep(0,3589)
##计算
for(i in 1:3589){
  cell <- all[all$label == i,]
  cyto_cell <- cyto[cyto$label == i,]
  nucl_cell <- nucl[nucl$label == i,]
  if(length(cyto_cell$geneID)<5 | length(nucl_cell$geneID)<5){##if(length(cyto_cell$geneID)+length(nucl_cell$geneID)<10){
    ##将细胞内剪接/未剪接RNA点数小于5的细胞过滤
    cyto_mean[i] <- NaN
    nucl_mean[i] <- NaN
  }
  else{
    center <- c(sum(nucl_cell$x * nucl_cell$MIDCount)/sum(nucl_cell$MIDCount),sum(nucl_cell$y * nucl_cell$MIDCount)/sum(nucl_cell$MIDCount))
    cell_mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount)/length(cell$x)
    cyto_mean[i] <- sum(sqrt((cyto_cell$x - center[1])^2 + (cyto_cell$y - center[2])^2)*cyto_cell$MIDCount)/length(cyto_cell$x)/cell_mean
    nucl_mean[i] <- sum(sqrt((nucl_cell$x - center[1])^2 + (nucl_cell$y - center[2])^2)*nucl_cell$MIDCount)/length(nucl_cell$x)/cell_mean
  }
}
UdS <- nucl_mean/cyto_mean 
UdS_sort <- cbind(UdS,seq(1,3589,1))
##计算每个细胞剪接/未剪接RNA距离细胞核中心点相对距离的差异
wilcox.test(cyto_mean,nucl_mean,paired = T)
##将UdS由小到大排序
UdS_sort <- na.omit(UdS_sort[order(UdS_sort[,1],decreasing = F),])