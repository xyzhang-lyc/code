##分离细胞核内/外点
cell <- type[which(type$label == 1837),]
U <- cell[which(cell$Type == "U"),]
S <- cell[which(cell$Type == "S"),]
center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount)/length(cell$x)
Smean <- sum(sqrt((S$x - center[1])^2 + (S$y - center[2])^2)*S$MIDCount)/length(S$x)/mean
Umean <- sum(sqrt((U$x - center[1])^2 + (U$y - center[2])^2)*U$MIDCount)/length(U$x)/mean
seg <- Umean
dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount/mean
nucl <- cell[which(dis<seg),]
cyto <- cell[which(dis>seg),]
for(i in ch[2:length(ch)]){
  cell <- type[which(type$label == i),]
  U <- cell[which(cell$Type == "U"),]
  S <- cell[which(cell$Type == "S"),]
  center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
  mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount)/length(cell$x)
  Smean <- sum(sqrt((S$x - center[1])^2 + (S$y - center[2])^2)*S$MIDCount)/length(S$x)/mean
  Umean <- sum(sqrt((U$x - center[1])^2 + (U$y - center[2])^2)*U$MIDCount)/length(U$x)/mean
  seg <- Umean
  dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount/mean
  cyto <- rbind(cyto,cell[which(dis>seg),])
  nucl <- rbind(nucl,cell[which(dis<seg),])
}
write.table(cyto,"D:/Graduate/空间转录组课题/数据表格/26_02hq_cyto.txt",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(nucl,"D:/Graduate/空间转录组课题/数据表格/26_02hq_nucl.txt",sep = "\t",row.names = F,col.names = T,quote = F)
##检验细胞核内外点基因表达差异
wilcox.test(temp$Count.x/107289,temp$Count.y/111478,paired = T)