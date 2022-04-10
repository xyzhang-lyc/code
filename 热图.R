library(ggplot2)
library(reshape2)
library(ggpubr)
ch <- UdS_sort[which(UdS_sort[,1] < 1),2]
for(i in ch){
  cell <- type[which(type$label == i),]
  U <- cell[which(cell$Type == "U"),]
  S <- cell[which(cell$Type == "S"),]
  center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
  mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount)/length(cell$x)
  Smean <- sum(sqrt((S$x - center[1])^2 + (S$y - center[2])^2)*S$MIDCount)/length(S$x)/mean
  Umean <- sum(sqrt((U$x - center[1])^2 + (U$y - center[2])^2)*U$MIDCount)/length(U$x)/mean
  seg <- Umean
  dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount/mean
  nucl <- cell[which(dis<seg),]
  cyto <- cell[which(dis>=seg),]
  a <- ggplot(data = cyto, mapping = aes(x = x, y = y))+ 
    geom_point(alpha = 0.2, size=cyto$MIDCount, color = "#8E8BFE") +
    geom_point(data = nucl, mapping = aes(x = x, y = y),alpha = 1, size=nucl$MIDCount, color = "#FEA3A2") +
    theme(aspect.ratio=(max(cell$x)-min(cell$x))/(max(cell$y)-min(cell$y)))
  S <- S %>% group_by(geneID) %>% 
    summarise(Count = sum(MIDCount))
  U <- U %>% group_by(geneID) %>% 
    summarise(Count = sum(MIDCount))
  nucl <- nucl %>% group_by(geneID) %>% 
    summarise(Count = sum(MIDCount))
  cyto <- cyto %>% group_by(geneID) %>% 
    summarise(Count = sum(MIDCount))
  all <- U %>% full_join(S, by = "geneID") %>% 
    full_join(nucl, by = "geneID") %>% 
    full_join(cyto, by = "geneID")
  colnames(all)[2] <- "U"
  colnames(all)[3] <- "S"
  colnames(all)[4] <- "nucl"
  colnames(all)[5] <- "cyto"
  all[is.na(all)] <- 0
  allhm <- cor(all[,2:5])
  allhm.m <-melt(allhm)
  b <- ggplot(allhm.m, aes(Var1, Var2)) + 
    geom_tile(aes(fill = value),colour = "white") + 
    scale_fill_gradient(name="Value", low = "#FEA3A2",high = "white") +
    labs(x = "Var1", y = "Var2", title = "correlation") +
    theme(plot.title = element_text(size = 13,hjust = 0.5))+
    theme(aspect.ratio=1) +
    geom_text(aes(Var1, Var2, label = format(round(value, 3))), color = "black", size = 3)
  c <- ggpubr::ggarrange(a,b,nrow = 1,ncol = 2,labels = c("A","B"))
  ggsave(paste0("D:/Graduate/空间转录组课题/图片观察/HQ_02/819个cell/cell_",i,".png"),c)
}