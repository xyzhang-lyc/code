##不同层未剪接/剪接RNA差异分析
cell <- type[which(type$label == 1837),]
U <- cell[which(cell$Type == "U"),]
center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount)/length(cell$x)
dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount/mean
for(j in 1:5){
  assign(paste("cont",j,sep = ""),cell[which((j-1)*max(dis)/5 < dis & dis < j*max(dis)/5),])
}
for(i in UdS_sort[2:100,2]){
  cell <- type[which(type$label == i),]
  U <- cell[which(cell$Type == "U"),]
  center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
  mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount)/length(cell$x)
  dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)*cell$MIDCount/mean
  for(j in 1:5){
    assign(paste("cont",j,sep = ""),rbind(get(paste("cont",j,sep = "")),cell[which((j-1)*max(dis)/5 < dis & dis < j*max(dis)/5),]))
  }
}
Scont <- c(length(which(cont1$Type == "S"))/length(cont1$geneID))##1-10
Scont <- append(Scont,length(which(cont5$Type == "S"))/length(cont5$geneID))
Ucont <- c(length(which(cont1$Type == "U"))/length(cont1$geneID))
Ucont <- append(Ucont,length(which(cont5$Type == "U"))/length(cont5$geneID))

##箱型图观察
boxmatrix1 <- matrix(rep(0,500),ncol = 5)
boxmatrix2 <- matrix(rep(0,500),ncol = 5)
n <- 1
cell <- type[which(type$label == 1837),]
U <- cell[which(cell$Type == "U"),]
center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)
for(j in 1:5){
  cont <- cell[which((j-1)*max(dis)/5 < dis & dis < j*max(dis)/5),]
  boxmatrix1[n,j] <- length(cont[which(cont$Type == "U"),]$geneID)/length(cont$geneID)
  boxmatrix2[n,j] <- length(cont[which(cont$Type == "S"),]$geneID)/length(cont$geneID)
}
n <- n + 1
for(i in UdS_sort[2:100,2]){
  cell <- type[which(type$label == i),]
  U <- cell[which(cell$Type == "U"),]
  center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
  dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)
  for(j in 1:5){
    cont <- cell[which((j-1)*max(dis)/5 < dis & dis < j*max(dis)/5),]
    boxmatrix1[n,j] <- length(cont[which(cont$Type == "U"),]$geneID)/length(cont$geneID)
    boxmatrix2[n,j] <- length(cont[which(cont$Type == "S"),]$geneID)/length(cont$geneID)
  }
  n <- n + 1
}