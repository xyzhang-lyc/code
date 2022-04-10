library(tidyverse)
library(GO.db)
godb <- select(GO.db, keys(GO.db), columns(GO.db))
cytoGO <- read.table("C:/Users/10145/Desktop/cyto.txt", header = F)##colnames(cytoGO) <- "GOID"
nuclGO <- read.table("C:/Users/10145/Desktop/nucl.txt", header = F)##colnames(nuclGO) <- "GOID"
cytoanno <- cytoGO %>% left_join(godb)
nuclanno <- nuclGO %>% left_join(godb)
celgo <- read.table("D:/Graduate/空间转录组课题/数据表格/17_FP200000440BR_A2_83_HQ03.celgo.txt", header = F)
m <- rep(0,1842)
nucl <- c()
cyto <- c()
for(i in 1:1842){
  m[i] <- mean(na.omit(celgo[,i]))
  if(is.na(m[i])){next}
  if(m[i] < 0.8){nucl <- append(nucl,i)}
  if(m[i] > 1.2){cyto <- append(cyto,i)}
}
hist(na.omit(m),breaks = seq(0,2,0.01),col = "#BEB8DC",border = "#BEB8DC")
for(i in 1:length(celgo[1,])){
  celgo[which(celgo[,i] == Inf),i] <- NA
}





library(tidyverse)
library(GO.db)
godb <- select(GO.db, keys(GO.db), columns(GO.db))
cytoGO <- read.table("C:/Users/10145/Desktop/cyto.txt", header = F)##colnames(cytoGO) <- "GOID"
nuclGO <- read.table("C:/Users/10145/Desktop/nucl.txt", header = F)##colnames(nuclGO) <- "GOID"
cytoanno <- cytoGO %>% left_join(godb)
nuclanno <- nuclGO %>% left_join(godb)
celgo <- read.table("D:/Graduate/空间转录组课题/数据表格/17_FP200000440BR_A2_83_HQ03.celgo.txt", header = F)
m <- rep(0,1842)
nucl <- c()
cyto <- c()
for(i in 1:1842){
  m[i] <- mean(na.omit(celgo[,i]))
  if(is.na(m[i])){next}
  if(m[i] < 0.8){nucl <- append(nucl,i)}
  if(m[i] > 1.2){cyto <- append(cyto,i)}
}
hist(na.omit(m),breaks = seq(0,2,0.01),col = "#BEB8DC",border = "#BEB8DC")
for(i in 1:length(celgo[1,])){
  celgo[which(celgo[,i] == Inf),i] <- NA
}
write.table(nuclanno,"D:/Graduate/空间转录组课题/数据表格/nuclanno.txt",sep = "\t",row.names = F,col.names = F,quote = F)








##WEGO注释
library(tidyverse)##管道符
library(GO.db)

godb <- select(GO.db, keys(GO.db), columns(GO.db))


interpro <- read.table("D:/Graduate/空间转录组课题/数据表格/soybean.interproscan.tsv", header = F, fill = T, sep = "\t")
gene2go <- interpro %>% dplyr::select(Gene = V1, GOID = V14) %>% na.omit() %>%
  separate(GOID, paste0("V", 1:(max(str_count(.$GOID,"\\|"))+1), seq = ""), sep = "\\|")  %>% 
  gather(key = "V", value = "GOID", -Gene) %>% dplyr::select(Gene, GOID) %>%
  na.omit() %>% base::unique()


go_annot <- gene2go %>% left_join(godb)

go2gene <- gene2go %>% group_by(GOID) %>%
  summarise(Gene = str_c(Gene, collapse = ",")) %>%
  mutate(Count = str_count(Gene, ",")+1) %>% arrange(desc(Count)) %>%
  left_join(godb)

wego <- gene2go %>% group_by(Gene) %>%
  summarise(GOID = str_c(GOID, collapse = ",")) %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,","))+1), seq = ""), sep = ",")
write_tsv(wego, "soybean.interproscan.wego.txt", col_names = FALSE, na = "")

ipr2gene <- interpro %>% group_by(V12) %>% 
  summarise(Gene = str_c(V1, collapse = ",")) %>%
  mutate(Count = str_count(Gene, ",")+1) %>% arrange(desc(Count)) %>%
  left_join(ipr) %>% base::unique()





##筛选符合条件的点
Gcount <- celabel %>% group_by(label) %>% 
  summarise(Gene = str_c(geneID, collapse = ",")) %>%
  mutate(Count = str_count(Gene, ",")+1) %>% arrange(desc(Count))
chlabel <- Gcount[which(Gcount$Count>=50),]##找出genecount>=50的细胞
chlabel <- chlabel$label[-c(1,2,3,5,8,15)]##排除不是细胞的部分

labcount <- celabel %>% group_by(geneID) %>% 
  summarise(label = str_c(unique(label), collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
chgene <- labcount$geneID[which(labcount$Count>3)]
chgene2go <- gene2go[gene2go$Gene %in% chgene,]

go2chgene <- chgene2go %>% group_by(GOID) %>%
  summarise(Gene = str_c(Gene, collapse = ",")) %>%
  mutate(Count = str_count(Gene, ",")+1) %>% arrange(desc(Count)) %>%
  left_join(godb)
chpoint <- celabel[celabel$geneID %in% chgene,]
chpoint <- chpoint[chpoint$label %in% chlabel,]
cel1 <- chpoint[which(chpoint$label == 4599),]
center1 <- c(sum(cel1$x)/length(cel1$x),sum(cel1$y)/length(cel1$y))
cel1mean <- sum(sqrt((cel1$x - center1[1])^2 + (cel1$y - center1[2])^2))/length(cel1$x)
f = seq(0,2*pi,0.001)
x = cel1mean*sin(f) + center1[1]
y = cel1mean*cos(f) + center1[2]
plot(x,y,type='l',xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)),asp = 1)
par(new=TRUE)
plot(cel1$x,cel1$y,xlab = "",ylab = "", pch = 20,col = rgb(45, 67,121, 100, maxColorValue = 255),xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)), cex = cel1$MIDCount)
par(new=TRUE)
plot(center1[1],center1[2],xlab = "",ylab = "",xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)))






##interpro数据库注释结果简单分析及可视化
interpro <- read_tsv("D:/Graduate/空间转录组课题/soybean.interproscan.type.tsv", na = "N/A", col_names = c("geneID","ipr_acc","ipr_desc")) %>% na.omit()

ipr <- interpro %>% select(geneID,ipr_acc,ipr_desc) %>% group_by(geneID,ipr_acc) %>%
  summarise(ipr_desc = ipr_desc) %>% group_by(ipr_acc,ipr_desc) %>% summarise(Count=n())%>%
  arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))

p <- ggplot(ipr) +
     geom_bar(aes(x = ipr_desc, y = Percent, fill = ipr_desc), stat = "identity") +
     scale_y_continuous(labels = scales::percent, limits = c(0, 0.08),name = "Percent of Domain") +
     scale_x_discrete(limits = ipr$ipr_desc[1:20], name = NULL) + scale_fill_discrete(guide = FALSE)+
     theme(axis.text.x=element_text(angle=60,vjust=1, hjust=1))





Coord <- read.table("D:/Graduate/空间转录组课题/17_FP200000440BR_A2_83_HQ03.sort.gem", header = F)
go <- strsplit(go2gene$Gene[1],split = ',')[[1]]
gocoord <- Coord[which(Coord$V1[1:500] %in% go),]
##找到go表中排名第一的go号




##library(png)
##imgpng <- readPNG("D:/Graduate/空间转录组课题/17_FP200000440BR_A2_83_HQ03.png")
##plot(c(0,ncol(imgpng)*4.75),c(0,nrow(imgpng)*4.75),type = "n",xlab = "",ylab = "",asp = 1)
##rasterImage(imgpng,0,0,ncol(imgpng)*4.75,nrow(imgpng)*4.75)
##par(new = TRUE) 
##plot(Coord$V2-11900,Coord$V3-5250,col = rgb(45, 67,121, 10, maxColorValue = 255),xlab = "",ylab = "",xlim = c(0,ncol(imgpng)*4.75),ylim = c(0,nrow(imgpng)*4.75))
##par(new = TRUE) 
##f=seq(0,2*pi,0.001)
##x=50*sin(f)+13475-11900
##y=50*cos(f)+6910-5250
##plot(x,y,type='l',xlim = c(0,ncol(imgpng)*4.75),ylim = c(0,nrow(imgpng)*4.75),asp=1)
##根据图片找到近似细胞的圆





celabel <- read.table("D:/Graduate/空间转录组课题/数据表格/Cell_GetExp_gene.txt",header = T)
cel1 <- celabel[which(celabel$label == 2655),]
center1 <- c(sum(cel1$x)/length(cel1$x),sum(cel1$y)/length(cel1$y))
cel1mean <- sum(sqrt((cel1$x - center1[1])^2 + (cel1$y - center1[2])^2))/length(cel1$x)
f = seq(0,2*pi,0.001)
x = cel1mean*sin(f) + center1[1]
y = cel1mean*cos(f) + center1[2]
plot(x,y,type='l',xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)),asp = 1)
par(new=TRUE)
plot(cel1$x,cel1$y,xlab = "",ylab = "", pch = 20,col = rgb(45, 67,121, 100, maxColorValue = 255),xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)))
par(new=TRUE)
plot(center1[1],center1[2],xlab = "",ylab = "",xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)))
go1 <- strsplit(go2gene$Gene[15],split = ',')[[1]]
gocoord1 <- cel1[which(cel1$geneID == "SoyZH13_10G190700"),]
par(new=TRUE) 
plot(gocoord1$x,gocoord1$y, col = '#FF7F50', xlim = c(min(cel1$x),max(cel1$x)), ylim = c(min(cel1$y),max(cel1$y)), xlab = "", ylab = "", pch = 20, cex = gocoord1$MIDCount, main = "")




##整个样本图
plot(celabel$x,celabel$y,col = rgb(45, 67,121, 100, maxColorValue = 255),xlim = c(min(celabel$x),max(celabel$x)), ylim = c(min(celabel$y),max(celabel$y)),xlab = "",ylab = "",pch = ".")






##找到点最多的cell
m <- 0
c <- 0
for(i in 1:3589){
  cel1 <- all[which(all$label == i),]
  if(length(cel1[,1])>m && c != 557){
    m <- length(cel1[,1])
    c <- i
  }
}
##2464,2655,2864,2819,4071,4621,3920,2776,2422,2263,3360,3926,2813,4770,2511,4645,4646,3290
##2649,3207,2914,4649不确定




cel1 <- celabel[which(celabel$label == 2464),]
center1 <- c(sum(cel1$x)/length(cel1$x),sum(cel1$y)/length(cel1$y))
plot(cel1$x,cel1$y,xlab = "",ylab = "", pch = 20,col = rgb(45, 67,121, 100, maxColorValue = 255),xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)))
par(new=TRUE)
plot(center1[1],center1[2],xlab = "",ylab = "",xlim = c(min(cel1$x),max(cel1$x)),ylim = c(min(cel1$y),max(cel1$y)))
ipr1 <- strsplit(ipr2gene$Gene[1],split = ',')[[1]]
iprcoord1 <- cel1[which(cel1$geneID %in% ipr1),]
par(new=TRUE) 
plot(iprcoord1$x,iprcoord1$y, col = '#FF7F50', xlim = c(min(cel1$x),max(cel1$x)), ylim = c(min(cel1$y),max(cel1$y)), xlab = "", ylab = "", pch = 20, cex = iprcoord1$MIDCount, main = "Pentatricopeptide repeat")





##找到每个细胞每种GO的到中心点的相对距离
celgo <- matrix(0,max(celabel$label),length(go2gene$GOID))
for(i in 1:length(celgo[1,])){
  go <- strsplit(go2gene$Gene[i],split = ",")[[1]]
  for(j in 1:length(celgo[,1])){
    cel <- celabel[which(celabel$label == j),]
    center <- c(sum(cel$x * cel$MIDCount)/sum(cel$MIDCount),sum(cel$y * cel$MIDCount)/sum(cel$MIDCount))
    gocoord <- cel[which(cel$geneID %in% go),]
    if(length(gocoord$geneID) == 0){
      celgo[j,i] <- NA
    }
    else{
      ##gomean <- sum(sqrt((gocoord$x - center[1])^2 + (gocoord$y - center[2])^2))##/length(gocoord$x)
      gomean <- sum(sqrt((gocoord$x - center[1])^2 + (gocoord$y - center[2])^2)/gocoord$MIDCount)/length(gocoord$x)
      ##celmean <- sum(sqrt((cel$x - center[1])^2 + (cel$y - center[2])^2))/length(cel$x)
      celmean <- sum(sqrt((cel$x - center[1])^2 + (cel$y - center[2])^2)/cel$MIDCount)/length(cel$x)
      celgo[j,i] <- gomean/celmean
    }
  }
}
celgo <- data.frame(celgo)
write_tsv(celgo, "D:/Graduate/空间转录组课题/17_FP200000440BR_A2_83_HQ03.celgo.txt", col_names = FALSE, na = "")
##celgo <- read.table("D:/Graduate/空间转录组课题/17_FP200000440BR_A2_83_HQ03.celgo.txt",header = F)
m <- rep(0,1772)
nucl <- c()
cyto <- c()
for(i in 1:1772){
  m[i] <- mean(na.omit(celgo[,i]))
  if(is.na(m[i])){next}
  if(m[i] < 0.8){nucl <- append(nucl,i)}
  if(m[i] > 1.2){cyto <- append(cyto,i)}
}
hist(na.omit(celgo[,16]),breaks = seq(0,3,0.05),col = "#BEB8DC",border = "#BEB8DC")
write.table(go2gene$GOID[cyto],"C:/Users/10145/Desktop/cyto.txt",sep = "\n",row.names = F,col.names = F,quote = F)


##找同一基因在不同细胞相对位置
gene1 <- rep(0,length(unique(chpoint$label)))
gene2 <- rep(0,length(unique(chpoint$label)))
gene3 <- rep(0,length(unique(chpoint$label)))

j <- 1
for(i in unique(chpoint$label)){
  cell <- chpoint[chpoint$label == i,]
  center <- c(sum(cell$x)/length(cell$x),sum(cell$y)/length(cell$y))
  gene <- cell[which(cell$geneID == "SoyZH13_19G006300"),]
  if(length(gene$geneID) == 0){
    gene3[j] <- NA
    j <- j+1
  }
  else{
    genemean <- sum(sqrt((gene$x - center[1])^2 + (gene$y - center[2])^2))/length(gene$x)
    celmean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2))/length(cell$x)
    gene3[j] <- genemean/celmean
    j <- j+1
  }
}






##找到每个细胞每种ipr的到中心点的相对距离
celipr <- matrix(0,max(celabel$label),length(ipr2gene$ipr_acc))
for(i in 1:length(celipr[1,])){
  ipracc <- strsplit(ipr2gene$Gene[i],split = ",")[[1]]
  for(j in 1:length(celipr[,1])){
    cel <- celabel[which(celabel$label == j),]
    center <- c(sum(cel$x)/length(cel$x),sum(cel$y)/length(cel$y))
    iprcoord <- cel[which(cel$geneID %in% ipracc),]
    if(length(iprcoord$geneID) == 0){
      celipr[j,i] <- 0
    }
    else{
      iprmean <- sum(sqrt((iprcoord$x - center[1])^2 + (iprcoord$y - center[2])^2))/length(iprcoord$x)
      celmean <- sum(sqrt((cel$x - center[1])^2 + (cel$y - center[2])^2))/length(cel$x)
      celipr[j,i] <- iprmean/celmean
    }
  }
}
celipr <- data.frame(celipr)





for(i in 1:length(celgo[1,])){
  celgo[which(celgo[,i] == Inf),i] <- NA
}
library(cluster)
clusbest <- clusGap(celgo[-c(373,417,3288,3307,771,1062,2649,3207,2914),c(7,16)],FUNcluster = pam,K.max = 10)


celk <- kmeans(celgo[-c(373,417,3288,3307,771,1062),], 6, iter.max = 20, nstart = 6, algorithm = "Hartigan-Wong")
##plot(celgo[-c(373,417,3288,3307,771,1062,2469,3207,2914),],col = celk$cluster,xlab = 1, ylab = 2, pch = 20)
##legend("topright",
##       c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8"), 
##       col = c(1,3,4,5,6,7,8),
##       text.col = "green4", 
##       pch = 20,
##       cex=0.7,
##       ncol=2,
##       bg = "gray90")


for(i in 1:6){
  plot(celabel$x[which(celabel$label %in% which(celk$cluster == i))],celabel$y[which(celabel$label %in% which(celk$cluster == i))]
       ,col = i,xlim = c(min(celabel$x),max(celabel$x)), ylim = c(min(celabel$y),max(celabel$y)),xlab = "",ylab = "",pch = ".")
  par(new = T)
}
legend("topright",
       c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8"), 
       col = c(1,2,3,4,5,6,7,8),
       text.col = "green4", 
       pch = 20,
       cex=0.7,
       ncol=2,
       bg = "gray90")







##箱型图看差异
library(ggsignif)
a <- celgo[,1]
##a <- celgo[which(celgo[,1] != Inf),1]
box <- cbind(a,rep(1,length(a)))
for(i in 2:20){
  a <- celgo[,i]
  ##a <- celgo[which(celgo[,2] != Inf),2]
  a <- cbind(a,rep(i,length(a)))
  box <- rbind(box,a)
}


box <- cbind(gene1[!is.na(gene1)],rep(1,length(gene1[!is.na(gene1)])))
box <- rbind(box,cbind(gene2[!is.na(gene2)],rep(2,length(gene2[!is.na(gene2)]))))
box <- rbind(box,cbind(gene3[!is.na(gene3)],rep(3,length(gene3[!is.na(gene3)]))))
box <- rbind(box,cbind(gene4[!is.na(gene4)],rep(4,length(gene4[!is.na(gene4)]))))
compaired <- list(c(1,2),c(1,3),c(2,3))
ggplot(data.frame(box),aes(box[,2],box[,1],group=box[,2]))+geom_boxplot()+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test)


a <- celgo[,2]
##a <- celgo[which(celgo[,2] != Inf),2]
a <- cbind(a,rep(2,length(a)))
box <- rbind(box,a)

a <- celgo[,3]
##a <- celgo[which(celgo[,3] != Inf),3]
a <- cbind(a,rep(3,length(a)))
box <- rbind(box,a)

a <- celgo[,7]
##a <- celgo[which(celgo[,7] != Inf),7]
a <- cbind(a,rep(4,length(a)))
box <- rbind(box,a)

a <- celgo[,9]
##a <- celgo[which(celgo[,9] != Inf),9]
a <- cbind(a,rep(5,length(a)))
box <- rbind(box,a)

a <- celgo[,11]
##a <- celgo[which(celgo[,11] != Inf),11]
a <- cbind(a,rep(6,length(a)))
box <- rbind(box,a)

a <- celgo[,16]
##a <- celgo[which(celgo[,16] != Inf),16]
a <- cbind(a,rep(7,length(a)))
box <- rbind(box,a)

a <- celgo[,20]
##a <- celgo[which(celgo[,20] != Inf),20]
a <- cbind(a,rep(8,length(a)))
box <- rbind(box,a)

compaired <- list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
ggplot(data.frame(box),aes(box[,2],box[,1],group=box[,2]))+geom_boxplot()+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test)

a <- celgo[,1]
box <- cbind(a,rep(1,length(a)))

a <- celgo[,5]
a <- cbind(a,rep(2,length(a)))
box <- rbind(box,a)

compaired <- list(c(1,2))
ggplot(data.frame(box),aes(V2,a,group=V2))+geom_boxplot()+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)




##chcor <- Coord[which(sqrt((Coord$V2-13475)^2+(Coord$V3-6910)^2) < 50),]##找到在范围内的坐标
##go <- strsplit(go2gene$Gene[1],split = ',')[[1]]
##gocoord <- chcor[which(chcor$V1 %in% go),]
##plot(chcor$V2,chcor$V3,xlim = c(13425,13525),ylim = c(6860,6960),xlab = "",ylab = "", pch = 20,col = rgb(45, 67,121, 100, maxColorValue = 255))
##par(new=TRUE) 
##plot(gocoord$V2,gocoord$V3,col = '#FF7F50',xlim = c(13650,13750),ylim = c(6860,6960),xlab = "",ylab = "",main='GO:0005515在细胞内的表达情况', pch = 20, cex = gocoord$V4)







##ks检验
hist(celgo[which(celgo[,8]!=Inf),8])
ks.test(celgo[which(celgo[,1]!=Inf),1],celgo[which(celgo[,10]!=Inf),10],exact = T,alternative = "two.sided")
##机器学习
set.seed(1)
go1test <- sample(1:4444,444,replace = FALSE)
set.seed(1)
go2test <- sample(1:3177,318,replace = FALSE)
x.test <- c(celgo[which(celgo[,1]!=Inf),1][go1test],celgo[which(celgo[,10]!=Inf),10][go2test])
y.test <- c(rep(1,444),rep(0,318))
x.train <- c(celgo[which(celgo[,1]!=Inf),1][-go1test],celgo[which(celgo[,10]!=Inf),10][-go2test])
y.train <- c(rep(1,4000),rep(0,2859))
##LOGISTIC
x <- x.train
glm.fit <- glm(y.train~x,family = binomial)
x <- x.test
glm.probs <- predict(glm.fit,newdata = data.frame(x),type = 'response')
glm.pred <- rep('0',length(y.test))
glm.pred[glm.probs>0.5] = '1'
table(predict=glm.pred,truth=y.test)
##LDA
library(MASS)
x <- x.train
lda.fit <- lda(y.train~x)
x <- x.test
lda.pred <- predict(lda.fit,data.frame(x))
table(predict=lda.pred$class,truth=y.test)






cluster1 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/Cell_GetExp_gene_01.txt",header = T)
cluster2 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/Cell_GetExp_gene_02.txt",header = T)
cluster3 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/Cell_GetExp_gene_03.txt",header = T)
cluster4 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/Cell_GetExp_gene_04.txt",header = T)
cluster5 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/Cell_GetExp_gene_05.txt",header = T)
cluster6 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/Cell_GetExp_gene_06.txt",header = T)

plot(cluster1$x,cluster1$y,col = rgb(192,192,192, 255, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(cluster2$x,cluster2$y,col = rgb(252,230,202, 255, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(cluster3$x,cluster3$y,col = rgb(202,235,216, 255, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(cluster4$x,cluster4$y,col = rgb(250,240,230, 255, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(cluster5$x,cluster5$y,col = rgb(255,235,205, 255, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(cluster6$x,cluster6$y,col = rgb(255,245,238, 255, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
legend("topright",
       c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6"), 
       col = c("#808A87","#FFE384","#FFC0CB","#40E0D0","#BC8F8F","#00C78C"),
       text.col = "green4", 
       pch = 20,
       cex=0.7,
       ncol=2,
       bg = "gray90")








c1_unq <- cluster1 %>% group_by(x,y) %>% summarise(count = length(geneID))
c1_unq$cluster <- rep(1,length(c1_unq$x))
c2_unq <- cluster2 %>% group_by(x,y) %>% summarise(count = length(geneID))
c2_unq$cluster <- rep(2,length(c2_unq$x))
c3_unq <- cluster3 %>% group_by(x,y) %>% summarise(count = length(geneID))
c3_unq$cluster <- rep(3,length(c3_unq$x))
c4_unq <- cluster4 %>% group_by(x,y) %>% summarise(count = length(geneID))
c4_unq$cluster <- rep(4,length(c4_unq$x))
c5_unq <- cluster5 %>% group_by(x,y) %>% summarise(count = length(geneID))
c5_unq$cluster <- rep(5,length(c5_unq$x))
c6_unq <- cluster6 %>% group_by(x,y) %>% summarise(count = length(geneID))
c6_unq$cluster <- rep(6,length(c6_unq$x))
c_unq <- rbind(c1_unq,c2_unq,c3_unq,c4_unq,c5_unq,c6_unq)
cluster <- type %>% left_join(c_unq) %>% na.omit()











c1len <- cluster1 %>% group_by(x,y) %>% summarise(count = length(label))
c1len <- length(c1len$count)
c1 <- cluster1 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c1len)
c2len <- cluster2 %>% group_by(x,y) %>% summarise(count = length(label))
c2len <- length(c2len$count)
c2 <- cluster2 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c2len)
c3len <- cluster3 %>% group_by(x,y) %>% summarise(count = length(label))
c3len <- length(c3len$count)
c3 <- cluster3 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c3len)
c4len <- cluster4 %>% group_by(x,y) %>% summarise(count = length(label))
c4len <- length(c4len$count)
c4 <- cluster4 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c4len)
c5len <- cluster5 %>% group_by(x,y) %>% summarise(count = length(label))
c5len <- length(c5len$count)
c5 <- cluster5 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c5len)
c6len <- cluster6 %>% group_by(x,y) %>% summarise(count = length(label))
c6len <- length(c6len$count)
c6 <- cluster6 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c6len)
cor_all <- c1 %>% full_join(c2, by = "geneID") %>% 
  full_join(c3, by = "geneID") %>% 
  full_join(c4, by = "geneID") %>% 
  full_join(c5, by = "geneID") %>% 
  full_join(c6, by = "geneID")
colnames(cor_all)[2] <- "c1"
colnames(cor_all)[3] <- "c2"
colnames(cor_all)[4] <- "c3"
colnames(cor_all)[5] <- "c4"
colnames(cor_all)[6] <- "c5"
colnames(cor_all)[7] <- "c6"
cor_all[is.na(cor_all)] <- 0
allhm <- cor(cor_all[,2:7])
allhm.m <-melt(allhm)
ggplot(allhm.m, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(name="Value", low = "#FEA3A2",high = "white") +
  labs(x = "Var1", y = "Var2", title = "correlation") +
  theme(plot.title = element_text(size = 13,hjust = 0.5))+
  theme(aspect.ratio=1)









cluster6 <- cluster[which(cluster$cluster == 6),]









c12cor <- rep(0,37)
c1 <- cluster1[which(18950 <= cluster1$y & cluster1$y <= 19050),]
c1len <- c1 %>% group_by(x,y) %>% summarise(count = length(label))
c1len <- length(c1len$count)
c1 <- c1 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c1len)
for(i in 1:37){
  c2 <- cluster2[which((15410+100*(i-1))<= cluster2$y & cluster2$y <=(15410+100*i)),]## & cluster2$x >= 10939
  c2len <- c2 %>% group_by(x,y) %>% summarise(count = length(label))
  c2len <- length(c2len$count)
  c2 <- c2 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c2len)
  c12 <- c1 %>% full_join(c2, by = "geneID")
  colnames(c12)[2] <- "c1"
  colnames(c12)[3] <- "c2"
  c12[is.na(c12)] <- 0
  c12cor[i] <- cor(c12[,2:3])[2,1]
}
plot(c1$x,c1$y,col = rgb(128,138,135, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(c2$x,c2$y,col = rgb(255,227,132, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")








c35cor <- rep(0,16)
c3 <- cluster3[which(16434 <= cluster3$y & cluster3$y <= 16454),]
c3len <- c3 %>% group_by(x,y) %>% summarise(count = length(label))
c3len <- length(c3len$count)
c3 <- c3 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c3len)
for(i in 1:16){
  c5 <- cluster5[which((16450+100*(i-1))<= cluster5$y & cluster5$y <= (16450+100*i)),]
  c5len <- c5 %>% group_by(x,y) %>% summarise(count = length(label))
  c5len <- length(c5len$count)
  c5 <- c5 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c5len)
  c35 <- c3 %>% full_join(c5, by = "geneID")
  colnames(c35)[2] <- "c3"
  colnames(c35)[3] <- "c5"
  c35[is.na(c35)] <- 0
  c35cor[i] <- cor(c35[,2:3])[2,1]
}
plot(c3$x,c3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(c5$x,c5$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")








c36cor <- rep(0,10)
c3 <- cluster3[which(16354 <= cluster3$y & cluster3$y <= 16454 & cluster3$x >= 11457 & cluster3$x <= 11557),]
c3len <- c3 %>% group_by(x,y) %>% summarise(count = length(label))
c3len <- length(c3len$count)
c3 <- c3 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c3len)
for(i in 1:10){
  c6 <- cluster6[which((16690+100*(i-1))<= cluster6$y & cluster6$y <=(17690+100*i)),]
  c6len <- c6 %>% group_by(x,y) %>% summarise(count = length(label))
  c6len <- length(c6len$count)
  c6 <- c6 %>% group_by(geneID) %>% summarise(Count = sum(MIDCount)/c6len)
  c36 <- c3 %>% full_join(c6, by = "geneID")
  colnames(c36)[2] <- "c1"
  colnames(c36)[3] <- "c2"
  c36[is.na(c36)] <- 0
  c36cor[i] <- cor(c36[,2:3])[2,1]
}










temp6 <- cluster6[which(16697<= cluster6$y & cluster6$y <=16797),]
temp5 <- cluster5[which(16667<= cluster5$y & cluster5$y <=16797 & cluster5$x >=11307 & cluster5$x <=11675),]
plot(temp5$x,temp5$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp6$x,temp6$y,col = rgb(0,199,140, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp5.1 <- temp5 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5.1[,c(1,3)]
temp6.1 <- temp6.1[,c(1,3)]
temp <- temp5.1 %>% full_join(temp6.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##接受原假设


temp5 <- cluster5[which(cluster5$y <= 16554),]
temp3 <- cluster3[which(cluster3$y <= 16554 & cluster3$y >= 16386 & cluster3$x >=11341 & cluster3$x <=11679),]
plot(temp3$x,temp3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp5$x,temp5$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp3.1 <- temp3 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5.1[,c(1,3)]
temp3.1 <- temp3.1[,c(1,3)]
temp <- temp5.1 %>% full_join(temp3.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##接受原假设


temp4 <- cluster4[which(cluster4$y < 17500 & cluster4$x >= 11050 & cluster4$x <=11300 & cluster4$y > 17000),]
temp6 <- cluster6[which(17250 <= cluster6$y & cluster6$x <=11250),]
plot(temp4$x,temp4$y,col = rgb(64,224,205, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp6$x,temp6$y,col = rgb(0,199,140, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp4.1 <- temp4 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6.1[,c(1,3)]
temp4.1 <- temp4.1[,c(1,3)]
temp <- temp4.1 %>% full_join(temp6.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##接受原假设






##背景观察
temp3 <- cluster3[which(cluster3$x <= 11600 & cluster3$x >= 11400 & cluster3$y <=18600 & cluster3$y >= 18400),]
temp6 <- cluster6[which(cluster6$x <= 11600 & cluster6$x >= 11400 & cluster6$y <=17200 & cluster6$y >= 17000),]
plot(temp3$x,temp3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp6$x,temp6$y,col = rgb(0,199,140, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp3.1 <- temp3 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6.1[,c(1,3)]
temp3.1 <- temp3.1[,c(1,3)]
temp <- temp3.1 %>% full_join(temp6.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##拒绝原假设





temp6 <- cluster6[which(cluster6$x <= 11600 & cluster6$x >= 11400 & cluster6$y <=17200 & cluster6$y >= 17000),]
temp4 <- cluster4[which(cluster4$y < 17570 & cluster4$x >= 11050 & cluster4$x <=11300 & cluster4$y > 17000),]
plot(temp4$x,temp4$y,col = rgb(64,224,205, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp6$x,temp6$y,col = rgb(0,199,140, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp4.1 <- temp4 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6.1[,c(1,3)]
temp4.1 <- temp4.1[,c(1,3)]
temp <- temp4.1 %>% full_join(temp6.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##拒绝原假设




temp5 <- cluster5[which(cluster5$x >= 10869 & cluster5$x <= 11069 & cluster5$y >= 17000 & cluster5$y <= 17205),]
temp6 <- cluster6[which(cluster6$x <= 11600 & cluster6$x >= 11400 & cluster6$y <=17200 & cluster6$y >= 17000),]
plot(temp5$x,temp5$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp6$x,temp6$y,col = rgb(0,199,140, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp5.1 <- temp5 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp6.1 <- temp6.1[,c(1,3)]
temp5.1 <- temp5.1[,c(1,3)]
temp <- temp5.1 %>% full_join(temp6.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##接受原假设





temp1.1 <- cluster1 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp2.1 <- cluster2 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp1.1 <- temp1.1[,c(1,3)]
temp2.1 <- temp2.1[,c(1,3)]
temp <- temp1.1 %>% full_join(temp2.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##拒绝原假设






temp3 <- cluster3[which(cluster3$x <= 11570 & cluster3$x >= 11400 & cluster3$y <=16200 & cluster3$y >= 16000),]
temp4 <- cluster4[which(cluster4$y < 17570 & cluster4$x >= 11050 & cluster4$x <=11300 & cluster4$y > 17000),]
plot(temp3$x,temp3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp4$x,temp4$y,col = rgb(64,224,205, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp4.1 <- temp4 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp3.1 <- temp3 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp4.1 <- temp4.1[,c(1,3)]
temp3.1 <- temp3.1[,c(1,3)]
temp <- temp4.1 %>% full_join(temp3.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##p-value = 0.007412
temp3 <- cluster3[which(cluster3$x <= 11570 & cluster3$x >= 11400 & cluster3$y <=18600 & cluster3$y >= 18400),]
temp4 <- cluster4[which(cluster4$y < 17570 & cluster4$x >= 11050 & cluster4$x <=11300 & cluster4$y > 17000),]
plot(temp3$x,temp3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp4$x,temp4$y,col = rgb(64,224,205, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp4.1 <- temp4 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp3.1 <- temp3 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp4.1 <- temp4.1[,c(1,3)]
temp3.1 <- temp3.1[,c(1,3)]
temp <- temp4.1 %>% full_join(temp3.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##接受




temp5 <- cluster5[which(cluster5$y <= 16585),]
temp3 <- cluster3[which(cluster3$x <= 11570 & cluster3$x >= 11400 & cluster3$y <=16200 & cluster3$y >= 16000),]
plot(temp3$x,temp3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp5$x,temp5$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp3.1 <- temp3 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5.1[,c(1,3)]
temp3.1 <- temp3.1[,c(1,3)]
temp <- temp5.1 %>% full_join(temp3.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##p-value = 0.005258
temp5 <- cluster5[which(cluster5$y <= 16585),]
temp3 <- cluster3[which(cluster3$x <= 11570 & cluster3$x >= 11400 & cluster3$y <=18600 & cluster3$y >= 18400),]
plot(temp3$x,temp3$y,col = rgb(255,192,203, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(temp5$x,temp5$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
temp3.1 <- temp3 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5 %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
temp5.1 <- temp5.1[,c(1,3)]
temp3.1 <- temp3.1[,c(1,3)]
temp <- temp5.1 %>% full_join(temp3.1, by = "geneID")
wilcox.test(temp$Count.x,temp$Count.y,paired = T)##接受




##新生RNA研究
type <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq.gem",header = T)
nucl <- type[which(type$label != -1 & type$Type == "U"),]
cyto <- type[which(type$label != -1 & type$Type == "S"),]
N <- type[which(type$label != -1 & type$Type == "N"),]
all <- rbind(nucl,cyto)
all <- rbind(all,N)
plot(cyto$x,cyto$y,col = rgb(64,224,205, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
par(new=TRUE)
plot(nucl$x,nucl$y,col = rgb(188,143,143, 100, maxColorValue = 255),xlim = c(10139,12746), ylim = c(15413,19050),xlab = "",ylab = "",pch = ".")
legend("topright",
       c("Cytoplasm","Nucleus"), 
       col = c("#40E0D0","#BC8F8F"),
       text.col = "green4", 
       pch = 20,
       cex=0.7,
       ncol=1,
       bg = "gray90")






cellC <- all[which(all$label == 3385),]
cytoC <- cyto[which(cyto$label == 3385),]
nuclC <- nucl[which(nucl$label == 3385),]
NC <- N[which(N$label == 3385),]
cel_center <- c(sum(cellC$x * cellC$MIDCount)/sum(cellC$MIDCount),sum(cellC$y * cellC$MIDCount)/sum(cellC$MIDCount))
##
plot(NC$x,NC$y,xlim = c(min(cellC$x),max(cellC$x)), ylim = c(min(cellC$y),max(cellC$y)),xlab = "",ylab = "",pch = 20,cex = NC$MIDCount)
par(new=TRUE)
plot(cel_center[1],cel_center[2],xlim = c(min(cellC$x),max(cellC$x)), ylim = c(min(cellC$y),max(cellC$y)),xlab = "",ylab = "")
par(new=TRUE)
plot(cytoC$x,cytoC$y,col = rgb(64,224,205, 255, maxColorValue = 255),xlim = c(min(cellC$x),max(cellC$x)), ylim = c(min(cellC$y),max(cellC$y)),xlab = "",ylab = "",pch = 20,cex = cytoC$MIDCount)
par(new=TRUE)
plot(nuclC$x,nuclC$y,col = rgb(188,143,143, 255, maxColorValue = 255),xlim = c(min(cellC$x),max(cellC$x)), ylim = c(min(cellC$y),max(cellC$y)),xlab = "",ylab = "",pch = 20,cex = nuclC$MIDCount)
legend("topright",
       c("Cytoplasm","Nucleus"), 
       col = c("#40E0D0","#BC8F8F"),
       text.col = "green4", 
       pch = 20,
       cex=0.7,
       ncol=1,
       bg = "gray90")











cyto_mean <- rep(0,3589)
nucl_mean <- rep(0,3589)
UdS <- rep(0,3589)
##cyto_cell$MIDCount  nucl_cell$MIDCount
for(i in 1:3589){
  cell <- all[all$label == i,]
  cyto_cell <- cyto[cyto$label == i,]
  nucl_cell <- nucl[nucl$label == i,]
  if(length(cyto_cell$geneID)<5 | length(nucl_cell$geneID)<5){##if(length(cyto_cell$geneID)+length(nucl_cell$geneID)<10){
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

wilcox.test(cyto_mean,nucl_mean,paired = T)

hist(UdS,breaks = seq(0,3,0.01),col = "#BEB8DC",border = "#BEB8DC")
abline(v=1,lwd=1,col="black",lty = 3)

(mean(nucl_mean[!is.na(nucl_mean)])+mean(cyto_mean[!is.na(cyto_mean)]))/2

UdS_sort <- na.omit(UdS_sort[order(UdS_sort[,1],decreasing = F),])



cell <- type[which(type$label == 1006),]##1006
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
plot(cyto$x,cyto$y,xlim = c(min(cell$x),max(cell$x)),col = rgb(142,139,254, 255, maxColorValue = 255),ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20,cex = cyto$MIDCount)
par(new=TRUE)
plot(nucl$x,nucl$y,xlim = c(min(cell$x),max(cell$x)),col = rgb(254,163,162, 255, maxColorValue = 255),ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20,cex = nucl$MIDCount)
par(new=TRUE)
plot(center[1],center[2],xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
legend("topright",
       c("Cyto","Nucl"), 
       col = c("#8E8BFE","#FEA3A2"),
       pch = 20,
       cex=1,
       ncol=2)
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
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)#换个好看的颜色
hM <- format(round(allhm, 3))#对数据保留2位小数

heatmap.2(allhm,
           trace="none",#不显示trace
           col = coul,
           density.info = "none",#图例取消density
           key.xlab ='Correlation',
           key.title = "",
           cexRow = 1,cexCol = 1,#修改横纵坐标字体
           Rowv = F,Colv = F, #去除聚类
           margins = c(6, 6),
           cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
)
library(ggplot2)
library(reshape2)
library(ggpubr)
a <- ggplot(data = cyto, mapping = aes(x = x, y = y))+ 
  geom_point(alpha = 0.2, size=cyto$MIDCount, color = "#8E8BFE") +
  geom_point(data = nucl, mapping = aes(x = x, y = y),alpha = 1, size=nucl$MIDCount, color = "#FEA3A2") + 
  theme(aspect.ratio=1)
allhm.m <-melt(allhm)
b <- ggplot(allhm.m, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(name="Value", low = "#FEA3A2",high = "white") +
  labs(x = "Var1", y = "Var2", title = "correlation") +
  theme(plot.title = element_text(size = 13,hjust = 0.5))+
  theme(aspect.ratio=1)
ggpubr::ggarrange(a,b,nrow = 1,ncol = 2,labels = c("A","B"))
##ggplot(data = type[which(type$label != -1),], mapping = aes(x = x, y = y)) + geom_point(alpha = type[which(type$label != -1),]$MIDCount, size=0.05, color = "green4")













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









cyto_mean <- rep(0,3589)
nucl_mean <- rep(0,3589)
UdS <- rep(0,3589)
##cyto_cell$MIDCount  nucl_cell$MIDCount
for(i in 1:3589){
  cell <- all[all$label == i,]
  cyto_cell <- cyto[cyto$label == i,]
  nucl_cell <- nucl[nucl$label == i,]
  if(length(cyto_cell$geneID)<5 | length(nucl_cell$geneID)<5){
    cyto_mean[i] <- NaN
    nucl_mean[i] <- NaN
  }
  else{
    center <- c(sum(cell$x * cell$MIDCount)/sum(cell$MIDCount),sum(cell$y * cell$MIDCount)/sum(cell$MIDCount))
    cell_mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)/cell$MIDCount)/length(cell$x)
    cyto_mean[i] <- sum(sqrt((cyto_cell$x - center[1])^2 + (cyto_cell$y - center[2])^2)/cyto_cell$MIDCount)/length(cyto_cell$x)/cell_mean
    nucl_mean[i] <- sum(sqrt((nucl_cell$x - center[1])^2 + (nucl_cell$y - center[2])^2)/nucl_cell$MIDCount)/length(nucl_cell$x)/cell_mean
  }
}
UdS <- nucl_mean/cyto_mean 
UdS_sort <- cbind(UdS,seq(1,3589,1))

wilcox.test(cyto_mean,nucl_mean,paired = T)

hist(UdS,breaks = seq(0,2.5,0.01),col = "#BEB8DC",border = "#BEB8DC")
abline(v=1,lwd=1,col="black",lty = 3)

(mean(nucl_mean[!is.na(nucl_mean)])+mean(cyto_mean[!is.na(cyto_mean)]))/2

UdS_sort <- na.omit(UdS_sort[order(UdS_sort[,1],decreasing = F),])



##原始定义
nucl <- type[which(type$label != -1 & type$Type == "U"),]
cyto <- type[which(type$label != -1 & type$Type == "S"),]
all <- rbind(nucl,cyto)
cyto_mean <- rep(0,3589)
nucl_mean <- rep(0,3589)
for(i in 1:3589){
  cell <- all[all$label == i,]
  cyto_cell <- cyto[cyto$label == i,]
  nucl_cell <- nucl[nucl$label == i,]
  center <- c(sum(cell$x)/length(cell$label),sum(cell$y)/length(cell$label))
  cell_mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)/cell$MIDCount)/length(cell$x)
  cyto_mean[i] <- sum(sqrt((cyto_cell$x - center[1])^2 + (cyto_cell$y - center[2])^2)/cyto_cell$MIDCount)/length(cyto_cell$x)/cell_mean
  nucl_mean[i] <- sum(sqrt((nucl_cell$x - center[1])^2 + (nucl_cell$y - center[2])^2)/nucl_cell$MIDCount)/length(nucl_cell$x)/cell_mean
}
UdS <- nucl_mean/cyto_mean
write.table(Uds_sort[order(Uds_sort[,1],decreasing = T),][,2][1:100],"D:/Graduate/空间转录组课题/数据表格/26_02hq_UdS_sort100.txt",sep = "\n",row.names = F,col.names = F,quote = F)




##修改中心点加入N
nucl <- type[which(type$label != -1 & type$Type == "U"),]
cyto <- type[which(type$label != -1 & type$Type == "S"),]
N <- type[which(type$label != -1 & type$Type == "N"),]
all <- rbind(nucl,cyto)
all <- rbind(all,N)
cyto_mean <- rep(0,3589)
nucl_mean <- rep(0,3589)
for(i in 1:3589){
  cell <- all[all$label == i,]
  cyto_cell <- cyto[cyto$label == i,]
  nucl_cell <- nucl[nucl$label == i,]
  center <- c(sum(cell$x)/length(cell$label),sum(cell$y)/length(cell$label))
  cell_mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)/cell$MIDCount)/length(cell$x)
  cyto_mean[i] <- sum(sqrt((cyto_cell$x - center[1])^2 + (cyto_cell$y - center[2])^2)/cyto_cell$MIDCount)/length(cyto_cell$x)/cell_mean
  nucl_mean[i] <- sum(sqrt((nucl_cell$x - center[1])^2 + (nucl_cell$y - center[2])^2)/nucl_cell$MIDCount)/length(nucl_cell$x)/cell_mean
}
UdS <- nucl_mean/cyto_mean
hist(UdS[c(1:2338,2340:3589)],breaks = seq(0,5.5,0.05))




##修改中心点定义（包括N）
nucl <- type[which(type$label != -1 & type$Type == "U"),]
cyto <- type[which(type$label != -1 & type$Type == "S"),]
N <- type[which(type$label != -1 & type$Type == "N"),]
all <- rbind(nucl,cyto)
all <- rbind(all,N)
cyto_mean <- rep(0,3589)
nucl_mean <- rep(0,3589)
for(i in 1:3589){
  cell <- all[all$label == i,]
  cyto_cell <- cyto[cyto$label == i,]
  nucl_cell <- nucl[nucl$label == i,]
  center <- c(sum(cell$x*cell$label)/sum(cell$label),sum(cell$y*cell$label)/sum(cell$label))
  cell_mean <- sum(sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)/cell$MIDCount)/length(cell$x)
  cyto_mean[i] <- sum(sqrt((cyto_cell$x - center[1])^2 + (cyto_cell$y - center[2])^2)/cyto_cell$MIDCount)/length(cyto_cell$x)/cell_mean
  nucl_mean[i] <- sum(sqrt((nucl_cell$x - center[1])^2 + (nucl_cell$y - center[2])^2)/nucl_cell$MIDCount)/length(nucl_cell$x)/cell_mean
}





ch <- UdS_sort[which(UdS_sort[,1] < 1),2]##
plot(1:length(ch),nucl_mean[sort(ch)],col = rgb(254,163,162, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20, ylim = c(0,3))
par(new=TRUE)
plot(1:length(ch),cyto_mean[sort(ch)],col = rgb(142,139,254, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20, ylim = c(0,3))
abline(h = c(0.5,1,1.5,2), v = c(0,500,1000,1500,2000,2500), col = "lightgray", lty = 3)
legend("topright",
       c("Splice","Unsplice"), 
       col = c("#8E8BFE","#FEA3A2"),
       pch = 20,
       cex=1,
       ncol=2)
abline(h=mean(nucl_mean[!is.na(nucl_mean)]),lwd=1,col="black")
abline(h=mean(cyto_mean[!is.na(cyto_mean)]),lwd=1,col="black",lty = 3)
h1 <- nucl_mean[ch]
abline(h=mean(h1[!is.na(h1)]),lwd=1,col="black")
h2 <- cyto_mean[ch]
abline(h=mean(h2[!is.na(h2)]),lwd=1,col="black",lty = 3)











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
cyto_ <- read.table("D:/Graduate/空间转录组课题/数据表格/26_02hq_cyto.txt",header = T)
nucl_ <- read.table("D:/Graduate/空间转录组课题/数据表格/26_02hq_nucl.txt",header = T)
anno <- read_tsv("D:/Graduate/空间转录组课题/数据表格/soybean.interproscan.type.tsv", na = "N/A", col_names = c("geneID","ipr_acc","ipr_desc")) %>% na.omit()
cyto_ <- cyto_ %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
nucl_ <- nucl_ %>% group_by(geneID) %>% 
  summarise(label = str_c(label, collapse = ",")) %>%
  mutate(Count = str_count(label, ",")+1) %>% arrange(desc(Count))
cyto_ <- cyto_[,c(1,3)]
nucl_ <- nucl_[,c(1,3)]
temp <- cyto_ %>% full_join(nucl_, by = "geneID")
wilcox.test(temp$Count.x/107289,temp$Count.y/111478,paired = T)
plot(temp$Count.y[1:100]/111478,type = "o",col = rgb(254,163,162, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20,ylim = c(0,0.1))
par(new=TRUE)
plot(temp$Count.x[1:100]/107289,type = "o",col = rgb(142,139,254, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20,ylim = c(0,0.1))
abline(h = 0:7, v = seq(0,100,10), col = "lightgray", lty = 3)
legend("topright",
       c("Splice","Unsplice"), 
       col = c("#8E8BFE","#FEA3A2"),
       pch = "——",
       cex=1,
       ncol=2)
table <- temp[which(temp$Count.x > 100 & temp$Count.y >100),]
Anno <- temp %>% left_join(anno, by = "geneID")
write.table(temp,"D:/Graduate/空间转录组课题/数据表格/26_02hq_C_N.txt",sep = "\t",row.names = F,col.names = T,quote = F)
blastp <- read.table("D:/Graduate/空间转录组课题/数据表格/soybean.protein.1.fasta.blastp.anno",header = F,sep = '\t',quote = "")
for(i in 1:length(blastp$V1)){
  blastp$V1[i] <- unlist(strsplit(blastp$V1[i],split = '[.]'))[1]
}

colnames(blastp)[1] <- "geneID"
colnames(blastp)[2] <- "ipr_acc"
colnames(blastp)[3] <- "ipr_desc"
write.table(blastp,"D:/Graduate/空间转录组课题/数据表格/blastp.anno.txt",sep = "\t",row.names = F,col.names = T,quote = F)








Anno <- temp %>% left_join(anno, by = "geneID")
nu <- Anno[which(Anno$Count.y > Anno$Count.x),]
cy <- Anno[which(Anno$Count.y < Anno$Count.x),]


NU <- nu %>% select(geneID,ipr_acc,ipr_desc,Count.y) %>% group_by(geneID,ipr_acc) %>%
  summarise(ipr_desc = ipr_desc,Count.y = Count.y) %>% group_by(ipr_acc,ipr_desc) %>% summarise(Count=sum(Count.y))%>%
  arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))
CY <- cy %>% select(geneID,ipr_acc,ipr_desc,Count.x) %>% group_by(geneID,ipr_acc) %>%
  summarise(ipr_desc = ipr_desc,Count.x = Count.x) %>% group_by(ipr_acc,ipr_desc) %>% summarise(Count=sum(Count.x))%>%
  arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))


NU <- NU[which(NU$ipr_desc %in% setdiff(NU$ipr_desc,CY$ipr_desc)),]
CY <- CY[which(CY$ipr_desc %in% setdiff(CY$ipr_desc,NU$ipr_desc)),]
ggplot(CY) +
  geom_bar(aes(x = ipr_desc, y = Percent, fill = ipr_desc), stat = "identity") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.009),name = "Percent of Domain") +
  scale_x_discrete(limits = CY$ipr_desc[1:20], name = NULL) + scale_fill_discrete(guide = FALSE)+
  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))






temp$p_value <- rep(0,length(temp$geneID))
temp <- na.omit(temp)
for(i in 1:length(temp$geneID)){
  table <- matrix(c(temp[i,]$Count.x,(temp[i,]$Count.x + temp[i,]$Count.y/111478*107289)/2,temp[i,]$Count.y/111478*107289,(temp[i,]$Count.x + temp[i,]$Count.y/111478*107289)/2),nrow = 2,ncol = 2)
  temp$p_value[i] <- chisq.test(table)$p.value
}
Anno <- temp %>% left_join(blastp, by = "geneID")
ipr <- Anno %>% group_by(geneID) %>% summarise(nucl = unique(Count.y),cyto = unique(Count.x),p_value = unique(p_value),ipr = str_c(unique(ipr_acc), collapse = ","),ipr_desc = str_c(unique(ipr_desc), collapse = ","))
write.table(Anno,"D:/Graduate/空间转录组课题/数据表格/26_02hq_anno.txt",sep = "\t",row.names = F,col.names = T,quote = F)






##cyto 107289 nucl 111478
cyto_gene <- temp[which(temp$Count.x>temp$Count.y/111478*107289),]
nucl_gene <- temp[which(temp$Count.x<temp$Count.y/111478*107289),]
cyto_gene <- cyto_gene %>% left_join(blastp, by = "geneID")
nucl_gene <- nucl_gene %>% left_join(blastp, by = "geneID")
write.table(cyto_gene,"D:/Graduate/空间转录组课题/数据表格/26_02hq_cyto_gene.txt",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(nucl_gene,"D:/Graduate/空间转录组课题/数据表格/26_02hq_nucl_gene.txt",sep = "\t",row.names = F,col.names = T,quote = F)










hq2602 <- type[which(type$label != -1),]
topcyto <- hq2602[which(hq2602$geneID %in% cyto_gene[order(cyto_gene[,4],decreasing = F),]$geneID[1:50]),]



cell <- type[which(type$label == 3201),]
center <- c(sum(cell$x)/length(cell$label),sum(cell$y)/length(cell$label))
f=seq(0,2*pi,0.001)
x=20*sin(f)+center[1]
y=20*cos(f)+center[2]
plot(cell$x,cell$y,col = rgb(0,0,0, 100, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
par(new=TRUE)
plot(x,y,type='l',col = rgb(142,139,254, 255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)),ylim = c(min(cell$y),max(cell$y)),asp=1)
par(new=TRUE)
x=8*sin(f)+center[1]
y=8*cos(f)+center[2]
plot(x,y,type='l',col = rgb(254,163,162, 255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)),ylim = c(min(cell$y),max(cell$y)),asp=1)
par(new=TRUE)
x=14*sin(f)+center[1]
y=14*cos(f)+center[2]
plot(x,y,type='l',xlim = c(min(cell$x),max(cell$x)),ylim = c(min(cell$y),max(cell$y)),asp=1)
write.table(topnucl,"D:/Graduate/空间转录组课题/数据表格/26_02hq_topnucl.gem",sep = "\t",row.names = F,col.names = T,quote = F)








##将细胞从内向外分成十份计算每一份U、S数量（比例）
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

boxplot(boxmatrix2)









cell <- type[which(type$label == 1006),]
U <- cell[which(cell$Type == "U"),]
center <- c(sum(U$x * U$MIDCount)/sum(U$MIDCount),sum(U$y * U$MIDCount)/sum(U$MIDCount))
dis <- sqrt((cell$x - center[1])^2 + (cell$y - center[2])^2)
for(j in 1:5){
  assign(paste("cont",j,sep = ""),cell[which((j-1)*max(dis)/5 < dis & dis < j*max(dis)/5),])
}
plot(cont1$x,cont1$y,col = rgb(254,163,162,255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
par(new = TRUE)
plot(cont2$x,cont2$y,col = rgb(3,168,158,255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
par(new = TRUE)
plot(cont3$x,cont3$y,col = rgb(153,51,250,255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
par(new = TRUE)
plot(cont4$x,cont4$y,col = rgb(218,112,214,255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
par(new = TRUE)
plot(cont5$x,cont5$y,col = rgb(142,139,254,255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)








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

name <- "SoyZH13_18G003700"
gene <- c(length(which(cont1$geneID == name))/length(cont1$geneID))
gene <- append(gene,length(which(cont2$geneID == name))/length(cont2$geneID))
gene <- append(gene,length(which(cont3$geneID == name))/length(cont3$geneID))
gene <- append(gene,length(which(cont4$geneID == name))/length(cont4$geneID))
gene <- append(gene,length(which(cont5$geneID == name))/length(cont5$geneID))
##gene <- append(gene,length(which(cont6$geneID == name))/length(cont6$geneID))
##gene <- append(gene,length(which(cont7$geneID == name))/length(cont7$geneID))
##gene <- append(gene,length(which(cont8$geneID == name))/length(cont8$geneID))
##gene <- append(gene,length(which(cont9$geneID == name))/length(cont9$geneID))
##gene <- append(gene,length(which(cont10$geneID == name))/length(cont10$geneID))
Midc<- sum(cont1[which(cont1$geneID == name),]$MIDCount)/sum(cont1$MIDCount)
Midc <- append(Midc,sum(cont2[which(cont2$geneID == name),]$MIDCount))/sum(cont2$MIDCount)
Midc <- append(Midc,sum(cont3[which(cont3$geneID == name),]$MIDCount))/sum(cont3$MIDCount)
Midc <- append(Midc,sum(cont4[which(cont4$geneID == name),]$MIDCount))/sum(cont4$MIDCount)
Midc <- append(Midc,sum(cont5[which(cont5$geneID == name),]$MIDCount))/sum(cont5$MIDCount)
##Midc <- append(Midc,sum(cont6[which(cont6$geneID == name),]$MIDCount))/sum(cont6$MIDCount)
##Midc <- append(Midc,sum(cont7[which(cont7$geneID == name),]$MIDCount))/sum(cont7$MIDCount)
##Midc <- append(Midc,sum(cont8[which(cont8$geneID == name),]$MIDCount))/sum(cont8$MIDCount)
##Midc <- append(Midc,sum(cont9[which(cont9$geneID == name),]$MIDCount))/sum(cont9$MIDCount)
##Midc <- append(Midc,sum(cont10[which(cont10$geneID == name),]$MIDCount))/sum(cont10$MIDCount)
plot(gene,type = "o",col = rgb(254,163,162, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20)##nucl
plot(Midc,type = "o",col = rgb(254,163,162, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20)

plot(gene,type = "o",col = rgb(142,139,254, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20)##cyto
plot(Midc,type = "o",col = rgb(142,139,254, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20)


plot(Ucont,type = "o",col = rgb(254,163,162, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20,ylim = c(0,0.2))
par(new=TRUE)
plot(Scont,type = "o",col = rgb(142,139,254, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20,ylim = c(0,0.2))



plot(Midc,type = "o",col = rgb(254,163,162, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20)##,ylim = c(0.06,0.15)
plot(Midc,type = "o",col = rgb(142,139,254, 255, maxColorValue = 255),xlab = "",ylab = "",pch = 20)









U <- type[which(type$Type == "U" & type$label != -1 & type$label %in% ch),]
NS <- type[which(type$Type != "U" & type$label != -1 & type$label %in% ch),]

plot(NS$x,NS$y,col = rgb(142,139,254, 10, maxColorValue = 255),xlim = c(min(type$x),max(type$x)), ylim = c(min(type$y),max(type$y)),xlab = "",ylab = "",pch = ".")
par(new = TRUE)
plot(U$x,U$y,col = rgb(0,0,0, U$MIDCount/max(U$MIDCount)*255, maxColorValue = 255),xlim = c(min(type$x),max(type$x)), ylim = c(min(type$y),max(type$y)),xlab = "",ylab = "",pch = ".")
par(new = TRUE)
plot(11546,16681,col = rgb(254,163,162, 255, maxColorValue = 255),xlim = c(min(type$x),max(type$x)), ylim = c(min(type$y),max(type$y)),xlab = "",ylab = "",pch = ".")





nucl <- cell[which(cell$Type == "U"),]
plot(cell$x,cell$y,col = rgb(255,0,0, cell$MIDCount/max(cell$MIDCount)*255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
par(new=TRUE)
plot(nucl$x,nucl$y,col = rgb(0,0,0, nucl$MIDCount/max(nucl$MIDCount)*255, maxColorValue = 255),xlim = c(min(cell$x),max(cell$x)), ylim = c(min(cell$y),max(cell$y)),xlab = "",ylab = "",pch = 20)
write.table(hq2602[,2:3],"D:/Graduate/空间转录组课题/数据表格/HQ_02/coord.txt",sep = "\t",row.names = F,col.names = T,quote = F)








##Pearson corelation
hq2602 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq.gem",header = T)
type <- hq2602[which(hq2602$label != -1),]
U <- hq2602[which(hq2602$Type == "U" & hq2602$label %in% ch[1:100]),c(1,4)]##
S <- hq2602[which(hq2602$Type == "S" & hq2602$label %in% ch[1:100]),c(1,4)]##
wall <- hq2602[which(hq2602$label == -1),]
wall <- wall[which(wall$x <= max(type$x) & wall$x >= min(type$x)),]
wall <- wall[which(wall$y <= max(type$y) & wall$y >= min(type$y)),]
wall <- wall[which(-2357*wall$x-500*wall$y+32555700 <= 0),]
wall <- wall[which(2041*wall$x-1233*wall$y+501025 >= 0),]
wall <- wall[which(-1416*wall$x-986*wall$y+35484700 >= 0),]
wall <- wall[which(2482*wall$x-615*wall$y-20682025 <= 0),]
wall <- wall[which(1657*wall$x-1450*wall$y+3316100 <= 0),]
wall <- wall[which(1157*wall$x-2000*wall$y+17906100 <= 0),]
##wall <- wall[,c(1,4)]
nucl <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq_nucl.txt",header = T,sep = '\t')
nucl <- nucl[which(nucl$Type != "N" & nucl$label %in% ch[1:100]),c(1,4)]##which(hq2602$label %in% ch)
cyto <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq_cyto.txt",header = T,sep = '\t')
cyto <- cyto[which(cyto$Type != "N" & cyto$label %in% ch[1:100]),c(1,4)]##which(hq2602$label %in% ch)
##nucl_gene <- c("SoyZH13_08G123800","SoyZH13_05G160200","SoyZH13_15G093000","SoyZH13_13G194100","SoyZH13_05G160200","SoyZH13_17G022400","SoyZH13_07G234800","SoyZH13_08G123800","SoyZH13_16G150400","SoyZH13_02G077000","SoyZH13_17G100400","SoyZH13_07G222600","SoyZH13_17G033500","SoyZH13_16G180400","SoyZH13_15G148800","SoyZH13_09G137600","SoyZH13_09G045100","SoyZH13_02G147000","SoyZH13_16G032300","SoyZH13_18G026200","SoyZH13_11G227800","SoyZH13_09G155600","SoyZH13_20G213200","SoyZH13_04G178100","SoyZH13_08G053200","SoyZH13_07G065100","SoyZH13_09G233700","SoyZH13_19G191400","SoyZH13_06G170200","SoyZH13_04G165800","SoyZH13_03G186600")
##nucl <- hq2602[which(hq2602$geneID %in% nucl_gene),]
##chlo_gene <- c("SoyZH13_06G001400","SoyZH13_04G001700","SoyZH13_14G184900","SoyZH13_02G213800","SoyZH13_05G054400","SoyZH13_11G055100","SoyZH13_17G166900","SoyZH13_20G030500")
##chlo <- hq2602[which(hq2602$geneID %in% chlo_gene),]
##统计基因类型及表达量
S <- S %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount)/1402)
U <- U %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount)/1108)
nucl <- nucl %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount)/1137)
cyto <- cyto %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount)/1373)
##wall <- wall %>% group_by(geneID) %>% 
##  summarise(Count = sum(MIDCount))

##合并所有
all <- U %>% full_join(S, by = "geneID") %>% 
  full_join(nucl, by = "geneID") %>% 
  full_join(cyto, by = "geneID") %>%
##  full_join(wall, by = "geneID") %>%
  na.omit()



colnames(all)[2] <- "U"
colnames(all)[3] <- "S"
colnames(all)[4] <- "nucl"
colnames(all)[5] <- "cyto"
##colnames(all)[6] <- "wall"

allhm <- cor(all[,2:5])

library(gplots)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)#换个好看的颜色
hM <- format(round(allhm, 3))#对数据保留2位小数

heatmap.2(allhm,
          trace="none",#不显示trace
          col = coul,
          density.info = "none",#图例取消density
          key.xlab ='Correlation',
          key.title = "",
          cexRow = 1,cexCol = 1,#修改横纵坐标字体
          Rowv = F,Colv = F, #去除聚类
          margins = c(6, 6),
          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
)





##分区pearson
library(gplots)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)#换个好看的颜色
type <- hq2602[which(hq2602$label != -1),]
U <- hq2602[which(hq2602$Type == "U"),]
S <- hq2602[which(hq2602$Type == "S"),]
wall <- hq2602[which(hq2602$label == -1),]
wall <- wall[which(wall$x <= max(type$x) & wall$x >= min(type$x)),]
wall <- wall[which(wall$y <= max(type$y) & wall$y >= min(type$y)),]
wall <- wall[which(-2357*wall$x-500*wall$y+32555700 <= 0),]
wall <- wall[which(2041*wall$x-1233*wall$y+501025 >= 0),]
wall <- wall[which(-1416*wall$x-986*wall$y+35484700 >= 0),]
wall <- wall[which(2482*wall$x-615*wall$y-20682025 <= 0),]
wall <- wall[which(1657*wall$x-1450*wall$y+3316100 <= 0),]
wall <- wall[which(1157*wall$x-2000*wall$y+17906100 <= 0),]
nucl <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq_nucl.txt",header = T,sep = '\t')
##nucl_gene <- c("SoyZH13_08G123800","SoyZH13_05G160200","SoyZH13_15G093000","SoyZH13_13G194100","SoyZH13_05G160200","SoyZH13_17G022400","SoyZH13_07G234800","SoyZH13_08G123800","SoyZH13_16G150400","SoyZH13_02G077000","SoyZH13_17G100400","SoyZH13_07G222600","SoyZH13_17G033500","SoyZH13_16G180400","SoyZH13_15G148800","SoyZH13_09G137600","SoyZH13_09G045100","SoyZH13_02G147000","SoyZH13_16G032300","SoyZH13_18G026200","SoyZH13_11G227800","SoyZH13_09G155600","SoyZH13_20G213200","SoyZH13_04G178100","SoyZH13_08G053200","SoyZH13_07G065100","SoyZH13_09G233700","SoyZH13_19G191400","SoyZH13_06G170200","SoyZH13_04G165800","SoyZH13_03G186600")
##nucl <- hq2602[which(hq2602$geneID %in% nucl_gene),]
nucl <- nucl[which(nucl$Type != "N"),]##
cyto <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq_cyto.txt",header = T,sep = '\t')
##chlo_gene <- c("SoyZH13_06G001400","SoyZH13_04G001700","SoyZH13_14G184900","SoyZH13_02G213800","SoyZH13_05G054400","SoyZH13_11G055100","SoyZH13_17G166900","SoyZH13_20G030500")
##chlo <- hq2602[which(hq2602$geneID %in% chlo_gene),]
cyto <- cyto[which(cyto$Type != "N"),]##

y1 <- floor(min(hq2602$y)/100)*100
x1 <- floor(min(hq2602$x)/100)*100
i <- 1
out <- file("D:/Graduate/空间转录组课题/数据表格/HQ_02/heatmap.txt","w")
while(TRUE){
  if(y1 > max(hq2602$y)){break}
  while(TRUE){
    if(x1 > max(hq2602$x)){break}
    wall1 <- wall[which(wall$x >= x1 & wall$x <= x1+100 & wall$y >= y1 & wall$y <= y1+100),c(1,4)]
    U1 <- U[which(U$x >= x1 & U$x <= x1+100 & U$y >= y1 & U$y <= y1+100),c(1,4)]
    S1 <- S[which(S$x >= x1 & S$x <= x1+100 & S$y >= y1 & S$y <= y1+100),c(1,4)]
    nucl1 <- nucl[which(nucl$x >= x1 & nucl$x <= x1+100 & nucl$y >= y1 & nucl$y <= y1+100),c(1,4)]
    chlo1 <- chlo[which(chlo$x >= x1 & chlo$x <= x1+100 & chlo$y >= y1 & chlo$y <= y1+100),c(1,4)]
    ##if(length(nucl1$geneID) == 0 | length(chlo1$geneID) == 0 | length(wall1$geneID) == 0 | length(U1$geneID) == 0 | length(S1$geneID) == 0){i <- i + 1;x1 <- x1 + 100;next;}
    chlo1 <- chlo1 %>% group_by(geneID) %>% 
      summarise(Count = sum(MIDCount))
    nucl1 <- nucl1 %>% group_by(geneID) %>% 
      summarise(Count = sum(MIDCount))
    S1 <- S1 %>% group_by(geneID) %>% 
      summarise(Count = sum(MIDCount))
    U1 <- U1 %>% group_by(geneID) %>% 
      summarise(Count = sum(MIDCount))
    wall1 <- wall1 %>% group_by(geneID) %>% 
      summarise(Count = sum(MIDCount))
    all1 <- U1 %>% full_join(S1, by = "geneID") %>% 
      full_join(nucl1, by = "geneID") %>% 
      full_join(chlo1, by = "geneID") %>%
      full_join(wall1, by = "geneID") ##%>%
      ##na.omit()
    for(i in 2:6){all1[is.na(all1[,i]),i] <- 0}
    ##if(length(all1$geneID) <= 2){i <- i + 1;x1 <- x1 + 100;next;}
    colnames(all1)[2] <- "U"
    colnames(all1)[3] <- "S"
    colnames(all1)[4] <- "nucl"
    colnames(all1)[5] <- "chlo"
    colnames(all1)[6] <- "wall"
    allhm1 <- cor(all1[,2:6])
    hM <- format(round(allhm1, 3))#对数据保留2位小数
    write(sprintf("X%dY%d",x1/100,y1/100),out,append = T)
    write(hM,out,append = T,ncolumns =5)
    write("\n",out,append = T)
    ##png(filename = paste0("D:/Graduate/空间转录组课题/图片观察/HQ_02/heatmap100/X",x1/100,"Y",y1/100,".png"))
    ##heatmap.2(allhm1,
    ##          trace="none",#不显示trace
    ##          col = coul,
    ##          density.info = "none",#图例取消density
    ##          key.xlab ='Correlation',
    ##          key.title = "",
    ##          cexRow = 1,cexCol = 1,#修改横纵坐标字体
    ##          Rowv = F,Colv = F, #去除聚类
    ##          margins = c(6, 6),
    ##          breaks = seq(-1,1,0.08), 
    ##          cellnote = hM,notecol='black'#添加相关系数的值及修改字体颜色
    ##)
    ##dev.off()
    x1 <- x1 + 100
    i <- i + 1
  }
  x1 <- floor(min(hq2602$x)/100)*100
  y1 <- y1 + 100
}
close(out)






wall1 <- wall[which(wall$x >= 10400 & wall$x <= 10400+200 & wall$y >= 15700 & wall$y <= 15700+200),c(1,4)]
U1 <- U[which(U$x >= 10400 & U$x <= 10400+200 & U$y >= 15700 & U$y <= 15700+200),c(1,4)]
S1 <- S[which(S$x >= 10400 & S$x <= 10400+200 & S$y >= 15700 & S$y <= 15700+200),c(1,4)]
nucl1 <- nucl[which(nucl$x >= 10400 & nucl$x <= 10400+200 & nucl$y >= 15700 & nucl$y <= 15700+200),c(1,4)]
cyto1 <- cyto[which(cyto$x >= 10400 & cyto$x <= 10400+200 & cyto$y >= 15700 & cyto$y <= 15700+200),c(1,4)]
cyto1 <- cyto1 %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount))
nucl1 <- nucl1 %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount))
S1 <- S1 %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount))
U1 <- U1 %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount))
wall1 <- wall1 %>% group_by(geneID) %>% 
  summarise(Count = sum(MIDCount))
all1 <- U1 %>% full_join(S1, by = "geneID") %>% 
  full_join(nucl1, by = "geneID") %>% 
  full_join(cyto1, by = "geneID") %>%
  full_join(wall1, by = "geneID") %>%
  na.omit()





hq2602 <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/hq2602_U_cluster.txt",header = F)
hq2602 <- hq2602[rep(1:nrow(hq2602),hq2602$MIDCount),]
write.table(test[,2:3],"D:/Graduate/空间转录组课题/数据表格/HQ_02/coord.txt",sep = "\t",row.names = F,col.names = T,quote = F)







test <- type[which(type$x > 10500 & type$x < 11500 & type$y > 15500 & type$y < 16500 & type$Type == "U"),]
test <- test[rep(1:nrow(test),test$MIDCount),]
plot(test$x,test$y,col = rgb(0,0,0, test$MIDCount/max(test$MIDCount)*255, maxColorValue = 255),xlim = c(min(hq2602$V1),max(hq2602$V1)), ylim = c(min(hq2602$V2),max(hq2602$V2)),xlab = "",ylab = "",pch = ".")







type <- read.table("D:/Graduate/空间转录组课题/数据表格/HQ_02/26_02hq.gem",header = T)
U <- type[which(type$Type == "U" & type$label != -1),]## & type$label %in% ch
NS <- type[which(type$Type != "U" & type$label != -1),]## & type$label %in% ch
cluster1 <- hq2602[which(hq2602$V3 == 1),]

 
plot(hq2602$V1,hq2602$V2,col = hq2602$V3,xlim = c(min(hq2602$V1),max(hq2602$V1)), ylim = c(min(hq2602$V2),max(hq2602$V2)),xlab = "",ylab = "",pch = ".")
par(new = TRUE)


plot(NS$x,NS$y,col = rgb(142,139,254, 10, maxColorValue = 255),xlim = c(min(hq2602$x),max(hq2602$x)), ylim = c(min(hq2602$y),max(hq2602$y)),xlab = "",ylab = "",pch = ".")
par(new = TRUE)
plot(U$x,U$y,col = rgb(0,0,0, U$MIDCount/max(U$MIDCount)*255, maxColorValue = 255),xlim = c(min(hq2602$x),max(hq2602$x)), ylim = c(min(hq2602$y),max(hq2602$y)),xlab = "",ylab = "",pch = ".")







wall500 <- wall[which(wall$x>=10500 & wall$x<=11000 & wall$y>=16500 & wall$y<=17000),]
plot(wall500$x,wall500$y,col = rgb(0,0,0,100, maxColorValue = 255),xlim = c(min(wall500$x),max(wall500$x)), ylim = c(min(wall500$y),max(wall500$y)),xlab = "",ylab = "",pch = ".")





plot(b$x,b$y,col = rgb(142,139,254, 10, maxColorValue = 255),xlim = c(min(type$x),max(type$x)), ylim = c(min(type$y),max(type$y)),xlab = "",ylab = "",pch = ".")
par(new = TRUE)
plot(a$x,a$y,col = rgb(0,0,0, a$MIDCount/max(a$MIDCount)*255, maxColorValue = 255),xlim = c(min(type$x),max(type$x)), ylim = c(min(type$y),max(type$y)),xlab = "",ylab = "",pch = ".")