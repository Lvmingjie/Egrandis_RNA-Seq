Include_reads = read.table("Include_reads.txt", header = TRUE, sep = "\t")
Skipped_reads = read.table("skipped_reads.txt", header = TRUE, sep = "\t")
Include_level = read.table("Include_level.txt", header = TRUE, sep = "\t")
dpsi = read.table("rMATS_dpsi.txt", header = TRUE, sep = "\t")


dpsi$SignificantlyDAS  = 0
sdeglist = dpsi$SignificantlyDAS

pvcf = 0.001   #FDR cutoff
dpcf = 0.2   # delta PSI  cutoff
colnum = length(dpsi[1,]) - 2
rownum = length(dpsi[,1])

for (g in 1:rownum) {
  sdeg = 0
  for (n in 1: (colnum/3)) {
    colpv = n*3 + 1
    coldp = n*3 - 1
    newpv = dpsi[g,colpv]
    newdp = abs(dpsi[g,coldp])
    if (  ! (is.na(newpv) | is.na(newdp)) ){
      if (newpv < pvcf & newdp > dpcf){
        sdeg = sdeg + 1
      }
    }
  }
  sdeglist[g] = sdeg
}

dpsi$SignificantlyDAS <- sdeglist
diffAS <- dpsi[dpsi$SignificantlyDAS >= 1, ]
length(dpsi[,1])
length(diffAS[,1])

diffASlong <- merge(diffAS, Include_reads, by = "EventID", all = FALSE)
diffASlong1 <- merge(diffASlong, Skipped_reads, by = "EventID", all = FALSE)
diffASlong2 <- merge(diffASlong1, Include_level, by = "EventID", all = FALSE)
diffASlong2 <- na.omit(diffASlong2)
nrow(diffASlong2)
outfile <- "different_AS.txt"
write.table( as.data.frame(diffASlong2), file= outfile, sep="\t",col.names=NA) 

##################cluster#########
raw <- read.table (outfile , header = TRUE, sep="\t") 
n = ncol(Include_level)-1
used <- raw[,(ncol(raw)-n+1):ncol(raw)]
samplenames = c('Mock',  # sample names
                'Bdef_6h',
                'Bdef_24h',
                'Bdef_48h',
                'Bdef_96h',
                'Bdef_21d'
)
colnames(used) <- samplenames

samplenum = n 
mydata <- used
length(mydata[,1])
fit <- kmeans(mydata, 50,iter.max = 10)
clumean <- aggregate(mydata,by=list(fit$cluster),FUN=mean) 
mydata <- data.frame(mydata, fit$cluster)
n = 1:length(clumean[,1])
for (i in 1:length(clumean[,1])){
  n[i] = length (mydata$fit.cluster[mydata$fit.cluster == i])
}
n   # number of clusters
##plot clusters

plot(-100,-100,xlim=c(0,samplenum),ylim=c(-0.5,1.5), xaxt = "n") 
axis(1,1:samplenum -0.5,samplenames) 
for (i in 1:length(clumean[,1])){
  lines(1:samplenum-0.5,clumean[i,2:(samplenum+1)],col = clumean$Group.1[i] ,lwd = 2)
  points(1:samplenum-0.5,clumean[i,2:(samplenum+1)],col = clumean$Group.1[i],pch = 19 ) 
  num = length (mydata$fit.cluster[mydata$fit.cluster == i])
  text(samplenum,clumean[i,samplenum], i,col =clumean$Group.1[i],font = 10,cex = 1.5 ) 
  text(samplenum,clumean[i,samplenum], num,col =clumean$Group.1[i] )
}
fcot <- data.frame(raw, fit$cluster)
write.table(fcot, file = "AS.clusters.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
lblue <- rgb(0, 205, 255,120, maxColorValue=255)
lred <- "coral";
pdf("AS_patterns.pdf",onefile = TRUE,width = 7.5, height = 5,pointsize = 10)
par(oma=c(8,3,2.6,3),mar=c(2,5,3,1))  #内外边界大小
for (i in 1:length(clumean[,1])){
  oneclu <- mydata[mydata$fit.cluster == i, ]
  ti = paste("Group",i,length(oneclu[,1]),sep = " ")
  plot(-1,-1,xlim=c(0,(samplenum)),ylim=c(0,1),xaxt = "n",yaxt = "n",main = ti,cex.main=2,font=2,xlab="",ylab="PSI",cex.lab=2,font.lab=2,font.main=2)  #X轴和Y轴，需要改，根据模式图来修改纵坐标
  axis(1,1:(samplenum)-0.5,samplenames,las=0,cex.axis=1,mgp=c(100,2,0),font=2,font.lab=2,las=2) 
  axis(2,las=2,cex.axis=2,at=c(-8:8),mgp=c(100,1,0),font.lab=2,font=2)
  abline(h=1,lty =2,lwd = 2, col ="darkgrey")
  abline(h=0,lty =2, lwd = 2,col ="darkgrey")
  for (j in 1:length(oneclu[,1])){
    lines(1:samplenum  -0.5,  oneclu[j,1:samplenum],col = lblue ,lwd = 1) 
  }
  lines(1:samplenum  -0.5,  clumean[i,2:(samplenum+1)  ],   col = "darkblue" ,lwd = 5)
  points( 1:samplenum  -0.5,clumean[i,2:(samplenum+1)  ],   col = "darkblue",pch = 19 ,cex = 2.6)
  points( 1:samplenum  -0.5,clumean[i,2:(samplenum+1)  ],   col = "white",pch = 19 ,cex = 1)
}

dev.off()
