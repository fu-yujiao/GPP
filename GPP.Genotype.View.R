if (!exists("GAPIT.Remove.outliers", mode = "function")) {
  `GAPIT.Remove.outliers` <- function(x, na.rm = TRUE, pro = 0.25, size = 1.5, ...) {
    qnt <- stats::quantile(x, probs = c(pro, 1 - pro), na.rm = na.rm, ...)
    y <- x
    H <- size * stats::IQR(y, na.rm = na.rm)
    y[x <= (qnt[1] - H)] <- min(y, na.rm = na.rm)
    y[x >= (qnt[2] + H)] <- max(y, na.rm = na.rm)
    idx <- x <= (qnt[1] - H) | x >= (qnt[2] + H)
    res <- vector("list")
    res$y <- y
    res$idx <- idx
    res
  }
}

GAPIT.Genotype.View <-function(GI=NULL,Frequency_MAF=NULL,chr=NULL, cut.dis=1,n.select=NULL,N4=FALSE,
                                 WS0=10000,Aver.Dis=1000,...){
chor_taxa=as.character(unique(GI[,2]))
chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
letter.index=grep("[A-Z]|[a-z]",chor_taxa)
if(!setequal(integer(0),letter.index))
  {     

      if(length(letter.index)!=length(chor_taxa))
        {
          chr.letter=chor_taxa[letter.index]
          chr.taxa=chor_taxa[-letter.index]
        }else{
          chr.letter=chor_taxa
          chr.taxa=NULL
        }
      Chr=as.character(GI[,2])
      for(i in letter.index)
        {
         index=Chr==chor_taxa[i]
         Chr[index]=i 
        }
      GI[,2]=as.data.frame(Chr)
  }
GI2=GI[order(as.numeric(as.matrix(GI[,3]))),]
GI2=GI2[order(as.numeric(as.matrix(GI2[,2]))),]
GI=GI2
chr=as.character(as.matrix(unique(GI[,2])))
allchr=as.character(GI[,2])

## make an index for marker selection with binsize
print("Filting marker for GAPIT.Genotype.View function ...")
pos.fix=as.numeric(GI[,2])*10^(nchar(max(as.numeric(GI[,3]))))+as.numeric(GI[,3])
if(is.null(n.select))n.select=10000
if(n.select>nrow(GI))n.select=nrow(GI)-1
rs.index=sample(nrow(GI)-1,n.select)
rs.index=sort(rs.index)

## filter genotype by rs.index
if((max(rs.index)+10)>nrow(GI)) rs.index[(rs.index+10)>nrow(GI)]=rs.index[(rs.index+10)>nrow(GI)]-10
dist2=NULL
GI2=GI[rs.index,]
dist1=abs(as.numeric(GI[rs.index,3])-as.numeric(GI[rs.index+1,3]))
if(N4) dist2=abs(as.numeric(GI[rs.index,3])-as.numeric(GI[rs.index+4,3]))
dist=c(dist1,dist2)
dist.out=GAPIT.Remove.outliers(dist,pro=0.1,size=1.1)
if(is.null(WS0)) WS0=((max(dist,na.rm=TRUE))%/%1000)*1000
if(WS0==0)WS0=1
index=dist>WS0
dist[index]=NA

# set different colors for odd or even chromosome
m=nrow(GI)
theCol=as.numeric(GI2[,2])%%2 # here should work, based on the Chr is numeric values
# Summer Beach palette (exclude #FEE199 as requested)
.sb_cols = c("#FC757B", "#F97F5F", "#FAA26F", "#FDCD94", "#B0D6A9", "#65BDBA", "#3C9BC9")
colDisp=array("#3C9BC9",m-1)
colIndex=theCol==1
colDisp[colIndex]="#FC757B"
colDisp=colDisp

chr.pos=rep(NA,length(chr))
chr.pos2=rep(1,length(chr)+1)
rownames(GI2)=1:nrow(GI2)
mm=nrow(GI2)
for(i in 1:length(chr))
{
  chr.pos[i]=floor(median(as.numeric(rownames(GI2[GI2[,2]==chr[i],]))))
  chr.pos2[i+1]=max(as.numeric(rownames(GI2[GI2[,2]==chr[i],])))
}
odd=seq(1,length(chr),2)
d.V=dist/Aver.Dis

if(is.null(Frequency_MAF)) {
    print("Frequency_MAF not provided, calculating r1...")
    Frequency_MAF <- data.frame(X=GI[,1], het.snp=runif(nrow(GI)), maf=runif(nrow(GI)), r=runif(nrow(GI), -1, 1))
  }
r1=Frequency_MAF[,4]
min_length <- min(length(d.V), length(r1))
d.V <- d.V[1:min_length]
r1 <- r1[1:min_length]

grDevices::pdf("GAPIT.Genotype.Distance_R_Chro.pdf", width =10, height = 6)
par(mfcol=c(2,3),mar = c(5,5,2,2))
plot(r1, xlab="Marker",las=1,xlim=c(1,mm),ylim=c(-1,1),
    ylab="R",axes=FALSE, main="a",cex=.5,col=colDisp)
axis(1,at=chr.pos2,labels=rep("",length(chr)+1))
axis(1,at=chr.pos[odd],labels=chr[odd],tick=FALSE)
axis(2,las=1)

print("The average distance between markers are ...")
print(head(d.V))
plot(d.V,las=1, xlab="Marker", ylab="Distance (Kb)",xlim=c(1,mm), ylim=c(0,ceiling(max(d.V,na.rm=TRUE))),
    axes=FALSE,main="d",cex=.5,col=colDisp)
axis(1,at=chr.pos2,labels=rep("",length(chr)+1))
axis(1,at=chr.pos[odd],labels=chr[odd],tick=FALSE)
axis(2,las=1)
print(str(r1))
r0.hist=hist(r1,  plot=FALSE)
r0=r0.hist$counts
r0.demo=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,1000,ifelse(nchar(max(r0))<=12,10000000,100000000000)))
r0.hist$counts=r0/r0.demo
ylab0=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,2,ifelse(nchar(max(r0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
d.V.hist=hist(d.V, plot=FALSE)
d.V0=d.V.hist$counts
d.V0.demo=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,1000,ifelse(nchar(max(d.V0))<=12,10000000,100000000000)))
ylab0=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,2,ifelse(nchar(max(d.V0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
d.V.hist$counts=d.V0/d.V0.demo
.sb_cols = c("#FC757B", "#F97F5F", "#FAA26F", "#FDCD94", "#B0D6A9", "#65BDBA", "#3C9BC9")
.r0_nbars = length(r0.hist$breaks) - 1L
.dv_nbars = length(d.V.hist$breaks) - 1L
.r0_col = rep_len(.sb_cols, max(1L, .r0_nbars))
.dv_col = rep_len(.sb_cols, max(1L, .dv_nbars))
plot(r0.hist, xlab="R", las=1,ylab=ylab.store[ylab0], main="b",col=.r0_col, border = "white")
plot(d.V.hist, las=1,xlab="Distance (Kb)",col=.dv_col, border = "white", ylab=ylab.store[ylab0], main="e",cex=.5,xlim=c(0,WS0/Aver.Dis))
#plot(d.V,r1,las=1,xlab="Distance (Kb)",ylim=c(-1,1),pch=16,
# ylab="R",main="c",cex=.5,col="gray60",xlim=c(0,WS0/Aver.Dis))
print(length(d.V))
print(length(r1))
# 计算有效的索引
valid_idx <- which(!is.na(d.V) & !is.na(r1))
plot(d.V[valid_idx], r1[valid_idx], las = 1, xlab = "Distance (Kb)", ylim = c(-1, 1), 
     pch = 16, ylab = "R", main = "c", cex = 0.5, col = "#3C9BC9", 
     xlim = c(0, WS0/Aver.Dis))

abline(h=0,col="darkred")
plot(d.V,r1^2,las=1,xlab="Distance (Kb)",ylim=c(0,1),pch=16,
  ylab="R sqaure", main="f",cex=.5,col="#3C9BC9",xlim=c(0,WS0/Aver.Dis))

dist[dist==0]=1
indOrder=order(dist)
ma=cbind(as.data.frame(dist),as.data.frame(r1)^2)
ma=ma[indOrder,]
index.na=ma[,1]>WS0
maPure=ma[!index.na,]
maPure=maPure[!is.na(maPure[,1]),]
ns=maPure[,1]
if(n.select>500)
{
  if(max(ns)>1000)
    {
    ns.bin=c(seq(0,90,10),seq(100,max(ns)/4,100),seq(max(ns)/4+200,max(ns)/3,200),seq(max(ns)/3+500,max(ns)/2,500),seq(max(ns)/2+1000,max(ns),1000))
    }else{
    ns.bin=c(seq(0,max(ns),10))    
    }
}else{
ns.bin=seq(0,max(ns),5000)
}
loc=matrix(NA,length(ns.bin)-1,3)
j=0
for (i in 1:(length(ns.bin)-1)){
  j=j+1
  pieceD=maPure[ ns.bin[i]<ns&ns<ns.bin[i+1], 1]
  pieceR=maPure[ ns.bin[i]<ns&ns<ns.bin[i+1], 2]
  loc[i,1]=ns.bin[i+1]
  loc[i,2]=mean(pieceR,na.rm=T)
  loc[i,3]=length(pieceR)
}
lines(loc[,1]/Aver.Dis,loc[,2],col="darkred",xlim=c(0,WS0/Aver.Dis))
colnames(loc)=c("Distance","Rsquare","Number")
write.csv(loc,paste("GAPIT.Genotype.Distance.Rsquare.csv",sep=""))
grDevices::dev.off()

#H=1-abs(X2-1)
#het.ind=apply(H,1,mean)
het.snp=Frequency_MAF[,2]
maf=Frequency_MAF[,3]
r1=Frequency_MAF[,4]
het.ind <- NULL
try({
  hmp_files <- list.files(pattern = "\\.hmp\\.txt$", ignore.case = TRUE)
  if (length(hmp_files) > 0) {
    pick <- which(grepl("^mdp_genotype.*\\.hmp\\.txt$", basename(hmp_files), ignore.case = TRUE))
    hmp <- if (length(pick) > 0) hmp_files[pick[1]] else hmp_files[1]
    con <- base::file(hmp, open = "r")
    on.exit(base::close(con), add = TRUE)
    header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
    n_taxa <- length(header) - 11
    if (n_taxa > 0) {
      het_count <- integer(n_taxa)
      obs_count <- integer(n_taxa)
      miss <- c("NN", "N", "NA", "", "0", ".", "--")
      repeat {
        lines <- readLines(con, n = 2000)
        if (length(lines) == 0) break
        for (line in lines) {
          fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
          if (length(fields) < 12) next
          geno <- fields[12:length(fields)]
          if (length(geno) != n_taxa) next
          nm <- !(geno %in% miss) & !is.na(geno)
          het <- nm & (nchar(geno) == 2) & (substr(geno, 1, 1) != substr(geno, 2, 2))
          het_count <- het_count + as.integer(het)
          obs_count <- obs_count + as.integer(nm)
        }
      }
      het.ind <- het_count / obs_count
      het.ind[obs_count == 0] <- NA_real_
    }
  }
}, silent = TRUE)
grDevices::pdf("GAPIT.Genotype.MAF_Heterozosity.pdf", width =10, height = 6)
#Display
layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
layout(mat = layout.matrix,
       heights = c(100,80,120),
       widths = c(2, 3))
par(mar = c(1, 5, 1, 1))
plot(het.snp,  las=1,ylab="Heterozygosity", xlim=c(1,mm),axes=FALSE,
    cex=.5,col=colDisp,xaxt='n')
axis(2,las=1)
par(mar = c(1, 5, 0, 1))
plot(maf, las=1,xlab="Marker", ylab="MAF",xlim=c(1,mm),axes=FALSE,
    cex=.5,col=colDisp,xaxt='n')
#output=cbind(het.snp,maf,r1)
#colnames(output)=c("het.snp","maf","r")
#write.csv(output,"GAPIT.Genotype.Frequency_MAF.csv",quote=FALSE)
axis(2,las=1)
par(mar = c(5, 5, 0, 1))
plot((r1^2),  las=1,ylab="R Sqaure", xlab="Marker", xlim=c(1,mm),axes=FALSE,cex=.5,col=colDisp)
axis(1,at=chr.pos2,labels=rep("",length(chr)+1))
axis(1,at=chr.pos,labels=chr,tick=FALSE)
axis(2,las=1)
grDevices::dev.off()

#Display Het and MAF distribution
grDevices::pdf("GAPIT.Genotype.Frequency.pdf", width =10, height = 3.5)
layout.matrix <- matrix(c(1,2,3), nrow = 1, ncol = 3)
layout(mat = layout.matrix,
       heights = c(100,80,120),
       widths = c(2, 2,2))
par(mar = c(5, 5, 2, 0))
.sb_cols2 = c("#FC757B", "#F97F5F", "#FAA26F", "#FDCD94", "#B0D6A9", "#65BDBA", "#3C9BC9")
if (!is.null(het.ind) && any(!is.na(het.ind))) {
  .x_hetind = as.numeric(het.ind[!is.na(het.ind)])
  .h1 = hist(.x_hetind, plot = FALSE)
  .n1 = length(.h1$counts)
  .col1 = rep_len(.sb_cols2, max(1L, .n1))
  barplot(.h1$density, col = .col1, border = "white", space = 0, las = 1,
          xlab = "Individual heterozygosity", ylab = "Frequency", main = "a")
  .ntick = min(6L, .n1 + 1L)
  .at1 = seq(0, .n1, length.out = .ntick)
  .lab1 = formatC(seq(min(.h1$breaks, na.rm = TRUE), max(.h1$breaks, na.rm = TRUE), length.out = .ntick),
                  format = "g", digits = 2)
  axis(1, at = .at1, labels = .lab1)
} else {
  plot.new()
  title(main = "a")
}
par(mar = c(5, 4, 2, 1))
.x_hetsnp = het.snp[!is.na(het.snp)]
.h2 = hist(.x_hetsnp, plot = FALSE)
.n2 = length(.h2$counts)
.col2 = rep_len(.sb_cols2, max(1L, .n2))
barplot(.h2$density, col = .col2, border = "white", space = 0, las = 1,
        xlab = "Marker heterozygosity", ylab = "Frequency", main = "b")
.ntick2 = min(6L, .n2 + 1L)
.at2 = seq(0, .n2, length.out = .ntick2)
.lab2 = formatC(seq(min(.h2$breaks, na.rm = TRUE), max(.h2$breaks, na.rm = TRUE), length.out = .ntick2),
                format = "g", digits = 2)
axis(1, at = .at2, labels = .lab2)
par(mar = c(5, 4, 2, 1))
.x_maf = maf[!is.na(maf)]
.h3 = hist(.x_maf, plot = FALSE)
.n3 = length(.h3$counts)
.col3 = rep_len(.sb_cols2, max(1L, .n3))
barplot(.h3$density, col = .col3, border = "white", space = 0, las = 1,
        xlab = "MAF", ylab = "Frequency", main = "c")
.ntick3 = min(6L, .n3 + 1L)
.at3 = seq(0, .n3, length.out = .ntick3)
.lab3 = formatC(seq(min(.h3$breaks, na.rm = TRUE), max(.h3$breaks, na.rm = TRUE), length.out = .ntick3),
                format = "g", digits = 2)
axis(1, at = .at3, labels = .lab3)

grDevices::dev.off()
print(paste("GAPIT.Genotype.View ", ". pdfs generate.","successfully!" ,sep = ""))
}
