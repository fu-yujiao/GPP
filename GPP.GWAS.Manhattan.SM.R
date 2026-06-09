`GAPIT.Manhattan` <-
function(model_store=NULL,GWAS_Results=NULL,GI.MP = NULL,GD=NULL,name.of.trait = "Trait",plot.type = "Genomewise",width0=18,height0=5.75,
DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,plot.style="Oceanic",CG=NULL,plot.bin=10^9,chor_taxa=NULL,byTraits=FALSE){#Manhattan-single model
environ_name=NULL
if(byTraits)
{
  for(i in 1:length(name.of.trait))
  {
     for(j in 1:length(model_store))
     {
        environ_name=c(environ_name,paste(model_store[j],".",name.of.trait[i],sep=""))
     }
  }
}else{
  for(i in 1:length(model_store))
  {
     for(j in 1:length(name.of.trait))
     {
        environ_name=c(environ_name,paste(model_store[i],".",name.of.trait[j],sep=""))
     }
  }
}
if(is.null(GWAS_Results))
  {
  for(i in 1:length(environ_name))
  {
  print(paste("Trying to read:", paste("GAPIT.Association.GWAS_Results.", environ_name[i], ".csv", sep = "")))
  file.exists(paste("GAPIT.Association.GWAS_Results", environ_name[i], ".csv", sep = ""))
  GWAS_Results=read.csv(paste("GAPIT.Association.GWAS_Results",environ_name[i],".csv",sep=""),head=T)

  GWAS_Results <- na.omit(GWAS_Results)
#str(GWAS_Results)
#head(GWAS_Results)
chm.to.analyze <- unique(GWAS_Results[,2])
#chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
# 获取所有唯一的染色体标签
chor_taxa = as.character(unique(GWAS_Results[,2]))

# 将数字染色体与字母染色体分开
num_chor_taxa = chor_taxa[!grepl("[A-Z]|[a-z]", chor_taxa)]
letter_chor_taxa = chor_taxa[grepl("[A-Z]|[a-z]", chor_taxa)]

# 对数字染色体进行排序，字母染色体保持原样
num_chor_taxa = num_chor_taxa[order(as.numeric(num_chor_taxa))]
# 合并数字和字母染色体
sorted_chor_taxa = c(num_chor_taxa, letter_chor_taxa)

GWAS_Results[,2] <- factor(GWAS_Results[,2], levels = sorted_chor_taxa)

numCHR= length(chm.to.analyze)
numMarker=nrow(GWAS_Results)
bonferroniCutOff=-log10(cutOff/numMarker)
sp=sort(GWAS_Results[,4])
#print("sp values:")
#print(head(sp))
spd=abs(cutOff-sp*numMarker/cutOff)
spd <- na.omit(spd)
index_fdr=grep(min(spd),spd)[1]
FDRcutoff=-log10(cutOff*index_fdr/numMarker)
if(is.null(GWAS_Results)) return
  #Handler of lable position only indicated by negatie position
    position.only=F
    if(!is.null(seqQTN)){
      if(seqQTN[1]<0){
        seqQTN=-seqQTN
        position.only=T
      }      
    }  
borrowSlot=9	
GWAS_Results[,borrowSlot]=0 
#Inicial as 0   
if(!is.null(seqQTN))GWAS_Results[seqQTN,borrowSlot]=1 	
if(plot.type == "Genomewise")
    {
        nchr=length(chm.to.analyze)

    #Set color schem            
        ncycle=ceiling(nchr/band)
        ncolor=band*ncycle
        #palette(rainbow(ncolor+1))
        cycle1=seq(1,nchr,by= ncycle)
        thecolor=cycle1
        for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
        col.Rainbow=rainbow(ncolor+1)[thecolor]         
          col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
          col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
            col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
            col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))        
            col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))       
            col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
            col.Oceanic=rep(c(  '#EC5f67',      '#FAC863',  '#99C794',      '#6699CC',  '#C594C5'),ceiling(numCHR/5))
            col.cougars=rep(c(  '#990000',      'dimgray'),ceiling(numCHR/2))
        
        if(plot.style=="Rainbow")plot.color= col.Rainbow
        if(plot.style =="FarmCPU")plot.color= col.Rainbow
        if(plot.style =="Rushville")plot.color= col.Rushville
        if(plot.style =="Congress")plot.color= col.Congress
        if(plot.style =="Ocean")plot.color= col.Ocean
        if(plot.style =="PLINK")plot.color= col.PLINK
            if(plot.style =="Beach")plot.color= col.Beach
            if(plot.style =="Oceanic")plot.color= col.Oceanic
            if(plot.style =="cougars")plot.color= col.cougars
        
        #FarmCPU uses filled dots
        mypch=1
        if(plot.style =="FarmCPU")mypch=20
                
        GWAS_Results <- GWAS_Results[order(GWAS_Results[,3]),]
        GWAS_Results <- GWAS_Results[order(GWAS_Results[,2]),]

        ticks=NULL
        lastbase=0
        
        
        
        #change base position to accumulatives (ticks)
        for (i in sorted_chor_taxa)
        {
            index=(GWAS_Results[,2]==i)
            ticks <- c(ticks, lastbase+mean(GWAS_Results[index,3], na.rm = TRUE))
            GWAS_Results[index,3]=GWAS_Results[index,3]+lastbase
            lastbase=max(GWAS_Results[index,3], na.rm = TRUE)
        }
        
        
        
        x0 <- as.numeric(GWAS_Results[,3])
        y0 <- -log10(as.numeric(GWAS_Results[,4]))
		#print(summary(x0))
		#print(summary(y0))
        z0 <- as.factor(GWAS_Results[,2])
        position=order(y0,decreasing = TRUE)
        index0=GAPIT.Pruning(y0[position],DPP=DPP)
        index=position[index0]
        
        x=x0[index]
        y=y0[index]
        z=z0[index]

        
        QTN=GWAS_Results[which(GWAS_Results[,borrowSlot]==1),]
       
        #Draw circles with same size and different thikness
        size=1 #1
        ratio=10 #5
        base=1 #1
        themax=ceiling(max(y))
        themin=floor(min(y))
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        s=size-wd/ratio/2
        
        
        if(plot.style =="FarmCPU"){
        pdf(paste("FarmCPU.", name.of.trait,".Manhattan.Plot.Genomewise.pdf" ,sep = ""), width = width0,height=height0)
        }else{
        pdf(paste("GAPIT.Association.Manhattan_Geno.", model_store, ".",name.of.trait,".pdf" ,sep = ""), width = width0,height=height0)
        }
            par(mar = c(3,6,5,1))
            plot(y~x,xlab="",ylab=expression(-log[10](italic(p))) ,las=1,
            cex.axis=1, cex.lab=1.3 ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=wd,cex=s+.3,main = paste(model_store,name.of.trait,sep = "."),cex.main=2.5)
        
        #Label QTN positions
        if(is.vector(QTN)){
		print("QTN is a vector")
          if(position.only){abline(v=QTN[3], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[3], QTN[4], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[3], QTN[4], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        }else{
		#print("QTN is a data.frame")
          if(position.only){abline(v=QTN[,3], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,3], QTN[,4], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[,3], QTN[,4], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        }
        
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,col="forestgreen")
        #Add FDR line
        abline(h=FDRcutoff,col="forestgreen",lty=2)
        
        if(length(chor_taxa)!=length(ticks))chor_taxa=NULL
        
        if(!is.null(sorted_chor_taxa))
        {axis(1, at=ticks,cex.axis=1,labels=sorted_chor_taxa,tick=T,gap.axis=0.25)
        }else{axis(1, at=ticks,cex.axis=1,labels=chm.to.analyze,tick=F)}
        axis(2, at=1:themax,cex.axis=1,las=1,labels=1:themax,gap.axis=3,tick=F)

        box()
        palette("default")
        dev.off()
        
    } }
  }else{
   GWAS_Results <- na.omit(GWAS_Results)
#str(GWAS_Results)
#head(GWAS_Results)
chm.to.analyze <- unique(GWAS_Results[,2])
#chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
# 获取所有唯一的染色体标签
chor_taxa = as.character(unique(GWAS_Results[,2]))

# 将数字染色体与字母染色体分开
num_chor_taxa = chor_taxa[!grepl("[A-Z]|[a-z]", chor_taxa)]
letter_chor_taxa = chor_taxa[grepl("[A-Z]|[a-z]", chor_taxa)]

# 对数字染色体进行排序，字母染色体保持原样
num_chor_taxa = num_chor_taxa[order(as.numeric(num_chor_taxa))]
# 合并数字和字母染色体
sorted_chor_taxa = c(num_chor_taxa, letter_chor_taxa)

GWAS_Results[,2] <- factor(GWAS_Results[,2], levels = sorted_chor_taxa)

numCHR= length(chm.to.analyze)
numMarker=nrow(GWAS_Results)
bonferroniCutOff=-log10(cutOff/numMarker)
sp=sort(GWAS_Results[,4])
#print("sp values:")
#print(head(sp))
spd=abs(cutOff-sp*numMarker/cutOff)
spd <- na.omit(spd)
index_fdr=grep(min(spd),spd)[1]
FDRcutoff=-log10(cutOff*index_fdr/numMarker)
if(is.null(GWAS_Results)) return
  #Handler of lable position only indicated by negatie position
    position.only=F
    if(!is.null(seqQTN)){
      if(seqQTN[1]<0){
        seqQTN=-seqQTN
        position.only=T
      }      
    }  
borrowSlot=9	
GWAS_Results[,borrowSlot]=0 
#Inicial as 0   
if(!is.null(seqQTN))GWAS_Results[seqQTN,borrowSlot]=1 	
if(plot.type == "Genomewise")
    {
        nchr=length(chm.to.analyze)

    #Set color schem            
        ncycle=ceiling(nchr/band)
        ncolor=band*ncycle
        #palette(rainbow(ncolor+1))
        cycle1=seq(1,nchr,by= ncycle)
        thecolor=cycle1
        for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
        col.Rainbow=rainbow(ncolor+1)[thecolor]         
          col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
          col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
            col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
            col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))        
            col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))       
            col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
            col.Oceanic=rep(c(  '#EC5f67',      '#FAC863',  '#99C794',      '#6699CC',  '#C594C5'),ceiling(numCHR/5))
            col.cougars=rep(c(  '#990000',      'dimgray'),ceiling(numCHR/2))
        
        if(plot.style=="Rainbow")plot.color= col.Rainbow
        if(plot.style =="FarmCPU")plot.color= col.Rainbow
        if(plot.style =="Rushville")plot.color= col.Rushville
        if(plot.style =="Congress")plot.color= col.Congress
        if(plot.style =="Ocean")plot.color= col.Ocean
        if(plot.style =="PLINK")plot.color= col.PLINK
            if(plot.style =="Beach")plot.color= col.Beach
            if(plot.style =="Oceanic")plot.color= col.Oceanic
            if(plot.style =="cougars")plot.color= col.cougars
        
        #FarmCPU uses filled dots
        mypch=1
        if(plot.style =="FarmCPU")mypch=20
                
        GWAS_Results <- GWAS_Results[order(GWAS_Results[,3]),]
        GWAS_Results <- GWAS_Results[order(GWAS_Results[,2]),]

        ticks=NULL
        lastbase=0
        
        
        
        #change base position to accumulatives (ticks)
        for (i in sorted_chor_taxa)
        {
            index=(GWAS_Results[,2]==i)
            ticks <- c(ticks, lastbase+mean(GWAS_Results[index,3], na.rm = TRUE))
            GWAS_Results[index,3]=GWAS_Results[index,3]+lastbase
            lastbase=max(GWAS_Results[index,3], na.rm = TRUE)
        }
        
        
        
        x0 <- as.numeric(GWAS_Results[,3])
        y0 <- -log10(as.numeric(GWAS_Results[,4]))
		#print(summary(x0))
		#print(summary(y0))
        z0 <- as.factor(GWAS_Results[,2])
        position=order(y0,decreasing = TRUE)
        index0=GAPIT.Pruning(y0[position],DPP=DPP)
        index=position[index0]
        
        x=x0[index]
        y=y0[index]
        z=z0[index]

        
        QTN=GWAS_Results[which(GWAS_Results[,borrowSlot]==1),]
       
        #Draw circles with same size and different thikness
        size=1 #1
        ratio=10 #5
        base=1 #1
        themax=ceiling(max(y))
        themin=floor(min(y))
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        s=size-wd/ratio/2
        
        
        if(plot.style =="FarmCPU"){
        pdf(paste("FarmCPU.", name.of.trait,".Manhattan.Plot.Genomewise.pdf" ,sep = ""), width = width0,height=height0)
        }else{
        pdf(paste("GAPIT.Association.Manhattan_Geno.", name.of.trait,".pdf" ,sep = ""), width = width0,height=height0)
        }
            par(mar = c(3,6,5,1))
            plot(y~x,xlab="",ylab=expression(-log[10](italic(p))) ,las=1,
            cex.axis=1, cex.lab=1.3 ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=wd,cex=s+.3,main = paste(name.of.trait,sep="             "),cex.main=2.5)
        
        #Label QTN positions
        if(is.vector(QTN)){
		print("QTN is a vector")
          if(position.only){abline(v=QTN[3], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[3], QTN[4], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[3], QTN[4], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        }else{
		print("QTN is a data.frame")
          if(position.only){abline(v=QTN[,3], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,3], QTN[,4], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[,3], QTN[,4], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        }
        
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,col="forestgreen")
        #Add FDR line
        abline(h=FDRcutoff,col="forestgreen",lty=2)
        
        if(length(chor_taxa)!=length(ticks))chor_taxa=NULL
        
        if(!is.null(sorted_chor_taxa))
        {axis(1, at=ticks,cex.axis=1,labels=sorted_chor_taxa,tick=T,gap.axis=0.25)
        }else{axis(1, at=ticks,cex.axis=1,labels=chm.to.analyze,tick=F)}
        axis(2, at=1:themax,cex.axis=1,las=1,labels=1:themax,gap.axis=3,tick=F)

        box()
        palette("default")
        dev.off()
        
    } }}
