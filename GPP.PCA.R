GAPIT.PCA<-function(X=NULL,PCs=NULL,eigenvalues=NULL,taxa=NULL, PC.number = NULL,radius=1,
  file.output=TRUE,PCA.total=0,PCA.col=NULL,
  PCA.3d=FALSE,PCA.legend=NULL,
  PCA.pch.2D=20,PCA.pch.3D=20,
  PCA.cex.lab=1.4,PCA.cex.symbols=1,
  PCA.lwd.2D=2,PCA.lwd.3D=3,
  PCA.bty.2D="o",
  PCA.2D.pairs=NULL,
  PCA.cex.2D=1,PCA.alpha=1,
  PCA.legend.cex=1,
  PCA.legend.outside.3D=FALSE,PCA.legend.inset.3D=0,
  PCA.mar.right.3D=12){
# Object: Conduct a principal component analysis, and output the prinicpal components into the workspace,
#         a text file of the principal components, and a pdf of the scree plot
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011  
############################################################################################## 
#Conduct the PCA 
print("Calling prcomp...")
if(!is.null(X))
  {
   PCA.X <- stats::prcomp(X)
   eigenvalues <- PCA.X$sdev^2
   evp=eigenvalues/sum(eigenvalues)
  }else{
   PCA.X=list()
   PCA.X$x=PCs[,-1]
   eigenvalues=as.numeric(as.matrix(eigenvalues))
   evp=eigenvalues/sum(eigenvalues)
   taxa=PCs[,1]
  }
if (is.null(PC.number)) {
  PC.number <- min(ncol(PCA.X$x), nrow(PCA.X$x))
}
nout=min(10,length(as.numeric(as.matrix(evp))))
xout=1:nout
if(is.null(PCA.col)) PCA.col="red"
# if(!is.null(PCA.legend)) PCA.col0=
PCA.col.plot <- PCA.col
if (!is.null(PCA.alpha) && is.numeric(PCA.alpha) && length(PCA.alpha) == 1 && PCA.alpha < 1) {
  PCA.col.plot <- grDevices::adjustcolor(PCA.col.plot, alpha.f = PCA.alpha)
}
print("Creating PCA graphs...")
#Create a Scree plot 

if(file.output & PC.number>1) 
  {
   grDevices::pdf("GAPIT.Genotype.PCA_eigenValue.pdf", width = 12, height = 12)
   graphics::par(mar=c(5,5,4,5)+.1,cex=2)
   #par(mar=c(10,9,9,10)+.1)
   plot(xout,eigenvalues[xout],type="b",col="blue",xlab="Principal components",ylab="Variance")
   graphics::par(new=TRUE)
   plot(xout,evp[xout]*100,type="n",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
   graphics::axis(4)
   graphics::mtext("Percentage (%)",side=4,line=3,cex=2)
   grDevices::dev.off()
  }
grDevices::pdf("GAPIT.Genotype.PCA_2D.pdf", width = 8, height = 8)
graphics::par(mar = c(5,5,5,5),xpd=TRUE)
maxPlot=min(as.numeric(PC.number[1]),3)

pairs_to_plot <- PCA.2D.pairs
if (!is.null(pairs_to_plot)) {
  if (is.matrix(pairs_to_plot) || is.data.frame(pairs_to_plot)) {
    pairs_to_plot <- split(pairs_to_plot, seq_len(nrow(pairs_to_plot)))
    pairs_to_plot <- lapply(pairs_to_plot, function(x) as.integer(x[1, 1:2]))
  } else if (is.numeric(pairs_to_plot) && length(pairs_to_plot) == 2) {
    pairs_to_plot <- list(as.integer(pairs_to_plot))
  }
}

if (is.null(pairs_to_plot)) {
  for(i in 1:(maxPlot-1))
  {
    for(j in (i+1):(maxPlot))
    {
      plot(PCA.X$x[,i],PCA.X$x[,j],xlab=paste("PC",i," (evp=",round(evp[i],4)*100,"%)",sep=""),ylab=paste("PC",j," (evp=",round(evp[j],4)*100,"%)",sep=""),pch=PCA.pch.2D,col=PCA.col.plot,cex=PCA.cex.2D,cex.axis=1.3,cex.lab=PCA.cex.lab, cex.axis=1.2, lwd=PCA.lwd.2D,las=1,bty=PCA.bty.2D)
      
	  if (!is.null(PCA.legend)) {
            legend_pch <- PCA.legend$pch
            if (is.null(legend_pch)) legend_pch <- PCA.pch.2D
            legend_inset <- PCA.legend$inset
            if (is.null(legend_inset)) legend_inset <- 0.02
            legend_bty <- PCA.legend$bty
            if (is.null(legend_bty)) legend_bty <- "o"
            legend_box_lwd <- PCA.legend$box.lwd
            if (is.null(legend_box_lwd)) legend_box_lwd <- 1
            legend_cex <- PCA.legend$cex
            if (is.null(legend_cex)) legend_cex <- PCA.legend.cex
            old_par <- graphics::par(xpd=FALSE)
            legend(PCA.legend$legend.pos, legend = PCA.legend$taxa, pch = legend_pch, col = PCA.legend$color,
                   ncol = PCA.legend$ncol, bty = legend_bty, inset = legend_inset, box.lwd = legend_box_lwd, cex = legend_cex)
            graphics::par(old_par)
        }
    }
  }
} else {
  for (pair in pairs_to_plot) {
    if (is.null(pair) || length(pair) < 2) next
    i <- as.integer(pair[1])
    j <- as.integer(pair[2])
    if (is.na(i) || is.na(j)) next
    if (i < 1 || j < 1) next
    if (i > ncol(PCA.X$x) || j > ncol(PCA.X$x)) next
    plot(PCA.X$x[,i],PCA.X$x[,j],xlab=paste("PC",i," (evp=",round(evp[i],4)*100,"%)",sep=""),ylab=paste("PC",j," (evp=",round(evp[j],4)*100,"%)",sep=""),pch=PCA.pch.2D,col=PCA.col.plot,cex=PCA.cex.2D,cex.axis=1.3,cex.lab=PCA.cex.lab, cex.axis=1.2, lwd=PCA.lwd.2D,las=1,bty=PCA.bty.2D)
    if (!is.null(PCA.legend)) {
      legend_pch <- PCA.legend$pch
      if (is.null(legend_pch)) legend_pch <- PCA.pch.2D
      legend_inset <- PCA.legend$inset
      if (is.null(legend_inset)) legend_inset <- 0.02
      legend_bty <- PCA.legend$bty
      if (is.null(legend_bty)) legend_bty <- "o"
      legend_box_lwd <- PCA.legend$box.lwd
      if (is.null(legend_box_lwd)) legend_box_lwd <- 1
      legend_cex <- PCA.legend$cex
      if (is.null(legend_cex)) legend_cex <- PCA.legend.cex
      old_par <- graphics::par(xpd=FALSE)
      legend(PCA.legend$legend.pos, legend = PCA.legend$taxa, pch = legend_pch, col = PCA.legend$color,
             ncol = PCA.legend$ncol, bty = legend_bty, inset = legend_inset, box.lwd = legend_box_lwd, cex = legend_cex)
      graphics::par(old_par)
    }
  }
}
grDevices::dev.off()


#output 3D plot
if(PCA.3d==TRUE)
{   
  if(1>2)
  {
#  if(!require(lattice)) install.packages("lattice")
#   library(lattice)
   pca=as.data.frame(PCA.X$x)
   
   grDevices::png(file="example%03d.png", width=500, heigh=500)
    for (i in seq(10, 80 , 1)){
        print(lattice::cloud(PC1~PC2*PC3,data=pca,screen=list(x=i,y=i-40),pch=20,color="red",
        col.axis="blue",cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3))
        }
    grDevices::dev.off()
    system("convert -delay 40 *.png GAPIT.PCA.3D.gif")
    
    # cleaning up
    file.remove(list.files(pattern=".png"))
    }

   if(!requireNamespace("rgl", quietly=TRUE)) install.packages("rgl")
   if(!requireNamespace("htmltools", quietly=TRUE)) install.packages("htmltools")
   if(!requireNamespace("manipulateWidget", quietly=TRUE)) install.packages("manipulateWidget")

    .gpp_rglwidget <- function(...) {
      if (exists("rglwidget", where = asNamespace("rgl"), inherits = FALSE)) {
        return(rgl::rglwidget(...))
      }
      if (requireNamespace("rglwidget", quietly = TRUE) && exists("rglwidget", where = asNamespace("rglwidget"), inherits = FALSE)) {
        return(rglwidget::rglwidget(...))
      }
      stop("rglwidget is not available.", call. = FALSE)
    }

    .gpp_toggle_widget <- function(widget, ids, label) {
      if (exists("toggleWidget", where = asNamespace("rgl"), inherits = FALSE)) {
        return(rgl::toggleWidget(widget, ids = ids, label = label))
      }
      if (requireNamespace("manipulateWidget", quietly = TRUE) && exists("toggleWidget", where = asNamespace("manipulateWidget"), inherits = FALSE)) {
        return(manipulateWidget::toggleWidget(widget, ids = ids, label = label))
      }
      widget
    }

    PCA1 <- PCA.X$x[,1]
    PCA2 <- PCA.X$x[,2]
    PCA3 <- PCA.X$x[,3]
    rgl::plot3d(min(PCA1), min(PCA2), min(PCA3),xlim=c(min(PCA1),max(PCA1)),
     ylim=c(min(PCA2),max(PCA2)),zlim=c(min(PCA3),max(PCA3)),
     xlab="PCA1",ylab="PCA2",zlab="PCA3",
     col = grDevices::rgb(255, 255, 255, 100, maxColorValue=255),radius=radius*0.01)
    group_labels <- NULL
    group_colors <- NULL
    if (!is.null(PCA.legend) && !is.null(PCA.legend$taxa) && !is.null(PCA.legend$color)) {
      df <- data.frame(label = as.character(PCA.legend$taxa), color = as.character(PCA.legend$color), stringsAsFactors = FALSE)
      df <- df[!is.na(df$label) & nzchar(df$label) & !is.na(df$color) & nzchar(df$color), , drop = FALSE]
      df <- df[!duplicated(df$label), , drop = FALSE]
      group_labels <- df$label
      group_colors <- df$color
    } else {
      group_colors <- unique(as.character(PCA.col))
      group_colors <- group_colors[!is.na(group_colors) & nzchar(group_colors)]
      group_labels <- paste0("Group ", seq_along(group_colors))
    }

    any_group <- FALSE
    group_ids <- list()
    for (k in seq_along(group_colors)) {
      idx <- as.character(PCA.col) == group_colors[k]
      if (!any(idx, na.rm = TRUE)) next
      any_group <- TRUE
      sids <- rgl::spheres3d(PCA1[idx], PCA2[idx], PCA3[idx], col = PCA.col.plot[idx], radius = radius)
      group_ids[[length(group_ids) + 1]] <- list(ids = sids, label = group_labels[k])
    }
    if (!any_group) {
      sids <- rgl::spheres3d(PCA1, PCA2, PCA3, col = PCA.col.plot, radius = radius)
      group_ids[[1]] <- list(ids = sids, label = "PCA")
    }
    widgets <- .gpp_rglwidget(width = 900, height = 900)
    for (g in group_ids) {
      widgets <- .gpp_toggle_widget(widgets, ids = g$ids, label = g$label)
    }
    if (interactive()) widgets
    htmltools::save_html(widgets, "Interactive.PCA.html")
}
#    if(!require(scatterplot3d)) install.packages("scatterplot3d")
#    library(scatterplot3d)

    if(!require(scatterplot3d)) install.packages("scatterplot3d")

    grDevices::pdf("GAPIT.Genotype.PCA_3D.pdf", width = 7, height = 7)
    legend_outside3d <- FALSE
    if (!is.null(PCA.legend)) {
      legend_outside3d <- PCA.legend$outside.3D
      if (is.null(legend_outside3d)) legend_outside3d <- PCA.legend.outside.3D
      legend_outside3d <- isTRUE(legend_outside3d)
    }
    if (legend_outside3d) {
      graphics::par(mar = c(5,5,5,PCA.mar.right.3D),xpd=TRUE)
    } else {
      graphics::par(mar = c(5,5,5,5),xpd=TRUE)
    }
    scatterplot3d::scatterplot3d(PCA.X$x[,1],
                  PCA.X$x[,2],
                  PCA.X$x[,3],
                  xlab = paste("PC",1," (evp=",round(evp[1],4)*100,"%)",sep=""),
                  ylab = paste("PC",2," (evp=",round(evp[2],4)*100,"%)",sep=""),
                  zlab = paste("PC",3," (evp=",round(evp[3],4)*100,"%)",sep=""),
                  pch = PCA.pch.3D,
                  color = PCA.col.plot,
                  col.axis = "blue",
                  cex.symbols = PCA.cex.symbols,
                  cex.lab = PCA.cex.lab,
                  cex.axis = 1.2,
                  lwd = PCA.lwd.3D,
                  angle = 55,
                  scale.y = 0.7)
    if (!is.null(PCA.legend)) {
    legend_pch <- PCA.legend$pch
    if (is.null(legend_pch)) legend_pch <- PCA.pch.3D
    legend_inset <- PCA.legend$inset.3D
    if (is.null(legend_inset)) legend_inset <- PCA.legend$inset
    if (is.null(legend_inset)) legend_inset <- PCA.legend.inset.3D
    legend_bty <- PCA.legend$bty
    if (is.null(legend_bty)) legend_bty <- "o"
    legend_box_lwd <- PCA.legend$box.lwd
    if (is.null(legend_box_lwd)) legend_box_lwd <- 1
    legend_cex <- PCA.legend$cex
    if (is.null(legend_cex)) legend_cex <- PCA.legend.cex
    old_par <- graphics::par(xpd=if (legend_outside3d) NA else FALSE)
    legend(PCA.legend$legend.pos, legend = PCA.legend$taxa, pch = legend_pch, col = PCA.legend$color,
           ncol = PCA.legend$ncol, bty = legend_bty, inset = legend_inset, box.lwd = legend_box_lwd, cex = legend_cex)
    graphics::par(old_par)
}
grDevices::dev.off()
print("Joining taxa...")
#Extract number of PCs needed
PCs <- cbind(taxa,as.data.frame(PCA.X$x))

#Remove duplicate (This is taken care by QC)
#PCs.unique <- unique(PCs[,1])
#PCs <-PCs[match(PCs.unique, PCs[,1], nomatch = 0), ]



print("Exporting PCs...")
#Write the PCs into a text file
if(file.output) utils::write.table(PCs[,1:(PCA.total+1)], "GAPIT.Genotype.PCA.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

# if(file.output) utils::write.table(PCA.X$rotation[,1:PC.number], "GAPIT.Genotype.PCA_loadings.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

if(file.output) utils::write.table(eigenvalues, "GAPIT.Genotype.PCA_eigenvalues.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#Return the PCs
return(list(PCs=PCs,EV=PCA.X$sdev^2,nPCs=NULL))}
