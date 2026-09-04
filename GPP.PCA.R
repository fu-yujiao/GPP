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
  PCA.mar.right.3D=12,
  PCA.dynamic3D.backend=c("plotly","rgl","python"),
  PCA.python.cmd=NULL,
  PCA.dynamic3D.title="Interactive 3D PCA",
  PCA.dynamic3D.marker.size=3,
  PCA.dynamic3D.opacity=0.7){

  PCA.dynamic3D.backend <- match.arg(PCA.dynamic3D.backend)

  make_group_from_legend <- function(PCA.col, taxa, PCA.legend) {
    n <- length(taxa)
    g <- rep_len("Other", n)
    if (!is.null(PCA.legend) && !is.null(PCA.legend$taxa) && !is.null(PCA.legend$color)) {
      labels <- as.character(PCA.legend$taxa)
      cols   <- as.character(PCA.legend$color)
      keep <- !is.na(labels) & nzchar(labels) & !is.na(cols) & nzchar(cols)
      labels <- labels[keep]; cols <- cols[keep]
      if (length(labels) > 0) {
        ac <- as.character(PCA.col)
        for (k in seq_along(labels)) {
          g[!is.na(ac) & ac == cols[k]] <- labels[k]
        }
      }
    }
    g
  }

  GAPIT.3D.PCA.plotly <- function(PCs_matrix, evp, colors, groups, taxa_vec,
                                  title, marker_size, marker_opacity,
                                  output_dir) {
    if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")
    stopifnot(requireNamespace("plotly", quietly = TRUE))
    stopifnot(requireNamespace("htmlwidgets", quietly = TRUE))

    pandoc_found <- tryCatch({
      if (requireNamespace("rmarkdown", quietly = TRUE) &&
          isTRUE(rmarkdown::pandoc_available())) {
        TRUE
      } else {
        cand <- c(
          Sys.which("pandoc"),
          Sys.which("pandoc.exe"),
          file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc", "pandoc.exe"),
          file.path(Sys.getenv("ProgramFiles"), "Pandoc", "pandoc.exe"),
          file.path(Sys.getenv("ProgramFiles(x86)"), "Pandoc", "pandoc.exe")
        )
        extra <- tryCatch({
          dir(file.path(Sys.getenv("ProgramFiles"), "RStudio", "resources"),
              pattern = "pandoc\\.exe$", recursive = TRUE, full.names = TRUE)
        }, error = function(e) character(0))
        cand <- c(cand, extra)
        cand <- cand[nzchar(cand) & file.exists(cand)]
        if (length(cand) > 0) {
          Sys.setenv(RSTUDIO_PANDOC = dirname(cand[1]))
          TRUE
        } else {
          FALSE
        }
      }
    }, error = function(e) FALSE)

    pc1_var <- evp[1] * 100; pc2_var <- evp[2] * 100; pc3_var <- evp[3] * 100
    total_var <- sum(c(pc1_var, pc2_var, pc3_var))

    x_label <- sprintf("PC1 (%.2f%%)", pc1_var)
    y_label <- sprintf("PC2 (%.2f%%)", pc2_var)
    z_label <- sprintf("PC3 (%.2f%%)", pc3_var)
    main_title <- sprintf("%s \u2014 GPP.Genotype.PCA", title)
    sub_title <- sprintf("Cumulative variance explained by first 3 PCs: %.2f%%", total_var)

    df <- data.frame(
      taxa = as.character(taxa_vec),
      PC1  = as.numeric(PCs_matrix[,1]),
      PC2  = as.numeric(PCs_matrix[,2]),
      PC3  = as.numeric(PCs_matrix[,3]),
      color = as.character(colors),
      group = as.character(groups),
      stringsAsFactors = FALSE
    )
    uniq_group <- unique(df$group)
    color_map <- stats::setNames(rep_len("#4C78A8", length(uniq_group)), uniq_group)
    if (!is.null(PCA.legend) && !is.null(PCA.legend$taxa) && !is.null(PCA.legend$color)) {
      labels <- as.character(PCA.legend$taxa)
      cols   <- as.character(PCA.legend$color)
      keep   <- !is.na(labels) & nzchar(labels) & labels %in% uniq_group
      labels <- labels[keep]; cols <- cols[keep]
      if (length(labels) > 0) color_map[labels] <- cols
    }

    p <- plotly::plot_ly(
      data = df,
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = ~group,
      colors = color_map,
      marker = list(size = marker_size, line = list(width = 0)),
      opacity = marker_opacity,
      type = "scatter3d",
      mode = "markers",
      hovertext = ~paste0("Taxa: ", taxa, "<br>Group: ", group),
      hoverinfo = "text+x+y+z",
      showlegend = FALSE
    )
    line_color <- "#d3d3d3"
    axis_style <- function(lab) {
      list(
        title = list(text = lab),
        showgrid = TRUE,
        gridcolor = line_color,
        gridwidth = 1,
        showzeroline = TRUE,
        zerolinecolor = line_color,
        zerolinewidth = 1,
        showline = TRUE,
        linecolor = line_color,
        linewidth = 1,
        showspikes = FALSE,
        showticklabels = TRUE
      )
    }
    p <- plotly::layout(p,
      paper_bgcolor = "#ffffff",
      plot_bgcolor = "#ffffff",
      title = list(text = paste0("<b>", main_title, "</b><br><sup>", sub_title, "</sup>"),
                   x = 0.02, xanchor = "left"),
      scene = list(
        xaxis = axis_style(x_label),
        yaxis = axis_style(y_label),
        zaxis = axis_style(z_label),
        bgcolor = "#ffffff",
        camera = list(eye = list(x = 1.4, y = 1.4, z = 1.2))
      ),
      showlegend = FALSE,
      margin = list(l = 0, r = 0, t = 90, b = 0)
    )

    out_html <- file.path(output_dir, "GAPIT.Genotype.PCA_3DPCA.html")
    if (pandoc_found) {
      htmlwidgets::saveWidget(p, file = out_html, selfcontained = TRUE, title = main_title)
    } else {
      savedir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
      wd0 <- setwd(savedir)
      on.exit({ tryCatch(setwd(wd0), error = function(e) invisible(NULL)) }, add = TRUE)
      htmlwidgets::saveWidget(p, file = basename(out_html), selfcontained = FALSE,
                              title = main_title,
                              libdir = paste0(tools::file_path_sans_ext(basename(out_html)), "_files"))
    }
    out_html
  }
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
  wd0 <- getwd()
  on.exit({
    tryCatch(setwd(wd0), error = function(e) invisible(NULL))
  }, add = TRUE)

  color_file <- NULL
  if (!is.null(PCA.col) && !is.null(taxa) && length(PCA.col) == length(taxa)) {
    mycols <- data.frame(taxa = as.character(taxa), color = as.character(PCA.col),
                         stringsAsFactors = FALSE)
    color_file <- file.path(getwd(), "color_file.csv")
    try(utils::write.csv(mycols, color_file, quote = FALSE, row.names = FALSE), silent = TRUE)
    if (!file.exists(color_file)) color_file <- NULL
  }

  legend_groups <- make_group_from_legend(PCA.col, taxa, PCA.legend)

  backend_used <- NA_character_
  use_plotly <- identical(tolower(PCA.dynamic3D.backend), "plotly")
  use_python <- identical(tolower(PCA.dynamic3D.backend), "python")
  use_rgl    <- identical(tolower(PCA.dynamic3D.backend), "rgl")

  if (use_plotly) {
    ok <- tryCatch({
      out <- GAPIT.3D.PCA.plotly(
        PCs_matrix = as.matrix(PCA.X$x[, 1:min(3, ncol(PCA.X$x))]),
        evp = evp,
        colors = PCA.col,
        groups = legend_groups,
        taxa_vec = taxa,
        title = PCA.dynamic3D.title,
        marker_size = PCA.dynamic3D.marker.size,
        marker_opacity = PCA.dynamic3D.opacity,
        output_dir = getwd()
      )
      is.character(out) && length(out) == 1 && nzchar(out) && file.exists(out)
    }, error = function(e) {
      warning("GAPIT.3D.PCA.plotly failed: ", e$message)
      FALSE
    })
    if (ok) backend_used <- "plotly"
  }

  if (is.na(backend_used) && (use_python || (!use_rgl && exists("GAPIT.3D.PCA.python", mode = "function")))) {
    ok <- tryCatch({
      fun <- get("GAPIT.3D.PCA.python", mode = "function")
      res <- do.call(fun, list(
        color_file = color_file,
        input_dir = getwd(),
        output_dir = getwd(),
        python_cmd = PCA.python.cmd,
        title = PCA.dynamic3D.title,
        marker_size = PCA.dynamic3D.marker.size,
        marker_opacity = PCA.dynamic3D.opacity
      ))
      isTRUE(res) || (length(dir(getwd(), pattern = "_3DPCA\\.html$", ignore.case = TRUE)) > 0)
    }, error = function(e) {
      warning("GAPIT.3D.PCA.python failed: ", e$message)
      FALSE
    })
    if (ok) backend_used <- "python"
  }

  if (is.na(backend_used)) {
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
    backend_used <- "rgl"
  }

  cat("Dynamic 3D backend used:", backend_used, "\n")
}
#    if(!require(scatterplot3d)) install.packages("scatterplot3d")
#    library(scatterplot3d)

    if(!require(scatterplot3d)) install.packages("scatterplot3d")

    grDevices::pdf("GAPIT.Genotype.PCA_3D.pdf", width = 7, height = 7)
    legend_outside3d <- isTRUE(PCA.legend.outside.3D)
    if (!is.null(PCA.legend)) {
      if (isTRUE(PCA.legend$outside.3D) || isTRUE(PCA.legend.outside.3D)) {
        legend_outside3d <- TRUE
      } else if (identical(PCA.legend$outside.3D, FALSE) || identical(PCA.legend.outside.3D, FALSE)) {
        legend_outside3d <- FALSE
      }
    }
    if (legend_outside3d) {
      graphics::par(mar = c(5,5,5,PCA.mar.right.3D),xpd=NA)
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
    if (is.null(legend_pch) || length(legend_pch) == 0) legend_pch <- PCA.pch.3D
    legend_inset <- PCA.legend$inset.3D
    if (is.null(legend_inset) || length(legend_inset) == 0) legend_inset <- PCA.legend.inset.3D
    if (is.null(legend_inset) || length(legend_inset) == 0) legend_inset <- PCA.legend$inset
    if (is.null(legend_inset) || length(legend_inset) == 0) legend_inset <- 0.02
    legend_bty <- PCA.legend$bty
    if (is.null(legend_bty) || !nzchar(legend_bty)) legend_bty <- "n"
    legend_box_lwd <- PCA.legend$box.lwd
    if (is.null(legend_box_lwd) || length(legend_box_lwd) == 0) legend_box_lwd <- 1
    legend_cex <- PCA.legend$cex
    if (is.null(legend_cex) || length(legend_cex) == 0) legend_cex <- PCA.legend.cex
    if (is.null(legend_cex) || length(legend_cex) == 0) legend_cex <- 1
    legend_pos <- PCA.legend$legend.pos
    if (is.null(legend_pos) || !nzchar(legend_pos)) legend_pos <- "topright"
    legend_ncol <- PCA.legend$ncol
    if (is.null(legend_ncol) || !is.finite(as.integer(legend_ncol))) legend_ncol <- 1
    old_par <- graphics::par(xpd=NA)
    legend(legend_pos, legend = PCA.legend$taxa, pch = legend_pch, col = PCA.legend$color,
           ncol = as.integer(legend_ncol), bty = legend_bty, inset = legend_inset,
           box.lwd = legend_box_lwd, cex = legend_cex)
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


############################################################################################
#  GAPIT.3D.PCA.python: from gapit_functions2.txt
############################################################################################
GAPIT.3D.PCA.python <- function(color_file = NULL, input_dir = getwd(), output_dir = getwd(),
                                  python_cmd = NULL,
                                  title = "Interactive 3D PCA",
                                  marker_size = 3,
                                  marker_opacity = 0.7,
                                  plotly_template = "plotly_white",
                                  offline = TRUE) {
  if (is.null(python_cmd)) {
    python_cmd <- Sys.which("python3")
    if (nchar(python_cmd) == 0) python_cmd <- Sys.which("python")
    if (nchar(python_cmd) == 0) python_cmd <- "python"
  }

  marker_size <- as.numeric(marker_size)
  if (is.na(marker_size) || length(marker_size) == 0) marker_size <- 3
  marker_opacity <- as.numeric(marker_opacity)
  if (is.na(marker_opacity) || length(marker_opacity) == 0) marker_opacity <- 0.7
  if (marker_opacity > 1) marker_opacity <- 1
  if (marker_opacity < 0) marker_opacity <- 0
  if (is.null(plotly_template) || !nzchar(plotly_template)) plotly_template <- "plotly_white"
  if (is.null(title) || !nzchar(title)) title <- "Interactive 3D PCA"
  offline <- isTRUE(offline)

  TITLE_PY <- deparse1(title)
  MARKER_SIZE_PY <- as.character(marker_size)
  MARKER_OPACITY_PY <- as.character(marker_opacity)
  PLOTLY_TEMPLATE_PY <- deparse1(plotly_template)
  OFFLINE_PY <- if (offline) "True" else "False"

  python_code <- paste(
    "import os",
    "import sys",
    "import glob",
    "import argparse",
    "import pandas as pd",
    "import plotly.express as px",
    "import numpy as np",
    "",
    paste0("TITLE = ", TITLE_PY),
    paste0("MARKER_SIZE = ", MARKER_SIZE_PY),
    paste0("MARKER_OPACITY = ", MARKER_OPACITY_PY),
    paste0("PLOTLY_TEMPLATE = ", PLOTLY_TEMPLATE_PY),
    paste0("OFFLINE = ", OFFLINE_PY),
    "",
    "def parse_arguments():",
    "    parser = argparse.ArgumentParser(description='Generate interactive 3D visualization for GAPIT PCA results')",
    "    parser.add_argument('-i', '--input', required=True)",
    "    parser.add_argument('-o', '--output', required=True)",
    "    parser.add_argument('-c', '--color', default=None)",
    "    return parser.parse_args()",
    "",
    "def ensure_dirs(input_dir, output_dir):",
    "    if not os.path.isdir(input_dir):",
    "        print('Error: Input directory does not exist:', input_dir)",
    "        sys.exit(1)",
    "    os.makedirs(output_dir, exist_ok=True)",
    "",
    "def load_color_file(color_file_path, pca_df):",
    "    try:",
    "        color_df = pd.read_csv(color_file_path)",
    "        if color_df.shape[1] >= 2:",
    "            color_df = color_df.iloc[:, :2]",
    "            color_df.columns = ['taxa', 'color']",
    "        else:",
    "            raise ValueError('Color file must have at least 2 columns')",
    "        color_dict = dict(zip(color_df['taxa'], color_df['color']))",
    "        pca_df['color'] = pca_df['taxa'].map(color_dict)",
    "        missing_colors = pca_df['color'].isna().sum()",
    "        if missing_colors > 0:",
    "            print('  Warning:', missing_colors, 'samples have no corresponding color, using default color')",
    "            pca_df.loc[pca_df['color'].isna(), 'color'] = 'black'",
    "        return pca_df",
    "    except Exception as e:",
    "        print('  Warning: Failed to read color file', color_file_path, ':', e)",
    "        pca_df['color'] = 'blue'",
    "        return pca_df",
    "",
    "def load_eigenvalues(eigenvalues_path):",
    "    try:",
    "        eigen_df = pd.read_csv(eigenvalues_path, header=0)",
    "        eigenvalues = eigen_df.iloc[:, 0].tolist()",
    "        return eigenvalues",
    "    except Exception as e:",
    "        print('  Warning: Failed to parse eigenvalue file', eigenvalues_path, ':', e)",
    "        return []",
    "",
    "def load_gapit_pca(pca_path, eigenvalues_path=None):",
    "    try:",
    "        df = pd.read_csv(pca_path)",
    "    except Exception:",
    "        try:",
    "            df = pd.read_csv(pca_path, sep='\\t')",
    "        except Exception:",
    "            df = pd.read_csv(pca_path, sep='\\s+')",
    "    if df.shape[1] < 4:",
    "        raise ValueError('Insufficient columns (%d), at least 4 needed: sample, PC1, PC2, PC3' % df.shape[1])",
    "    if df.shape[1] == 4:",
    "        df.columns = ['taxa', 'PC1', 'PC2', 'PC3']",
    "    else:",
    "        df = df.iloc[:, :4]",
    "        df.columns = ['taxa', 'PC1', 'PC2', 'PC3']",
    "    if eigenvalues_path and os.path.exists(eigenvalues_path):",
    "        try:",
    "            eigenvalues = load_eigenvalues(eigenvalues_path)",
    "            total_variance = np.sum(eigenvalues)",
    "            pc1_var = eigenvalues[0] / total_variance * 100 if len(eigenvalues) > 0 else 0",
    "            pc2_var = eigenvalues[1] / total_variance * 100 if len(eigenvalues) > 1 else 0",
    "            pc3_var = eigenvalues[2] / total_variance * 100 if len(eigenvalues) > 2 else 0",
    "            df.attrs['variance_explained'] = {",
    "                'PC1': pc1_var,",
    "                'PC2': pc2_var,",
    "                'PC3': pc3_var",
    "            }",
    "        except Exception as e:",
    "            print('  Warning: Failed to read eigenvalue file', eigenvalues_path, ':', e)",
    "    return df",
    "",
    "def make_figure(pca_df, title):",
    "    variance_explained = pca_df.attrs.get('variance_explained', {})",
    "    if variance_explained:",
    "        x_label = 'PC1 (%.2f%%)' % variance_explained.get('PC1', 0)",
    "        y_label = 'PC2 (%.2f%%)' % variance_explained.get('PC2', 0)",
    "        z_label = 'PC3 (%.2f%%)' % variance_explained.get('PC3', 0)",
    "    else:",
    "        x_label, y_label, z_label = 'PC1', 'PC2', 'PC3'",
    "    if 'color' in pca_df.columns:",
    "        fig = px.scatter_3d(",
    "            data_frame=pca_df,",
    "            x='PC1',",
    "            y='PC2',",
    "            z='PC3',",
    "            hover_name='taxa',",
    "            title=title,",
    "            labels={'PC1': x_label, 'PC2': y_label, 'PC3': z_label},",
    "            template=PLOTLY_TEMPLATE,",
    "            color='color',",
    "            color_discrete_map='identity'",
    "        )",
    "        fig.update_traces(marker=dict(size=MARKER_SIZE, opacity=MARKER_OPACITY))",
    "    else:",
    "        fig = px.scatter_3d(",
    "            data_frame=pca_df,",
    "            x='PC1',",
    "            y='PC2',",
    "            z='PC3',",
    "            hover_name='taxa',",
    "            title=title,",
    "            labels={'PC1': x_label, 'PC2': y_label, 'PC3': z_label},",
    "            template=PLOTLY_TEMPLATE",
    "        )",
    "        fig.update_traces(marker=dict(size=MARKER_SIZE, opacity=MARKER_OPACITY, color='blue'))",
    "    if variance_explained:",
    "        total_var = sum([variance_explained.get('PC%d' % i, 0) for i in range(1, 4)])",
    "        fig.update_layout(",
    "            title_text='%s<br><sup>Cumulative variance explained by first 3 PCs: %.2f%%</sup>' % (title, total_var)",
    "        )",
    "    fig.update_layout(",
    "        scene_camera=dict(eye=dict(x=1.4, y=1.4, z=1.2)),",
    "        margin=dict(l=0, r=0, t=80, b=0)",
    "    )",
    "    return fig",
    "",
    "def export_html(fig, out_path, offline=True):",
    "    include = 'inline' if offline else 'cdn'",
    "    fig.write_html(out_path, include_plotlyjs=include, full_html=True)",
    "",
    "def main():",
    "    args = parse_arguments()",
    "    input_dir = args.input",
    "    output_dir = args.output",
    "    color_file = args.color",
    "    ensure_dirs(input_dir, output_dir)",
    "    pca_path = os.path.join(input_dir, 'GAPIT.Genotype.PCA.csv')",
    "    eigenvalues_path = os.path.join(input_dir, 'GAPIT.Genotype.PCA_eigenvalues.csv')",
    "    if not os.path.exists(pca_path):",
    "        print('Error: PCA file not found:', pca_path)",
    "        return",
    "    if not os.path.exists(eigenvalues_path):",
    "        print('  No eigenvalue file found, using default axis labels')",
    "        eigenvalues_path = None",
    "    else:",
    "        print('  Eigenvalue file found')",
    "    ok, skip, fail = 0, 0, 0",
    "    base = os.path.basename(pca_path)",
    "    stem, ext = os.path.splitext(base)",
    "    out_html = os.path.join(output_dir, '%s_3DPCA.html' % stem)",
    "    print('Processing:', base)",
    "    try:",
    "        df = load_gapit_pca(pca_path, eigenvalues_path)",
    "        print('  Successfully read data:', df.shape[0], 'samples,', df.shape[1], 'columns')",
    "        if color_file and os.path.exists(color_file):",
    "            df = load_color_file(color_file, df)",
    "            print('  Color information loaded')",
    "        else:",
    "            print('  No color file provided or file not found, using default color')",
    "            df['color'] = 'red'",
    "    except ValueError as ve:",
    "        print('  Skipped:', ve)",
    "        skip += 1",
    "    except Exception as e:",
    "        print('  Failed to read:', e)",
    "        fail += 1",
    "    else:",
    "        try:",
    "            fig = make_figure(df, TITLE)",
    "            export_html(fig, out_html, offline=OFFLINE)",
    "            print('  Saved:', out_html)",
    "            ok += 1",
    "        except Exception as e:",
    "            print('  Failed to export:', e)",
    "            fail += 1",
    "    print('Processing complete!')",
    "",
    "if __name__ == '__main__':",
    "    main()",
    sep = "\n"
  )

  python_file <- "temp_GAPIT_3DPCA.py"
  writeLines(python_code, python_file, useBytes = FALSE)
  on.exit({
    if (file.exists(python_file)) {
      try(file.remove(python_file), silent = TRUE)
    }
  }, add = TRUE)

  out_html_expected <- file.path(output_dir, paste0(
    sub("\\.csv$", "", basename(Sys.glob(file.path(input_dir, "GAPIT.Genotype.PCA.csv"))[1])),
    "_3DPCA.html"
  ))
  if (is.na(out_html_expected) || !nzchar(out_html_expected)) {
    out_html_expected <- file.path(output_dir, "GAPIT.Genotype.PCA_3DPCA.html")
  }

  if (.Platform$OS.type == "windows") {
    python_file_win <- normalizePath(python_file, winslash = "\\", mustWork = FALSE)
    input_win    <- normalizePath(input_dir,   winslash = "\\", mustWork = FALSE)
    output_win   <- normalizePath(output_dir,  winslash = "\\", mustWork = FALSE)
    cmd <- paste(paste0("\"", python_cmd, "\""),
                 paste0("\"", python_file_win, "\""),
                 "-i", paste0("\"", input_win, "\""),
                 "-o", paste0("\"", output_win, "\""))
    if (!is.null(color_file) && file.exists(color_file)) {
      color_win <- normalizePath(color_file, winslash = "\\", mustWork = FALSE)
      cmd <- paste(cmd, "-c", paste0("\"", color_win, "\""))
    }
  } else {
    cmd <- paste(shQuote(python_cmd), shQuote(python_file),
                "-i", shQuote(normalizePath(input_dir, mustWork = FALSE)),
                "-o", shQuote(normalizePath(output_dir, mustWork = FALSE)))
    if (!is.null(color_file) && file.exists(color_file)) {
      cmd <- paste(cmd, "-c", shQuote(normalizePath(color_file)))
    }
  }
  if (!is.null(color_file) && file.exists(color_file)) {
    cat("Using color file:", color_file, "\n")
  } else if (!is.null(color_file)) {
    warning("Color file not found: ", color_file)
  }
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Running Python script with:", python_cmd, "\n")

  if (file.exists(out_html_expected)) {
    try(file.remove(out_html_expected), silent = TRUE)
  }

  py_stdout <- tempfile(fileext = ".txt")
  py_stderr <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(py_stdout)) try(file.remove(py_stdout), silent = TRUE)
    if (file.exists(py_stderr)) try(file.remove(py_stderr), silent = TRUE)
  }, add = TRUE)

  status <- tryCatch({
    if (.Platform$OS.type == "windows") {
      .cmd <- paste(cmd, "1>", paste0("\"", py_stdout, "\""),
                          "2>", paste0("\"", py_stderr, "\""))
      system(.cmd)
    } else {
      system2(python_cmd, args = c(
        python_file,
        "-i", normalizePath(input_dir, mustWork = FALSE),
        "-o", normalizePath(output_dir, mustWork = FALSE),
        if (!is.null(color_file) && file.exists(color_file)) c("-c", normalizePath(color_file)) else character(0)
      ), stdout = py_stdout, stderr = py_stderr)
    }
  }, error = function(e) {
    warning("Python script failed to run: ", e$message)
    999L
  })

  if (file.exists(py_stdout)) {
    out_txt <- try(readLines(py_stdout, warn = FALSE), silent = TRUE)
    if (inherits(out_txt, "character") && length(out_txt) > 0) {
      cat(paste0("  [py-out] ", out_txt, collapse = "\n"), "\n")
    }
  }
  if (file.exists(py_stderr)) {
    err_txt <- try(readLines(py_stderr, warn = FALSE), silent = TRUE)
    if (inherits(err_txt, "character") && length(err_txt) > 0 && any(nzchar(err_txt))) {
      message(paste0("  [py-err] ", err_txt, collapse = "\n"))
    }
  }

  html_found <- file.exists(out_html_expected)
  if (!html_found) {
    cand <- dir(output_dir, pattern = "_3DPCA\\.html$", ignore.case = TRUE, full.names = TRUE)
    if (length(cand) > 0) {
      html_found <- TRUE
      out_html_expected <- cand[1]
    }
  }
  if (html_found) {
    cat("Python script finished! -> Saved:", out_html_expected, "\n")
  } else {
    cat("Python script finished, but NO _3DPCA.html output was produced (see [py-out]/[py-err] above)\n")
  }
  invisible(html_found)
}

