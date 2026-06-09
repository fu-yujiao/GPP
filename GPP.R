GPP <- function(file_path,
                graph = c("PCA", "genotype", "GWAS"),
                PCA.col = NULL,
                PCA.3d = NULL,
                PCA.legend = NULL,
                PCA.pch.2D = 20,
                PCA.pch.3D = 20,
                PCA.cex.lab = 1.4,
                PCA.cex.symbols = 1,
                PCA.lwd.2D = 2,
                PCA.lwd.3D = 3,
                PCA.bty.2D = "o",
                PCA.2D.pairs = NULL,
                PCA.cex.2D = 1,
                PCA.alpha = 1,
                PCA.legend.cex = 1,
                PCA.legend.outside.3D = FALSE,
                PCA.legend.inset.3D = 0,
                PCA.mar.right.3D = 12,
                PCA.meta_path = NULL,
                PCA.meta_taxa_col = 1,
                PCA.meta_group_col = 2,
                PCA.meta_color_col = NULL,
                name.of.trait = NULL,
                model_store = NULL,
                Y.names = NULL,
                ...) {
  .gpp_do <- function(fun, args) {
    f <- get(fun, mode = "function")
    keep <- intersect(names(args), names(formals(f)))
    do.call(f, args[keep])
  }

  .gpp_is_color <- function(x) {
    if (is.null(x)) return(FALSE)
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) return(FALSE)
    tolower(x) %in% tolower(grDevices::colors())
  }

  .gpp_read_meta <- function(meta_path) {
    cls <- tryCatch(
      utils::read.csv(meta_path, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(cls) || ncol(cls) < 2) {
      cls <- utils::read.table(meta_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
    }
    cls
  }

  .gpp_pick_col_by_name <- function(df, patterns) {
    if (is.null(df) || ncol(df) == 0) return(NA_integer_)
    nms <- names(df)
    if (is.null(nms)) return(NA_integer_)
    low <- tolower(trimws(nms))
    pats <- tolower(patterns)
    for (p in pats) {
      i <- which(low == p)
      if (length(i) > 0) return(i[1])
    }
    for (p in pats) {
      i <- which(grepl(p, low, fixed = TRUE))
      if (length(i) > 0) return(i[1])
    }
    NA_integer_
  }

  .gpp_pca_col_legend_from_meta <- function(taxa, meta_path, taxa_col, group_col, color_col) {
    cls <- .gpp_read_meta(meta_path)
    if (is.null(cls) || ncol(cls) < 2) return(list(col = NULL, legend = NULL))

    taxa_col <- as.integer(taxa_col)
    group_col <- as.integer(group_col)

    if (is.na(taxa_col) || taxa_col < 1 || taxa_col > ncol(cls)) {
      taxa_col <- .gpp_pick_col_by_name(cls, c("taxa", "id", "sample", "line", "name"))
      if (is.na(taxa_col)) taxa_col <- 1L
    }

    if (is.na(group_col) || group_col < 1 || group_col > ncol(cls)) {
      group_col <- .gpp_pick_col_by_name(cls, c("group", "pop", "population", "species", "type", "class"))
    }

    if (is.null(color_col)) {
      color_col <- NA_integer_
      color_col_guess <- .gpp_pick_col_by_name(cls, c("color", "colour", "col"))
      if (!is.na(color_col_guess)) {
        color_col <- color_col_guess
      } else if (ncol(cls) >= 3 && all(.gpp_is_color(cls[[3]]) | is.na(cls[[3]]) | !nzchar(as.character(cls[[3]])))) {
        color_col <- 3L
      } else if (ncol(cls) >= 2 && all(.gpp_is_color(cls[[2]]) | is.na(cls[[2]]) | !nzchar(as.character(cls[[2]])))) {
        color_col <- 2L
      }
    } else {
      color_col <- as.integer(color_col)
      if (is.na(color_col) || color_col < 1 || color_col > ncol(cls)) color_col <- NA_integer_
    }

    cls_taxa <- trimws(as.character(cls[[taxa_col]]))
    cls_group <- if (!is.na(group_col) && group_col >= 1 && group_col <= ncol(cls)) trimws(as.character(cls[[group_col]])) else rep(NA_character_, nrow(cls))
    cls_color <- if (!is.na(color_col)) trimws(as.character(cls[[color_col]])) else rep(NA_character_, nrow(cls))

    taxa_to_group <- stats::setNames(cls_group, cls_taxa)
    taxa_to_color <- stats::setNames(cls_color, cls_taxa)

    group_for_taxa <- unname(taxa_to_group[taxa])
    color_for_taxa <- unname(taxa_to_color[taxa])

    if (!is.na(color_col) && any(.gpp_is_color(color_for_taxa), na.rm = TRUE)) {
      col_out <- rep("grey70", length(taxa))
      ok <- !is.na(color_for_taxa) & nzchar(color_for_taxa) & .gpp_is_color(color_for_taxa)
      col_out[ok] <- color_for_taxa[ok]

      legend_label <- group_for_taxa
      legend_label[is.na(legend_label) | !nzchar(legend_label)] <- col_out[is.na(legend_label) | !nzchar(legend_label)]
      keep <- !is.na(legend_label) & nzchar(legend_label) & col_out != "grey70"
      if (!any(keep)) return(list(col = col_out, legend = NULL))

      df <- data.frame(label = legend_label[keep], color = col_out[keep], stringsAsFactors = FALSE)
      idx <- !duplicated(df$label)
      legend <- list(
        taxa = df$label[idx],
        color = df$color[idx],
        pch = 20,
        legend.pos = "topright",
        ncol = 1,
        inset = 0.02,
        inset.3D = 0,
        bty = "n"
      )
      return(list(col = col_out, legend = legend))
    }

    if (any(.gpp_is_color(group_for_taxa), na.rm = TRUE)) {
      col_out <- rep("grey70", length(taxa))
      ok <- !is.na(group_for_taxa) & nzchar(group_for_taxa) & .gpp_is_color(group_for_taxa)
      col_out[ok] <- group_for_taxa[ok]
      df <- data.frame(label = col_out[ok], color = col_out[ok], stringsAsFactors = FALSE)
      idx <- !duplicated(df$label)
      legend <- list(
        taxa = df$label[idx],
        color = df$color[idx],
        pch = 20,
        legend.pos = "topright",
        ncol = 1,
        inset = 0.02,
        inset.3D = 0,
        bty = "n"
      )
      return(list(col = col_out, legend = legend))
    }

    groups <- group_for_taxa
    groups[is.na(groups) | !nzchar(groups)] <- NA_character_
    uniq <- unique(groups[!is.na(groups)])
    if (length(uniq) == 0) return(list(col = NULL, legend = NULL))

    pal_fun <- NULL
    if (exists("hcl.colors", where = asNamespace("grDevices"), mode = "function")) {
      pal_fun <- function(n) grDevices::hcl.colors(n, palette = "Dark 3")
    } else {
      pal_fun <- function(n) grDevices::rainbow(n)
    }
    cols <- pal_fun(length(uniq))
    names(cols) <- uniq

    col_out <- rep("grey70", length(taxa))
    ok <- !is.na(groups)
    col_out[ok] <- cols[groups[ok]]

    legend <- list(
      taxa = uniq,
      color = unname(cols[uniq]),
      pch = 20,
      legend.pos = "topright",
      ncol = 1,
      inset = 0.02,
      inset.3D = 0,
      bty = "n"
    )
    list(col = col_out, legend = legend)
  }

  graph_vec <- as.character(graph)
  graph_vec <- graph_vec[!is.na(graph_vec) & nzchar(graph_vec)]
  graph_up <- toupper(graph_vec)

  map <- c(
    GENOTYPE = "GENOTYPE",
    G = "GENOTYPE",
    PCA = "PCA",
    GWAS = "GWAS",
    LD = "LD",
    LDPLOT = "LD",
    LDDECAY = "LD",
    MULTISOFTWARE = "GWAS_MULTISOFTWARE",
    GWAS_MULTISOFTWARE = "GWAS_MULTISOFTWARE",
    GWAS.MULTISOFTWARE = "GWAS_MULTISOFTWARE"
  )
  graph_key <- unique(unname(ifelse(graph_up %in% names(map), map[graph_up], graph_up)))

  dots <- list(...)
  result_list <- list()

  if ("PCA" %in% graph_key) {
    PCs <- read.csv(file.path(file_path, "GAPIT.Genotype.PCA.csv"), header = TRUE, stringsAsFactors = FALSE)
    eigenvalues <- read.csv(file.path(file_path, "GAPIT.Genotype.PCA_eigenvalues.csv"), header = TRUE)

    taxa <- trimws(as.character(PCs[[1]]))
    meta_path_use <- PCA.meta_path
    if (is.null(meta_path_use) || !nzchar(meta_path_use)) {
      candidates <- c("file.csv", "PCA.meta.csv", "PCA.group.csv", "PCA.groups.csv")
      found <- candidates[file.exists(file.path(file_path, candidates))]
      if (length(found) > 0) meta_path_use <- file.path(file_path, found[1])
    } else {
      if (!grepl("^[A-Za-z]:[/\\\\]", meta_path_use)) meta_path_use <- file.path(file_path, meta_path_use)
    }

    if ((is.null(PCA.col) || is.null(PCA.legend)) && !is.null(meta_path_use) && file.exists(meta_path_use)) {
      auto <- .gpp_pca_col_legend_from_meta(
        taxa = taxa,
        meta_path = meta_path_use,
        taxa_col = PCA.meta_taxa_col,
        group_col = PCA.meta_group_col,
        color_col = PCA.meta_color_col
      )
      if (is.null(PCA.col) && !is.null(auto$col)) PCA.col <- auto$col
      if (is.null(PCA.legend) && !is.null(auto$legend)) PCA.legend <- auto$legend
    }

    if (!is.null(PCA.col)) {
      PCA.col <- rep(PCA.col, length.out = nrow(PCs))
    } else {
      message("PCA.col is NULL. Proceeding without color mapping.")
    }

    if (!is.null(PCA.legend)) {
      ok_legend <- is.list(PCA.legend) && !is.null(PCA.legend$taxa) && !is.null(PCA.legend$color)
      if (ok_legend) {
        message("Using provided PCA.legend data to generate legend.")
      } else {
        message("Invalid PCA.legend format. Skipping legend generation.")
        PCA.legend <- NULL
      }
    }

    args <- c(
      list(
        PCs = PCs,
        eigenvalues = eigenvalues,
        PCA.total = min(ncol(PCs) - 1, nrow(PCs) - 1),
        PCA.col = PCA.col,
        PCA.3d = PCA.3d,
        PCA.legend = PCA.legend,
        PCA.pch.2D = PCA.pch.2D,
        PCA.pch.3D = PCA.pch.3D,
        PCA.cex.lab = PCA.cex.lab,
        PCA.cex.symbols = PCA.cex.symbols,
        PCA.lwd.2D = PCA.lwd.2D,
        PCA.lwd.3D = PCA.lwd.3D,
        PCA.bty.2D = PCA.bty.2D,
        PCA.2D.pairs = PCA.2D.pairs,
        PCA.cex.2D = PCA.cex.2D,
        PCA.alpha = PCA.alpha,
        PCA.legend.cex = PCA.legend.cex,
        PCA.legend.outside.3D = PCA.legend.outside.3D,
        PCA.legend.inset.3D = PCA.legend.inset.3D,
        PCA.mar.right.3D = PCA.mar.right.3D
      ),
      dots
    )
    result_list$PCA <- .gpp_do("GAPIT.PCA", args)
  }

  if ("GENOTYPE" %in% graph_key) {
    GI <- read.table(file.path(file_path, "mdp_SNP_information.txt"), header = TRUE)
    Frequency_MAF <- read.csv(file.path(file_path, "GAPIT.Genotype.Frequency_MAF.csv"), header = TRUE)
    args <- c(list(GI = GI, Frequency_MAF = Frequency_MAF), dots)
    result_list$genotype <- .gpp_do("GAPIT.Genotype.View", args)
  }

  if ("GWAS" %in% graph_key) {
    if (length(name.of.trait) == 1 && length(model_store) == 1) {
      args <- c(list(model_store = model_store, name.of.trait = name.of.trait), dots)
      result_list$GWAS <- .gpp_do("GAPIT.Manhattan", args)
    } else {
      GM0 <- NULL
      if (!("GM" %in% names(dots)) || is.null(dots$GM)) {
        environ_name <- paste(rep(model_store, each = length(Y.names)), rep(Y.names, times = length(model_store)), sep = ".")
        fn <- file.path(file_path, paste0("GAPIT.Association.GWAS_Results.", environ_name, ".csv"))
        fn <- fn[file.exists(fn)]
        if (length(fn) > 0) {
          tmp <- read.csv(fn[1], header = TRUE, stringsAsFactors = FALSE)
          if (ncol(tmp) >= 3) GM0 <- tmp[, c(1:3)]
        }
      }

      GM_use <- if (("GM" %in% names(dots)) && !is.null(dots$GM)) dots$GM else GM0
      if (!is.null(GM_use)) {
        args_mm <- c(list(model_store = model_store, Y.names = Y.names, GM = GM_use), dots)
        result_list$GWAS_multiple <- .gpp_do("GAPIT.Multiple.Manhattan", args_mm)
      }

      args_circ <- c(list(model_store = model_store, Y.names = Y.names), dots)
      GMM <- .gpp_do("GAPIT.Circle.Manhattan.Plot", args_circ)
      result_list$Circle <- GMM
      if (is.list(GMM)) {
        if (!is.null(GMM$multip_mapP)) result_list$multip_mapP <- GMM$multip_mapP
        if (!is.null(GMM$xz)) result_list$xz <- GMM$xz
      }
    }
  }

  if ("GWAS_MULTISOFTWARE" %in% graph_key) {
    result_list$GWAS_multisoftware <- .gpp_do("GPP.GWAS.Manhattan.Horizontal.MultiSoftware", dots)
  }

  if ("LD" %in% graph_key) {
    if (exists("lddecay_plot", mode = "function")) {
      result_list$LD <- .gpp_do("lddecay_plot", dots)
    } else if (exists("LDplot", mode = "function")) {
      result_list$LD <- .gpp_do("LDplot", dots)
    } else {
      stop("LD module not loaded: lddecay_plot()/LDplot() not found.", call. = FALSE)
    }
  }

  invisible(result_list)
}
