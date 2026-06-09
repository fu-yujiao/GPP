GPP.GWAS.Read.Manhattan <- function(path, trait = NULL, software = NULL, model = NULL) {
  if (is.null(path) || length(path) != 1) stop("path must be a single file path.")
  path <- normalizePath(as.character(path), winslash = "\\", mustWork = TRUE)
  ext <- tolower(tools::file_ext(path))
  df <- tryCatch(
    {
      if (ext %in% c("csv")) {
        utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        utils::read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      }
    },
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) < 1) stop("Failed to read file: ", path)
  cn <- names(df)
  cn_norm <- gsub("[^a-z0-9]+", "_", tolower(cn))
  infer_software <- function(df, path) {
    bn <- basename(path)
    if (grepl("^gapit\\.", bn, ignore.case = TRUE) || any(cn_norm %in% c("h_b_p_value", "h_b_p_value_"))) return("GAPIT")
    if (any(cn_norm %in% c("test", "beta", "stat", "nmiss")) && any(cn_norm %in% c("p"))) return("PLINK")
    if (any(grepl("^trait_", cn_norm)) || grepl("^trait\\.", bn, ignore.case = TRUE)) return("RMVP")
    if (any(cn_norm %in% c("se", "effect")) && any(cn_norm %in% c("maf"))) return("RMVP")
    "Unknown"
  }
  infer_model_trait <- function(path, df, software, trait_in) {
    bn <- basename(path)
    ext0 <- tolower(tools::file_ext(bn))
    bn0 <- if (ext0 %in% c("csv", "txt", "tsv")) tools::file_path_sans_ext(bn) else bn
    out <- list(model = NULL, trait = trait_in)
    if (software == "GAPIT") {
      m <- regexec("^GAPIT\\.Association\\.GWAS_Results\\.([^.]+)\\.([^.]+)$", bn0, ignore.case = TRUE)
      mm <- regmatches(bn0, m)[[1]]
      if (length(mm) == 3) {
        out$model <- mm[2]
        out$trait <- mm[3]
      }
      if (is.null(out$model) || is.null(out$trait)) {
        m2 <- regexec("^GAPIT\\.([^.]+)\\.([^.]+)", bn0, ignore.case = TRUE)
        mm2 <- regmatches(bn0, m2)[[1]]
        if (length(mm2) >= 3) {
          out$model <- mm2[2]
          out$trait <- mm2[3]
        }
      }
    } else if (software == "PLINK") {
      m <- regexec("^GWAS_([^\\.]+)\\.assoc\\.(.+)$", bn0, ignore.case = TRUE)
      mm <- regmatches(bn0, m)[[1]]
      if (length(mm) == 3) {
        out$trait <- mm[2]
        out$model <- paste0("assoc.", mm[3])
      } else {
        if (grepl("\\.assoc\\.", bn0, ignore.case = TRUE)) {
          out$model <- sub("^.*\\.assoc\\.", "assoc.", bn0, ignore.case = TRUE)
        }
      }
      if (is.null(out$trait) || is.null(out$model)) {
        m2 <- regexec("^PLINK\\.([^.]+)\\.([^.]+)", bn0, ignore.case = TRUE)
        mm2 <- regmatches(bn0, m2)[[1]]
        if (length(mm2) >= 3) {
          out$model <- mm2[2]
          out$trait <- mm2[3]
        } else {
          parts <- strsplit(bn0, "\\.")[[1]]
          if (length(parts) >= 3 && toupper(parts[1]) == "PLINK") {
            out$model <- parts[2]
            out$trait <- parts[3]
          }
        }
      }
    } else if (software == "RMVP") {
      parts <- strsplit(bn0, "\\.")[[1]]
      if (length(parts) >= 2) {
        if (!is.null(parts[1]) && !tolower(parts[1]) %in% c("trait", "gwas")) out$trait <- parts[1]
        out$model <- parts[length(parts)]
      }
      if (is.null(out$trait)) {
        pcols <- cn[grepl("^trait\\.", cn, ignore.case = TRUE)]
        if (length(pcols) >= 1) {
          pparts <- strsplit(pcols[1], "\\.")[[1]]
          if (length(pparts) >= 2 && !tolower(pparts[1]) %in% c("trait")) out$trait <- pparts[1]
          if (length(pparts) >= 2 && is.null(out$model)) out$model <- pparts[length(pparts)]
        }
      }
      m2 <- regexec("^RMVP\\.([^.]+)\\.([^.]+)", bn0, ignore.case = TRUE)
      mm2 <- regmatches(bn0, m2)[[1]]
      if (length(mm2) >= 3) {
        out$model <- mm2[2]
        out$trait <- mm2[3]
      }
    }
    if (is.null(out$trait) && !is.null(trait_in)) out$trait <- trait_in
    out
  }
  if (is.null(software)) software <- infer_software(df, path)
  mt <- infer_model_trait(path, df, software, trait)
  if (is.null(model)) model <- mt$model
  trait <- mt$trait
  pick_col <- function(norm_names, candidates) {
    idx <- which(norm_names %in% candidates)
    if (length(idx) >= 1) return(idx[1])
    integer(0)
  }
  to_num <- function(x) suppressWarnings(as.numeric(x))
  snp_idx <- pick_col(cn_norm, c("snp", "rs", "rsid", "marker", "id"))
  chr_idx <- pick_col(cn_norm, c("chr", "chrom", "chromosome"))
  bp_idx <- pick_col(cn_norm, c("bp", "pos", "position", "bp_position"))
  p_candidates <- which(cn_norm %in% c("p", "p_value", "pvalue", "pval", "p_value_", "p.value", "p_value_wald", "p_wald"))
  p_candidates <- unique(c(p_candidates, which(grepl("^trait_", cn_norm))))
  pick_p <- function(df, idxs) {
    if (length(idxs) < 1) return(integer(0))
    scores <- rep(-Inf, length(idxs))
    for (i in seq_along(idxs)) {
      v <- to_num(df[[idxs[i]]])
      ok <- is.finite(v) & v > 0 & v <= 1
      scores[i] <- sum(ok)
    }
    if (all(scores <= 0)) return(integer(0))
    idxs[which.max(scores)]
  }
  p_idx <- pick_p(df, p_candidates)
  if (length(p_idx) < 1) {
    num_cols <- which(vapply(df, function(x) is.numeric(suppressWarnings(as.numeric(x))), logical(1)))
    p_idx <- pick_p(df, num_cols)
  }
  if (length(chr_idx) < 1 || length(bp_idx) < 1 || length(p_idx) < 1) {
    stop("Cannot find required columns (CHR/BP/P) in file: ", path)
  }
  out <- data.frame(
    SNP = if (length(snp_idx) >= 1) as.character(df[[snp_idx]]) else NA_character_,
    CHR = df[[chr_idx]],
    BP = df[[bp_idx]],
    P = df[[p_idx]],
    stringsAsFactors = FALSE
  )
  out$CHR <- as.character(out$CHR)
  out$BP <- to_num(out$BP)
  out$P <- to_num(out$P)
  out$Software <- software
  out$Model <- if (is.null(model) || is.na(model) || model == "") NA_character_ else as.character(model)
  out$Trait <- if (is.null(trait) || is.na(trait) || trait == "") NA_character_ else as.character(trait)
  out$File <- basename(path)
  out <- out[is.finite(out$BP) & is.finite(out$P), , drop = FALSE]
  out <- out[!is.na(out$CHR) & out$CHR != "", , drop = FALSE]
  out
}

GPP.GWAS.Significant.Table <- function(dat,
                                      cut_off = 0.05,
                                      max_per_chr = 1,
                                      label_priority = c("Gene", "SNP"),
                                      annot = NULL,
                                      window_kb = 10,
                                      dedup_by_position = TRUE,
                                      dedup_pos_bin = 1) {
  if (is.null(dat) || nrow(dat) < 1) return(NULL)
  label_priority <- match.arg(label_priority, several.ok = TRUE)
  d <- dat
  if (!is.null(annot)) {
    ann <- annot
    if (is.character(ann) && length(ann) == 1) {
      ext <- tolower(tools::file_ext(ann))
      if (ext %in% c("gtf", "gff", "gff3")) {
        ann <- tryCatch({
          if (requireNamespace("data.table", quietly = TRUE)) {
            tmp <- data.table::fread(ann, sep = "\t", header = FALSE, stringsAsFactors = FALSE, select = c(1, 3, 4, 5, 9))
            tmp <- tmp[tmp[[2]] == "gene", ]
            g_id <- sapply(tmp[[5]], function(x) {
              m <- regexpr("gene_name\\s+\"([^\"]+)\"", x)
              if (m == -1) m <- regexpr("gene_id\\s+\"([^\"]+)\"", x)
              if (m != -1) sub(".*\"([^\"]+)\".*", "\\1", regmatches(x, m)) else NA_character_
            })
            data.frame(CHR = as.character(tmp[[1]]), START = tmp[[3]], END = tmp[[4]], Gene = g_id, stringsAsFactors = FALSE)
          } else {
            tmp <- utils::read.table(ann, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
            tmp <- tmp[tmp[, 3] == "gene", ]
            g_id <- sapply(tmp[, 9], function(x) {
              m <- regexpr("gene_name\\s+\"([^\"]+)\"", x)
              if (m == -1) m <- regexpr("gene_id\\s+\"([^\"]+)\"", x)
              if (m != -1) sub(".*\"([^\"]+)\".*", "\\1", regmatches(x, m)) else NA_character_
            })
            data.frame(CHR = as.character(tmp[, 1]), START = as.numeric(tmp[, 4]), END = as.numeric(tmp[, 5]), Gene = g_id, stringsAsFactors = FALSE)
          }
        }, error = function(e) NULL)
        
        if (!is.null(ann) && nrow(ann) > 0) {
          ann <- ann[!is.na(ann$Gene) & ann$Gene != "", ]
          d$SNP_Start <- d$BP - window_kb * 1000
          d$SNP_End <- d$BP + window_kb * 1000
          
          if (requireNamespace("data.table", quietly = TRUE)) {
            dt_d <- data.table::as.data.table(d)
            dt_ann <- data.table::as.data.table(ann)
            data.table::setkey(dt_ann, CHR, START, END)
            overlaps <- data.table::foverlaps(dt_d, dt_ann, by.x = c("CHR", "SNP_Start", "SNP_End"), type = "any", nomatch = NULL)
            if (nrow(overlaps) > 0) {
              overlaps <- overlaps[!is.na(overlaps$Gene), ]
              # Take first overlapping gene per SNP
              overlaps <- overlaps[!duplicated(overlaps$SNP), ]
              d <- merge(d, overlaps[, c("SNP", "Gene"), with = FALSE], by = "SNP", all.x = TRUE, sort = FALSE)
            }
          } else {
            d$Gene <- NA_character_
            for (i in seq_len(nrow(d))) {
              c_chr <- as.character(d$CHR[i])
              sub_ann <- ann[ann$CHR == c_chr, , drop = FALSE]
              if (nrow(sub_ann) > 0) {
                ov <- which(d$SNP_Start[i] <= sub_ann$END & d$SNP_End[i] >= sub_ann$START)
                if (length(ov) > 0) d$Gene[i] <- sub_ann$Gene[ov[1]]
              }
            }
          }
          d$SNP_Start <- NULL
          d$SNP_End <- NULL
        }
      } else {
        ann <- tryCatch(
          {
            if (ext %in% c("csv")) utils::read.csv(ann, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) else utils::read.table(ann, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
          },
          error = function(e) NULL
        )
        if (!is.null(ann) && is.data.frame(ann) && nrow(ann) > 0) {
          cn <- names(ann)
          cn_norm <- gsub("[^a-z0-9]+", "_", tolower(cn))
          snp_col <- which(cn_norm %in% c("snp", "rs", "rsid", "marker", "id"))
          gene_col <- which(cn_norm %in% c("gene", "geneid", "gene_name", "genename", "symbol"))
          if (length(snp_col) >= 1) {
            snp_col <- snp_col[1]
            gene_col <- if (length(gene_col) >= 1) gene_col[1] else integer(0)
            ann2 <- data.frame(SNP = as.character(ann[[snp_col]]), stringsAsFactors = FALSE)
            if (length(gene_col) >= 1) ann2$Gene <- as.character(ann[[gene_col]]) else ann2$Gene <- NA_character_
            ann2 <- ann2[!is.na(ann2$SNP) & ann2$SNP != "", , drop = FALSE]
            if (nrow(ann2) > 0) {
              d <- merge(d, ann2, by = "SNP", all.x = TRUE, sort = FALSE)
            }
          }
        }
      }
    }
  }
  d <- d[is.finite(d$BP) & is.finite(d$P), , drop = FALSE]
  d <- d[d$P > 0 & d$P <= 1, , drop = FALSE]
  if (nrow(d) < 1) return(NULL)
  d$LOG10P <- -log10(d$P)
  num_marker <- nrow(d)
  bonf <- -log10(cut_off / num_marker)
  d_sig <- d[d$LOG10P >= bonf, , drop = FALSE]
  if (nrow(d_sig) < 1) return(NULL)
  if (!("SNP" %in% names(d_sig))) d_sig$SNP <- NA_character_
  d_sig$CHR <- as.character(d_sig$CHR)
  d_sig <- d_sig[order(d_sig$CHR, -d_sig$LOG10P, d_sig$BP), , drop = FALSE]
  if (any(!is.na(d_sig$SNP))) {
    key <- paste(d_sig$CHR, d_sig$SNP, sep = "\t")
    keep <- !duplicated(key)
    d_sig <- d_sig[keep, , drop = FALSE]
  }
  if (isTRUE(dedup_by_position)) {
    bin <- suppressWarnings(as.numeric(dedup_pos_bin))
    if (!is.finite(bin) || bin <= 0) bin <- 1
    key2 <- paste(d_sig$CHR, floor(d_sig$BP / bin), sep = "\t")
    keep2 <- !duplicated(key2)
    d_sig <- d_sig[keep2, , drop = FALSE]
  }
  max_per_chr2 <- suppressWarnings(as.numeric(max_per_chr))
  if (!is.null(max_per_chr) && is.finite(max_per_chr2) && max_per_chr2 > 0) {
    keep2 <- unlist(tapply(seq_len(nrow(d_sig)), d_sig$CHR, function(idx) head(idx, as.integer(max_per_chr2))))
    keep2 <- keep2[is.finite(keep2)]
    d_sig <- d_sig[keep2, , drop = FALSE]
  }
  d_sig$Label <- NA_character_
  for (k in label_priority) {
    if (k %in% names(d_sig)) {
      empty <- is.na(d_sig$Label) | d_sig$Label == ""
      d_sig$Label[empty] <- as.character(d_sig[[k]])[empty]
    } else if (k == "SNP") {
      empty <- is.na(d_sig$Label) | d_sig$Label == ""
      d_sig$Label[empty] <- as.character(d_sig$SNP)[empty]
    }
  }
  d_sig$Label[is.na(d_sig$Label) | d_sig$Label == ""] <- as.character(d_sig$SNP)[is.na(d_sig$Label) | d_sig$Label == ""]
  d_sig$Bonferroni <- bonf
  d_sig
}

GPP.GWAS.Manhattan.Stacked <- function(dat_list,
                                      main = NULL,
                                      cut_off = 0.05,
                                      dpp = 50000,
                                      pch = 1,
                                      cex_points = 0.5,
                                      plot_style = c("Oceanic", "Gray", "PLINK"),
                                      col_alt = NULL,
                                      alpha = 0.9,
                                      point_lwd = 0.9,
                                      sig_cex = NULL,
                                      sig_lwd = 1.6,
                                      show_fdr = TRUE,
                                      line_col = "forestgreen",
                                      label_side = c("right", "none"),
                                      label_cex = 1.3,
                                      annotate_sig = TRUE,
                                      max_label_per_chr = Inf,
                                      annot = NULL,
                                      label_priority = c("Gene", "SNP"),
                                      label_cex_sig = 1.2,
                                      label_font = 3,
                                      arrow_col = "grey70",
                                      arrow_lwd = 1.2,
                                      arrow_len = 0.08,
                                      arrow_style = c("arrow", "segment"),
                                      label_offset = 0.85,
                                      min_arrow_len = 1.05,
                                      prefer_diagonal = TRUE,
                                      label_rows = 2,
                                      label_row_gap = 0.35,
                                      label_pad_x = 0.0,
                                      label_top_gap = 0.15,
                                      label_dx_frac = 0.018,
                                      label_direction = c("left", "right", "alternate"),
                                      label_placement = c("around", "near_up"),
                                      marker_line = FALSE,
                                      marker_line_col = "red",
                                      marker_line_lwd = 0.8,
                                      marker_line_alpha = 0.35,
                                      highlight_label_points = TRUE,
                                      highlight_label_col = "red3",
                                      highlight_label_pch = 19,
                                      y_expand = 1.1,
                                      out_sig = NULL,
                                      cex_axis = 1.35,
                                      cex_lab = 1.5,
                                      cex_main = 1.45,
                                      axis_tick_len = -0.01,
                                      title_line = 0.2,
                                      show_main = TRUE,
                                      mar_left = 6.5,
                                      mar_right = 6.5,
                                      mar_top = 0.45,
                                      mar_bottom = 0.4,
                                      mar_bottom_last = 2.4,
                                      draw_chr_lines = FALSE,
                                      repel_iter = 60,
                                      repel_step = 0.22,
                                      window_kb = 10,
                                      dedup_pos_bin = 1,
                                      cluster_x_frac = 0.01) {
  plot_style <- match.arg(plot_style)
  label_side <- match.arg(label_side)
  label_priority <- match.arg(label_priority, several.ok = TRUE)
  label_direction <- match.arg(label_direction)
  arrow_style <- match.arg(arrow_style)
  label_placement <- match.arg(label_placement)
  if (is.null(dat_list) || length(dat_list) < 1) stop("dat_list is empty.")
  dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
  if (length(dat_list) < 1) stop("dat_list is empty.")
  order_chr <- function(chr_vec) {
    u <- unique(as.character(chr_vec))
    u_num <- suppressWarnings(as.numeric(u))
    is_num <- is.finite(u_num)
    c(u[is_num][order(u_num[is_num])], u[!is_num])
  }
  all_chr <- unlist(lapply(dat_list, function(d) as.character(d$CHR)))
  chr_levels <- order_chr(all_chr)
  chr_max <- setNames(rep(0, length(chr_levels)), chr_levels)
  for (d in dat_list) {
    d2 <- d
    d2$CHR <- factor(as.character(d2$CHR), levels = chr_levels)
    d2 <- d2[is.finite(d2$BP) & is.finite(d2$P) & d2$P > 0 & d2$P <= 1 & !is.na(d2$CHR), , drop = FALSE]
    if (nrow(d2) < 1) next
    mx <- tapply(d2$BP, as.character(d2$CHR), max, na.rm = TRUE)
    for (k in names(mx)) chr_max[k] <- max(chr_max[k], mx[[k]])
  }
  chr_start <- setNames(rep(0, length(chr_levels)), chr_levels)
  ticks <- rep(NA_real_, length(chr_levels))
  chr_end <- rep(NA_real_, length(chr_levels))
  lastbase <- 0
  for (i in seq_along(chr_levels)) {
    chr_start[chr_levels[i]] <- lastbase
    ticks[i] <- lastbase + chr_max[chr_levels[i]] / 2
    lastbase <- lastbase + chr_max[chr_levels[i]]
    chr_end[i] <- lastbase
  }
  xlim <- c(0, lastbase)
  if (is.null(col_alt)) {
    if (plot_style == "Oceanic") {
      col_alt <- c("#EC5f67", "#FAC863", "#99C794", "#6699CC", "#C594C5")
    } else {
      col_alt <- c("gray10", "gray70")
    }
  }
  chr_cols <- rep(col_alt, length.out = length(chr_levels))
  if (is.null(sig_cex)) sig_cex <- cex_points
  arrow_tail_from_box <- function(x0, y0, xlab, ylab, box, h0) {
    dx <- x0 - xlab
    dy <- y0 - ylab
    if (!is.finite(dx) || !is.finite(dy) || (dx == 0 && dy == 0)) return(c(xlab, ylab))
    tt <- numeric(0)
    if (dx > 0) tt <- c(tt, (box[2] - xlab) / dx) else if (dx < 0) tt <- c(tt, (box[1] - xlab) / dx)
    if (dy > 0) tt <- c(tt, (box[4] - ylab) / dy) else if (dy < 0) tt <- c(tt, (box[3] - ylab) / dy)
    tt <- tt[is.finite(tt) & tt > 0]
    if (length(tt) < 1) return(c(xlab, ylab))
    t0 <- min(tt)
    xt <- xlab + t0 * dx
    yt <- ylab + t0 * dy
    len <- sqrt(dx ^ 2 + dy ^ 2)
    eps <- 0.22 * if (is.finite(h0) && h0 > 0) h0 else 0.1
    if (is.finite(len) && len > 0 && is.finite(eps) && eps > 0) {
      xt <- xt + eps * dx / len
      yt <- yt + eps * dy / len
    }
    c(xt, yt)
  }
  label_layout <- function(x0, labels, xlim, cex, rows = 2, pad = 0, dx = 1, direction = "alternate", iter = 60) {
    if (length(x0) < 1) return(list(x = numeric(0), row = integer(0), w = numeric(0), h = numeric(0)))
    rows <- suppressWarnings(as.integer(rows))
    if (!is.finite(rows) || rows < 1) rows <- 1L
    w <- graphics::strwidth(labels, units = "user", cex = cex)
    h <- graphics::strheight(labels, units = "user", cex = cex)
    row_id <- ((seq_along(x0) - 1) %% rows) + 1L
    dir_vec <- rep(1, length(x0))
    if (direction == "left") dir_vec <- rep(-1, length(x0))
    if (direction == "right") dir_vec <- rep(1, length(x0))
    if (direction == "alternate") dir_vec <- ifelse((seq_along(x0) %% 2) == 0, 1, -1)
    x <- x0 + dir_vec * dx
    x <- pmin(pmax(x, xlim[1] + w / 2), xlim[2] - w / 2)
    for (r in seq_len(rows)) {
      idx <- which(row_id == r)
      if (length(idx) < 2) next
      o <- idx[order(x[idx])]
      for (k in 2:length(o)) {
        i <- o[k]
        j <- o[k - 1]
        need <- (w[j] / 2 + w[i] / 2 + pad)
        if (x[i] - x[j] < need) x[i] <- x[j] + need
      }
      o2 <- rev(o)
      for (k in 2:length(o2)) {
        i <- o2[k]
        j <- o2[k - 1]
        need <- (w[i] / 2 + w[j] / 2 + pad)
        if (x[j] - x[i] < need) x[i] <- x[j] - need
      }
      x[idx] <- pmin(pmax(x[idx], xlim[1] + w[idx] / 2), xlim[2] - w[idx] / 2)
      o3 <- idx[order(x[idx])]
      for (k in 2:length(o3)) {
        i <- o3[k]
        j <- o3[k - 1]
        need <- (w[j] / 2 + w[i] / 2 + pad)
        if (x[i] - x[j] < need) x[i] <- x[j] + need
      }
      x[idx] <- pmin(pmax(x[idx], xlim[1] + w[idx] / 2), xlim[2] - w[idx] / 2)
    }
    list(x = x, row = row_id, w = w, h = h)
  }
  prune_points <- function(values, dpp) {
    if (length(values) <= dpp) return(seq_along(values))
    v <- sqrt(values)
    rv <- stats::runif(length(v))
    v <- v + rv
    v <- v[order(v, decreasing = TRUE)]
    the_min <- min(v)
    the_max <- max(v)
    rng <- the_max - the_min
    if (!is.finite(rng) || rng <= 0) return(seq_len(min(length(values), dpp)))
    interval <- rng / dpp
    ladder <- round(v / interval)
    ladder2 <- c(ladder[-1], 0)
    keep <- ladder - ladder2
    which(keep > 0)
  }
  make_title <- function(d) {
    meta <- unique(d[, c("Software", "Model", "Trait"), drop = FALSE])
    s0 <- if ("Software" %in% names(meta)) meta$Software[1] else NA_character_
    m0 <- if ("Model" %in% names(meta)) meta$Model[1] else NA_character_
    t0 <- if ("Trait" %in% names(meta)) meta$Trait[1] else NA_character_
    paste(na.omit(c(s0, m0, t0)), collapse = ".")
  }
  n_panel <- length(dat_list)
  graphics::par(mfrow = c(n_panel, 1), oma = c(0, 0, 2.2, 0))
  sig_all <- list()
  for (i in seq_along(dat_list)) {
    if (i == n_panel) {
      graphics::par(mar = c(mar_bottom_last, mar_left, mar_top, mar_right))
    } else {
      graphics::par(mar = c(mar_bottom, mar_left, mar_top, mar_right))
    }
    d <- dat_list[[i]]
    d$CHR <- factor(as.character(d$CHR), levels = chr_levels)
    d <- d[is.finite(d$BP) & is.finite(d$P) & d$P > 0 & d$P <= 1 & !is.na(d$CHR), , drop = FALSE]
    if (nrow(d) < 1) next
    p0 <- d$P
    p0[p0 <= 0] <- min(p0[p0 > 0], na.rm = TRUE)
    p0[p0 <= 0] <- 1e-300
    y <- -log10(p0)
    o <- order(y, decreasing = TRUE)
    idx0 <- prune_points(y[o], dpp = dpp)
    idx <- o[idx0]
    d2 <- d[idx, , drop = FALSE]
    y2 <- y[idx]
    x <- d2$BP + chr_start[as.character(d2$CHR)]
    chr_id <- as.integer(d2$CHR)
    col_raw <- chr_cols[chr_id]
    title0 <- make_title(d)
    if (i == 1 && is.null(main)) {
      if (isTRUE(show_main)) {
        t0 <- unique(unlist(lapply(dat_list, function(dd) unique(dd$Trait))))
        t0 <- t0[!is.na(t0) & t0 != ""]
        main <- if (length(t0) >= 1) paste0("Trait: ", t0[1]) else ""
      } else {
        main <- ""
      }
    }
    y_max <- max(y2, finite = TRUE)
    y_max <- max(1, y_max)
    d_sig0 <- NULL
    if (annotate_sig) d_sig0 <- GPP.GWAS.Significant.Table(d,
                                                           cut_off = cut_off,
                                                           max_per_chr = max_label_per_chr,
                                                           label_priority = label_priority,
                                                           annot = annot,
                                                           window_kb = window_kb,
                                                           dedup_by_position = TRUE,
                                                           dedup_pos_bin = dedup_pos_bin)

    extra_h <- if (!is.null(d_sig0) && nrow(d_sig0) > 0) {
      label_offset + min(10, suppressWarnings(as.integer(repel_iter))) * repel_step + 0.8
    } else {
      0
    }
    y_lim <- y_max * y_expand + extra_h
    if (y_lim > y_max + 5) y_lim <- y_max + 5
    if (y_lim < y_max + 2) y_lim <- y_max + 2
    graphics::plot.default(NA,
                           xlim = xlim,
                           ylim = c(0, y_lim),
                           xlab = "",
                           ylab = expression(-log[10](italic(p))),
                           axes = FALSE,
                           main = "",
                           cex.axis = cex_axis,
                           cex.lab = cex_lab,
                           xaxs = "i",
                           yaxs = "i")
    usr <- graphics::par("usr")
    pin <- graphics::par("pin")
    cm_in <- 0.5 / 2.54
    dx_1cm <- if (is.finite(pin[1]) && pin[1] > 0) diff(usr[1:2]) / pin[1] * cm_in else NA_real_
    dy_1cm <- if (is.finite(pin[2]) && pin[2] > 0) diff(usr[3:4]) / pin[2] * cm_in else NA_real_
    if (isTRUE(draw_chr_lines)) {
      for (b in chr_end) {
        if (is.finite(b)) graphics::abline(v = b, col = grDevices::adjustcolor("grey70", alpha.f = 0.6), lwd = 1)
      }
    }
    col_point <- grDevices::adjustcolor(col_raw, alpha.f = alpha)
    graphics::points(x, y2, pch = pch, cex = cex_points, col = col_point, lwd = point_lwd)
    num_marker <- nrow(d)
    bonf <- -log10(cut_off / num_marker)
    graphics::abline(h = bonf, col = line_col, lwd = 2)
    if (show_fdr) {
      p_sorted <- sort(d$P[is.finite(d$P) & d$P > 0 & d$P <= 1])
      spd <- abs(cut_off - p_sorted * num_marker / cut_off)
      spd <- spd[is.finite(spd)]
      index_fdr <- if (length(spd) >= 1) which.min(spd) else NA_integer_
      fdr_cut <- if (is.finite(index_fdr)) -log10(cut_off * index_fdr / num_marker) else NA_real_
      if (is.finite(fdr_cut)) graphics::abline(h = fdr_cut, col = line_col, lwd = 2, lty = 2)
    }
    sig_idx <- which(y2 >= bonf)
    if (length(sig_idx) >= 1) {
      graphics::points(x[sig_idx], y2[sig_idx],
                       pch = pch,
                       cex = sig_cex,
                       col = grDevices::adjustcolor(col_raw[sig_idx], alpha.f = 1),
                       lwd = sig_lwd)
    }
    if (i == n_panel) {
      graphics::axis(1, at = ticks, labels = chr_levels, las = 1, cex.axis = cex_axis, tck = axis_tick_len)
    } else {
      graphics::axis(1, at = ticks, labels = FALSE, tick = FALSE, tck = axis_tick_len)
    }
    y_tick_min <- 0
    y_rng <- y_lim - y_tick_min
    y_step <- if (y_rng > 20) 5 else if (y_rng > 10) 2 else 1
    y0 <- ceiling(y_tick_min / y_step) * y_step
    y1 <- floor(y_lim / y_step) * y_step
    y_ticks <- if (is.finite(y0) && is.finite(y1) && y1 >= y0) seq(y0, y1, by = y_step) else numeric(0)
    graphics::axis(2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = cex_axis, tck = axis_tick_len)
    if (label_side == "right") graphics::mtext(title0, side = 4, line = 0.6, cex = label_cex)
    if (annotate_sig) {
      d_sig <- d_sig0
      if (!is.null(d_sig) && nrow(d_sig) > 0) {
        d_sig$Panel <- title0
        d_sig$X <- d_sig$BP + chr_start[d_sig$CHR]
        d_sig$Y <- d_sig$LOG10P
        d_sig <- d_sig[order(-d_sig$Y, d_sig$X), , drop = FALSE]
        
        # --- Visual Deduplication (Keep highest SNP in same vertical column) ---
        keep_idx <- rep(TRUE, nrow(d_sig))
        x_threshold <- diff(xlim) * suppressWarnings(as.numeric(cluster_x_frac))
        if (!is.finite(x_threshold) || x_threshold <= 0) x_threshold <- diff(xlim) * 0.01
        
        for (j in seq_len(nrow(d_sig))) {
          if (!keep_idx[j]) next
          close_idx <- which(abs(d_sig$X - d_sig$X[j]) < x_threshold & d_sig$CHR == d_sig$CHR[j])
          close_idx <- close_idx[close_idx > j]
          if (length(close_idx) > 0) {
            keep_idx[close_idx] <- FALSE
          }
        }
        d_sig <- d_sig[keep_idx, , drop = FALSE]
        # -----------------------------------------------------------------------
        
        labs <- as.character(d_sig$Label)
        labs[is.na(labs) | labs == ""] <- as.character(d_sig$SNP)[is.na(labs) | labs == ""]
        dx0 <- diff(xlim) * suppressWarnings(as.numeric(label_dx_frac))
        if (!is.finite(dx0) || dx0 <= 0) dx0 <- diff(xlim) * 0.015
        local_xpd <- graphics::par("xpd")
        graphics::par(xpd = NA)
        placed <- list()
        for (j in seq_len(nrow(d_sig))) {
          x0 <- d_sig$X[j]
          y0 <- d_sig$Y[j]
          lab <- labs[j]
          w0 <- graphics::strwidth(lab, units = "user", cex = label_cex_sig)
          h0 <- graphics::strheight(lab, units = "user", cex = label_cex_sig)
          y_base <- y0 + max(label_offset, min_arrow_len)
          if (!is.finite(y_base)) y_base <- y0
          min_dy <- max(0.6, 0.55 * max(label_offset, min_arrow_len))
          dir0 <- 1
          if (label_direction == "left") dir0 <- -1
          if (label_direction == "right") dir0 <- 1
          if (label_direction == "alternate") dir0 <- if ((j %% 2) == 0) 1 else -1
          best <- NULL
          base_off <- c(dx0)
          off_right <- base_off
          off_left <- -base_off
          if (dir0 < 0) {
            pref <- off_left
            other <- off_right
          } else {
            pref <- off_right
            other <- off_left
          }
          
          # 优先尝试：若附近（约 0.5cm）有空位，优先在附近标注；否则再逐步拉长箭头避开重叠
          if (is.null(best) && is.finite(dx_1cm) && is.finite(dy_1cm) && dx_1cm > 0 && dy_1cm > 0) {
            crowd <- 0L
            if (length(placed) >= 1) {
              x_low <- x0 - 2 * dx_1cm
              x_high <- x0 + 2 * dx_1cm
              for (bb in placed) {
                if (!(bb[2] < x_low || bb[1] > x_high)) crowd <- crowd + 1L
              }
            }
            
            # 如果当前区域已经拥挤（附近已有 >=2 个标签），优先把新标签放到右侧 0.5cm 附近；
            # 若还是放不下，再逐步“延伸”距离（箭头变长），直到找到空位
            if (crowd >= 2L) {
              dx_r <- abs(dx_1cm)
              dy_r <- max(dy_1cm, min_dy)
              cand_lane <- list(
                c(+1.0 * dx_r, 1.0 * dy_r),
                c(+1.0 * dx_r, 1.8 * dy_r),
                c(+1.8 * dx_r, 1.0 * dy_r),
                c(+1.8 * dx_r, 1.8 * dy_r),
                c(+3.0 * dx_r, 2.2 * dy_r),
                c(+4.2 * dx_r, 2.8 * dy_r)
              )
              for (xy in cand_lane) {
                x_try <- x0 + xy[1]
                y_try <- y0 + xy[2]
                x_try <- min(max(x_try, xlim[1] + w0 / 2), xlim[2] - w0 / 2)
                if (y_try > y_lim - h0) y_try <- y_lim - h0
                if (!is.finite(y_try) || y_try < y0 + min_dy) next
                if (y_try < h0 * 0.6) next
                box_try <- c(x_try - w0 / 2, x_try + w0 / 2, y_try - h0 / 2, y_try + h0 / 2)
                overlap <- FALSE
                if (length(placed) >= 1) {
                  for (bb in placed) {
                    if (!(box_try[2] < bb[1] || box_try[1] > bb[2] || box_try[4] < bb[3] || box_try[3] > bb[4])) {
                      overlap <- TRUE
                      break
                    }
                  }
                }
                if (!overlap) {
                  best <- list(x = x_try, y = y_try, w = w0, h = h0, box = box_try)
                  break
                }
              }
            }
          }
          
          if (is.null(best) && is.finite(dx_1cm) && is.finite(dy_1cm) && dx_1cm > 0 && dy_1cm > 0) {
            x_near <- if (dir0 < 0) -abs(dx_1cm) else abs(dx_1cm)
            x_near_other <- -x_near
            y_near <- max(dy_1cm, min_dy)
            cand_near <- list(
              c(x_near, y_near),
              c(x_near, 1.2 * y_near),
              c(x_near_other, y_near),
              c(x_near_other, 1.2 * y_near)
            )
            for (xy in cand_near) {
              x_try <- x0 + xy[1]
              y_try <- y0 + xy[2]
              x_try <- min(max(x_try, xlim[1] + w0 / 2), xlim[2] - w0 / 2)
              if (y_try > y_lim - h0) y_try <- y_lim - h0
              if (!is.finite(y_try) || y_try < y0 + min_dy) next
              if (y_try < h0 * 0.6) next
              box_try <- c(x_try - w0 / 2, x_try + w0 / 2, y_try - h0 / 2, y_try + h0 / 2)
              overlap <- FALSE
              if (length(placed) >= 1) {
                for (bb in placed) {
                  if (!(box_try[2] < bb[1] || box_try[1] > bb[2] || box_try[4] < bb[3] || box_try[3] > bb[4])) {
                    overlap <- TRUE
                    break
                  }
                }
              }
              if (!overlap) {
                best <- list(x = x_try, y = y_try, w = w0, h = h0, box = box_try)
                break
              }
            }
          }
          
          for (k in seq_len(max(1, suppressWarnings(as.integer(repel_iter))))) {
            kk <- k - 1
            if (label_placement == "near_up") {
              y_try0 <- y_base + kk * repel_step
              if (y_try0 > y_lim - h0) y_try0 <- y_lim - h0
              cand <- list(c(pref[1], y_try0 - y0), c(other[1], y_try0 - y0))
            } else {
              dy0 <- max(label_offset, h0 * 0.9, min_arrow_len)
              aspect <- w0 / max(abs(dx0), .Machine$double.eps)
              if (is.finite(aspect) && aspect > 3) dy0 <- dy0 * min(2.2, aspect / 3)
              rscale <- 1 + kk * 0.35
              dyA <- dy0 * rscale
              cand <- list(
                c(pref[1], dyA),
                c(pref[1], 0.65 * dyA),
                c(other[1], dyA),
                c(other[1], 0.65 * dyA)
              )
              if (j == 1) cand[[length(cand) + 1]] <- c(other[1], -0.55 * dy0)
            }
            if (!is.null(best)) break
            for (xy in cand) {
              x_try <- x0 + xy[1]
              y_try <- y0 + xy[2]
              if (!is.finite(y_try)) next
              if (xy[2] >= 0 && y_try < y0 + min_dy) next
              x_try <- min(max(x_try, xlim[1] + w0 / 2), xlim[2] - w0 / 2)
              if (y_try > y_lim - h0) y_try <- y_lim - h0
              if (!is.finite(y_try)) next
              if (xy[2] >= 0 && y_try < y0 + min_dy) next
              if (y_try < h0 * 0.6) next
              box_try <- c(x_try - w0 / 2, x_try + w0 / 2, y_try - h0 / 2, y_try + h0 / 2)
              overlap <- FALSE
              if (length(placed) >= 1) {
                for (bb in placed) {
                  if (!(box_try[2] < bb[1] || box_try[1] > bb[2] || box_try[4] < bb[3] || box_try[3] > bb[4])) {
                    overlap <- TRUE
                    break
                  }
                }
              }
              if (!overlap) {
                best <- list(x = x_try, y = y_try, w = w0, h = h0, box = box_try)
                break
              }
            }
            if (!is.null(best)) break
          }
          if (is.null(best)) {
            x_try <- min(max(x0 + dir0 * dx0, xlim[1] + w0 / 2), xlim[2] - w0 / 2)
            y_try <- min(y_lim - h0, y_base)
            best <- list(x = x_try, y = y_try, w = w0, h = h0, box = c(x_try - w0 / 2, x_try + w0 / 2, y_try - h0 / 2, y_try + h0 / 2))
          }
          placed[[length(placed) + 1]] <- best$box
          d_sig$X_text[j] <- best$x
          d_sig$Y_text[j] <- best$y
          if (isTRUE(highlight_label_points)) {
            graphics::points(x0, y0, pch = highlight_label_pch, cex = cex_points, col = highlight_label_col)
          }
          tail_xy <- arrow_tail_from_box(x0, y0, best$x, best$y, best$box, best$h)
          if (arrow_style == "segment") {
            graphics::segments(tail_xy[1], tail_xy[2], x0, y0, col = arrow_col, lwd = arrow_lwd)
          } else {
            graphics::arrows(tail_xy[1], tail_xy[2], x0, y0, length = arrow_len, col = arrow_col, lwd = arrow_lwd)
          }
          graphics::text(best$x, best$y, labels = lab, cex = label_cex_sig, font = label_font, adj = c(0.5, 0.5))
        }
        graphics::par(xpd = local_xpd)
        sig_all[[length(sig_all) + 1]] <- d_sig
      }
    }
  }
  if (isTRUE(show_main) && !is.null(main) && !is.na(main) && main != "") {
    graphics::mtext(main, side = 3, outer = TRUE, line = title_line, cex = cex_main)
  }
  sig_df <- if (length(sig_all) >= 1) do.call(rbind, sig_all) else NULL
  if (!is.null(out_sig)) {
    out_sig <- normalizePath(as.character(out_sig), winslash = "\\", mustWork = FALSE)
    ext <- tolower(tools::file_ext(out_sig))
    if (ext %in% c("xlsx")) {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Significant")
        openxlsx::writeData(wb, "Significant", sig_df)
        openxlsx::saveWorkbook(wb, out_sig, overwrite = TRUE)
      } else {
        utils::write.csv(sig_df, sub("\\.xlsx$", ".csv", out_sig, ignore.case = TRUE), row.names = FALSE, quote = TRUE)
      }
    } else {
      utils::write.csv(sig_df, out_sig, row.names = FALSE, quote = TRUE)
    }
  }
  invisible(list(sig = sig_df, chr_levels = chr_levels, ticks = ticks))
}

GPP.GWAS.Manhattan.Horizontal <- function(dat,
                                         main = NULL,
                                         cut_off = 0.05,
                                         dpp = 50000,
                                         pch = 21,
                                         cex_points = 0.35,
                                         plot_style = c("Oceanic", "Gray", "PLINK"),
                                         col_alt = NULL,
                                         alpha = 0.35,
                                         border_alpha = 0.9,
                                         highlight = TRUE,
                                         cex_sig = 0.65,
                                         line_col = "forestgreen",
                                         show_fdr = TRUE,
                                         label_side = c("right", "top", "none"),
                                         label_cex = 0.95,
                                         show_x_axis = TRUE) {
  plot_style <- match.arg(plot_style)
  label_side <- match.arg(label_side)
  if (is.null(dat) || nrow(dat) < 1) stop("dat is empty.")
  dat <- dat[is.finite(dat$BP) & is.finite(dat$P), , drop = FALSE]
  dat <- dat[dat$P > 0 & dat$P <= 1, , drop = FALSE]
  if (nrow(dat) < 1) stop("No valid P values in dat.")
  order_chr <- function(chr_vec) {
    u <- unique(as.character(chr_vec))
    u_num <- suppressWarnings(as.numeric(u))
    is_num <- is.finite(u_num)
    c(u[is_num][order(u_num[is_num])], u[!is_num])
  }
  chr_levels <- order_chr(dat$CHR)
  dat$CHR <- factor(as.character(dat$CHR), levels = chr_levels)
  dat <- dat[order(dat$CHR, dat$BP), , drop = FALSE]
  n <- nrow(dat)
  cumpos <- numeric(n)
  ticks <- numeric(length(chr_levels))
  chr_end <- numeric(length(chr_levels))
  lastbase <- 0
  for (i in seq_along(chr_levels)) {
    idx <- which(dat$CHR == chr_levels[i])
    if (length(idx) < 1) next
    bp0 <- dat$BP[idx]
    cumpos[idx] <- bp0 + lastbase
    ticks[i] <- lastbase + mean(bp0, na.rm = TRUE)
    lastbase <- lastbase + max(bp0, na.rm = TRUE)
    chr_end[i] <- lastbase
  }
  p0 <- dat$P
  p0[p0 <= 0] <- min(p0[p0 > 0], na.rm = TRUE)
  p0[p0 <= 0] <- 1e-300
  logp <- -log10(p0)
  prune_points <- function(values, dpp) {
    if (length(values) <= dpp) return(seq_along(values))
    v <- sqrt(values)
    rv <- stats::runif(length(v))
    v <- v + rv
    v <- v[order(v, decreasing = TRUE)]
    the_min <- min(v)
    the_max <- max(v)
    rng <- the_max - the_min
    if (!is.finite(rng) || rng <= 0) return(seq_len(min(length(values), dpp)))
    interval <- rng / dpp
    ladder <- round(v / interval)
    ladder2 <- c(ladder[-1], 0)
    keep <- ladder - ladder2
    which(keep > 0)
  }
  o <- order(logp, decreasing = TRUE)
  idx_keep0 <- prune_points(logp[o], dpp = dpp)
  idx_keep <- o[idx_keep0]
  logp <- logp[idx_keep]
  cumpos <- cumpos[idx_keep]
  chr_id <- as.integer(dat$CHR[idx_keep])
  num_marker <- sum(is.finite(dat$P) & dat$P > 0 & dat$P <= 1)
  bonf <- -log10(cut_off / num_marker)
  p_sorted <- sort(dat$P[is.finite(dat$P) & dat$P > 0 & dat$P <= 1])
  spd <- abs(cut_off - p_sorted * num_marker / cut_off)
  spd <- spd[is.finite(spd)]
  index_fdr <- if (length(spd) >= 1) which.min(spd) else NA_integer_
  fdr_cut <- if (is.finite(index_fdr)) -log10(cut_off * index_fdr / num_marker) else NA_real_
  xlim <- c(0, max(logp, finite = TRUE) * 1.05)
  ylim <- range(cumpos, finite = TRUE)
  if (is.null(col_alt)) {
    if (plot_style == "Oceanic") {
      col_alt <- c("#EC5f67", "#FAC863", "#99C794", "#6699CC", "#C594C5")
    } else {
      col_alt <- c("gray10", "gray70")
    }
  }
  col_raw <- rep(col_alt, length.out = length(chr_levels))[chr_id]
  col_border <- grDevices::adjustcolor(col_raw, alpha.f = border_alpha)
  col_fill <- grDevices::adjustcolor(col_raw, alpha.f = alpha)
  if (is.null(main)) {
    meta <- unique(dat[, c("Software", "Model", "Trait"), drop = FALSE])
    s0 <- if ("Software" %in% names(meta)) meta$Software[1] else NA_character_
    m0 <- if ("Model" %in% names(meta)) meta$Model[1] else NA_character_
    t0 <- if ("Trait" %in% names(meta)) meta$Trait[1] else NA_character_
    main <- paste(na.omit(c(s0, m0, t0)), collapse = " | ")
    if (is.na(main) || main == "") main <- basename(dat$File[1])
  }
  main0 <- if (label_side == "top") main else ""
  xlab0 <- if (show_x_axis) expression(-log[10](italic(p))) else ""
  graphics::plot(cumpos ~ logp,
                 type = "n",
                 xlab = xlab0,
                 ylab = "",
                 axes = FALSE,
                 main = main0,
                 xlim = xlim,
                 ylim = ylim,
                 xaxs = "i",
                 yaxs = "i")
  graphics::points(cumpos ~ logp, pch = pch, cex = cex_points, col = col_border, bg = col_fill)
  if (highlight) {
    sig <- logp >= bonf
    if (any(sig)) {
      graphics::points(cumpos[sig] ~ logp[sig],
                       pch = pch,
                       cex = cex_sig,
                       col = grDevices::adjustcolor(col_raw[sig], alpha.f = 1),
                       bg = grDevices::adjustcolor(col_raw[sig], alpha.f = min(1, alpha * 2.2)))
    }
  }
  if (show_x_axis) graphics::axis(1, las = 1)
  graphics::axis(2, at = ticks, labels = chr_levels, las = 1, tick = FALSE)
  for (i in seq_along(chr_end)) {
    if (is.finite(chr_end[i]) && chr_end[i] > 0) graphics::abline(h = chr_end[i], col = grDevices::adjustcolor("grey70", alpha.f = 0.6), lwd = 1)
  }
  graphics::abline(v = bonf, col = line_col, lwd = 2)
  if (show_fdr && is.finite(fdr_cut)) graphics::abline(v = fdr_cut, col = line_col, lwd = 2, lty = 2)
  if (label_side == "right") graphics::mtext(main, side = 4, line = 0.6, cex = label_cex)
  graphics::box()
  invisible(list(bonferroni = bonf, fdr = fdr_cut, markers = num_marker))
}

GPP.GWAS.Manhattan.Layered <- function(dat_list,
                                       main = NULL,
                                       cut_off = 0.05,
                                       dpp = 50000,
                                       pch = 21,
                                       cex_points = 0.35,
                                       plot_style = c("Oceanic", "Gray", "PLINK"),
                                       col_alt = NULL,
                                       alpha = 0.35,
                                       border_alpha = 0.9,
                                       track_gap = 1.5,
                                       show_fdr = TRUE,
                                       line_col = "forestgreen",
                                       label_side = c("right", "none"),
                                       label_cex = 0.95,
                                       axis_cex = 0.95,
                                       axis_tick_len = -0.015) {
  plot_style <- match.arg(plot_style)
  label_side <- match.arg(label_side)
  if (is.null(dat_list) || length(dat_list) < 1) stop("dat_list is empty.")
  dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
  if (length(dat_list) < 1) stop("dat_list is empty.")
  order_chr <- function(chr_vec) {
    u <- unique(as.character(chr_vec))
    u_num <- suppressWarnings(as.numeric(u))
    is_num <- is.finite(u_num)
    c(u[is_num][order(u_num[is_num])], u[!is_num])
  }
  all_chr <- unlist(lapply(dat_list, function(d) as.character(d$CHR)))
  chr_levels <- order_chr(all_chr)
  chr_max <- setNames(rep(0, length(chr_levels)), chr_levels)
  for (d in dat_list) {
    d$CHR <- factor(as.character(d$CHR), levels = chr_levels)
    d <- d[is.finite(d$BP) & is.finite(d$P) & d$P > 0 & d$P <= 1 & !is.na(d$CHR), , drop = FALSE]
    if (nrow(d) < 1) next
    mx <- tapply(d$BP, as.character(d$CHR), max, na.rm = TRUE)
    for (k in names(mx)) chr_max[k] <- max(chr_max[k], mx[[k]])
  }
  chr_start <- setNames(rep(0, length(chr_levels)), chr_levels)
  ticks <- rep(NA_real_, length(chr_levels))
  lastbase <- 0
  for (i in seq_along(chr_levels)) {
    chr_start[chr_levels[i]] <- lastbase
    ticks[i] <- lastbase + chr_max[chr_levels[i]] / 2
    lastbase <- lastbase + chr_max[chr_levels[i]]
  }
  xlim <- c(0, lastbase)
  if (is.null(col_alt)) {
    if (plot_style == "Oceanic") {
      col_alt <- c("#EC5f67", "#FAC863", "#99C794", "#6699CC", "#C594C5")
    } else {
      col_alt <- c("gray10", "gray70")
    }
  }
  chr_cols <- rep(col_alt, length.out = length(chr_levels))
  prune_points <- function(values, dpp) {
    if (length(values) <= dpp) return(seq_along(values))
    v <- sqrt(values)
    rv <- stats::runif(length(v))
    v <- v + rv
    v <- v[order(v, decreasing = TRUE)]
    the_min <- min(v)
    the_max <- max(v)
    rng <- the_max - the_min
    if (!is.finite(rng) || rng <= 0) return(seq_len(min(length(values), dpp)))
    interval <- rng / dpp
    ladder <- round(v / interval)
    ladder2 <- c(ladder[-1], 0)
    keep <- ladder - ladder2
    which(keep > 0)
  }
  tracks <- vector("list", length(dat_list))
  track_max <- rep(0, length(dat_list))
  track_titles <- character(length(dat_list))
  track_bonf <- rep(NA_real_, length(dat_list))
  track_fdr <- rep(NA_real_, length(dat_list))
  for (i in seq_along(dat_list)) {
    d <- dat_list[[i]]
    d$CHR <- factor(as.character(d$CHR), levels = chr_levels)
    d <- d[is.finite(d$BP) & is.finite(d$P) & d$P > 0 & d$P <= 1 & !is.na(d$CHR), , drop = FALSE]
    if (nrow(d) < 1) {
      tracks[[i]] <- NULL
      track_max[i] <- 0
      track_titles[i] <- ""
      next
    }
    p0 <- d$P
    p0[p0 <= 0] <- min(p0[p0 > 0], na.rm = TRUE)
    p0[p0 <= 0] <- 1e-300
    logp <- -log10(p0)
    o <- order(logp, decreasing = TRUE)
    idx0 <- prune_points(logp[o], dpp = dpp)
    idx <- o[idx0]
    d2 <- d[idx, , drop = FALSE]
    logp2 <- logp[idx]
    x <- d2$BP + chr_start[as.character(d2$CHR)]
    chr_id <- as.integer(d2$CHR)
    col_raw <- chr_cols[chr_id]
    tracks[[i]] <- list(x = x, y = logp2, col = col_raw, data = d)
    track_max[i] <- max(logp2, finite = TRUE)
    meta <- unique(d[, c("Software", "Model", "Trait"), drop = FALSE])
    s0 <- if ("Software" %in% names(meta)) meta$Software[1] else NA_character_
    m0 <- if ("Model" %in% names(meta)) meta$Model[1] else NA_character_
    t0 <- if ("Trait" %in% names(meta)) meta$Trait[1] else NA_character_
    track_titles[i] <- paste(na.omit(c(s0, m0, t0)), collapse = " | ")
    num_marker <- sum(is.finite(d$P) & d$P > 0 & d$P <= 1)
    track_bonf[i] <- -log10(cut_off / num_marker)
    p_sorted <- sort(d$P[is.finite(d$P) & d$P > 0 & d$P <= 1])
    spd <- abs(cut_off - p_sorted * num_marker / cut_off)
    spd <- spd[is.finite(spd)]
    index_fdr <- if (length(spd) >= 1) which.min(spd) else NA_integer_
    track_fdr[i] <- if (is.finite(index_fdr)) -log10(cut_off * index_fdr / num_marker) else NA_real_
  }
  track_height <- ceiling(pmax(track_max, 1))
  offsets <- numeric(length(track_height))
  for (i in seq_along(track_height)) offsets[i] <- sum(track_height[seq_len(i - 1)], na.rm = TRUE) + track_gap * (i - 1)
  ylim <- c(0, max(offsets + track_height, na.rm = TRUE) + 0.5)
  if (is.null(main)) {
    t0 <- unique(unlist(lapply(dat_list, function(d) unique(d$Trait))))
    t0 <- t0[!is.na(t0) & t0 != ""]
    main <- if (length(t0) >= 1) paste0("Trait: ", t0[1]) else ""
  }
  graphics::plot.default(NA,
                         xlim = xlim,
                         ylim = ylim,
                         xlab = "Chromosome",
                         ylab = expression(-log[10](italic(p))),
                         axes = FALSE,
                         main = main,
                         xaxs = "i",
                         yaxs = "i")
  graphics::axis(1, at = ticks, labels = chr_levels, las = 1, cex.axis = axis_cex)
  for (i in seq_along(track_height)) {
    off <- offsets[i]
    h <- track_height[i]
    graphics::abline(h = off, col = grDevices::adjustcolor("grey70", alpha.f = 0.7), lwd = 1)
    if (!is.null(tracks[[i]])) {
      tt <- tracks[[i]]
      col_border <- grDevices::adjustcolor(tt$col, alpha.f = border_alpha)
      col_fill <- grDevices::adjustcolor(tt$col, alpha.f = alpha)
      graphics::points(tt$x, tt$y + off, pch = pch, cex = cex_points, col = col_border, bg = col_fill)
      graphics::abline(h = off + track_bonf[i], col = line_col, lwd = 2)
      if (show_fdr && is.finite(track_fdr[i])) graphics::abline(h = off + track_fdr[i], col = line_col, lwd = 2, lty = 2)
      yt <- pretty(c(0, h), n = 4)
      yt <- yt[yt >= 0 & yt <= h]
      if (length(yt) >= 1) graphics::axis(2, at = off + yt, labels = yt, las = 1, cex.axis = axis_cex, tck = axis_tick_len)
      if (label_side == "right") graphics::mtext(track_titles[i], side = 4, at = off + h / 2, line = 0.6, cex = label_cex)
    }
  }
  graphics::box()
  invisible(list(chr_levels = chr_levels, ticks = ticks, offsets = offsets, track_height = track_height))
}

GPP.GWAS.Manhattan.Horizontal.MultiSoftware <- function(dir = NULL,
                                                       files = NULL,
                                                       trait = NULL,
                                                       out = NULL,
                                                       width = 14,
                                                       height = 16,
                                                       cut_off = 0.05,
                                                       dpp = 50000,
                                                       cex_points = 0.35,
                                                       plot_style = c("Oceanic", "Gray", "PLINK"),
                                                       alpha = 0.35,
                                                       border_alpha = 0.9,
                                                       highlight = TRUE,
                                                       cex_sig = 0.65,
                                                       label_side = c("right", "top", "none"),
                                                       label_cex = 0.95,
                                                       plot_mode = c("stacked", "layered", "sideways"),
                                                       show_main = TRUE,
                                                       annotate_sig = TRUE,
                                                       max_label_per_chr = Inf,
                                                       annot = NULL,
                                                       label_priority = c("Gene", "SNP"),
                                                       out_sig = NULL,
                                                       label_dx_frac = 0.018,
                                                       label_offset = 0.85,
                                                       min_arrow_len = 1.05,
                                                       prefer_diagonal = TRUE,
                                                       label_rows = 2,
                                                       label_top_gap = 0.15,
                                                       label_placement = c("around", "near_up"),
                                                       window_kb = 10,
                                                       dedup_pos_bin = 1,
                                                       cluster_x_frac = 0.01,
                                                       axis_tick_len = -0.01) {
  plot_style <- match.arg(plot_style)
  label_side <- match.arg(label_side)
  plot_mode <- match.arg(plot_mode)
  label_priority <- match.arg(label_priority, several.ok = TRUE)
  label_placement <- match.arg(label_placement)
  if (is.null(files)) {
    if (is.null(dir)) stop("Either dir or files must be provided.")
    dir <- normalizePath(as.character(dir), winslash = "\\", mustWork = TRUE)
    cand <- list.files(dir, full.names = TRUE)
    cand <- cand[file.info(cand)$isdir == FALSE]
    cand <- cand[grepl("\\.csv$|\\.assoc\\.|\\.txt$|\\.tsv$|\\.linear$", cand, ignore.case = TRUE)]
    if (!is.null(trait)) {
      cand_trait <- cand[grepl(trait, basename(cand), ignore.case = TRUE)]
      if (length(cand_trait) >= 1) cand <- unique(c(cand_trait, cand))
    }
    files <- cand
  }
  files <- unique(normalizePath(as.character(files), winslash = "\\", mustWork = TRUE))
  if (length(files) < 1) stop("No input files.")
  if (is.null(out)) {
    base_out <- if (!is.null(trait)) trait else "Trait"
    base_out <- gsub("[^A-Za-z0-9_\\-]+", "_", base_out)
    out <- file.path(if (!is.null(dir)) dir else dirname(files[1]), paste0(base_out, ".MultiSoftware.Horizontal.Manhattan.pdf"))
  }
  out <- normalizePath(as.character(out), winslash = "\\", mustWork = FALSE)
  dat_list <- lapply(files, function(f) {
    tryCatch(GPP.GWAS.Read.Manhattan(f, trait = trait), error = function(e) NULL)
  })
  dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
  if (length(dat_list) < 1) stop("No readable GWAS files.")
  has_trait <- vapply(dat_list, function(d) any(!is.na(d$Trait) & d$Trait != ""), logical(1))
  if (!is.null(trait)) {
    keep <- vapply(dat_list, function(d) any(tolower(d$Trait) == tolower(trait) | grepl(trait, d$File, ignore.case = TRUE)), logical(1))
    if (any(keep)) dat_list <- dat_list[keep]
  } else if (any(has_trait)) {
    t0 <- unique(unlist(lapply(dat_list[has_trait], function(d) unique(d$Trait))))
    t0 <- t0[!is.na(t0) & t0 != ""]
    if (length(t0) >= 1) trait <- t0[1]
  }
  dat_list <- dat_list[order(vapply(dat_list, function(d) paste0(d$Software[1], "_", d$Model[1]), character(1)))]
  grDevices::pdf(out, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  if (is.null(out_sig)) {
    out_sig <- sub("\\.pdf$", ".Significant.csv", out, ignore.case = TRUE)
  }
  if (plot_mode == "stacked") {
    GPP.GWAS.Manhattan.Stacked(dat_list,
                              cut_off = cut_off,
                              dpp = dpp,
                              cex_points = cex_points,
                              plot_style = plot_style,
                              alpha = alpha,
                              show_fdr = TRUE,
                              show_main = show_main,
                              label_side = if (label_side == "none") "none" else "right",
                              label_cex = label_cex,
                              annotate_sig = annotate_sig,
                              max_label_per_chr = max_label_per_chr,
                              annot = annot,
                              label_priority = label_priority,
                              label_dx_frac = label_dx_frac,
                              label_offset = label_offset,
                              min_arrow_len = min_arrow_len,
                              prefer_diagonal = prefer_diagonal,
                              label_rows = label_rows,
                              label_top_gap = label_top_gap,
                              label_placement = label_placement,
                              window_kb = window_kb,
                              dedup_pos_bin = dedup_pos_bin,
                              cluster_x_frac = cluster_x_frac,
                              axis_tick_len = axis_tick_len,
                              out_sig = out_sig)
  } else if (plot_mode == "layered") {
    graphics::par(mar = c(4.5, 6.5, 3.2, 6.5))
    GPP.GWAS.Manhattan.Layered(dat_list,
                              cut_off = cut_off,
                              dpp = dpp,
                              cex_points = cex_points,
                              plot_style = plot_style,
                              alpha = alpha,
                              border_alpha = border_alpha,
                              show_fdr = TRUE,
                              label_side = if (label_side == "none") "none" else "right",
                              label_cex = label_cex)
  } else {
    n_panel <- length(dat_list)
    graphics::par(mfrow = c(n_panel, 1), mar = c(4.2, 6.5, 3.2, 2.0))
    for (i in seq_along(dat_list)) {
      GPP.GWAS.Manhattan.Horizontal(dat_list[[i]],
                                   cut_off = cut_off,
                                   dpp = dpp,
                                   cex_points = cex_points,
                                   plot_style = plot_style,
                                   alpha = alpha,
                                   border_alpha = border_alpha,
                                   highlight = highlight,
                                   cex_sig = cex_sig,
                                   label_side = label_side,
                                   label_cex = label_cex,
                                   show_x_axis = (i == n_panel))
    }
  }
  invisible(list(out = out, trait = trait, data = dat_list))
}
