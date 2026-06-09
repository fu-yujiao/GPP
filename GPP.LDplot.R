lddecay_read_stat <- function(file) {
  d <- read.table(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  n0 <- names(d)
  n1 <- tolower(gsub("[^a-z0-9]+", "_", n0))
  pick <- function(rx, fallback) {
    i <- which(grepl(rx, n1))
    if (length(i)) i[1] else fallback
  }
  i_dist <- pick("^dist$|^distance$|dist_bp|distbase|distbp", 1)
  i_sumr2 <- pick("^sum_r_2$|^sum_r2$|^sum_r_\\^2$|sum_r", 4)
  i_npairs <- pick("^numberpairs$|^npairs$|^n_pairs$|pairs", 6)
  i_meanr2 <- pick("^mean_r_2$|^mean_r2$|^mean_r_\\^2$|^r2$|rsquare|r_square", 2)
  list(
    dist_bp = suppressWarnings(as.numeric(d[[i_dist]])),
    sum_r2 = suppressWarnings(as.numeric(d[[i_sumr2]])),
    npairs = suppressWarnings(as.numeric(d[[i_npairs]])),
    mean_r2 = suppressWarnings(as.numeric(d[[i_meanr2]]))
  )
}

lddecay_bin <- function(file, xmax_kb = 500, bin_kb = 5, min_pairs = 50, min_bp = 15, r2_source = c("sum", "mean")) {
  r2_source <- match.arg(r2_source)
  x <- lddecay_read_stat(file)
  dist_bp <- x$dist_bp
  dis_kb <- dist_bp / 1000
  if (r2_source == "sum") {
    sum_r2 <- x$sum_r2
    npairs <- x$npairs
    ok <- is.finite(dis_kb) & is.finite(sum_r2) & is.finite(npairs) & npairs > 0
    ok <- ok & is.finite(dist_bp) & dist_bp >= min_bp & dis_kb >= 0 & dis_kb <= xmax_kb
    dis_kb <- dis_kb[ok]
    sum_r2 <- sum_r2[ok]
    npairs <- npairs[ok]
    br <- seq(0, xmax_kb, by = bin_kb)
    if (tail(br, 1) < xmax_kb) br <- c(br, xmax_kb)
    g <- cut(dis_kb, breaks = br, include.lowest = TRUE, right = FALSE)
    sumr2_bin <- tapply(sum_r2, g, sum, na.rm = TRUE)
    n_bin <- tapply(npairs, g, sum, na.rm = TRUE)
    r2_bin <- sumr2_bin / n_bin
    centers <- br[-1] - bin_kb / 2
    df <- data.frame(dis = centers, r2 = as.numeric(r2_bin), n = as.numeric(n_bin))
    df <- df[is.finite(df$r2) & is.finite(df$n) & df$n >= min_pairs, , drop = FALSE]
    df[order(df$dis), , drop = FALSE]
  } else {
    r2 <- x$mean_r2
    ok <- is.finite(dis_kb) & is.finite(r2)
    ok <- ok & is.finite(dist_bp) & dist_bp >= min_bp & dis_kb >= 0 & dis_kb <= xmax_kb
    dis_kb <- dis_kb[ok]
    r2 <- r2[ok]
    br <- seq(0, xmax_kb, by = bin_kb)
    if (tail(br, 1) < xmax_kb) br <- c(br, xmax_kb)
    g <- cut(dis_kb, breaks = br, include.lowest = TRUE, right = FALSE)
    r2_bin <- tapply(r2, g, mean, na.rm = TRUE)
    n_bin <- tapply(r2, g, length)
    centers <- br[-1] - bin_kb / 2
    df <- data.frame(dis = centers, r2 = as.numeric(r2_bin), n = as.numeric(n_bin))
    df <- df[is.finite(df$r2) & is.finite(df$n) & df$n >= min_pairs, , drop = FALSE]
    df[order(df$dis), , drop = FALSE]
  }
}

lddecay_half <- function(df, half_method = c("interpolate", "nearest")) {
  half_method <- match.arg(half_method)
  df <- df[is.finite(df$dis) & is.finite(df$r2), , drop = FALSE]
  df <- df[order(df$dis), , drop = FALSE]
  if (!nrow(df)) return(c(x = NA_real_, y = NA_real_))
  max_r2 <- max(df$r2, na.rm = TRUE)
  half_y <- max_r2 / 2
  if (half_method == "nearest") {
    idx <- which.min(abs(df$r2 - half_y))
    return(c(x = df$dis[idx], y = half_y))
  }
  below <- which(df$r2 <= half_y)
  if (!length(below)) {
    idx <- which.min(abs(df$r2 - half_y))
    return(c(x = df$dis[idx], y = half_y))
  }
  i2 <- below[1]
  if (i2 == 1) return(c(x = df$dis[1], y = half_y))
  i1 <- i2 - 1
  x1 <- df$dis[i1]
  y1 <- df$r2[i1]
  x2 <- df$dis[i2]
  y2 <- df$r2[i2]
  if (!is.finite(y1) || !is.finite(y2) || y2 == y1) return(c(x = df$dis[i2], y = half_y))
  xh <- x1 + (half_y - y1) * (x2 - x1) / (y2 - y1)
  c(x = xh, y = half_y)
}

lddecay_plot <- function(
  files,
  labels = NULL,
  out = NULL,
  device = c("pdf", "png", "tiff", "none"),
  width = 7.5,
  height = 6.2,
  res = 300,
  useDingbats = FALSE,
  xmax_kb = 500,
  bin_kb = 5,
  min_pairs = 50,
  min_bp = 15,
  r2_source = c("sum", "mean"),
  geom = c("points", "lines", "both"),
  col = c("springgreen4", "tomato3", "steelblue3"),
  point_alpha = 0.85,
  pch = c(15, 16, 17),
  lty = 1,
  lwd = 1.8,
  cex = 0.9,
  xlab = "Distance (Kb)",
  ylab = "R square",
  xlim = NULL,
  ylim = NULL,
  x_ticks_by = 100,
  y_ticks_by = 0.1,
  mar = c(5.5, 5.2, 3.5, 2.2),
  legend_pos = "topright",
  legend_cex = 1,
  half = TRUE,
  half_method = c("interpolate", "nearest"),
  half_label_style = c("legend", "arrow", "both", "none"),
  half_label_arrows = NULL,
  half_label_digits_kb = 1,
  half_label_digits_r2 = 2,
  half_label_cex = 0.85,
  half_label_gap = 0.04,
  half_label_xpad = 0.03,
  half_label_ypad = 0.06
) {
  files <- files[file.exists(files)]
  if (!length(files)) stop("no input files found")
  if (is.null(labels)) labels <- sub("\\.Out\\.LDdecay\\.stat$", "", basename(files))
  device <- match.arg(device)
  geom <- match.arg(geom)
  r2_source <- match.arg(r2_source)
  half_method <- match.arg(half_method)
  half_label_style <- match.arg(half_label_style)
  show_half_arrows <- if (is.null(half_label_arrows)) {
    (half_label_style == "arrow") || (half_label_style == "both" && legend_pos == "none")
  } else {
    isTRUE(half_label_arrows)
  }

  df_list <- lapply(files, lddecay_bin, xmax_kb = xmax_kb, bin_kb = bin_kb, min_pairs = min_pairs, min_bp = min_bp, r2_source = r2_source)
  nonempty <- vapply(df_list, nrow, integer(1)) > 0
  if (!any(nonempty)) stop("all inputs are empty after filtering")
  df_list <- df_list[nonempty]
  files <- files[nonempty]
  labels <- labels[nonempty]

  n <- length(df_list)
  col <- rep(col, length.out = n)
  col_pt <- grDevices::adjustcolor(col, alpha.f = point_alpha)
  pch <- rep(pch, length.out = n)
  lty <- rep(lty, length.out = n)
  lwd <- rep(lwd, length.out = n)

  max_x_raw <- max(vapply(df_list, function(z) max(z$dis, na.rm = TRUE), numeric(1)), na.rm = TRUE)
  xmax <- max(300, x_ticks_by * ceiling(max_x_raw / x_ticks_by))
  if (is.null(xlim)) xlim <- c(0, xmax)

  max_y_raw <- max(vapply(df_list, function(z) max(z$r2, na.rm = TRUE), numeric(1)), na.rm = TRUE)
  ymax <- max(0.6, ceiling(max_y_raw / y_ticks_by) * y_ticks_by)
  if (is.null(ylim)) ylim <- c(0, ymax)

  open_device <- function() {
    if (device == "none") return(invisible(NULL))
    if (is.null(out) || !nzchar(out)) stop("out must be provided unless device='none'")
    if (device == "pdf") return(pdf(out, width = width, height = height, useDingbats = useDingbats))
    if (device == "png") return(png(out, width = width, height = height, units = "in", res = res))
    if (device == "tiff") return(tiff(out, width = width, height = height, units = "in", res = res, compression = "lzw"))
  }
  close_device <- function() {
    if (device != "none") dev.off()
  }

  open_device()
  on.exit(close_device(), add = TRUE)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = mar, las = 1)

  plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)

  if (geom %in% c("points", "both")) {
    for (i in seq_len(n)) {
      z <- df_list[[i]]
      points(z$dis, z$r2, pch = pch[i], cex = cex, col = col_pt[i])
    }
  }
  if (geom %in% c("lines", "both")) {
    for (i in seq_len(n)) {
      z <- df_list[[i]]
      lines(z$dis, z$r2, col = col[i], lty = lty[i], lwd = lwd[i])
    }
  }

  x_ticks <- seq(xlim[1], xlim[2], by = x_ticks_by)
  y_ticks <- seq(ylim[1], ylim[2], by = y_ticks_by)
  axis(1, at = x_ticks, tck = -0.01, cex.axis = 1)
  axis(2, at = y_ticks, tck = -0.01, las = 1, cex.axis = 1)
  mtext(xlab, side = 1, line = 3.5, cex = 1.2, font = 2)
  mtext(ylab, side = 2, line = 3.5, cex = 1.2, font = 2, las = 0)
  title(main = "Linkage Disequilibrium", font.main = 2, cex.main = 1.5)


  half_list <- vector("list", n)
  if (isTRUE(half)) {
    for (i in seq_len(n)) {
      half_list[[i]] <- lddecay_half(df_list[[i]], half_method = half_method)
      hx <- half_list[[i]][["x"]]
      hy <- half_list[[i]][["y"]]
      if (!is.finite(hx) || !is.finite(hy)) next
      segments(x0 = xlim[1], y0 = hy, x1 = hx, y1 = hy, lty = 2, col = col[i])
      segments(x0 = hx, y0 = ylim[1], x1 = hx, y1 = hy, lty = 2, col = col[i])
    }
  }

  make_label <- function(i) {
    hx <- half_list[[i]][["x"]]
    hy <- half_list[[i]][["y"]]
    kb <- formatC(hx, format = "f", digits = half_label_digits_kb)
    r2 <- formatC(hy, format = "f", digits = half_label_digits_r2)
    paste0(labels[i], ": ", kb, " kb (r2 = ", r2, ")")
  }

  if (isTRUE(half) && isTRUE(show_half_arrows)) {
    y_top <- ylim[2]
    x_rng <- diff(xlim)
    y_rng <- diff(ylim)
    x_pad <- x_rng * half_label_xpad
    y_pad <- y_rng * half_label_ypad
    used_y <- numeric(0)

    ord <- order(vapply(half_list, function(h) h[["y"]], numeric(1)), decreasing = TRUE)
    for (k in seq_along(ord)) {
      i <- ord[k]
      hx <- half_list[[i]][["x"]]
      hy <- half_list[[i]][["y"]]
      if (!is.finite(hx) || !is.finite(hy)) next
      to_right <- hx < (xlim[1] + 0.68 * x_rng)
      tx <- if (to_right) hx + x_pad else hx - x_pad
      tx <- min(max(tx, xlim[1] + x_pad), xlim[2] - x_pad)
      ty <- min(hy + y_pad, y_top - y_pad)
      if (length(used_y)) {
        mid <- (ylim[1] + ylim[2]) / 2
        for (iter in seq_len(200)) {
          if (all(abs(ty - used_y) >= half_label_gap)) break
          ty <- if (ty >= mid) ty - half_label_gap else ty + half_label_gap
          ty <- min(max(ty, ylim[1] + y_pad), y_top - y_pad)
        }
      }
      used_y <- c(used_y, ty)
      lab <- make_label(i)
      text(tx, ty, labels = lab, col = col[i], cex = half_label_cex, adj = if (to_right) c(0, 0.5) else c(1, 0.5))
      arrows(x0 = tx, y0 = ty, x1 = hx, y1 = hy, length = 0.08, col = col[i], angle = 18)
    }
  }

  legend_labels <- labels
  if (isTRUE(half) && half_label_style %in% c("legend", "both")) {
    legend_labels <- vapply(seq_len(n), make_label, character(1))
  }
  if (legend_pos != "none") {
    legend(legend_pos, legend = legend_labels, pch = if (geom == "lines") NA_integer_ else pch, col = col, lty = if (geom == "points") NA_integer_ else lty, lwd = if (geom == "points") NA_real_ else lwd, bty = "n", cex = legend_cex)
  }

  invisible(list(files = files, labels = labels, data = df_list, half = half_list, xlim = xlim, ylim = ylim))
}

LDplot <- function(files, out = "LDdecay.pdf", ...) {
  lddecay_plot(files = files, out = out, device = "pdf", ...)
}
