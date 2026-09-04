`GAPIT.Multiple.Manhattan` <-
function(model_store,DPP=50000,chor_taxa=NULL,cutOff=0.01,band=5,seqQTN=NULL,byTraits=FALSE,
    Y.names=NULL,GM=NULL,interQTN=NULL,WS=10e5,outpch=NULL,inpch=NULL,
    plot.style="MarineBreeze",plot.line=TRUE,plot.type=c("h","s","w")){
    #Object: Make a Manhattan Plot
    #Output: pdfs of the Multiple Manhattan Plot
    #Authors: Zhiwu Zhang and Jiabo Wang
    # Last update: AUG 24, 2022
    ##############################################################################################
  .gpp_extract_gwas_result <- function(environ_result) {
    df <- environ_result
    cn <- colnames(df)
    cn <- trimws(cn)
    cn <- sub("^\ufeff", "", cn)
    cn <- sub("^ï\\.\\.", "", cn)
    cn <- gsub("[[:space:][:cntrl:]]+", "", cn)
    cn_norm <- gsub("[^a-z0-9]+", "_", tolower(cn))
    cn_norm <- gsub("_+", "_", cn_norm)
    cn_norm <- gsub("^_+|_+$", "", cn_norm)
    pick_col <- function(norm_names, candidates) {
      idx <- which(norm_names %in% candidates)
      if (length(idx) >= 1) return(idx[1])
      integer(0)
    }
    to_num <- function(x) suppressWarnings(as.numeric(x))
    snp_idx <- pick_col(cn_norm, c("snp", "rs", "rsid", "marker", "id"))
    chr_idx <- pick_col(cn_norm, c("chr", "chrom", "chromosome"))
    pos_idx <- pick_col(cn_norm, c("pos", "position", "bp", "bp_position"))
    maf_idx <- pick_col(cn_norm, c("maf", "minor_allele_frequency", "minorallelefrequency"))
    p_candidates <- unique(c(which(grepl("^trait_", cn_norm)), which(cn_norm %in% c("p", "p_value", "pvalue", "pval", "p_wald", "p_value_wald"))))
    p_candidates <- setdiff(p_candidates, which(cn_norm %in% c("maf", "effect", "se", "stderr", "beta", "nobs")))
    pick_p <- function(df, idxs) {
      if (length(idxs) < 1) return(integer(0))
      scores01 <- rep(-Inf, length(idxs))
      scorespos <- rep(-Inf, length(idxs))
      for (i in seq_along(idxs)) {
        v <- to_num(df[[idxs[i]]])
        ok01 <- is.finite(v) & v > 0 & v <= 1
        okpos <- is.finite(v) & v > 0
        scores01[i] <- sum(ok01)
        scorespos[i] <- sum(okpos)
      }
      if (max(scores01, na.rm = TRUE) > 0) return(idxs[which.max(scores01)])
      if (max(scorespos, na.rm = TRUE) > 0) return(idxs[which.max(scorespos)])
      integer(0)
    }
    p_idx <- pick_p(df, p_candidates)
    if (length(p_idx) < 1) {
      p2 <- which(grepl("^p($|_)", cn_norm))
      p_idx <- pick_p(df, p2)
    }
    if (length(chr_idx) < 1 || length(pos_idx) < 1 || length(p_idx) < 1) {
      stop("Cannot find required columns (Chr/Pos/P) in GWAS result. Columns: ", paste(cn, collapse = ", "))
    }
    out <- data.frame(
      SNP = if (length(snp_idx) >= 1) as.character(df[[snp_idx]]) else df[[1]],
      Chr = df[[chr_idx]],
      Pos = df[[pos_idx]],
      P.value = df[[p_idx]],
      MAF = if (length(maf_idx) >= 1) df[[maf_idx]] else NA_real_,
      stringsAsFactors = FALSE
    )
    .gpp_clean_chr <- function(ch) {
      ch <- as.character(ch)
      ch <- trimws(ch)
      ch <- sub("^\ufeff", "", ch)
      ch <- sub("^ï\\.\\.", "", ch)
      ch <- gsub("[[:space:][:cntrl:]]+", "", ch)
      num_ch <- suppressWarnings(as.numeric(ch))
      not_missing = !is.na(ch) & nzchar(ch)
      has_non_numeric_label = any(not_missing & is.na(num_ch))
      if (has_non_numeric_label) {
        return(ch)
      }
      all_int = if (any(not_missing)) {
        all(is.finite(num_ch[not_missing]) & (num_ch[not_missing] == round(num_ch[not_missing])))
      } else TRUE
      if (all_int) {
        return(num_ch)
      }
      ch
    }
    out$Chr <- .gpp_clean_chr(out$Chr)
    out$Pos <- to_num(out$Pos)
    p <- to_num(out$P.value)
    if (sum(is.finite(p) & p > 0 & p <= 1) == 0 && sum(is.finite(p) & p > 0) > 0 && stats::median(p, na.rm = TRUE) > 1) {
      p <- 10^(-p)
    }
    out$P.value <- p
    out$MAF <- to_num(out$MAF)
    out
  }
  .gpp_clean_chr <- function(ch) {
    ch <- as.character(ch)
    ch <- trimws(ch)
    ch <- sub("^\ufeff", "", ch)
    ch <- sub("^ï\\.\\.", "", ch)
    ch <- gsub("[[:space:][:cntrl:]]+", "", ch)
    num_ch <- suppressWarnings(as.numeric(ch))
    not_missing = !is.na(ch) & nzchar(ch)
    has_non_numeric_label = any(not_missing & is.na(num_ch))
    if (has_non_numeric_label) {
      return(ch)
    }
    all_int = if (any(not_missing)) {
      all(is.finite(num_ch[not_missing]) & (num_ch[not_missing] == round(num_ch[not_missing])))
    } else TRUE
    if (all_int) {
      return(num_ch)
    }
    ch
  }
  .gpp_map_xy <- function(chr_vec) {
    chr_vec_c <- as.character(chr_vec)
    options(warn = -1)
    numeric.chr <- suppressWarnings(as.numeric(chr_vec_c))
    options(warn = 0)
    max.chr <- suppressWarnings(max(numeric.chr, na.rm = TRUE))
    if (!is.finite(max.chr)) max.chr <- 0
    map.xy.index <- which(!(numeric.chr %in% c(0:max.chr)))
    chr.xy <- character(0)
    if (length(map.xy.index) != 0) {
      chr.xy <- unique(chr_vec_c[map.xy.index])
      for (k in seq_along(chr.xy)) {
        chr_vec_c[chr_vec_c == chr.xy[k]] <- max.chr + k
      }
    }
    list(
      mapped = suppressWarnings(as.numeric(chr_vec_c)),
      chr.xy = chr.xy,
      max.chr = max.chr
    )
  }
  if (!is.null(GM)) {
    if (ncol(GM) >= 2) {
      gm_chr <- .gpp_clean_chr(GM[, 2])
      GM[, 2] <- gm_chr
    }
    if (ncol(GM) >= 3) {
      GM[, 3] <- suppressWarnings(as.numeric(GM[, 3]))
    }
  }
  Nenviron=length(model_store)*length(Y.names)
  environ_name=NULL
  new_xz=NULL
  if(byTraits)
  {
    for(i in 1:length(Y.names))
    {
       for(j in 1:length(model_store))
       {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
          environ_name=c(environ_name,paste(model_store[j],".",Y.names[i],sep=""))
       }
    }
  }else{
    for(i in 1:length(model_store))
    {
       for(j in 1:length(Y.names))
       {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
          environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
       }
    }
    }
# ---- Global chromosome set (over ALL traits + GM) ----
.global_combined_chrs = character(0)
if (!is.null(GM) && nrow(GM) >= 1 && ncol(GM) >= 2) {
  .global_combined_chrs = c(.global_combined_chrs, as.character(GM[, 2]))
}
for (.tmp_i in seq_along(environ_name)) {
  .tmp_f = paste("GAPIT.Association.GWAS_Results.", environ_name[.tmp_i], ".csv", sep = "")
  .tmp_df = tryCatch(read.csv(.tmp_f, head = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(.tmp_df) && nrow(.tmp_df) >= 1) {
    .tmp_ex = tryCatch(.gpp_extract_gwas_result(.tmp_df), error = function(e) NULL)
    if (!is.null(.tmp_ex) && nrow(.tmp_ex) >= 1) {
      .global_combined_chrs = c(.global_combined_chrs, as.character(.tmp_ex[, 2]))
    }
  }
}
.global_combined_chrs = .global_combined_chrs[!is.na(.global_combined_chrs) & nzchar(.global_combined_chrs)]
.global_map = .gpp_map_xy(.global_combined_chrs)
.global_unique_chr = as.character(unique(.global_combined_chrs))
.global_unique_ord = order(.global_map$mapped[match(.global_unique_chr, as.character(.global_combined_chrs))])
.chm_to_analyze = .global_unique_chr[.global_unique_ord]
.chm_to_analyze = .chm_to_analyze[!is.na(.chm_to_analyze) & nzchar(.chm_to_analyze)]
.nchr = length(.chm_to_analyze)
.chm_to_pos = setNames(seq_along(.chm_to_analyze), .chm_to_analyze)
print(paste("[Info] Global chromosomes detected (", .nchr, "): ",
            paste(head(.chm_to_analyze, 20), collapse = ", "),
            if (.nchr > 20) paste0(" ...(+", .nchr - 20, ")") else "", sep = ""))
non_num_chrs = as.character(.global_map$chr.xy)
non_num_chrs = non_num_chrs[!is.na(non_num_chrs) & nzchar(non_num_chrs)]
if (length(non_num_chrs) >= 1) {
  print(paste("[Info] Non-numeric chromosomes (X/Y/MT/...) preserved: ",
              paste(non_num_chrs, collapse = ", "), sep = ""))
} else {
  print("[Info] No non-numeric chromosomes found; all Chr were numeric.")
}
# ---- END Global chromosome set ----
sig_pos_snp=NULL
sig_pos_chro=NULL
sig_pos_newx=NULL
simulation=FALSE
if(!is.null(seqQTN)){
	simulation=TRUE
}
themax.y0=NULL
store.x=NULL
y_filter0=NULL
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[i],".csv",sep=""),head=T)
  environ_result=.gpp_extract_gwas_result(environ_result)
  num.markers=nrow(environ_result)
  .gpp_xy <- .gpp_map_xy(environ_result[, 2])
  keep_xy_idx <- which(is.finite(.gpp_xy$mapped) &
                         !is.na(environ_result[, 3]) & is.finite(environ_result[, 3]) &
                         !is.na(environ_result[, 4]) & is.finite(environ_result[, 4]))
  if (length(keep_xy_idx) < nrow(environ_result)) {
    environ_result <- environ_result[keep_xy_idx, , drop=FALSE]
    .gpp_xy$mapped <- .gpp_xy$mapped[keep_xy_idx]
  }
  ord1 <- order(environ_result[, 3])
  environ_result <- environ_result[ord1, , drop=FALSE]
  .gpp_xy$mapped <- .gpp_xy$mapped[ord1]
  ord2 <- order(.gpp_xy$mapped)
  environ_result <- environ_result[ord2, , drop=FALSE]
  .gpp_xy$mapped <- .gpp_xy$mapped[ord2]
  environ_filter=environ_result[!is.na(environ_result[,4]),]
  themax.y=round(max(-log10(environ_filter[,4])),0)
  themax.y0=round(max(c(themax.y,themax.y0)),0)

  y_filter=environ_filter[environ_filter[,4]<(cutOff/(num.markers)),,drop=FALSE]
  traits=environ_name[i]
  if(nrow(y_filter)>0)y_filter=cbind(as.matrix(y_filter[,1:5]),traits)
  y_filter0=rbind(y_filter0,y_filter)

  result=environ_result
  result=result[match(as.character(GM[,1]),as.character(result[,1])),]
  rownames(result)=seq_len(nrow(result))
  if(i==1){
    result0=result
    # Keep only core columns: SNP, Chr, Pos, <trait1>
    colnames(result0)[4]=environ_name[i]
    if (ncol(result0) > 4) result0 <- result0[, 1:4, drop=FALSE]
    }
  if(i!=1){
    result_join <- result[, c(1, 4), drop=FALSE]
    # Merge with explicit suffixes; we will rename the NEW trait column (always last)
    result0=merge(result0, result_join,
                  by.x = colnames(result0)[1],
                  by.y = colnames(result_join)[1],
                  suffixes = c(".x", ".y"),
                  all.x = TRUE)
    # The newly added trait p-value column is ALWAYS the LAST column after merge
    colnames(result0)[ncol(result0)] = environ_name[i]
    # Safety: clean any left-over .x / .y artifact names in case merge left them
    tmp_cn <- colnames(result0)
    if (any(grepl("\\.x$|\\.y$", tmp_cn[-c(1:3)]))) {
      tmp_cn[grepl("^P\\.value(\\.x|\\.y)?$", tmp_cn)] <- environ_name[i]
      colnames(result0) <- tmp_cn
    }
    }
  rownames(result)=seq_len(nrow(result))
  result[is.na(result[,4]),4]=1
  # Use SNP name based significant index instead of rownumber
  m_threshold = cutOff / nrow(result)
  p_col <- suppressWarnings(as.numeric(result[, 4]))
  sig_mask <- !is.na(p_col) & is.finite(p_col) & (p_col < m_threshold)
  if (any(sig_mask)) {
    sig_pos_snp <- c(sig_pos_snp, as.character(result[sig_mask, 1]))
    sig_pos_chro <- c(sig_pos_chro, as.character(result[sig_mask, 2]))
    sig_pos_newx <- c(sig_pos_newx, rep_len(2, sum(sig_mask)))
  }
}
  write.csv(y_filter0,paste("GAPIT.Association.Filter_GWAS_results.csv",sep=""),quote=FALSE)

# Build map_store globally from GM (sorted by mapped to keep X/Y at end)
gm_map_chrs = if (!is.null(GM) && nrow(GM) >= 1 && ncol(GM) >= 2) as.character(GM[, 2]) else character(0)
.gm_global_map = .gpp_map_xy(gm_map_chrs)
.ms_v1 = gm_map_chrs
.ms_v2 = if (!is.null(GM) && nrow(GM) >= 1 && ncol(GM) >= 3) suppressWarnings(as.numeric(GM[, 3])) else numeric(0)
.ms_ok = !is.na(.ms_v1) & nzchar(.ms_v1) & is.finite(.ms_v2)
if (length(.ms_ok) >= 1) {
  .ms_v1 = .ms_v1[.ms_ok]
  .ms_v2 = .ms_v2[.ms_ok]
  .gm_global_map$mapped = .gm_global_map$mapped[.ms_ok]
  if (!is.null(GM) && nrow(GM) >= 1) GM = GM[.ms_ok, , drop=FALSE]
}
map_store = if (length(.ms_v1) >= 1) {
  data.frame(V1 = as.character(.ms_v1), V2 = as.numeric(.ms_v2),
             stringsAsFactors = FALSE)
} else {
  data.frame(V1 = character(0), V2 = numeric(0), stringsAsFactors = FALSE)
}
if (nrow(map_store) >= 1) {
  ms_ord = order(.gm_global_map$mapped, as.numeric(map_store[, 2]), na.last = NA)
  if (length(ms_ord) >= 1) {
    map_store = map_store[ms_ord, , drop=FALSE]
    .gm_global_map$mapped = .gm_global_map$mapped[ms_ord]
    if (!is.null(GM) && nrow(GM) >= 1) GM = GM[ms_ord, , drop=FALSE]
  }
}
.ms_chrs = if (nrow(map_store) >= 1) as.character(map_store[, 1]) else character(0)
.ms_unique = if (length(.ms_chrs) >= 1) {
  ord_tmp = order(.gpp_map_xy(.ms_chrs)$mapped[match(as.character(unique(.ms_chrs)),
                                                      as.character(.ms_chrs))])
  as.character(unique(.ms_chrs))[ord_tmp]
} else character(0)
.ms_unique = .ms_unique[!is.na(.ms_unique) & nzchar(.ms_unique)]
if (length(.ms_unique) >= 1) .chm_to_analyze = unique(c(.ms_unique, .chm_to_analyze))
.chm_global_ord = order(.gpp_map_xy(.chm_to_analyze)$mapped[match(as.character(.chm_to_analyze),
                                                                   as.character(.chm_to_analyze))])
.chm_to_analyze = .chm_to_analyze[.chm_global_ord]
.chm_to_analyze = .chm_to_analyze[!is.na(.chm_to_analyze) & nzchar(.chm_to_analyze)]
.nchr = length(.chm_to_analyze)
.chm_to_pos = setNames(seq_along(.chm_to_analyze), .chm_to_analyze)
# Build cumulative positions + ticks (global, including X/Y)
lastbase=0
ticks=NULL
max.x=NULL
.ms_len = if (nrow(map_store) >= 1) seq_len(nrow(map_store)) else integer(0)
if (length(.ms_unique) >= 1 && nrow(map_store) >= 1) {
  .ms_allchrs = as.character(map_store[, 1])
  for (j in .ms_unique)
  {
    idx = which(.ms_allchrs == j)
    if (length(idx) == 0) next
    posv = as.numeric(map_store[idx, 2])
    posv[!is.finite(posv)] = 0
    mn = mean(posv, na.rm = TRUE)
    if (!is.finite(mn)) mn = 0
    ticks = c(ticks, lastbase + mn)
    posv2 = posv + lastbase
    posv2[!is.finite(posv2)] = lastbase
    map_store[idx, 2] = posv2
    mx = max(posv2, na.rm = TRUE)
    if (!is.finite(mx)) mx = lastbase
    lastbase = mx
    max.x = c(max.x, mx)
  }
  minv = suppressWarnings(min(as.numeric(map_store[, 2]), na.rm = TRUE))
  if (!is.finite(minv)) minv = 0
  max.x = c(minv, max.x)
} else {
  max.x = 0
}
store.x = as.numeric(map_store[, 2])
if (length(store.x) == 0) store.x = 0
# Build new_xz using SNP name (match on GM)
if (length(sig_pos_snp) >= 2) {
  sig_pos_snp_uniq = unique(as.character(sig_pos_snp))
  gm_snps = if (!is.null(GM) && nrow(GM) >= 1) as.character(GM[, 1]) else character(0)
  ms_rows = match(sig_pos_snp_uniq, gm_snps)
  ok_rows = which(is.finite(ms_rows) & ms_rows >= 1 & ms_rows <= nrow(map_store))
  if (length(ok_rows) >= 2) {
    sig_pos_snp_uniq = sig_pos_snp_uniq[ok_rows]
    ms_rows = ms_rows[ok_rows]
    new_xz0_df = data.frame(
      pos = as.character(sig_pos_snp_uniq),
      times = as.integer(2),
      chro = as.character(map_store[ms_rows, 1]),
      xlab = as.numeric(map_store[ms_rows, 2]),
      stringsAsFactors = FALSE
    )
    # Collapse duplicates (same new_xz) similar to old dayu logic
    xlab_vals = as.numeric(new_xz0_df$xlab)
    scom = sort(xlab_vals)
    if (length(scom) >= 2) {
      de.sc = scom[-1] - scom[-length(scom)]
      dayu1.index = duplicated(scom) | c(abs(de.sc) < WS, FALSE)
      if (any(dayu1.index)) {
        scom2 = scom[dayu1.index]
        sc.index = as.character(new_xz0_df$xlab) %in% as.character(scom2)
        new_xz_df = new_xz0_df[sc.index, , drop=FALSE]
      } else {
        new_xz_df = new_xz0_df
      }
    } else {
      new_xz_df = new_xz0_df
    }
    if (nrow(new_xz_df) >= 1) {
      dup_idx = duplicated(new_xz_df$xlab)
      if (any(dup_idx)) new_xz_df$times[dup_idx] = 1
      new_xz_df = new_xz_df[!duplicated(new_xz_df), , drop=FALSE]
      new_xz_df = new_xz_df[as.character(new_xz_df$times) != "0", , drop=FALSE]
      # Convert chromosome ids to mapped integers for new_xz[,3]
      chro_to_int = suppressWarnings(as.numeric(as.character(new_xz_df$chro)))
      map_int = .gpp_map_xy(as.character(new_xz_df$chro))
      chro_to_int = as.numeric(map_int$mapped)
      new_xz = cbind(
        pos = match(as.character(new_xz_df$pos), gm_snps),
        times = as.integer(new_xz_df$times),
        chro = chro_to_int,
        xlab = as.numeric(new_xz_df$xlab)
      )
      new_xz = matrix(as.numeric(new_xz), nrow = nrow(new_xz), ncol = 4)
    }
  }
}
if (!exists("new_xz", inherits = FALSE) || !is.matrix(new_xz) || nrow(new_xz) < 1) new_xz = NULL

# print(new_xz)
# setup colors
# print(head(result))
# chm.to.analyze <- unique(result[,2])
chm.to.analyze <- .chm_to_analyze
nchr <- .nchr
size=1 #1
ratio=10 #5
base=1 #1
numCHR=nchr
numMarker=nrow(GM)
bonferroniCutOff=-log10(cutOff/numMarker)

ncycle=ceiling(nchr/5)
ncolor=band*ncycle
thecolor=seq(1,nchr,by= ncycle)
col.Rainbow=rainbow(ncolor+1)     
col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))    
col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))     
col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(numCHR/5))
col.MarineBreeze=rep(c("#BFDFD2","#51999F","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A"),ceiling(numCHR/8))
col.cougars=rep(c(  '#990000',    'dimgray'),ceiling(numCHR/2))  
if(plot.style=="Rainbow")plot.color= col.Rainbow
if(plot.style =="FarmCPU")plot.color= col.Rainbow
if(plot.style =="Rushville")plot.color= col.Rushville
if(plot.style =="Congress")plot.color= col.Congress
if(plot.style =="Ocean")plot.color= col.Ocean
if(plot.style =="PLINK")plot.color= col.PLINK
if(plot.style =="Beach")plot.color= col.Beach
if(plot.style=="Oceanic")plot.color= col.Oceanic
if(plot.style=="MarineBreeze")plot.color= col.MarineBreeze
if(plot.style =="cougars")plot.color= col.cougars  

if("h"%in%plot.type)
{
    Max.high=6*Nenviron
    if(Max.high>8)Max.high=40
    pdf(paste("GAPIT.Association.Manhattans_High",".pdf" ,sep = ""), width = 20,height=6*Nenviron)
    par(mfrow=c(Nenviron,1))
    mypch=1
    for(k in 1:Nenviron)
    { 
       if(k==Nenviron)
        {
        par(mar = c(3.5,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(1.5,8,0.5,8))    
        }
       environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
       environ_result=.gpp_extract_gwas_result(environ_result)
       result=environ_result
       .gpp_h0_map <- .gpp_map_xy(result[, 2])
       keep_h0 <- which(is.finite(.gpp_h0_map$mapped) & !is.na(result[,3]) & !is.na(result[,4]))
       if (length(keep_h0) < nrow(result)) {
         result <- result[keep_h0, , drop=FALSE]
         .gpp_h0_map$mapped <- .gpp_h0_map$mapped[keep_h0]
       }
       result_h0_ord <- order(result[, 3])
       result <- result[result_h0_ord, , drop=FALSE]
       .gpp_h0_map$mapped <- .gpp_h0_map$mapped[result_h0_ord]
       result_h0_ord2 <- order(.gpp_h0_map$mapped)
       result <- result[result_h0_ord2, , drop=FALSE]
       .gpp_h0_map$mapped <- .gpp_h0_map$mapped[result_h0_ord2]
       result=result[match(as.character(GM[,1]),as.character(result[,1])),]
       rownames(result)=seq_len(nrow(result))
       GI.MP=result[,c(2:4)]
       borrowSlot=4
       GI.MP[,borrowSlot]=0 #Inicial as 0
       GI.MP[,5]=1:(nrow(GI.MP))
       GI.MP[,6]=1:(nrow(GI.MP)) 
       GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
       GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
       GI.MP[is.na(GI.MP[,3]),3]=1
    #Retain SNPs that have P values between 0 and 1 (not na etc)
       GI.MP <- GI.MP[GI.MP[,3]>0,]
       GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 only (preserve X/Y/99 etc non-zero non-standard chromosomes)
       zero_h1_chr <- suppressWarnings(as.numeric(as.character(GI.MP[,1])))
       zero_h1_mask <- is.finite(zero_h1_chr) & zero_h1_chr == 0
       GI.MP <- GI.MP[!zero_h1_mask,,drop=FALSE]
       total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
       if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
       numMarker=nrow(GI.MP)
       GI.MP[,3] <-  -log10(GI.MP[,3])
       GI.MP[,5]=1:numMarker
       y.lim <- ceiling(max(GI.MP[,3]))  
       .gpp_h1_map <- .gpp_map_xy(GI.MP[,1])
       .local_chr = as.character(unique(GI.MP[,1]))
       .merged_chr = unique(c(as.character(.chm_to_analyze), .local_chr))
       .map_merged = .gpp_map_xy(.merged_chr)
       chm.to.analyze = .merged_chr[order(.map_merged$mapped)]
       nchr = length(chm.to.analyze)
       GI.MP_order_chr <- .gpp_h1_map$mapped
       GI.MP <- GI.MP[order(GI.MP_order_chr, GI.MP[,2]), ]
       GI.MP[,6]=1:(nrow(GI.MP))
       MP_store=GI.MP
       index_GI=MP_store[,3]>=0
       MP_store <- MP_store[index_GI,]
       ticks=NULL
       lastbase=0
       # print(head(MP_store))
       for(i in chm.to.analyze)
          {
           index=(MP_store[,1]==i)
           ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
           MP_store[index,2]=MP_store[index,2]+lastbase
           lastbase=max(MP_store[index,2])
          }
       x0 <- as.numeric(MP_store[,2])
       y0 <- as.numeric(MP_store[,3])
       chor_taxa <- as.character(chm.to.analyze)
       chor_to_z <- setNames(seq_along(chor_taxa), chor_taxa)
       z0 <- as.numeric(chor_to_z[as.character(MP_store[,1])])
       if(anyNA(z0)) z0[is.na(z0)] <- 1L
       # print(ticks)
       x1=sort(x0)
       position=order(y0,decreasing = TRUE)
       values=y0[position]
       if(length(values)<=DPP)
         {
         index=position[c(1:length(values))]
         }else{       
          # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
        index=position[which(values>cut0)]
         }        
       x=x0[index]
       y=y0[index]
       z=z0[index]
        #Extract QTN
       QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
       
       themax=ceiling(max(y))
       themax2=ceiling((ceiling(themax/4))*4)

       themin=floor(min(y))
       wd=((y-themin+base)/(themax-themin+base))*size*ratio
       s=size-wd/ratio/2
       plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s+2.5,cex.main=2)
       mtext(side=2,expression(-log[10](italic(p))),line=3.5, cex=2.5)
       if(!simulation)
         {
          abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], pch=21, cex=2.8,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], pch=20, cex=1.8,lwd=1.5,col="dimgrey")
         }        
       if(plot.line)
         {
          if(!is.null(nrow(new_xz)))  
            {
             abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
            }else{
             abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
            }
         }
        #Add a horizontal line for bonferroniCutOff
       # print(ticks)
       abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
       axis(2, yaxp=c(0,themax2,4),cex.axis=2.3,tick=T,las=1,lwd=2.5)
       if(k==Nenviron)axis(1, at=max.x,cex.axis=2.5,labels=rep("",length(max.x)),tick=T,lwd=2.5)
       if(k==Nenviron)axis(1, at=ticks,cex.axis=2.5,labels=chm.to.analyze,tick=F,line=1)
       mtext(side=4,paste(environ_name[k],sep=""),line=3.2,cex=2)
    }#end of environ_name
       dev.off()
}#end of plot.type

if("w"%in%plot.type)
{
 pdf(paste("GAPIT.Association.Manhattans_Wide",".pdf" ,sep = ""), width = 16,height=8.5)
 par(mfrow=c(Nenviron,1))
 mtext.h=0.5
 size=2
 ratio=5
 for(k in 1:Nenviron)
 { 
  if(k==Nenviron)
        {#par(mfrow=c(Nenviron,1))
          # print(par())
        par(mar = c(2.5,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(2,8,0.5,8))    
         # par(pin=c(10,(8-mtext.h)/Nenviron))

        }
  environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
  environ_result=.gpp_extract_gwas_result(environ_result)
  #print(environ_result[as.numeric(new_xz[,1]),])
  result=environ_result
    .gpp_w0_map <- .gpp_map_xy(result[, 2])
    keep_w0 <- which(is.finite(.gpp_w0_map$mapped) & !is.na(result[,3]) & !is.na(result[,4]))
    if (length(keep_w0) < nrow(result)) {
      result <- result[keep_w0, , drop=FALSE]
      .gpp_w0_map$mapped <- .gpp_w0_map$mapped[keep_w0]
    }
    result_w0_ord <- order(result[, 3])
    result <- result[result_w0_ord, , drop=FALSE]
    .gpp_w0_map$mapped <- .gpp_w0_map$mapped[result_w0_ord]
    result_w0_ord2 <- order(.gpp_w0_map$mapped)
    result <- result[result_w0_ord2, , drop=FALSE]
    .gpp_w0_map$mapped <- .gpp_w0_map$mapped[result_w0_ord2]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=seq_len(nrow(result))
    GI.MP=result[,c(2:4)]
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP))
    
    
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 only (preserve X/Y/99 etc non-zero non-standard chromosomes)
    zero_w1_chr <- suppressWarnings(as.numeric(as.character(GI.MP[,1])))
    zero_w1_mask <- is.finite(zero_w1_chr) & zero_w1_chr == 0
    GI.MP <- GI.MP[!zero_w1_mask,,drop=FALSE]
    total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    .gpp_w1_map <- .gpp_map_xy(GI.MP[,1])
    .local_chr = as.character(unique(GI.MP[,1]))
    .merged_chr = unique(c(as.character(.chm_to_analyze), .local_chr))
    .map_merged = .gpp_map_xy(.merged_chr)
    chm.to.analyze = .merged_chr[order(.map_merged$mapped)]
    nchr = length(chm.to.analyze)
    GI.MP_order_chr <- .gpp_w1_map$mapped
    GI.MP <- GI.MP[order(GI.MP_order_chr, GI.MP[,2]), ]
    GI.MP[,6]=1:(nrow(GI.MP))
    MP_store=GI.MP
        index_GI=MP_store[,3]>=0
        MP_store <- MP_store[index_GI,]
        ticks=NULL
        lastbase=0
        for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
        x0 <- as.numeric(MP_store[,2])
        y0 <- as.numeric(MP_store[,3])
        chor_taxa <- as.character(chm.to.analyze)
        chor_to_z <- setNames(seq_along(chor_taxa), chor_taxa)
        z0 <- as.numeric(chor_to_z[as.character(MP_store[,1])])
       if(anyNA(z0)) z0[is.na(z0)] <- 1L
       x1=sort(x0)

       position=order(y0,decreasing = TRUE)
       values=y0[position]
        if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
            }else{      
        # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
      
        index=position[which(values>cut0)]
        }     
        x=x0[index]
        y=y0[index]
        z=z0[index]
        # print(length(x))
        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
        QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
        themax=ceiling(max(y))
        themax2=ceiling((ceiling(themax/4))*4)
        themin=floor(min(y))
        # size=5
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        # wd=0.5
        s=size-wd/ratio/2
        mypch=1
        bamboo=4
        # plot(y~x,xlab="",ylab="" ,ylim=c(0,themax),
        #     cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=0.5,cex=0.7,cex.main=2)
        plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s,cex.main=2)
        mtext(side=2,expression(-log[10](italic(p))),line=3, cex=1)
        if(plot.line)
        {
          if(!is.null(nrow(new_xz)))  {abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=2)
             }else{abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=2)
             }
        }
        if(!simulation){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          # print("$$$")
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=2.8,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], type="p",pch=20, cex=1.5,lwd=1.5,col="dimgrey")
          }
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,lty=1,untf=T,lwd=1,col="forestgreen")
        axis(2, yaxp=c(0,themax2,bamboo),cex.axis=1.5,las=1,tick=F)
        if(k==Nenviron)axis(1, at=ticks,cex.axis=1.5,line=0.001,labels=chm.to.analyze,tick=F)
        mtext(side=4,paste(environ_name[k],sep=""),line=2,cex=1)
 box()
 }#end of environ_name
 dev.off()
}#end of plot.type

if("s"%in%plot.type)
{
    # wd=((y-themin+base)/(themax-themin+base))*size*ratio
 wd=2
 if(is.null(outpch))
 {
   allpch0=c(0,1,2,5,6)
 }else{
   allpch0=outpch
 }
 if(is.null(inpch))
 {
   add.pch=c("+","*","-","#","<",">","^","$","=","|","?",as.character(1:9),letters[1:26],LETTERS[1:26]) 
 }else{
   add.pch=inpch
 }
 n.vals=ceiling(Nenviron/length(allpch0))-1
 s=size-wd/ratio/2
 DPP=500
 


 pdf(paste("GAPIT.Association.Manhattans_Symphysic",".pdf" ,sep = ""), width = 30,height=18)
 par(mfrow=c(1,1))
 par(mar = c(5,8,5,1))
 themax.y02=ceiling((ceiling(themax.y0/4))*4)

 plot(1~1,col="white",xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(store.x,na.rm=TRUE),max(store.x,na.rm=TRUE)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,axes=FALSE,
    pch=1,cex.main=4)
 # print(ticks)   
        #Add a horizontal line for bonferroniCutOff
 axis(1, at=max.x,cex.axis=2,labels=rep("",length(max.x)),tick=T,lwd=2.5)
 axis(1, at=ticks,cex.axis=2,labels=chm.to.analyze,tick=F,line=1)
 axis(2, yaxp=c(0,themax.y02,4),cex.axis=2,tick=T,las=1,lwd=2.5)
 if(!is.null(cutOff))abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
 if(plot.line)
    {
    if(!is.null(nrow(new_xz)))  
        {
        abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
        }else{
        abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
        }
    }
 mtext(side=2,expression(-log[10](italic(p))),line=4, cex=2.5)
 # legend("top",legend=paste(environ_name,sep=""),ncol=length(environ_name),
 #       col="black",pch=allpch[1:Nenviron],lty=0,lwd=1,cex=2,
 #       bty = "o", bg = "white",box.col="white")
 # step.vals=0
 # if(1>2){
 for(k in 1:Nenviron)
  { 
    step.vals=ceiling(k/length(allpch0))-1

    environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
    environ_result=.gpp_extract_gwas_result(environ_result)
    result=environ_result
    .gpp_s0_map <- .gpp_map_xy(result[, 2])
    keep_s0 <- which(is.finite(.gpp_s0_map$mapped) & !is.na(result[,3]) & !is.na(result[,4]))
    if (length(keep_s0) < nrow(result)) {
      result <- result[keep_s0, , drop=FALSE]
      .gpp_s0_map$mapped <- .gpp_s0_map$mapped[keep_s0]
    }
    result_s0_ord <- order(result[, 3])
    result <- result[result_s0_ord, , drop=FALSE]
    .gpp_s0_map$mapped <- .gpp_s0_map$mapped[result_s0_ord]
    result_s0_ord2 <- order(.gpp_s0_map$mapped)
    result <- result[result_s0_ord2, , drop=FALSE]
    .gpp_s0_map$mapped <- .gpp_s0_map$mapped[result_s0_ord2]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=seq_len(nrow(result))
    GI.MP=result[,c(2:4)]
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP)) 
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 only (preserve X/Y/99 etc non-zero non-standard chromosomes)
    zero_s1_chr <- suppressWarnings(as.numeric(as.character(GI.MP[,1])))
    zero_s1_mask <- is.finite(zero_s1_chr) & zero_s1_chr == 0
    GI.MP <- GI.MP[!zero_s1_mask,,drop=FALSE]
    total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    .gpp_s1_map <- .gpp_map_xy(GI.MP[,1])
    .local_chr = as.character(unique(GI.MP[,1]))
    .merged_chr = unique(c(as.character(.chm_to_analyze), .local_chr))
    .map_merged = .gpp_map_xy(.merged_chr)
    chm.to.analyze = .merged_chr[order(.map_merged$mapped)]
    nchr = length(chm.to.analyze)
    GI.MP_order_chr <- .gpp_s1_map$mapped
    GI.MP <- GI.MP[order(GI.MP_order_chr, GI.MP[,2]), ]
    GI.MP[,6]=1:(nrow(GI.MP))
    MP_store=GI.MP
    index_GI=MP_store[,3]>=0
    MP_store <- MP_store[index_GI,]
    ticks=NULL
    lastbase=0
    for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
    x0 <- as.numeric(MP_store[,2])
    y0 <- as.numeric(MP_store[,3])
    chor_taxa <- as.character(chm.to.analyze)
    chor_to_z <- setNames(seq_along(chor_taxa), chor_taxa)
    z0 <- as.numeric(chor_to_z[as.character(MP_store[,1])])
    if(anyNA(z0)) z0[is.na(z0)] <- 1L
    max.x=NULL
    for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            max.x=c(max.x,max(x0[index]))
        }
    max.x=c(min(x0),max.x)
    x1=sort(x0)

    position=order(y0,decreasing = TRUE)
    values=y0[position]
    if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
        }else{       
          # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)

        index=position[which(values>cut0)]
        }        
    x=x0[index]
    y=y0[index]
    z=z0[index]
        # print(length(x))

        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
    QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
    themax=ceiling(max(y))
    # themax.y02=ceiling((ceiling(themax.y0/4)+1)*4)
    # print(themax.y02)
    themin=floor(min(y))
    mypch=allpch0[k-step.vals*length(allpch0)]
    
   # if(k!=1) par(new=T)
    
    par(new=T)
    plot(y~x,xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(x),max(x)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,col=plot.color[z],axes=FALSE,
    pch=mypch,lwd=1,cex=s+2.5,cex.main=4)
    if(step.vals!=0)
    {
      points(y~x,pch=add.pch[step.vals],col=plot.color[z],cex=s+0.5,cex.main=4)
    }
    if(!simulation)
       {
        abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")
        }else{
        points(QTN[,2], QTN[,3], pch=20, cex=2,lwd=2.5,col="dimgrey")
       }        
    
 }#end of environ_name
 dev.off()
 # }
 ## Plot legend
 nchar.traits=1.5
 # environ_name=paste(environ_name,"1234",sep="")
 nchar0=max(nchar(environ_name))
 # print(nchar0)
 if(Nenviron>5)
 {
  yourpch=c(rep(allpch0,n.vals),allpch0[1:(Nenviron-length(allpch0)*n.vals)])
 }else{
  yourpch=allpch0[1:Nenviron]
 }
  yourpch2=NULL
  for(pp in 1:n.vals)
  {
    yourpch2=c(yourpch2,rep(add.pch[pp],length(allpch0)))
  }
  # yourpch2=c(rep(allpch0,n.vals),allpch0[Nenviron-length(allpch0)*n.vals])
  if(Nenviron>5){
  yourpch2=yourpch2[1:(Nenviron-length(allpch0))]
  yourpch2=c(rep(NA,length(allpch0)),yourpch2)
  }
 
 max.row=25
 max.pch=ifelse(Nenviron<max.row,Nenviron,max.row)
 n.col.pch=ceiling(Nenviron/max.row)
 ratio.cex=ceiling(Nenviron/5)
  if(ratio.cex>5)ratio.cex=5
  cex.Ne=ratio.cex*(0.05*ratio.cex+0.35) #the size of cex
  if(ratio.cex<3)cex.Ne=1
  c.t.d=c(0.5,1,1,1,1.5)[ratio.cex] # the different between size of cex and text
  cex.di=0.3*ratio.cex # the different size between cex and signal
  text.di=c(.02,0.01,0.01,0.02,0.02)[ratio.cex] #the different distance between cex and text

  high.Ne=2*ratio.cex # the total highth of figure
  cex.betw=c(.38,.5,.7,0.9,1)[ratio.cex] # distance between cexes
  x.di=c(1.12,1.09,0.5,0.8,1.3)[ratio.cex]/2 # distance between markers in x axis
  # print(nchar0)
  # x.di0=c(0)
  if(n.col.pch>1){
    text.di=(0.01*nchar0)/3+0.02
    # x.di=0.52*n.col.pch
    x.di=(0.1*(ceiling(nchar0/5)-1)+(n.col.pch-1)*0.1)*ceiling(nchar0/5)#*n.col.pch
  }
  # print(x.di)
 # if(Nenviron>5){
 #  cex.Ne=3
 #  cex.di=1.5
 #  text.di=.02
 #  high.Ne=10
 #  cex.betw=0.9
 #  }else{
 #  cex.Ne=1
 #  cex.di=0.3
 #  high.Ne=Nenviron/2 
 #  text.di=.02
 #  cex.betw=0.9
 #  }


 write.csv(environ_name,"GAPIT.Association.Manhattans_Symphysic_Traitsnames.csv",quote=FALSE)
 pdf(paste("GAPIT.Association.Manhattans_Symphysic_Legend",".pdf" ,sep = ""), width = 4+(x.di*(n.col.pch+1)),height=high.Ne)
 par(mfrow=c(1,1))
 par(mar = c(cex.Ne+1,2,cex.Ne+1,2))
 # print(length(yourpch))
 # print(length(yourpch2))
 plot(0,0,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),col="white")
 for(kk in 1:n.col.pch)
 {
 par(new=T)

 if(kk==n.col.pch)
 {
  if(n.col.pch==1)
  {
  # print(kk)
  max.pch2=Nenviron-(n.col.pch-1)*max.row
  
  plot(rep(0,max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,xlab="",ylab="" ,axes=FALSE,col="black",
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):Nenviron])
  if(Nenviron>5) points(rep(0,max.pch2),((max.pch):(max.pch-max.pch2+1))*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):Nenviron])
  text(rep((0+text.di),max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,labels=environ_name[((kk-1)*max.row+1):Nenviron],pos=4,cex=cex.Ne-c.t.d)
  }else{
  # print(kk)
  max.pch2=Nenviron-(n.col.pch-1)*max.row
  
  plot(rep((kk-1)*x.di,max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,xlab="",ylab="" ,axes=FALSE,col="black",
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):Nenviron])
  if(Nenviron>5) points(rep((kk-1)*x.di,max.pch2),((max.pch):(max.pch-max.pch2+1))*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):Nenviron])
  text(rep(((kk-1)*x.di+text.di),max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,labels=environ_name[((kk-1)*max.row+1):Nenviron],pos=4,cex=cex.Ne-c.t.d)
    
  }
 }else{
  if(kk==1)
  {
  print(kk)
  plot(rep(0,max.pch),(max.pch:1)*cex.betw,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  if(Nenviron>5) points(rep(0,max.pch),((max.pch):1)*cex.betw,
  xlim=c(1,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  text(rep((0+text.di),max.pch),(max.pch:1)*cex.betw,labels=environ_name[((kk-1)*max.row+1):((kk-1)*max.row+max.row)],pos=4,cex=cex.Ne-c.t.d)
  }else{
  print(kk)
  plot(rep((kk-1)*x.di,max.pch),(max.pch:1)*cex.betw,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  if(Nenviron>5) points(rep((kk-1)*x.di,max.pch),((max.pch):1)*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  text(rep(((kk-1)*x.di+text.di),max.pch),(max.pch:1)*cex.betw,labels=environ_name[((kk-1)*max.row+1):((kk-1)*max.row+max.row)],pos=4,cex=cex.Ne-c.t.d)
   
  }

 }
}#end of plot.type
 dev.off()

}
print("GAPIT.Association.Manhattans has done !!!")
return(list(multip_mapP=result0,xz=new_xz))
} #end of GAPIT.Manhattan
#=============================================================================================
GPP.Circle.Manhattan.Plot <-function(model_store,Y.names=NULL,environ_name=NULL,QTN.position=NULL,cutOff=0.05,byTraits=FALSE,plot.style="MarineBreeze")
    {
if(byTraits)
{
  for(i in 1:length(Y.names))
  {
     for(j in 1:length(model_store))
     {
        environ_name=c(environ_name,paste(model_store[j],".",Y.names[i],sep=""))
     }
  }
}else{
  for(i in 1:length(model_store))
  {
     for(j in 1:length(Y.names))
     {
        environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
     }
  }
}
.local_extract <- function(path) {
  df <- tryCatch(read.csv(path, head = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  cn <- colnames(df)
  cn <- trimws(cn)
  cn <- sub("^\ufeff", "", cn)
  cn <- sub("^ï\\.\\.", "", cn)
  cn <- gsub("[[:space:][:cntrl:]]+", "", cn)
  cn_norm <- gsub("[^a-z0-9]+", "_", tolower(cn))
  cn_norm <- gsub("_+", "_", cn_norm)
  cn_norm <- gsub("^_+|_+$", "", cn_norm)
  pick_col <- function(norm_names, candidates) {
    idx <- which(norm_names %in% candidates)
    if (length(idx) >= 1) return(idx[1])
    integer(0)
  }
  to_num <- function(x) suppressWarnings(as.numeric(x))
  clean_chr <- function(ch) {
    ch <- as.character(ch)
    ch <- trimws(ch)
    ch <- sub("^\ufeff", "", ch)
    ch <- sub("^ï\\.\\.", "", ch)
    ch <- gsub("[[:space:][:cntrl:]]+", "", ch)
    num_ch <- suppressWarnings(as.numeric(ch))
    not_missing = !is.na(ch) & nzchar(ch)
    has_non_numeric_label = any(not_missing & is.na(num_ch))
    if (has_non_numeric_label) {
      return(ch)
    }
    all_int = if (any(not_missing)) {
      all(is.finite(num_ch[not_missing]) & (num_ch[not_missing] == round(num_ch[not_missing])))
    } else TRUE
    if (all_int) {
      return(num_ch)
    }
    ch
  }
  snp_idx <- pick_col(cn_norm, c("snp", "rs", "rsid", "marker", "id"))
  chr_idx <- pick_col(cn_norm, c("chr", "chrom", "chromosome"))
  pos_idx <- pick_col(cn_norm, c("pos", "position", "bp", "bp_position"))
  if (length(snp_idx) < 1) snp_idx <- 1
  if (length(chr_idx) < 1) chr_idx <- 2
  if (length(pos_idx) < 1) pos_idx <- 3
  data.frame(
    SNP = as.character(df[[snp_idx]]),
    Chr = clean_chr(df[[chr_idx]]),
    Pos = to_num(df[[pos_idx]]),
    stringsAsFactors = FALSE
  )
}
gm_list <- list()
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  fpath <- paste("GAPIT.Association.GWAS_Results.",environ_name[i],".csv",sep="")
  df <- .local_extract(fpath)
  if (!is.null(df) && nrow(df) > 0) gm_list[[length(gm_list)+1]] <- df
}
if (length(gm_list) == 0) {
  stop("No valid GWAS result files found for GM construction.")
}
all_gm <- do.call(rbind, gm_list)
all_gm <- all_gm[!duplicated(all_gm$SNP), , drop=FALSE]
all_gm <- all_gm[!is.na(all_gm$SNP) & all_gm$SNP != "", , drop=FALSE]
all_gm <- all_gm[!is.na(all_gm$Pos), , drop=FALSE]
all_gm <- all_gm[order(all_gm$Pos), , drop=FALSE]
.gm_map <- (function() {
  chr_vec_c <- as.character(all_gm$Chr)
  options(warn = -1)
  numeric.chr <- suppressWarnings(as.numeric(chr_vec_c))
  options(warn = 0)
  max.chr <- suppressWarnings(max(numeric.chr, na.rm = TRUE))
  if (!is.finite(max.chr)) max.chr <- 0
  map.xy.index <- which(!(numeric.chr %in% c(0:max.chr)))
  chr.xy <- character(0)
  if (length(map.xy.index) != 0) {
    chr.xy <- unique(chr_vec_c[map.xy.index])
    for (k in seq_along(chr.xy)) {
      chr_vec_c[chr_vec_c == chr.xy[k]] <- max.chr + k
    }
  }
  list(
    mapped = suppressWarnings(as.numeric(chr_vec_c)),
    chr.xy = chr.xy,
    max.chr = max.chr
  )
})()
ord1 <- order(.gm_map$mapped, all_gm$Pos, na.last = NA)
all_gm <- all_gm[ord1, , drop=FALSE]
GM <- all_gm
GMM=GAPIT.Multiple.Manhattan(model_store=model_store,Y.names=Y.names,GM=GM,seqQTN=QTN.position,cutOff=cutOff,plot.type=c("w","h","s"))
GPP.Circle.Manhattan.Plot.Core(Pmap=GMM$multip_mapP,band=1,r=3,plot.type=c("c","q"),signal.line=1,xz=GMM$xz,threshold=cutOff,plot.style=plot.style)
return(GMM)
}

if (!exists("GAPIT.Circle.Manhattan.Plot", mode = "function") && exists("GPP.Circle.Manhattan.Plot", mode = "function")) {
  GAPIT.Circle.Manhattan.Plot <- GPP.Circle.Manhattan.Plot
}


circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
	{
		graphics::curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
		graphics::curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
	}
Densitplot <- function(
		map,
		col=c("darkblue", "white", "red"),
		main="SNP Density",
		bin=1e6,
		band=3,
		width=5,
		legend.len=10,
		legend.max=NULL,
		legend.pt.cex=3,
		legend.cex=1,
		legend.y.intersp=1,
		legend.x.intersp=1,
		plot=TRUE
	)
	{   #print(head(map))
		map <- as.matrix(map)
		map <- map[!is.na(map[, 2]), ]
		map <- map[!is.na(map[, 3]), ]
		map <- map[map[, 2] != 0, ]
		#map <- map[map[, 3] != 0, ]
		options(warn = -1)
		max.chr <- max(as.numeric(map[, 2]), na.rm=TRUE)
		if(is.infinite(max.chr))	max.chr <- 0
		map.xy.index <- which(!as.numeric(map[, 2]) %in% c(0 : max.chr))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(map[map.xy.index, 2])
			for(i in 1:length(chr.xy)){
				map[map[, 2] == chr.xy[i], 2] <- max.chr + i
			}
		}
		map <- map[order(as.numeric(map[, 2]), as.numeric(map[, 3])), ]
		chr <- as.numeric(map[, 2])
		pos <- as.numeric(map[, 3])
		chr.num <- unique(chr)
		#print(chr.num)
		chorm.maxlen <- max(pos)
		if(plot)	plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main,axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
		pos.x <- list()
		col.index <- list()
		maxbin.num <- NULL
		#print(chr.num)
		for(i in 1 : length(chr.num)){
			pos.x[[i]] <- pos[which(chr == chr.num[i])]
			cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
			if(cut.len <= 1){
				col.index[[i]] = 1
			}else{
				cut.r <- cut(pos.x[[i]], cut.len, labels=FALSE)
				eachbin.num <- table(cut.r)
		        #print(eachbin.num)

				maxbin.num <- c(maxbin.num, max(eachbin.num))
				col.index[[i]] <- rep(eachbin.num, eachbin.num)
			}
		}

		Maxbin.num <- max(maxbin.num)
		maxbin.num <- Maxbin.num
		if(!is.null(legend.max)){
			maxbin.num <- legend.max
		}
		#print(col)
		#print(maxbin.num)
		col = grDevices::colorRampPalette(col)(maxbin.num)
		col.seg=NULL
		for(i in 1 : length(chr.num)){
			if(plot)	graphics::polygon(c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])), 
				c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
				width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey", border="grey")
			if(!is.null(legend.max)){
				if(legend.max < Maxbin.num){
					col.index[[i]][col.index[[i]] > legend.max] <- legend.max
				}
			}
			col.seg <- c(col.seg, col[round(col.index[[i]] * length(col) / maxbin.num)])
			if(plot)	graphics::segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
			col=col[round(col.index[[i]] * length(col) / maxbin.num)], lwd=1)
		}
		if(length(map.xy.index) != 0){
			for(i in 1:length(chr.xy)){
				chr.num[chr.num == max.chr + i] <- chr.xy[i]
			}
		}
		chr.num <- rev(chr.num)
		if(plot)	graphics::mtext(at=seq(band, length(chr.num) * band, band),text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=0.6, line=0.2)
		if(plot)	graphics::axis(3, at=seq(0, chorm.maxlen, length=10), labels=c(NA, paste(round((seq(0, chorm.maxlen, length=10))[-1] / 1e6, 0), "Mb", sep="")),
			font=1, cex.axis=0.8, tck=0.01, lwd=2, padj=1.2)
		# image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
		# round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
		# t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
		legend.y <- round(seq(0, maxbin.num, length=legend.len))
		len <- legend.y[2]
		legend.y <- seq(0, maxbin.num, len)
		if(!is.null(legend.max)){
			if(legend.max < Maxbin.num){
				if(!maxbin.num %in% legend.y){
					legend.y <- c(legend.y, paste(">=", maxbin.num, sep=""))
					legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
				}else{
					legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
					legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
				}
			}else{
				if(!maxbin.num %in% legend.y){
					legend.y <- c(legend.y, maxbin.num)
				}
				legend.y.col <- c(legend.y[-1])
			}
		}else{
			if(!maxbin.num %in% legend.y){
				legend.y <- c(legend.y, paste(">", max(legend.y), sep=""))
				legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
			}else{
				legend.y.col <- c(legend.y[-1])
			}
		}
		legend.y.col <- as.numeric(legend.y.col)
		legend.col <- c("grey", col[round(legend.y.col * length(col) / maxbin.num)])
		if(plot)	graphics::legend(x=(chorm.maxlen + chorm.maxlen/100), y=( -width/2.5 - band * (length(chr.num) - length(chr.num) - 1)), title="", legend=legend.y, pch=15, pt.cex = legend.pt.cex, col=legend.col,
			cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
		if(!plot)	return(list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y))
	}

GPP.Circle.Manhattan.Plot.Core <- function(
	Pmap,
	col=c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
	#col=c("darkgreen", "darkblue", "darkyellow", "darkred"),
	plot.style = c("Oceanic", "Rainbow", "FarmCPU", "Rushville", "Congress", "Ocean", "PLINK", "Beach", "MarineBreeze", "cougars"),
	
	bin.size=1e6,
	bin.max=NULL,
	pch=19,
	band=1,
	cir.band=0.5,
	H=1.5,
	ylim=NULL,
	cex.axis=1,
	plot.type="c",
	multracks=TRUE,
	cex=c(0.5,0.8,1),
	r=0.3,
	xlab="Chromosome",
	ylab=expression(-log[10](italic(p))),
	xaxs="i",
	yaxs="r",
	outward=TRUE,
	threshold = 0.01, 
	threshold.col="red",
	threshold.lwd=1,
	threshold.lty=2,
	amplify= TRUE,     # is that available for remark signal pch col
	chr.labels=NULL,
	signal.cex = 2,
	signal.pch = 8,
	signal.col="red",
	signal.line=NULL,
	cir.chr=TRUE,
	cir.chr.h=1.3,
	chr.den.col=c("darkgreen", "yellow", "red"),
	#chr.den.col=c(126,177,153),
	cir.legend=TRUE,
	cir.legend.cex=0.8,
	cir.legend.col="grey45",
	LOG10=TRUE,
	box=FALSE,
	conf.int.col="grey",
	file.output=TRUE,
	file="pdf",
	dpi=300,
	xz=NULL,
	memo="",
	dummy.chr.den.as.bg=TRUE
)
{		#print("Starting Circular-Manhattan plot!",quote=F)
	taxa=colnames(Pmap)[-c(1:3)]
	# --- Clean taxa labels: remove .x/.y/P.value/P.value.x/P.value.y artifacts ---
	if (length(taxa) > 0) {
	  bad_pattern <- "^(P\\.value|pvalue|p_value|pval)(\\.x|\\.y)?$"
	  # Remove any name that looks like a pure p-value column
	  keep_idx <- which(!grepl(bad_pattern, taxa, ignore.case = TRUE))
	  if (length(keep_idx) > 0) {
	    taxa <- taxa[keep_idx]
	    # Also strip trailing .x/.y just in case
	    taxa <- sub("\\.x$", "", taxa, ignore.case = TRUE)
	    taxa <- sub("\\.y$", "", taxa, ignore.case = TRUE)
	  }
	  # If after filtering there are no good names, fall back
	  if (length(taxa) == 0) taxa <- paste0("Trait_", seq_len(length(colnames(Pmap)[-c(1:3)])))
	}
	if(!is.null(memo) && memo != "")	memo <- paste("_", memo, sep="")
	if(length(taxa) == 0)	taxa <- "Index"
	taxa <- paste(taxa, memo, sep="")

	# --- Apply plot.style chromosome colors (MarineBreeze etc) ---
	if (missing(plot.style) || is.null(plot.style) || length(plot.style) > 1) {
		plot.style <- plot.style[1]
	}
    {
        options(warn = -1)
        numeric.chr_p <- suppressWarnings(as.numeric(as.character(Pmap[, 1])))
        options(warn = 0)
        max.chr_p <- suppressWarnings(max(numeric.chr_p, na.rm = TRUE))
        if (!is.finite(max.chr_p)) max.chr_p <- 0
        map.xy.index_p <- which(!(numeric.chr_p %in% c(0:max.chr_p)))
        chr_p <- as.character(Pmap[, 1])
        if (length(map.xy.index_p) != 0) {
            chr.xy_p <- unique(chr_p[map.xy.index_p])
            for (i in seq_along(chr.xy_p)) {
                chr_p[chr_p == chr.xy_p[i]] <- max.chr_p + i
            }
        }
        chr_p_numeric <- suppressWarnings(as.numeric(chr_p))
    }
	nchr_p <- length(unique(chr_p_numeric))
	col_Rainbow <- grDevices::rainbow(max(2, nchr_p + 1))
	col_FarmCPU <- rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"), ceiling(nchr_p/5))
	col_Rushville <- rep(c("orangered","navyblue"), ceiling(nchr_p/2))
	col_Congress <- rep(c("deepskyblue3","firebrick"), ceiling(nchr_p/2))
	col_Ocean <- rep(c("steelblue4","cyan3"), ceiling(nchr_p/2))
	col_PLINK <- rep(c("gray10","gray70"), ceiling(nchr_p/2))
	col_Beach <- rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"), ceiling(nchr_p/5))
	col_Oceanic <- rep(c('#EC5f67','#FAC863','#99C794','#6699CC','#C594C5'), ceiling(nchr_p/5))
	col_MarineBreeze <- rep(c("#BFDFD2","#51999F","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A"), ceiling(nchr_p/8))
	col_cougars <- rep(c('#990000','dimgray'), ceiling(nchr_p/2))
	plot_color <- NULL
	if (plot.style == "Rainbow") plot_color <- col_Rainbow
	if (plot.style == "FarmCPU") plot_color <- col_FarmCPU
	if (plot.style == "Rushville") plot_color <- col_Rushville
	if (plot.style == "Congress") plot_color <- col_Congress
	if (plot.style == "Ocean") plot_color <- col_Ocean
	if (plot.style == "PLINK") plot_color <- col_PLINK
	if (plot.style == "Beach") plot_color <- col_Beach
	if (plot.style == "Oceanic") plot_color <- col_Oceanic
	if (plot.style == "MarineBreeze") plot_color <- col_MarineBreeze
	if (plot.style == "cougars") plot_color <- col_cougars
	if (!is.null(plot_color)) {
		col <- plot_color
	} else {
		# Default legacy colors (one track)
		col <- rep(c('#FF6A6A','#FAC863','#99C794','#6699CC','#C594C5'), ceiling(length(taxa)/5))
	}
	# --- End plot.style colors ---

    legend.bit=round(nrow(Pmap)/30)

    numeric.chr <- as.numeric(Pmap[, 1])
	options(warn = 0)
	max.chr <- max(numeric.chr, na.rm=TRUE)
    aa=Pmap[1:legend.bit,]
    aa[,2]=max.chr+1
    #print(aa[,3])
    aa[,3]=sample(1:10^7.5,legend.bit)
    aa[,-c(1:3)]=0
    Pmap=rbind(Pmap,aa)
    #print(unique(Pmap[,2]))
	#SNP-Density plot
	if("d" %in% plot.type){
		print("SNP_Density Plotting...")
		if(file.output){
			if(file=="jpg")	grDevices::jpeg(paste("SNP_Density.",paste(taxa,collapse="."),".jpg",sep=""), width = 9*dpi,height=7*dpi,res=dpi,quality = 100)
			if(file=="pdf")	grDevices::pdf(paste("GAPIT.Association.SNP_Density", taxa,".pdf" ,sep=""), width = 9,height=7)
			if(file=="tiff")	grDevices::tiff(paste("SNP_Density.",paste(taxa,collapse="."),".tiff",sep=""), width = 9*dpi,height=7*dpi,res=dpi)
			graphics::par(xpd=TRUE)
		}else{
			if(is.null(grDevices::dev.list()))	grDevices::dev.new(width = 9,height=7)
			graphics::par(xpd=TRUE)
		}

		Densitplot(map=Pmap[,c(1:3)], col=col, bin=bin.size, legend.max=bin.max, main=paste("The number of SNPs within ", bin.size/1e6, "Mb window size", sep=""))
		if(file.output)	grDevices::dev.off()
	}

	if(length(plot.type) !=1 | (!"d" %in% plot.type)){
	
		#order Pmap by the name of SNP
		#Pmap=Pmap[order(Pmap[,1]),]
		Pmap <- as.matrix(Pmap)

		#delete the column of SNPs names
		Pmap <- Pmap[,-1]
		Pmap[is.na(Pmap)]=1
		#print(dim(Pmap))

		#scale and adjust the parameters
		cir.chr.h <- cir.chr.h/5
		cir.band <- cir.band/5
		threshold=threshold/nrow(Pmap)
		if(!is.null(threshold)){
			threshold.col <- rep(threshold.col,length(threshold))
			threshold.lwd <- rep(threshold.lwd,length(threshold))
			threshold.lty <- rep(threshold.lty,length(threshold))
			signal.col <- rep(signal.col,length(threshold))
			signal.pch <- rep(signal.pch,length(threshold))
			signal.cex <- rep(signal.cex,length(threshold))
		}
		if(length(cex)!=3) cex <- rep(cex,3)
		if(!is.null(ylim)){
			if(length(ylim)==1) ylim <- c(0,ylim)
		}
		
		if(is.null(conf.int.col))	conf.int.col <- NA
		if(is.na(conf.int.col)){
			conf.int=FALSE
		}else{
			conf.int=TRUE
		}

		#get the number of traits
		R=ncol(Pmap)-2

		#replace the non-euchromosome
		options(warn = -1)
		numeric.chr <- as.numeric(Pmap[, 1])
		options(warn = 0)
		max.chr <- max(numeric.chr, na.rm=TRUE)
		if(is.infinite(max.chr))	max.chr <- 0
		map.xy.index <- which(!numeric.chr %in% c(0:max.chr))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(Pmap[map.xy.index, 1])
			for(i in 1:length(chr.xy)){
				Pmap[Pmap[, 1] == chr.xy[i], 1] <- max.chr + i
			}
		}

		Pmap <- matrix(as.numeric(Pmap), nrow(Pmap))

		#order the GWAS results by chromosome and position
		Pmap <- Pmap[order(Pmap[, 1], Pmap[,2]), ]

		dummy.row.index <- integer(0)
		dummy.chr.value <- NA
		if(dummy.chr.den.as.bg && ncol(Pmap) > 2){
			zero.row.index <- which(rowSums(Pmap[,-c(1:2)] == 0) == (ncol(Pmap)-2))
			if(length(zero.row.index) > 0){
				dummy.chr.value <- as.numeric(names(sort(table(Pmap[zero.row.index, 1]), decreasing=TRUE))[1])
				dummy.row.index <- zero.row.index[Pmap[zero.row.index, 1] == dummy.chr.value]
			}
		}

		#get the index of chromosome
		chr <- unique(Pmap[,1])
		chr.ori <- chr
		if(length(map.xy.index) != 0){
			for(i in 1:length(chr.xy)){
				chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
			}
		}

		pvalueT <- as.matrix(Pmap[,-c(1:2)])
		#print(dim(pvalueT))
		pvalue.pos <- Pmap[, 2]
		p0.index <- Pmap[, 1] == 0
		if(sum(p0.index) != 0){
			pvalue.pos[p0.index] <- 1:sum(p0.index)
		}
		pvalue.pos.list <- tapply(pvalue.pos, Pmap[, 1], list)
		
		#scale the space parameter between chromosomes
		if(!missing(band)){
			band <- floor(band*(sum(sapply(pvalue.pos.list, max))/100))
		}else{
			band <- floor((sum(sapply(pvalue.pos.list, max))/100))
		}
		if(band==0)	band=1
		
		if(LOG10){
			pvalueT[pvalueT <= 0] <- 1
			pvalueT[pvalueT > 1] <- 1
		}

		#set the colors for the plot
		#palette(heat.colors(1024)) #(heatmap)
		#T=floor(1024/max(pvalue))
		#plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
		
		#print(col)
		if(is.vector(col)){
			col <- matrix(col,R,length(col),byrow=TRUE)
		}
		if(is.matrix(col)){
			#try to transform the colors into matrix for all traits
			col <- matrix(as.vector(t(col)),R,dim(col)[2],byrow=TRUE)
		}

		Num <- as.numeric(table(Pmap[,1]))
		Nchr <- length(Num)
		N <- NULL
		#print(Nchr)
		#set the colors for each traits
		for(i in 1:R){
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			N[i] <- ceiling(Nchr/length(colx))
		}
		
		#insert the space into chromosomes and return the midpoint of each chromosome
		ticks <- NULL
		pvalue.posN <- NULL
		#pvalue <- pvalueT[,j]
		for(i in 0:(Nchr-1)){
			if (i==0){
				#pvalue <- append(pvalue,rep(Inf,band),after=0)
				pvalue.posN <- pvalue.pos.list[[i+1]] + band
				ticks[i+1] <- max(pvalue.posN)-floor(max(pvalue.pos.list[[i+1]])/2)
			}else{
				#pvalue <- append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
				pvalue.posN <- c(pvalue.posN, max(pvalue.posN) + band + pvalue.pos.list[[i+1]])
				ticks[i+1] <- max(pvalue.posN)-floor(max(pvalue.pos.list[[i+1]])/2)
			}
		}
		pvalue.posN.list <- tapply(pvalue.posN, Pmap[, 1], list)
		#NewP[[j]] <- pvalue
		
		#merge the pvalues of traits by column
		if(LOG10){
			logpvalueT <- -log10(pvalueT)
		}else{
			pvalueT <- abs(pvalueT)
			logpvalueT <- pvalueT
		}

		add <- list()
		for(i in 1:R){
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			add[[i]] <- c(Num,rep(0,N[i]*length(colx)-Nchr))
		}

		TotalN <- max(pvalue.posN)

		if(length(chr.den.col) > 1){
			cir.density=TRUE
			den.fold <- 20
			density.list <- Densitplot(map=Pmap[,c(1,1,2)], col=chr.den.col, plot=FALSE, bin=bin.size, legend.max=bin.max)
			if(dummy.chr.den.as.bg && length(dummy.row.index) > 0 && length(density.list$den.col) > 0){
				bg_col <- graphics::par("bg")
				if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
				dummy.row.index <- dummy.row.index[dummy.row.index >= 1 & dummy.row.index <= length(density.list$den.col)]
				if(length(dummy.row.index) > 0)	density.list$den.col[dummy.row.index] <- bg_col
			}
			#list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y)
		}else{
			cir.density=FALSE
		}


        #print(dim(pvalueT))

		if(is.null(xz)){
		signal.line.index <- NULL
		if(!is.null(threshold)){
			if(!is.null(signal.line)){
				for(l in 1:R){
					if(LOG10){
						signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min(threshold)))
					}else{
						signal.line.index <- c(signal.line.index,which(pvalueT[,l] > max(threshold)))
					}
				}
				signal.line.index <- unique(signal.line.index)
			}
		}
		signal.lty=rep(2,length(signal.line.index))
	    }else{
        signal.line.index=as.numeric(as.vector(xz[,1]))
        signal.lty=as.numeric(as.vector(xz[,2]))
	    }#end is.null(xz)
        
		signal.line.index <- pvalue.posN[signal.line.index]
	}
	    


    if("c" %in% plot.type)
    {
		if(file.output){
			if(file=="jpg")	grDevices::jpeg(paste("GAPIT.Manhattan.Multiple.Plot.circular.jpg",sep=""), width = 8*dpi,height=8*dpi,res=dpi,quality = 100)
			if(file=="pdf")	grDevices::pdf(paste("GAPIT.Association.Manhattans_Circular.pdf" ,sep=""), width = 10,height=10)
			if(file=="tiff")	grDevices::tiff(paste("GAPIT.Manhattan.Multiple.Plot.circular.tiff",sep=""), width = 8*dpi,height=8*dpi,res=dpi)
		}
		if(!file.output){
			if(!is.null(grDevices::dev.list()))	grDevices::dev.new(width=8, height=8)
			graphics::par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
		}
		graphics::par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
		RR <- r+H*R+cir.band*R
		if(cir.density){
			plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
		}else{
			plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
		}
		if(!is.null(signal.line)){
			if(!is.null(signal.line.index)){
				X1chr <- (RR)*sin(2*pi*(signal.line.index-round(band/2))/TotalN)
				Y1chr <- (RR)*cos(2*pi*(signal.line.index-round(band/2))/TotalN)
				X2chr <- (r)*sin(2*pi*(signal.line.index-round(band/2))/TotalN)
				Y2chr <- (r)*cos(2*pi*(signal.line.index-round(band/2))/TotalN)
				#print(signal.line)

				#print(dim(pvalueT))
				#print(head(pvalueT))
				#print(dim(xz))
				#print(xz)
				#print(head(pvalue.posN))
				graphics::segments(X1chr,Y1chr,X2chr,Y2chr,lty=signal.lty,lwd=signal.line,col="grey")
			}
		}
		for(i in 1:R){
		
			#get the colors for each trait
			colx <- col[i,]
			colx <- colx[!is.na(colx)]
			
			#debug
			#print(colx)
			
			#print(paste("Circular_Manhattan Plotting ",taxa[i],"...",sep=""))
			pvalue <- pvalueT[,i]
			logpvalue <- logpvalueT[,i]
			if(is.null(ylim)){
				if(LOG10){
					Max <- ceiling(-log10(min(pvalue[pvalue!=0])))
				}else{
					Max <- ceiling(max(pvalue[pvalue!=Inf]))
					if(Max<=1)
					Max <- max(pvalue[pvalue!=Inf])
				}
			}else{
				Max <- ylim[2]
			}
			Cpvalue <- (H*logpvalue/Max)

			if(outward==TRUE){
				if(cir.chr==TRUE){
					
					#plot the boundary which represents the chromosomes
					polygon.num <- 1000
					#print(length(chr))
					for(k in 1:length(chr)){
						if(k==1){
							polygon.index <- seq(round(band/2)+1,-round(band/2)+max(pvalue.posN.list[[1]]), length=polygon.num)
							#change the axis from right angle into circle format
							X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)

							#print(length(X1chr))
							if(is.null(chr.den.col)){
								poly_col <- rep(colx,ceiling(length(chr)/length(colx)))[k]
								poly_border <- poly_col
							}else{
								if(cir.density){
									poly_col <- "grey"
									poly_border <- "grey"
								}else{
									poly_col <- chr.den.col
									poly_border <- chr.den.col
								}
							}
							if(dummy.chr.den.as.bg && !is.na(dummy.chr.value) && chr[k] == dummy.chr.value){
								bg_col <- graphics::par("bg")
								if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
								poly_col <- bg_col
								poly_border <- bg_col
							}
							graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=poly_col,border=poly_border)
						}else{
							polygon.index <- seq(1+round(band/2)+max(pvalue.posN.list[[k-1]]),-round(band/2)+max(pvalue.posN.list[[k]]), length=polygon.num)
							X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
							if(is.null(chr.den.col)){
								poly_col <- rep(colx,ceiling(length(chr)/length(colx)))[k]
								poly_border <- poly_col
							}else{
								if(cir.density){
									poly_col <- "grey"
									poly_border <- "grey"
								}else{
									poly_col <- chr.den.col
									poly_border <- chr.den.col
								}
							}
							if(dummy.chr.den.as.bg && !is.na(dummy.chr.value) && chr[k] == dummy.chr.value){
								bg_col <- graphics::par("bg")
								if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
								poly_col <- bg_col
								poly_border <- bg_col
							}
							graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=poly_col,border=poly_border)
						}
					}
					
					if(cir.density){
						if(dummy.chr.den.as.bg && length(dummy.row.index) > 0 && length(density.list$den.col) > 0){
							bg_col <- graphics::par("bg")
							if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
							dummy.row.index <- dummy.row.index[dummy.row.index >= 1 & dummy.row.index <= length(density.list$den.col)]
							if(length(dummy.row.index) > 0)	density.list$den.col[dummy.row.index] <- bg_col
						}

						graphics::segments(
							(RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							col=density.list$den.col, lwd=0.1
						)
						graphics::legend(
							x=RR+4*cir.chr.h,
							y=(RR+4*cir.chr.h)/2,
							horiz=F,
							title="Density", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
							cex=1, bty="n",
							y.intersp=1,
							x.intersp=1,
							yjust=0.5, xjust=0, xpd=TRUE
						)
						
					}
					
					# XLine=(RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)
					if(cir.density){
						circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
						circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
					}else{
						circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
						circle.plot(myr=RR,lwd=1.5,add=TRUE)
					}

				}
				
				#plot the y axis of legend for each trait
				if(cir.legend==TRUE){
					#try to get the number after radix point
					if(Max<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=1
					}
					graphics::segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					graphics::segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					#text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					graphics::text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					graphics::text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					#text(-r/15,r+H*(i-0)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					#text(r/5,0.4*(i-1),taxa[i],adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				    
				}
				X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
				Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
				# plot point in figure
				graphics::points(X[1:(length(X)-legend.bit)],Y[1:(length(Y)-legend.bit)],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				
				# plot significant line
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							significantline1=ifelse(LOG10, H*(-log10(threshold[thr]))/Max, H*(threshold[thr])/Max)
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							# plot significant line
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								#if(thr==length(threshold))circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col="black",lwd=threshold.lwd[thr],lty=threshold.lty[thr])
								#print("!!!!!")
								circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
								#circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col="black",lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
					}
				}
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(amplify==TRUE){
							if(LOG10){
								threshold <- sort(threshold)
								significantline1=H*(-log10(max(threshold)))/Max
							}else{
								threshold <- sort(threshold, decreasing=TRUE)
								significantline1=H*(min(threshold))/Max
							}
							
							p_amp.index <- which(Cpvalue>=significantline1)
							HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							graphics::points(HX1,HY1,pch=19,cex=cex[1],col="white")
							
								for(ll in 1:length(threshold)){
									if(ll == 1){
										if(LOG10){
											significantline1=H*(-log10(threshold[ll]))/Max
										}else{
											significantline1=H*(threshold[ll])/Max
										}
										p_amp.index <- which(Cpvalue>=significantline1)
										HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									}else{
										if(LOG10){
											significantline0=H*(-log10(threshold[ll-1]))/Max
											significantline1=H*(-log10(threshold[ll]))/Max
										}else{
											significantline0=H*(threshold[ll-1])/Max
											significantline1=H*(threshold[ll])/Max
										}
										p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
										HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									}
								
									if(is.null(signal.col)){
										# print(signal.pch)
										graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
									}else{
										# print(signal.pch)
										graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=signal.col[ll])
									}
								}
						}
					}
				}
				if(cir.chr==TRUE){
					ticks1=1.07*(RR+cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.07*(RR+cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)
					label.index <- seq_along(ticks)
					if(dummy.chr.den.as.bg && !is.na(dummy.chr.value))	label.index <- label.index[chr[label.index] != dummy.chr.value]
					if(is.null(chr.labels)){
						#print(length(ticks))
						for(i in label.index){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in label.index){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}else{
					ticks1=(0.9*r)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=(0.9*r)*cos(2*pi*(ticks-round(band/2))/TotalN)
					label.index <- seq_along(ticks)
					if(dummy.chr.den.as.bg && !is.na(dummy.chr.value))	label.index <- label.index[chr[label.index] != dummy.chr.value]
					if(is.null(chr.labels)){
						for(i in label.index){
						angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in label.index){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}
			}
			if(outward==FALSE){
				if(cir.chr==TRUE){
					# XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)

					polygon.num <- 1000
					for(k in 1:length(chr)){
						if(k==1){
							polygon.index <- seq(round(band/2)+1,-round(band/2)+max(pvalue.posN.list[[1]]), length=polygon.num)
							X1chr=(2*cir.band+RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(2*cir.band+RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
							if(is.null(chr.den.col)){
								poly_col <- rep(colx,ceiling(length(chr)/length(colx)))[k]
								poly_border <- poly_col
							}else{
								if(cir.density){
									poly_col <- "grey"
									poly_border <- "grey"
								}else{
									poly_col <- chr.den.col
									poly_border <- chr.den.col
								}
							}
							if(dummy.chr.den.as.bg && !is.na(dummy.chr.value) && chr[k] == dummy.chr.value){
								bg_col <- graphics::par("bg")
								if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
								poly_col <- bg_col
								poly_border <- bg_col
							}
							graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=poly_col,border=poly_border)
						}else{
							polygon.index <- seq(1+round(band/2)+max(pvalue.posN.list[[k-1]]),-round(band/2)+max(pvalue.posN.list[[k]]), length=polygon.num)
							X1chr=(2*cir.band+RR)*sin(2*pi*(polygon.index)/TotalN)
							Y1chr=(2*cir.band+RR)*cos(2*pi*(polygon.index)/TotalN)
							X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
							Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
							if(is.null(chr.den.col)){
								poly_col <- rep(colx,ceiling(length(chr)/length(colx)))[k]
								poly_border <- poly_col
							}else{
								if(cir.density){
									poly_col <- "grey"
									poly_border <- "grey"
								}else{
									poly_col <- chr.den.col
									poly_border <- chr.den.col
								}
							}
							if(dummy.chr.den.as.bg && !is.na(dummy.chr.value) && chr[k] == dummy.chr.value){
								bg_col <- graphics::par("bg")
								if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
								poly_col <- bg_col
								poly_border <- bg_col
							}
							graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=poly_col,border=poly_border)
						}
					}
					if(cir.density){
						if(dummy.chr.den.as.bg && length(dummy.row.index) > 0 && length(density.list$den.col) > 0){
							bg_col <- graphics::par("bg")
							if(is.null(bg_col) || is.na(bg_col) || bg_col == "transparent")	bg_col <- "white"
							dummy.row.index <- dummy.row.index[dummy.row.index >= 1 & dummy.row.index <= length(density.list$den.col)]
							if(length(dummy.row.index) > 0)	density.list$den.col[dummy.row.index] <- bg_col
						}

						graphics::segments(
							(2*cir.band+RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(2*cir.band+RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(2*cir.band+RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
							(2*cir.band+RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
							col=density.list$den.col, lwd=0.1
						)
						graphics::legend(
							x=RR+4*cir.chr.h,
							y=(RR+4*cir.chr.h)/2,
							title="Density", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
							cex=1, bty="n",
							y.intersp=1,
							x.intersp=1,
							yjust=0.5, xjust=0, xpd=TRUE
						)
						
					}
					
					if(cir.density){
						circle.plot(myr=2*cir.band+RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
						circle.plot(myr=2*cir.band+RR,lwd=1.5,add=TRUE,col='grey')
					}else{
						circle.plot(myr=2*cir.band+RR+cir.chr.h,lwd=1.5,add=TRUE)
						circle.plot(myr=2*cir.band+RR,lwd=1.5,add=TRUE)
					}

				}


				if(cir.legend==TRUE){
					
					#try to get the number after radix point
					if(Max<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=2
					}
					graphics::segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					graphics::segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
					graphics::text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					#text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					graphics::text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					#text(-r/15,r+H*(i-1)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				    #text(r,0.4*(i-1),taxa[i],adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)

				}
				
				X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
				Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
				#points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				graphics::points(X[1:(length(X)-legend.bit)],Y[1:(length(Y)-legend.bit)],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
					
						for(thr in 1:length(threshold)){
							significantline1=ifelse(LOG10, H*(-log10(threshold[thr]))/Max, H*(threshold[thr])/Max)
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
						if(amplify==TRUE){
							if(LOG10){
								threshold <- sort(threshold)
								significantline1=H*(-log10(max(threshold)))/Max
							}else{
								threshold <- sort(threshold, decreasing=TRUE)
								significantline1=H*(min(threshold))/Max
							}
							p_amp.index <- which(Cpvalue>=significantline1)
							HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							graphics::points(HX1,HY1,pch=19,cex=cex[1],col="white")
							
								for(ll in 1:length(threshold)){
									if(ll == 1){
										if(LOG10){
											significantline1=H*(-log10(threshold[ll]))/Max
										}else{
											significantline1=H*(threshold[ll])/Max
										}
										p_amp.index <- which(Cpvalue>=significantline1)
										HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									}else{
										if(LOG10){
											significantline0=H*(-log10(threshold[ll-1]))/Max
											significantline1=H*(-log10(threshold[ll]))/Max
										}else{
											significantline0=H*(threshold[ll-1])/Max
											significantline1=H*(threshold[ll])/Max
										}
										p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
										HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
										HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
									
									}
								
									if(is.null(signal.col)){
										graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
									}else{
										graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=signal.col[ll])
									}
								}
						}
					}
				}
				
				if(cir.chr==TRUE){
					ticks1=1.1*(2*cir.band+RR)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.1*(2*cir.band+RR)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:(length(ticks)-1)){
						  angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						  graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}else{
					ticks1=1.0*(RR+cir.band)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.0*(RR+cir.band)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(chr.labels)){
						for(i in 1:length(ticks)){
						
							#adjust the angle of labels of circle plot
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}	
				}
			}
		}
		taxa=append("Centre",taxa,)
		taxa_col=rep("black",R)
		taxa_col=append("red",taxa_col)
		for(j in 1:(R+1)){
            graphics::text(r/5,0.4*(j-1),taxa[j],adj=1,col=taxa_col[j],cex=cir.legend.cex,font=2)
				    
		}
		taxa=taxa[-1]
		if(file.output) grDevices::dev.off()
	}

	if("q" %in% plot.type){
		#print("Starting QQ-plot!",quote=F)
		amplify=FALSE
		if(multracks){
			if(file.output){
				if(file=="jpg")	grDevices::jpeg(paste("GAPIT.Multracks.QQ.plot.jpg",sep=""), width = R*2.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
				if(file=="pdf")	grDevices::pdf(paste("GAPIT.Association.QQs_Tracks.pdf",sep=""), width = R*2.5,height=5.5)
				if(file=="tiff")	grDevices::tiff(paste("GAPIT.Multracks.QQ.plot.tiff",sep=""), width = R*2.5*dpi,height=5.5*dpi,res=dpi)
				graphics::par(mfcol=c(1,R),mar = c(0,1,4,1.5),oma=c(3,5,0,0),xpd=TRUE)
			}else{
				if(is.null(grDevices::dev.list()))	grDevices::dev.new(width = 2.5*R, height = 5.5)
				graphics::par(xpd=TRUE)
			}
			for(i in 1:R){
				print(paste("Multracks_QQ Plotting ",taxa[i],"...",sep=""))		
				P.values=as.numeric(Pmap[,i+2])
				# --- Robust cleanup for QQ (outcome217 style artifacts) ---
				P.values <- P.values[is.finite(P.values)]
				P.values <- P.values[!is.na(P.values)]
				P.values <- P.values[P.values > 0 & P.values < 1]
				if (length(P.values) < 10) {
				  P.values <- as.numeric(Pmap[,i+2])
				  P.values <- P.values[!is.na(P.values) & is.finite(P.values)]
				  P.values <- pmax(P.values, 1e-300)
				  P.values <- pmin(P.values, 1 - 1e-16)
				}
				if(LOG10){
					P.values=P.values[order(P.values)]
					N=length(P.values)
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
				if(LOG10){
					log.P.values <- -log10(P.values)
				}else{
					log.P.values <- P.values
				}
				# --- Official GAPIT lambda (gapit_functions.txt line 15595):
				#     lambda.estimated = median(P.values) / median(p_value_quantiles)
				# Do NOT de-duplicate (matches classic GAPIT GAPIT.QQ exactly)
				lambda <- NA_real_
				if(length(P.values) > 0){
				  med_p <- stats::median(as.numeric(P.values), na.rm = TRUE)
				  med_q <- stats::median(as.numeric(p_value_quantiles), na.rm = TRUE)
				  if (is.finite(med_p) && is.finite(med_q) && med_q > 0) {
				    lambda <- med_p / med_q
				  }
				}
				
				#calculate the confidence interval of QQ-plot
				if(conf.int){
					N1=length(log.Quantiles)
					c95 <- rep(NA,N1)
					c05 <- rep(NA,N1)
					for(j in 1:N1){
						xi=ceiling((10^-log.Quantiles[j])*N)
						if(xi==0)xi=1
						c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
						c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
					}
					index=length(c95):1
				}else{
					c05 <- 1
					c95 <- 1
				}
				
				# Robust Y limit (avoid top-whisker blowing up)
				ylow_ci <- suppressWarnings(max(-log10(c05), -log10(c95), na.rm = TRUE))
				if (!is.finite(ylow_ci)) ylow_ci <- 0
				yobs_p999 <- stats::quantile(log.P.values, probs = 0.999, na.rm = TRUE)
				yobs_max <- max(log.P.values, na.rm = TRUE)
				YlimMax <- max(floor(ylow_ci + 1), floor(max(yobs_p999, yobs_max) + 1))
				if (!is.finite(YlimMax) || YlimMax < 1) YlimMax <- 10
				plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,YlimMax),xlab ="", ylab="", main = taxa[i])
				if(!is.na(lambda)){
					usr <- graphics::par("usr")
					graphics::text(
						x = usr[1] + 0.05 * (usr[2] - usr[1]),
						y = usr[4] - 0.05 * (usr[4] - usr[3]),
						labels = bquote(lambda == .(format(lambda, digits=6))),
						adj = c(0, 1)
					)
				}
				graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
				graphics::axis(2, at=seq(0,YlimMax,ceiling(YlimMax/10)), labels=seq(0,YlimMax,ceiling(YlimMax/10)), cex.axis=cex.axis)
				
				#plot the confidence interval of QQ-plot
				
				if(conf.int)	graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
				
				if(!is.null(threshold.col)){
				    graphics::par(xpd=FALSE);
				    graphics::abline(a = 0, b = 1, col = threshold.col[1],lwd=2);
				    graphics::par(xpd=TRUE)
				}
				graphics::points(log.Quantiles, log.P.values, col = col[1],pch=1,cex=cex[3])
				#print(max(log.Quantiles))
				#	print(length(log.Quantiles))
				#	print(length(log.P.values))
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						thre.line=-log10(min(threshold))
						if(amplify==TRUE){
							thre.index=which(log.P.values>=thre.line)
							if(length(thre.index)!=0){
							
								#cover the points that exceed the threshold with the color "white"
								graphics::points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
								if(is.null(signal.col)){
									graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch[1],cex=signal.cex[1])
								}else{
									graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
								}
							}
						}
					}
				}
			}
			if(box)	box()
			if(file.output) grDevices::dev.off()
			if(R > 1){
				#qq_col=rainbow(R)
                qq_col=rep(c( '#FF6A6A',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(R/5))

				signal.col <- NULL
				if(file.output){
					if(file=="jpg")	grDevices::jpeg(paste("GAPIT.Multiple.QQ.plot.symphysic.jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
					if(file=="pdf")	grDevices::pdf(paste("GAPIT.Association.QQs_Symphysic.pdf",sep=""), width = 5.5,height=5.5)
					if(file=="tiff")	grDevices::tiff(paste("GAPIT.Multiple.QQ.plot.symphysic.tiff",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi)
					graphics::par(mar = c(5,5,4,2),xpd=TRUE)
				}else{
					grDevices::dev.new(width = 5.5, height = 5.5)
					graphics::par(xpd=TRUE)
				}
				P.values=as.numeric(Pmap[,i+2])
				# --- Robust cleanup for QQ (outcome217 style artifacts) ---
				P.values <- P.values[is.finite(P.values)]
				P.values <- P.values[!is.na(P.values)]
				P.values <- P.values[P.values > 0 & P.values < 1]
				if (length(P.values) < 10) {
				  P.values <- as.numeric(Pmap[,i+2])
				  P.values <- P.values[!is.na(P.values) & is.finite(P.values)]
				  P.values <- pmax(P.values, 1e-300)
				  P.values <- pmin(P.values, 1 - 1e-16)
				}
				if(LOG10){
					N=length(P.values)
					P.values=P.values[order(P.values)]
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
											
				# calculate the confidence interval of QQ-plot
				if(conf.int){
					N1=length(log.Quantiles)
					c95 <- rep(NA,N1)
					c05 <- rep(NA,N1)
					for(j in 1:N1){
						xi=ceiling((10^-log.Quantiles[j])*N)
						if(xi==0)xi=1
						c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
						c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
					}
					index=length(c95):1
				}
				
				if(!conf.int){c05 <- 1; c95 <- 1}
				
				Pmap.min <- Pmap[,3:(R+2)]

				ylow_ci <- suppressWarnings(max(-log10(c05), -log10(c95), na.rm = TRUE))
				if (!is.finite(ylow_ci)) ylow_ci <- 0
				pos_vals <- suppressWarnings(-log10(Pmap.min[Pmap.min > 0 & Pmap.min < 1]))
				if (length(pos_vals) > 0 && is.finite(max(pos_vals, na.rm = TRUE))) {
				  yobs_p999 <- stats::quantile(pos_vals, probs = 0.999, na.rm = TRUE)
				  YlimMax <- max(floor(ylow_ci + 1), floor(max(yobs_p999, max(pos_vals, na.rm = TRUE)) + 1))
				} else {
				  YlimMax <- max(floor(ylow_ci + 1), 10)
				}
				if (!is.finite(YlimMax) || YlimMax < 1) YlimMax <- 10
				plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0, floor(YlimMax+1)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = "QQ plot")
				#legend("topleft",taxa,col=t(col)[1:R],pch=1,pt.lwd=2,text.font=6,box.col=NA)			
				graphics::legend("topleft",taxa,col=qq_col[1:R],pch=1,pt.lwd=3,text.font=6,box.col=NA)
				graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
				graphics::axis(2, at=seq(0,floor(YlimMax+1),ceiling((YlimMax+1)/10)), labels=seq(0,floor((YlimMax+1)),ceiling((YlimMax+1)/10)), cex.axis=cex.axis)
				#print(log.Quantiles[index])
				#print(index)
				#print(length(log.Quantiles))

				# plot the confidence interval of QQ-plot
				if(conf.int)	graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
				
				for(i in 1:R){
					#print(paste("Multraits_QQ Plotting ",taxa[i],"...",sep=""))
					P.values=as.numeric(Pmap[,i+2])
					# --- Robust cleanup for QQ (outcome217 style artifacts) ---
					P.values <- P.values[is.finite(P.values)]
					P.values <- P.values[!is.na(P.values)]
					P.values <- P.values[P.values > 0 & P.values < 1]
					if (length(P.values) < 10) {
					  P.values <- as.numeric(Pmap[,i+2])
					  P.values <- P.values[!is.na(P.values) & is.finite(P.values)]
					  P.values <- pmax(P.values, 1e-300)
					  P.values <- pmin(P.values, 1 - 1e-16)
					}
				    if(LOG10){
					N=length(P.values)
					P.values=P.values[order(P.values)]
				    }else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				    }
				    p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				    log.Quantiles <- -log10(p_value_quantiles)
				    if(LOG10){
					log.P.values <- -log10(P.values)
				    }else{
					log.P.values <- P.values
				    }
				
						
					if((i == 1) & !is.null(threshold.col)){
					    graphics::par(xpd=FALSE);
					    graphics::abline(a = 0, b = 1, col = threshold.col[1],lwd=2);
					    graphics::par(xpd=TRUE)}
					#print(length(log.Quantiles))
				    #print("!!!!!") 
					#points(log.Quantiles, log.P.values, col = t(col)[i],pch=1,lwd=3,cex=cex[3])
					graphics::points(log.Quantiles, log.P.values, col = qq_col[i],pch=1,lwd=3,cex=cex[3])
					
					#print(max(log.Quantiles))
					#
	
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							thre.line=-log10(min(threshold))
							if(amplify==TRUE){
								thre.index=which(log.P.values>=thre.line)
								if(length(thre.index)!=0){
								
									# cover the points that exceed the threshold with the color "white"
									graphics::points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,lwd=3,cex=cex[3])
									if(is.null(signal.col)){
										graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = t(col)[i],pch=signal.pch[1],cex=signal.cex[1])
									}else{
										graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
									}
								}
							}
						}
					}
				}
					box()
				if(file.output) grDevices::dev.off()
			}

			for(iqq in 1:R){
				P.values=as.numeric(Pmap[,iqq+2])
				# --- Robust cleanup for QQ (outcome217 style artifacts) ---
				P.values <- P.values[is.finite(P.values)]
				P.values <- P.values[!is.na(P.values)]
				P.values <- P.values[P.values > 0 & P.values < 1]
				if (length(P.values) < 10) {
				  P.values <- as.numeric(Pmap[,iqq+2])
				  P.values <- P.values[!is.na(P.values) & is.finite(P.values)]
				  P.values <- pmax(P.values, 1e-300)
				  P.values <- pmin(P.values, 1 - 1e-16)
				}
				if(LOG10){
					N=length(P.values)
					P.values=P.values[order(P.values)]
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
				if(LOG10){
					log.P.values <- -log10(P.values)
				}else{
					log.P.values <- P.values
				}

				# --- Official GAPIT lambda (gapit_functions.txt line 15595):
				#     lambda.estimated = median(P.values) / median(p_value_quantiles)
				# Do NOT de-duplicate (matches classic GAPIT GAPIT.QQ exactly)
				lambda <- NA_real_
				if(length(P.values) > 0){
				  med_p <- stats::median(as.numeric(P.values), na.rm = TRUE)
				  med_q <- stats::median(as.numeric(p_value_quantiles), na.rm = TRUE)
				  if (is.finite(med_p) && is.finite(med_q) && med_q > 0) {
				    lambda <- med_p / med_q
				  }
				}

				if(conf.int){
					N1=length(log.Quantiles)
					c95 <- rep(NA,N1)
					c05 <- rep(NA,N1)
					for(j in 1:N1){
						xi=ceiling((10^-log.Quantiles[j])*N)
						if(xi==0)xi=1
						c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
						c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
					}
					index=length(c95):1
				}else{
					c05 <- 1
					c95 <- 1
				}

				ylow_ci <- suppressWarnings(max(-log10(c05), -log10(c95), na.rm = TRUE))
				if (!is.finite(ylow_ci)) ylow_ci <- 0
				yobs_p999 <- stats::quantile(log.P.values, probs = 0.999, na.rm = TRUE)
				yobs_max <- max(log.P.values, na.rm = TRUE)
				YlimMax <- max(floor(ylow_ci + 1), floor(max(yobs_p999, yobs_max) + 1))
				if (!is.finite(YlimMax) || YlimMax < 1) YlimMax <- 10
				if(file.output){
					if(file=="jpg")	grDevices::jpeg(paste0("GAPIT.Association.QQ.",taxa[iqq],".jpg"), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
					if(file=="pdf"){
						.qq_fn <- paste0("GAPIT.Association.QQ.",taxa[iqq],".pdf")
						.qq_ok <- FALSE
						tryCatch({grDevices::pdf(.qq_fn, width = 5.5,height=5.5); .qq_ok <- TRUE}, error=function(e){})
						if(!.qq_ok){
							for(.qq_k in 1:50){
								.qq_fn2 <- sub("\\.pdf$", paste0("-", .qq_k, ".pdf"), .qq_fn)
								.qq_ok2 <- FALSE
								tryCatch({grDevices::pdf(.qq_fn2, width = 5.5,height=5.5); .qq_ok2 <- TRUE}, error=function(e){})
								if(.qq_ok2){.qq_ok <- TRUE; break}
							}
						}
						if(!.qq_ok)	stop("Cannot open QQ pdf device.", call. = FALSE)
					}
					if(file=="tiff")	grDevices::tiff(paste0("GAPIT.Association.QQ.",taxa[iqq],".tiff"), width = 5.5*dpi,height=5.5*dpi,res=dpi)
					graphics::par(mar = c(5,5,4,2),xpd=TRUE)
				}else{
					if(is.null(grDevices::dev.list()))	grDevices::dev.new(width = 5.5, height = 5.5)
					graphics::par(xpd=TRUE)
				}
				plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,YlimMax),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[iqq]))
				graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
				graphics::axis(2, at=seq(0,YlimMax,ceiling(YlimMax/10)), labels=seq(0,YlimMax,ceiling(YlimMax/10)), cex.axis=cex.axis)
				if(conf.int)	graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col="grey85",border="grey85")
				if(!is.null(threshold.col)){
					graphics::par(xpd=FALSE)
					graphics::abline(a = 0, b = 1, col = "red",lwd=2)
					graphics::par(xpd=TRUE)
				}
				if(!is.na(lambda)){
					usr <- graphics::par("usr")
					graphics::text(
						x = usr[1] + 0.05 * (usr[2] - usr[1]),
						y = usr[4] - 0.05 * (usr[4] - usr[3]),
						labels = bquote(lambda == .(format(lambda, digits=6))),
						adj = c(0, 1)
					)
				}
				graphics::points(log.Quantiles, log.P.values, col = "blue",pch=1,cex=cex[3])
				box()
				if(file.output) grDevices::dev.off()
			}
		}else{
			for(i in 1:R){
				print(paste("Q_Q Plotting ",taxa[i],"...",sep=""))
				if(file.output){
					if(file=="jpg")	grDevices::jpeg(paste("QQplot.",taxa[i],".jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
					if(file=="pdf")	grDevices::pdf(paste("GAPIT.Association.QQ.",taxa[i],".pdf",sep=""), width = 5.5,height=5.5)
					if(file=="tiff")	grDevices::tiff(paste("QQplot.",taxa[i],".tiff",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi)
					graphics::par(mar = c(5,5,4,2),xpd=TRUE)
				}else{
					if(is.null(grDevices::dev.list()))	grDevices::dev.new(width = 5.5, height = 5.5)
					graphics::par(xpd=TRUE)
				}
				P.values=as.numeric(Pmap[,i+2])
				# --- Robust cleanup for QQ (outcome217 style artifacts) ---
				P.values <- P.values[is.finite(P.values)]
				P.values <- P.values[!is.na(P.values)]
				P.values <- P.values[P.values > 0 & P.values < 1]
				if (length(P.values) < 10) {
				  P.values <- as.numeric(Pmap[,i+2])
				  P.values <- P.values[!is.na(P.values) & is.finite(P.values)]
				  P.values <- pmax(P.values, 1e-300)
				  P.values <- pmin(P.values, 1 - 1e-16)
				}
				if(LOG10){
					N=length(P.values)
					P.values=P.values[order(P.values)]
				}else{
					N=length(P.values)
					P.values=P.values[order(P.values,decreasing=TRUE)]
				}
				p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
				log.Quantiles <- -log10(p_value_quantiles)
				if(LOG10){
					log.P.values <- -log10(P.values)
				}else{
					log.P.values <- P.values
				}
				# --- Official GAPIT lambda (gapit_functions.txt line 15595):
				#     lambda.estimated = median(P.values) / median(p_value_quantiles)
				# Do NOT de-duplicate (matches classic GAPIT GAPIT.QQ exactly)
				lambda <- NA_real_
				if(length(P.values) > 0){
				  med_p <- stats::median(as.numeric(P.values), na.rm = TRUE)
				  med_q <- stats::median(as.numeric(p_value_quantiles), na.rm = TRUE)
				  if (is.finite(med_p) && is.finite(med_q) && med_q > 0) {
				    lambda <- med_p / med_q
				  }
				}
				
				#calculate the confidence interval of QQ-plot
				if(conf.int){
					N1=length(log.Quantiles)
					c95 <- rep(NA,N1)
					c05 <- rep(NA,N1)
					for(j in 1:N1){
						xi=ceiling((10^-log.Quantiles[j])*N)
						if(xi==0)xi=1
						c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
						c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
					}
					index=length(c95):1
				}else{
					c05 <- 1
					c95 <- 1
				}
				# Robust Y limit
				ylow_ci <- suppressWarnings(max(-log10(c05), -log10(c95), na.rm = TRUE))
				if (!is.finite(ylow_ci)) ylow_ci <- 0
				yobs_p999 <- stats::quantile(log.P.values, probs = 0.999, na.rm = TRUE)
				yobs_max <- max(log.P.values, na.rm = TRUE)
				YlimMax <- max(floor(ylow_ci + 1), floor(max(yobs_p999, yobs_max) + 1))
				if (!is.finite(YlimMax) || YlimMax < 1) YlimMax <- 10
				plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,YlimMax),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[i]))
				if(!is.na(lambda)){
					usr <- graphics::par("usr")
					graphics::text(
						x = usr[1] + 0.05 * (usr[2] - usr[1]),
						y = usr[4] - 0.05 * (usr[4] - usr[3]),
						labels = bquote(lambda == .(format(lambda, digits=6))),
						adj = c(0, 1)
					)
				}
				graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
				graphics::axis(2, at=seq(0,YlimMax,ceiling(YlimMax/10)), labels=seq(0,YlimMax,ceiling(YlimMax/10)), cex.axis=cex.axis)
				
				#plot the confidence interval of QQ-plot
				#print(log.Quantiles[index])
				qq_col = grDevices::rainbow(R)
				#if(conf.int)	polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
				if(conf.int)	graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=qq_col[i],border=conf.int.col)
				
				if(!is.null(threshold.col)){
				    graphics::par(xpd=FALSE);
				    graphics::abline(a = 0, b = 1, col = threshold.col[1],lwd=2);
				    graphics::par(xpd=TRUE)
				    }
				 
				graphics::points(log.Quantiles, log.P.values, col = col[1],pch=19,cex=2)
				
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						thre.line=-log10(min(threshold))
						if(amplify==TRUE){
							thre.index=which(log.P.values>=thre.line)
							if(length(thre.index)!=0){
							    #print("!!!!")
								#cover the points that exceed the threshold with the color "white"
								graphics::points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,lwd=3,cex=cex[3])
								if(is.null(signal.col)){
									graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch[1],cex=signal.cex[1])
								}else{
									graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
								}
							}
						}
					}
				}
				box()
				if(file.output) grDevices::dev.off()
			}
		}
		print("Multiple QQ plot has been finished!")
	}
	}#End of Whole function 
