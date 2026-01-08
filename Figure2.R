## Figure 2a

library(data.table)
library(GenomicRanges)
library(IRanges)
library(CMplot)

NC_DH <- fread("NC.xpclr.pi.fst.merge.tsv", sep = "\t", header = TRUE)
SC_DH <- fread("SC.xpclr.pi.fst.merge.tsv", sep = "\t", header = TRUE)

make_snp_gr <- function(pts) {
  GRanges(
    seqnames = pts$chr,
    ranges   = IRanges(start = pts$pos, width = 1),
    SNP      = paste0(pts$chr, ":", pts$pos)
  )
}

gr_NC  <- make_snp_gr(pts_NC)
gr_SC  <- make_snp_gr(pts_SC)
gr_all <- unique(c(gr_NC, gr_SC))

gr_all$category <- "none"

gr_all[queryHits(findOverlaps(gr_all, NCSCboth))]$category <- "both"
gr_all[queryHits(findOverlaps(gr_all, NCunique))]$category <- "NCunique"
gr_all[queryHits(findOverlaps(gr_all, SCunique))]$category <- "SCunique"

out_cat <- data.table(
  chr = as.character(seqnames(gr_all)),
  pos = start(gr_all),
  SNP = mcols(gr_all)$SNP,
  cat = mcols(gr_all)$category
)

plot_manhattan_simple <- function(df,
                                  metric = c("pi_ratio", "MEAN_FST"),
                                  cohort = c("NC", "SC"),
                                  both, NCunique, SCunique,
                                  thr_q = 0.95,
                                  base_cols = c("grey40","grey70"),
                                  cex = 0.08,
                                  highlight_cex = 0.12) {
  metric <- match.arg(metric)
  cohort <- match.arg(cohort)

  df <- as.data.table(df)

  df[, chrom := as.character(chrom)]
  df[, end   := as.numeric(end)]
  df[, val   := as.numeric(get(metric))]
  df <- df[is.finite(end) & is.finite(val)]

  dat <- data.frame(
    SNP        = paste0(df$chrom, ":", df$end),
    Chromosome = df$chrom,
    Position   = df$end,
    P          = df$val
  )

  thr <- as.numeric(quantile(dat$P, probs = thr_q, na.rm = TRUE))
  dat_thr <- dat$SNP[dat$P > thr]

  snp_gr <- GRanges(dat$Chromosome, IRanges(dat$Position, width = 1))

  in_both <- dat$SNP[queryHits(findOverlaps(snp_gr, both))]
  in_NC   <- dat$SNP[queryHits(findOverlaps(snp_gr, NCunique))]
  in_SC   <- dat$SNP[queryHits(findOverlaps(snp_gr, SCunique))]

  if (cohort == "NC") {
    hl_both <- intersect(in_both, dat_thr)
    hl_uniq <- intersect(in_NC,   dat_thr)

    hl_snp <- c(hl_uniq, hl_both) 
    hl_col <- c(rep("#498DCB", length(hl_uniq)),
                rep("#F41520", length(hl_both)))
  } else {
    hl_both <- intersect(in_both, dat_thr)
    hl_uniq <- intersect(in_SC,   dat_thr)

    hl_snp <- c(hl_uniq, hl_both)
    hl_col <- c(rep("#255372", length(hl_uniq)),
                rep("#F41520", length(hl_both)))
  }

  # ---- plot ----
  CMplot(
    dat,
    type = "p",
    plot.type = "m",
    LOG10 = FALSE,
    col = base_cols,
    ylim = c(0, max(dat$P, na.rm = TRUE) * 1.05),
    cex  = cex,
    chr.labels.angle = 0,
    highlight = hl_snp,
    highlight.col = hl_col,
    highlight.cex = highlight_cex,   
    ylab = metric,
    main = paste0(cohort, " | ", metric, " | top", (1 - thr_q) * 100, "%"),
    file.output = FALSE
  )
}

plot_manhattan_simple(
  NC_DH, metric = "pi_ratio", cohort = "NC",
  both = NCSCboth, NCunique = NCunique, SCunique = SCunique
)

plot_manhattan_simple(
  SC_DH, metric = "pi_ratio", cohort = "SC",
  both = NCSCboth, NCunique = NCunique, SCunique = SCunique
)


## Figure 2c
library(data.table)

plot_region_sc_nc <- function(
  chrom, v1_bp, v2_bp,
  pad_mb = 0.3,
  file_sc = "SC_DH.xpclr.pi.fst.merge.tsv",
  file_nc = "NC_DH.xpclr.pi.fst.merge.tsv",
  col_sc  = "#234882",
  col_nc  = "#90bae0",
  pch_sc  = 17,
  pch_nc  = 20,
  gene_name = NULL
){
  stopifnot(is.numeric(v1_bp), is.numeric(v2_bp), is.numeric(pad_mb))

  pad_bp <- pad_mb * 1e6
  x_left_bp  <- min(v1_bp, v2_bp) - pad_bp
  x_right_bp <- max(v1_bp, v2_bp) + pad_bp


  res_sc <- extract_region(file = file_sc, chr = as.character(chrom),
                           start = x_left_bp, end = x_right_bp, units = "bp")
  res_nc <- extract_region(file = file_nc, chr = as.character(chrom),
                           start = x_left_bp, end = x_right_bp, units = "bp")


  prep_df <- function(d){
    d <- within(d, {
      end      <- as.numeric(end)
      pos_mb   <- end / 1e6
      pi_ratio <- as.numeric(pi_ratio)
      MEAN_FST <- as.numeric(MEAN_FST)
      xpclr    <- as.numeric(xpclr)
    })
    d <- d[is.finite(d$end) & is.finite(d$pos_mb), ]
    d <- d[order(d$end), ]
    d
  }

  d_sc <- prep_df(res_sc)
  d_nc <- prep_df(res_nc)


  d2_sc <- d_sc[is.finite(d_sc$xpclr), ]
  d2_nc <- d_nc[is.finite(d_nc$xpclr), ]


  d1_sc <- d_sc[is.finite(d_sc$pi_ratio) & d_sc$pi_ratio > 0 & is.finite(d_sc$MEAN_FST), ]
  d1_nc <- d_nc[is.finite(d_nc$pi_ratio) & d_nc$pi_ratio > 0 & is.finite(d_nc$MEAN_FST), ]


  colpick <- function(nm, alts) {
    hit <- which(nm %in% alts)
    if (length(hit) == 0) stop("Cannot find columns: ", paste(alts, collapse = "/"))
    nm[hit[1]]
  }

  fread_cols <- function(path){
    DT <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
    nm <- names(DT)
    c_end   <- colpick(nm, c("end","pos","position"))
    c_pi    <- colpick(nm, c("pi_ratio","piRatio","Pi_ratio"))
    c_fst   <- colpick(nm, c("MEAN_FST","mean_fst","FST","fst"))
    c_xpclr <- colpick(nm, c("xpclr","XPCLR","XP-CLR"))
    DT[, .(
      end      = as.numeric(get(c_end)),
      pi_ratio = as.numeric(get(c_pi)),
      MEAN_FST = as.numeric(get(c_fst)),
      xpclr    = as.numeric(get(c_xpclr))
    )]
  }

  q95 <- function(x) stats::quantile(x[is.finite(x)], 0.95, na.rm = TRUE, type = 7)

  DT_sc <- fread_cols(file_sc)
  DT_nc <- fread_cols(file_nc)

  thr <- list(
    sc = list(
      xpclr  = q95(DT_sc$xpclr),
      neg2pi = q95(-log2(DT_sc$pi_ratio[DT_sc$pi_ratio > 0])),
      fst    = q95(DT_sc$MEAN_FST)
    ),
    nc = list(
      xpclr  = q95(DT_nc$xpclr),
      neg2pi = q95(-log2(DT_nc$pi_ratio[DT_nc$pi_ratio > 0])),
      fst    = q95(DT_nc$MEAN_FST)
    )
  )


  xlim_shared <- range(c(d1_sc$pos_mb, d1_nc$pos_mb, d2_sc$pos_mb, d2_nc$pos_mb), na.rm = TRUE)
  xticks <- pretty(xlim_shared, n = 6)

  vlines_mb <- c(v1_bp, v2_bp) / 1e6
  vline_col <- "#44444480"
  vline_lty <- 3
  vline_lwd <- 1.2


  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
  op <- par(xaxs = "i")


  y_all_top <- c(d2_sc$xpclr, d2_nc$xpclr)
  y_all_top <- y_all_top[is.finite(y_all_top) & y_all_top >= 0]
  ylim_top  <- c(0, max(pretty(c(0, max(y_all_top, na.rm = TRUE)))))

  par(mar = c(4, 4, 1, 5.5))
  plot(d2_sc$pos_mb, d2_sc$xpclr,
       type = "b", pch = pch_sc, lwd = 2, col = col_sc,
       xlim = xlim_shared, ylim = ylim_top * 0.85,
       xlab = "", ylab = "XP-CLR", xaxt = "n")
  lines(d2_nc$pos_mb, d2_nc$xpclr, type = "b", pch = pch_nc, lwd = 1.8, col = col_nc)

  axis(1, at = xticks, labels = format(xticks, scientific = FALSE))
  legend("topleft", c("SC","NC"), col = c(col_sc, col_nc), pch = c(pch_sc, pch_nc), lwd = 2, bty = "n")

  abline(v = vlines_mb, lty = vline_lty, lwd = vline_lwd, col = vline_col)
  abline(h = thr$sc$xpclr, col = col_sc, lty = 2, lwd = 1.2)
  abline(h = thr$nc$xpclr, col = col_nc, lty = 4, lwd = 1.2)


  y_left_all <- c(-log2(d1_sc$pi_ratio), -log2(d1_nc$pi_ratio))
  y_left_all <- y_left_all[is.finite(y_left_all)]
  ylim_left  <- range(pretty(y_left_all))

  par(mar = c(4, 4, 1, 5.5))
  plot(d1_sc$pos_mb, -log2(d1_sc$pi_ratio),
       type = "l", lty = 2, lwd = 1.8, col = col_sc,
       xlim = xlim_shared, ylim = ylim_left,
       xaxt = "n",
       xlab = if (is.null(gene_name)) "Position (Mb)" else bquote(italic(.(gene_name))),
       ylab = expression(-log[2](pi~ratio)))
  lines(d1_nc$pos_mb, -log2(d1_nc$pi_ratio), lty = "C8", lwd = 1.8, col = col_nc)

  axis(1, at = xticks, labels = format(xticks, scientific = FALSE))
  abline(v = vlines_mb, lty = vline_lty, lwd = vline_lwd, col = vline_col)

  legend("topleft",
         legend = c(expression(-log[2](pi~ratio)), expression(italic(F)[ST])),
         lty = c(2, 1), lwd = 1.8, col = "black", bty = "n", cex = 0.9)

  y_right_all <- c(d1_sc$MEAN_FST, d1_nc$MEAN_FST)
  y_right_all <- y_right_all[is.finite(y_right_all) & y_right_all >= 0]
  ylim_right  <- c(0, max(pretty(c(0, max(y_right_all, na.rm = TRUE)))))

  par(new = TRUE)
  plot(d1_sc$pos_mb, d1_sc$MEAN_FST,
       type = "l", lty = 1, lwd = 1.8, col = col_sc,
       axes = FALSE, xlab = NA, ylab = NA,
       xlim = xlim_shared, ylim = ylim_right)
  lines(d1_nc$pos_mb, d1_nc$MEAN_FST, lty = 1, lwd = 1.8, col = col_nc)

  axis(4, at = pretty(ylim_right))
  mtext(expression("Mean " * italic(F)[ST]), side = 4, line = 3)

  par(op); layout(1)

  invisible(list(
    region_bp  = c(start_bp = x_left_bp, end_bp = x_right_bp),
    vlines_mb  = vlines_mb,
    thresholds = thr
  ))
}

plot_region_sc_nc(
  chrom = 10,
  v1_bp = 5900001,
  v2_bp = 5924000,
  gene_name = "APBA2"
)


## Figure 2d
library(clusterProfiler)
library(org.Ggallus.eg.db)
library(dplyr)
library(forcats)
library(ggplot2)
library(paletteer)


run_enrich_go <- function(gene_ids, group_name,
                          OrgDb = org.Ggallus.eg.db,
                          keyType = "GID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1) {
  ego <- enrichGO(
    gene          = gene_ids,
    OrgDb         = OrgDb,
    keyType       = keyType,
    ont           = ont,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff
  )
  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(df)   
  df$Group <- group_name
  df
}


df_NC   <- run_enrich_go(ann_NCunique$gene_id,  "NC")
df_SC   <- run_enrich_go(ann_SCunique$gene_id,  "SC")
df_Both <- run_enrich_go(ann_NCSCboth$gene_id,  "Both")

df_all <- bind_rows(df_NC, df_SC, df_Both)

df_sel <- df_all %>%
  filter(ID %in% go_ids) %>%
  mutate(
    Group = factor(Group, levels = c("NC", "SC", "Both")),
    Description = fct_inorder(as.factor(Description))
  )


p <- ggplot(df_sel, aes(x = Group, y = Description)) +
  geom_point(
    aes(fill = -log10(p.adjust), size = Count),
    shape = 21, colour = "black", stroke = 0.4
  ) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Blue") +
  scale_size(range = c(2, 8), name = "Gene Count") +
  labs(x = NULL, y = NULL, fill = "-log10(padj)") +
  theme_bw(base_size = 6) +
  theme(
    axis.text     = element_text(size = 5),
    legend.text   = element_text(size = 5),
    legend.title  = element_text(size = 5)
  )