## Figure 1a

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(terra)       
library(tidyterra)  
library(paletteer)

r <- rast("NE1_HR_LC_SR_W_DR.tif")

sample_location <- read.table("sample_location.txt",sep = '\t', header = TRUE, check.names = FALSE)

make_map <- function(xlim, ylim) {
  rast_win <- crop(r, ext(c(xlim, ylim)))
  df_win <- sample_location %>%
    filter(
      lon_off >= xlim[1], lon_off <= xlim[2],
      lat_off >= ylim[1], lat_off <= ylim[2]
    )
  
  ggplot() +
    geom_spatraster_rgb(data = rast_win, inherit.aes = FALSE) +
    geom_point(
      data = df_win,
      aes(x = lon_off, y = lat_off, fill = species, size = n, shape = shape),
      color = "black", stroke = 0.25 / ggplot2::.pt, alpha = 0.95,
      show.legend = c(size = TRUE, fill = FALSE, shape = FALSE)
    ) +
    scale_shape_identity() +
    scale_fill_manual(values = species_cols, drop = FALSE) +
    scale_size_area(
      limits   = size_limits,
      max_size = 3,
      trans    = "log10",
      breaks   = size_breaks,
      labels   = as.character(size_breaks),
      name     = "Sample count"
    ) +
    geom_text_repel(
      data = df_win,
      aes(x = lon_off, y = lat_off, label = label, color = species),
      parse = TRUE,
      size = 5 / ggplot2::.pt,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.15,
      segment.colour = NA,  
      show.legend = FALSE
    ) +
    scale_color_manual(values = species_cols, drop = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void() +
    guides(
      size = guide_legend(
        override.aes = list(shape = 21, fill = "grey70", color = "black", alpha = 1)
      )
    )
}

size_breaks <- c(1, 100, 200)
size_breaks <- size_breaks[size_breaks >= size_limits[1] & size_breaks <= size_limits[2]]

p_asia <- make_map(c( 60, 140), c(-10, 54))


## Figure 1b

# =======================
# Libraries
# =======================
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)

# =======================
# Parameters (filenames only)
# =======================
tree_file  <- "prune.in.fmiss0.1.maf0.05.hwee3.min4.nexus.varsites.phy.treefile"
info_file  <- "Sample.info"
color_file <- "species_color.txt"

BRANCH_COLOR_BY_SPECIES <- TRUE
outgroup_vec <- c("147209", "ypt4001")

OPEN_ANGLE <- 180
TREE_SIZE  <- 0.25

species_color <- read.table(
  color_file, sep = "\t", header = FALSE, comment.char = "",
  stringsAsFactors = FALSE
)
species_cols <- setNames(species_color$V2, species_color$V1)

tree_new <- read.tree(tree_file)

if (length(outgroup_vec) > 0) {
  tree_new <- root(tree_new, outgroup = outgroup_vec, resolve.root = TRUE)
}

# =======================
# Base plot & tip order
# =======================
p_base <- ggtree(
  tree_new,
  layout     = "fan",
  open.angle = OPEN_ANGLE,
  ladderize  = TRUE,
  size       = TREE_SIZE,
  color      = "grey50"
)

tip_order <- p_base$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

info <- read.delim(info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

info2 <- info %>%
  transmute(
    label   = tip_order,
    species = trimws(species)
  ) %>%
  filter(label %in% tree_new$tip.label)

anno_df <- tibble(label = tip_order) %>%
  left_join(info2, by = "label") %>%
  mutate(species = replace_na(species, "Unknown"))

# =======================
# Tree coloring
# =======================
if (BRANCH_COLOR_BY_SPECIES) {
  group_list <- split(anno_df$label, anno_df$species)

  p_tree <- ggtree(
    tree_new,
    layout     = "fan",
    open.angle = OPEN_ANGLE,
    ladderize  = TRUE,
    size       = TREE_SIZE,
    color      = NA
  ) %>%
    groupOTU(group_list) +
    geom_tree(aes(color = group), size = TREE_SIZE, lineend = "round") +
    scale_color_manual(values = species_cols, guide = "none")

} else {
  p_tree <- p_base
}

p_tree

## Figure 1c

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

eigenvec_file <- "prune_in.3070_sample.rmoutgroup.fmissing.maf0.05.hwee3.eigenvec"
eigenval_file <- "prune_in.3070_sample.rmoutgroup.fmissing.maf0.05.hwee3.eigenval"
info_file     <- "Sample.info"
color_file    <- "species_color.txt"


eigenvec <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
eigenval <- scan(eigenval_file)
colnames(eigenvec)[1:2] <- c("ind", "sample")
colnames(eigenvec)[3:ncol(eigenvec)] <- paste0("PC", seq_len(ncol(eigenvec) - 2))



var_exp <- eigenval / sum(eigenval) * 100
info <- read.delim(info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

meta <- info %>%
  transmute(
    sample  = .data[[id_col]],
    region  = ifelse(is.na(china_region) | china_region == "", "Other_region", china_region),
    species = trimws(species)
  )

eigenvec <- eigenvec %>%
  left_join(meta, by = "sample") %>%
  mutate(
    region  = replace_na(region, "Other_region"),
    species = replace_na(species, "Unknown")
  )

region_cols_custom <- c(
  "East China"      = "#EDC948",
  "North China"     = "#76B7B2",
  "South China"     = "#E15759",
  "Northwest China" = "#59A14F",
  "Northeast China" = "#4E79A7",
  "Southwest China" = "#F28E2B",
  "Central China"   = "#B07AA1",
  "Other_region"    = "grey50"
)

species_color <- read.table(color_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
species_cols  <- setNames(species_color$V2, species_color$V1)

plot_pca <- function(df, pc_x, pc_y, color_by, palette, show_legend = TRUE, label_mode = "plain",
                     point_size = 0.2, axis_text_size = 6, axis_title_size = 6, label_size = 1.2) {

  ix <- as.integer(sub("PC", "", pc_x))
  iy <- as.integer(sub("PC", "", pc_y))
  x_lab <- paste0(pc_x, " (", sprintf("%.2f", var_exp[ix]), "%)")
  y_lab <- paste0(pc_y, " (", sprintf("%.2f", var_exp[iy]), "%)")

  p <- ggplot(df, aes(x = .data[[pc_x]], y = .data[[pc_y]], color = .data[[color_by]])) +
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_manual(values = palette, na.value = "grey70", drop = FALSE) +
    theme_classic(base_size = axis_text_size) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      axis.title = element_text(size = axis_title_size)
    ) +
    labs(x = x_lab, y = y_lab, color = color_by) +
    coord_equal()

  if (!is.null(label_mode)) {
    lab_df <- df %>%
      mutate(.grp = .data[[color_by]]) %>%
      group_by(.grp) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(.lab = as.character(.grp))

    if (label_mode == "repel") {
      p <- p + geom_text_repel(data = lab_df, aes(label = .lab),
                               size = label_size, show.legend = FALSE,
                               max.overlaps = Inf, segment.color = NA)
    } else if (label_mode == "plain") {
      p <- p + geom_text(data = lab_df, aes(label = .lab),
                         size = label_size, show.legend = FALSE, vjust = -0.5)
    }
  }

  p
}

p1 <- plot_pca(eigenvec, "PC1", "PC2", "species", species_cols, show_legend = FALSE, label_mode = "plain",
               point_size = 0.1, axis_text_size = 5, axis_title_size = 5, label_size = 1.2)

p2 <- plot_pca(eigenvec, "PC1", "PC2", "region", region_cols_custom, show_legend = TRUE, label_mode = "plain",
               point_size = 0.1, axis_text_size = 5, axis_title_size = 5, label_size = 1.2)

## Figure 1d

## From ADMIXTURE

## Figure 1e

make_pie <- function(df, group_name, color_map) {
  df_pie <- df %>%
    filter(group == group_name) %>%
    count(china_region) %>%
    mutate(
      colour = ifelse(china_region %in% names(color_map),
                      color_map[china_region], "#D5D5D5"),
      label  = paste0(china_region, "\n(n=", n, ")")
    ) %>%
    arrange(desc(china_region)) %>%
    mutate(prop = n / sum(n),
           ypos = cumsum(prop) - 0.5 * prop)
  
  ggplot(df_pie, aes(x = "", y = n, fill = china_region)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y", start = 270) +
    scale_fill_manual(values = color_map) +
    geom_text(aes(y = ypos * sum(n), label = label),
              color = "black", size = 3) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(group_name)
}

groups <- unique(group_df1$group)  # 找到所有 group
plots  <- lapply(groups, make_pie, df = group_df1, color_map = color_map)

## Figure 1f

raster_background_points <- function(...) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(..., raster.dpi = 600)
  } else {
    ggplot2::geom_point(...)
  }
}

make_track <- function(dt, cohort_name, metric_name = c("FST", "PI"), show_x = FALSE, show_title = FALSE) {
  metric_name <- match.arg(metric_name)

  value_col <- if (metric_name == "FST") "WEIGHTED_FST" else "pi_ratio"
  pass_col  <- if (metric_name == "FST") "pass_fst" else "pass_pi"

  plot_dt <- copy(dt[is.finite(get(value_col)) & is.finite(genome_pos)])
  plot_dt[, y := as.numeric(get(value_col))]
  if (metric_name == "PI") plot_dt[, y := -y]

  plot_dt[, point_color := background_color]
  plot_dt[get(pass_col) == TRUE & category == "NC",     point_color := COLORS[["NC"]]]
  plot_dt[get(pass_col) == TRUE & category == "SC",     point_color := COLORS[["SC"]]]
  plot_dt[get(pass_col) == TRUE & category == "Shared", point_color := COLORS[["BOTH"]]]

  background_dt <- plot_dt[point_color %in% c(COLORS[["Grey1"]], COLORS[["Grey2"]])]
  highlight_dt  <- plot_dt[!(point_color %in% c(COLORS[["Grey1"]], COLORS[["Grey2"]]))]

  max_fst <- max(dt$WEIGHTED_FST, na.rm = TRUE)
  max_pi  <- max(dt$pi_ratio, na.rm = TRUE)

  y_title <- if (metric_name == "FST") expression(italic(F)[ST]) else expression(pi~ratio)

  p <- ggplot() +
    raster_background_points(
      data = background_dt,
      aes(x = genome_pos, y = y, colour = point_color),
      size = POINT_SIZE_BACKGROUND, alpha = 0.75, shape = 16
    ) +
    geom_point(
      data = highlight_dt,
      aes(x = genome_pos, y = y, colour = point_color),
      size = POINT_SIZE_HIGHLIGHT, alpha = 0.98, shape = 16
    ) +
    geom_rect(
      data = gap_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "white", colour = NA
    ) +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = genome_layout$chr_center,
      labels = genome_layout$chrom,
      limits = c(PLOT_X_MIN, PLOT_X_MAX),
      expand = expansion(mult = c(0, 0)),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    labs(
      x = NULL,
      y = y_title,
      title = if (show_title) paste0(cohort_name, "/DH") else NULL
    ) +
    coord_cartesian(xlim = c(PLOT_X_MIN, PLOT_X_MAX), clip = "off", expand = FALSE) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(
        size = 10, face = "bold", hjust = 0,
        colour = if (cohort_name == "NC") COLORS[["NC"]] else COLORS[["SC"]],
        margin = margin(b = 1)
      ),
      axis.title.y = element_text(size = 8.5, margin = margin(r = 3)),
      axis.text.y = element_text(size = 7, colour = "black"),
      axis.text.x = if (show_x) {
        element_text(size = 6, face = "italic", colour = "black", margin = margin(t = 1))
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      axis.line.x  = if (show_x) element_line(linewidth = 0.30) else element_blank(),
      axis.line.y  = element_line(linewidth = 0.30),
      axis.ticks.y = element_line(linewidth = 0.25),
      plot.margin = margin(t = 1.2, r = 4, b = if (show_x) 1.0 else 0.1, l = 2, unit = "mm")
    )

  if (metric_name == "FST") {
    p <- p + scale_y_continuous(limits = c(0, max_fst * 1.08), expand = expansion(mult = c(0, 0.02)))
  } else {
    p <- p + scale_y_continuous(
      limits = c(-max_pi * 1.08, 0),
      breaks = pretty(c(0, max_pi)),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0.02, 0.02))
    )
  }

  labels_this_track <- label_points[cohort == cohort_name & metric == metric_name]

  if (nrow(labels_this_track) > 0L) {
    p <- p +
      geom_segment(
        data = labels_this_track,
        aes(x = x, xend = x, y = y, yend = label_y),
        inherit.aes = FALSE,
        linewidth = 0.18,
        colour = "grey55"
      ) +
      geom_text(
        data = labels_this_track,
        aes(x = x, y = label_y, label = gene_name, colour = label_color),
        inherit.aes = FALSE,
        size = 2.35,
        fontface = "italic",
        vjust = labels_this_track$vjust,
        show.legend = FALSE
      ) +
      scale_colour_identity()
  }

  p
}

p_NC_FST <- make_track(NC_scan, "NC", "FST", show_x = FALSE, show_title = TRUE)
p_NC_PI  <- make_track(NC_scan, "NC", "PI",  show_x = TRUE,  show_title = FALSE)
p_SC_FST <- make_track(SC_scan, "SC", "FST", show_x = FALSE, show_title = TRUE)
p_SC_PI  <- make_track(SC_scan, "SC", "PI",  show_x = TRUE,  show_title = FALSE)

patchwork::wrap_plots(
  p_NC_FST,
  p_NC_PI,
  p_SC_FST,
  p_SC_PI,
  ncol = 1,
  heights = c(1, 1, 1, 1)
)


## Figure 1g

## Drawn directly in Adobe Illustrator

## Figure 1h

library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(ggplot2)
library(dplyr)
library(grid)

traits <- read.delim('trait.meta.txt', header = TRUE, row.names = 1)

pheno <- read.table("chicken.plink.fam", head=FALSE)
rownames(pheno) <- pheno[,2]
out <- pheno[, -c(3:5,7)]
colnames(out) <- c("FID", "IID", sprintf("T%s", 1:(ncol(out)-2)))

pheno <- out[,-(1:2)] 
colnames(pheno) <- 1:ncol(pheno)

p.data <- pheno 
colnames(p.data) <- traits[colnames(pheno), 1]
p.cor <- cor(p.data, method='s', use="pairwise.complete.obs") 

k <- 6
phm <- pheatmap(p.cor, cutree_rows = k, cutree_cols = k, 
                cellwidth = 10, cellheight = 10, 
                border_color = 'white')
group <- cutree(phm$tree_row, k = k) %>% 
  as.data.frame() %>% rename(group = '.') %>%
  mutate(group = sprintf("P%s", group))

pord <- phm$tree_row$labels[phm$tree_row$order]

gcol <- setNames(brewer.pal(n = k, name = "Set1"), 
                 sort(unique(group$group)))

library(paletteer)
cols <- setNames(colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(length(pord)),
                 pord)

gcolx <- cols %>% as.data.frame() %>% rename(col='.') %>% mutate(id=rownames(.), group=group[id,'group']) %>% arrange(group) %>%  group_by(group) %>%
  arrange(col, .by_group = TRUE) %>%  # 按亮度排序
  summarise(
    median_col = col[ceiling(n() / 2)]  # 取中位数颜色
  )
gcol <- setNames(gcolx$median_col, gcolx$group)

group <- group %>% mutate(variable = rownames(.)) %>%
  mutate(group = factor(group, levels = sort(unique(group))))

pcor <- p.cor[pord, pord]
pmat <- cor.mtest(p.data, conf.level = .95)$p
pmat <- matrix(p.adjust(pmat[pord, pord], method = 'BY'), nrow = nrow(pmat))
rownames(pmat) <- colnames(pmat) <- pord

variable_color <- gcol[group[pord, 'group']]

p1 <- ggplot(group[pord,], aes(x = 10, y = seq_along(variable))) +
  geom_rect(aes(xmin=0, xmax=0.4, 
                ymin=seq_along(variable)-0.5, ymax=seq_along(variable)+0.5, 
                fill = variable, color = variable)) +
  geom_text(aes(x = 0.5, y = seq_along(variable), label = variable), 
            hjust = 0, size=2.5, color="black") +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  xlab(NULL) + ylab(NULL) + 
  scale_y_reverse() +
  theme_void() +
  theme(legend.position = "none")

p2 <- ggplot(group[pord,], aes(x = 10, y = seq_along(variable))) +
  geom_rect(aes(xmin=0, xmax=0.4, 
                ymin=seq_along(variable)-0.5, ymax=seq_along(variable)+0.5, 
                fill = group, color = group)) +
  geom_text(aes(x = 0.5, y = seq_along(variable), label = variable), 
            hjust = 0, size=2.5, color="black") +
  scale_fill_manual(values=gcol) +
  scale_color_manual(values=gcol) +
  xlab(NULL) + ylab(NULL) + 
  scale_y_reverse() +
  theme_void() +
  theme(legend.position = "none")

genetic.corr <- read.table("genetic.corr.out", row.names=1, head=TRUE, check.names=FALSE)
genetic.corr <- as.matrix(genetic.corr)
diag(genetic.corr) <- 1
genetic.pvls <- read.table("genetic.corr.pvl.out", row.names=1, head=TRUE, check.names=FALSE)
genetic.pvls <- as.matrix(genetic.pvls)
diag(genetic.pvls) <- 0

g.cor <- genetic.corr
g.val <- genetic.pvls
for(i in 1:nrow(genetic.corr)){
  for(j in 1:ncol(genetic.corr)){
    g.cor[i,j] <- mean(c(genetic.corr[i,j], genetic.corr[j,i]), na.rm=TRUE)
    g.val[i,j] <- min(c(genetic.pvls[i,j], genetic.pvls[j,i]), na.rm=TRUE)
  }
}
rownames(g.cor) <- colnames(g.cor) <- rownames(g.val) <- colnames(g.val) <- colnames(pheno)
g.cor[is.na(g.cor)] <- 0
g.val[is.infinite(g.val)] <- 1
rownames(g.cor) <- traits[rownames(g.cor), 1]
colnames(g.cor) <- traits[colnames(g.cor), 1]
rownames(g.val) <- traits[rownames(g.val), 1]
colnames(g.val) <- traits[colnames(g.val), 1]

gcor <- g.cor[pord, pord]
gmat <- matrix(p.adjust(g.val[pord, pord], method = 'BY'), 
               nrow = length(pord))
rownames(gmat) <- colnames(gmat) <- pord


dev.off()

pdf("correlation_heatmap.pdf", width = 8.27, height = 8.27)
# color bar
print(cowplot::plot_grid(p1, p2))
col <- brewer.pal(n = 9, name = "RdGy") %>% rev()
corrplot(
    pcor,
    method = "square",
    type = "lower",
    order = "original",
    p.mat = pmat,
    insig = 'label_sig',
    sig.level = c(0.0001, 0.001, 0.01),
    col = col,
    cl.pos = "b",
    tl.col = variable_color,
    tl.cex = 0.7,
    pch.cex = 0.8,
    mar = c(0, 0, 0, 1)
)

vie1 <- viewport(
    width = 0.38,
    height = 0.585,
    x = 1,
    y = 0.63
)
print(p1, vp = vie1)

col <- brewer.pal(n = 9, name = "PiYG") %>% rev()
corrplot(
    gcor,
    method = "square",
    type = "upper",
    order = "original",
    p.mat = gmat,
    insig = 'label_sig',
    sig.level = c(0.001, 0.01, 0.1),
    col = col,
    cl.pos = "b",
    tl.col = variable_color,
    tl.cex = 0.7,
    pch.cex = 0.8,
    mar = c(0, 0, 0, 1)
)

dev.off()



phm1 <- pheatmap(pcor)
phm2 <- pheatmap(gcor)

library(vegan)
mantel.tst <- mantel(pcor, gcor, permutations=9999)

pdf("correlation_tree.pdf", width=8.27, heigh=4, pointsize=10)
hc1 <- phm1$tree_row
hc2 <- phm2$tree_row
l <- length(hc1$order)
# The matrix to draw the arrows:
ord_arrow <- cbind((1:l)[order(hc1$order)], (1:l)[order(hc2$order)])
# The two vectors of ordered leave labels:
leaves1 <- hc1$labels[hc1$order]
leaves1.col <- gcol[group[leaves1, 'group']]
names(leaves1.col) <- leaves1
leaves2 <- hc2$labels[hc2$order]
leaves2.col <- gcol[group[leaves2, 'group']]
names(leaves2.col) <- leaves2
# And the plot:
layout(matrix(1:5, nrow = 1), width = c(5, 3, 5.5, 3, 5))
# The first dendrogram:
op <- par(mar = c(3, 3, 3, 0))
plot(
    as.dendrogram(hc1),
    horiz = TRUE,
    leaflab = "none",
    ylim = c(0, l),
    main = "Phenotypic correlation"
)
# The first serie of labels (i draw them separately because, for the second serie, I didn't find a simple way to draw them nicely on the cluster):
par(op)
op <- par(mar = c(3, 0, 3, 0))
plot(
    NA,
    bty = "n",
    axes = FALSE,
    xlim = c(0, 1),
    ylim = c(0, l),
    ylab = "",
    xlab = ""
)
sapply(1:l, function(x) {
    points(0, x, pch = 15, col = leaves1.col[x])
    text(
        x = 0,
        y = x,
        col = leaves1.col[x],
        labels = leaves1[x],
        pos = 4,
        cex = 0.8
    )
})
# The arrows:
par(op)
op <- par(mar = c(3, 0, 3, 0))
plot(
    NA,
    bty = "n",
    axes = FALSE,
    xlim = c(0, 1),
    ylim = c(0, l),
    ylab = "",
    xlab = "",
    main = mantel.tst$signif
)
apply(ord_arrow, 1, function(x) {
    arrows(0,
           x[1],
           1,
           x[2],
           code = 3,
           length = 0.05,
           col = leaves1.col[x])
})
# The second serie of labels:
par(op)
op <- par(mar = c(3, 0, 3, 0))
plot(
    NA,
    bty = "n",
    axes = FALSE,
    xlim = c(0, 1),
    ylim = c(0, l),
    ylab = "",
    xlab = ""
)
sapply(1:l, function(x) {
    points(1, x, pch = 15, col = leaves2.col[x])
    text(
        x = 1,
        y = x,
        col = leaves2.col[x],
        labels = leaves2[x],
        pos = 2,
        cex = 0.8
    )
})
# And the second dendrogram (to reverse it I reversed the xlim vector:
par(op)
op <- par(mar = c(3, 0, 3, 3))
plot(
    as.dendrogram(hc2),
    horiz = TRUE,
    xlim = c(0, max(hc2$height)),
    leaflab = "none",
    ylim = c(0, l),
    main = "Genetic correlation"
)
par(op)
par(mfrow = c(1, 1))

dev.off()

## Figure 1i
library(CMplot)

pval <- 1e-3

sigsnp <- rowSums(data[,-(1:3)] < pval) > 0

CMplot(
    data[sigsnp,],
    col = mycols,
    type = "p",
    plot.type = "c",
    r = 0.4,
    cex = 0.05,
    cir.axis = TRUE,
    outward = FALSE,
    cir.axis.col = "black",
    cir.chr.h = 2.5,
    ylim = c(3, 10),
    chr.den.col = "black",
    file = "pdf",
    file.name = pval,
    file.output = TRUE,
    verbose = TRUE,
    width = 10,
    height = 10
) 

CMplot(
    data[sigsnp,],
    col = mycols,
    type = "p",
    plot.type = "c",
    r = 0.4,
    cex = 0.05,
    cir.chr.h = 2.5,
    ylim = c(3, 10),
    chr.den.col = c("darkgreen", "yellow", "red"),
    bin.size = 1e6,
    outward = FALSE,
    file = "pdf",
    file.name = sprintf("density-%s", pval),
    file.output = TRUE,
    verbose = TRUE,
    width = 10,
    height = 10
)

## Figure 1j

CMplot(
    data[, c("SNP", "chr", "pos", trait)],
    col = mycols,
    cex = 0.1,
    signal.cex = 0.75,
    plot.type = "m",
    LOG10 = TRUE,
    axis.cex = 0.5,
    axis.lwd = 0.5,
    lab.cex = 0.75,
    highlight = snps,
    highlight.text = genes,
    highlight.col = rep("red", length(snps)),
    highlight.cex = 0.25,
    highlight.text.cex = 0.5,
    highlight.text.col = rep("red", length(snps)),
    ylim = c(2, ymax),
    threshold = 10 ^ -limy,
    threshold.lty = 2,
    threshold.col = "red",
    file = "pdf",
    file.name = sprintf('highlight-%s', trait),
    file.output = TRUE,
    verbose = TRUE,
    width = 8.27,
    height = 4.5,
    chr.labels.angle = 0,
    main.cex = 0.75,
    main = trait
)

## Figure 1k

# Load required libraries for data manipulation and visualization
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# Load metadata and set row names for easy indexing
meta <- readRDS('trait.meta.rds')
rownames(meta) <- meta$id

# Load and preprocess GWAS annotation data
sankdat <- read.table('gwas.assoc.sig.full.annotation.bed')
sankdat <- sankdat %>% 
    # Rename columns for clarity
    rename(chr=V1, start=V2, end=V3, SNP=V4, gene=V7, dist=V8, category=V9, OCR=V10) %>%
    # Reformat OCR distances into discrete categorical bins (e.g., inPeak, <1K, >3K)
    mutate(OCR=OCR/1000, OCR=ifelse(OCR==0, "inPeak", ifelse(OCR>3, ">3K", sprintf("<%sK", ceiling(OCR))))) %>%
    # Map trait names from metadata using ID column (V5)
    mutate(trait=meta[as.character(V5),'variable'])

# Define node labels and order for the Sankey diagram
showlabs <- c('CDS','UTR5','UTR3','intron','proximal','distal','inPeak','<1K','<2K','<3K','>3K')
nodex <- c(meta$variable, showlabs)

# Assign specific colors to each node using custom palettes
ncols <- c(meta$tcol, setNames(c(paletteer::paletteer_d("MetBrewer::Juarez"), 
                                 paletteer::paletteer_d("lisa::GustavKlimt")), showlabs))

# Generate the Sankey diagram
library(ggsankey)
pdf('sankey.pdf', width = 8.27, height = 5)

sankdat %>% 
    # Select columns representing the flow layers: Trait -> Category -> OCR
    dplyr::select(trait, category, OCR) %>%
    # Transform data into long format required by ggsankey
    make_long(trait, category, OCR) %>% 
    # Filter labels to only show specific categories/distances
    mutate(label=ifelse(node %in% showlabs, node, NA)) %>%
    # Set factor levels to control the vertical stacking order
    mutate(node = factor(node, levels=rev(nodex))) %>%
    # Initialize ggplot with ggsankey aesthetics
    ggplot(aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node, 
               label = label,
               fill = node)) +
    # Draw flows and nodes
    geom_sankey(flow.alpha = 0.9, width = 0.1, alpha = 1) +  
    # Add text labels to nodes
    geom_sankey_label(size = 3, color = 'white', fill = NA) +
    # Apply theme and color scales
    theme_sankey(base_size = 18) +
    scale_fill_manual(values = ncols) + 
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5)) +
    ggtitle("GWAS significant SNPs")

dev.off()