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

## Figure 1c/d

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

p3 <- plot_pca(eigenvec, "PC1", "PC3", "species", species_cols, show_legend = FALSE, label_mode = "plain",
               point_size = 0.1, axis_text_size = 5, axis_title_size = 5, label_size = 1.2)

p4 <- plot_pca(eigenvec, "PC1", "PC3", "region", region_cols_custom, show_legend = TRUE, label_mode = "plain",
               point_size = 0.1, axis_text_size = 5, axis_title_size = 5, label_size = 1.2)

(p1 + p2) / (p3 + p4)