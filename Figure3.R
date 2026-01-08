## Figure 3b
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

pdf("Fig3b.correlation_heatmap.pdf", width = 8.27, height = 8.27)
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

pdf("Fig3b.correlation_tree.pdf", width=8.27, heigh=4, pointsize=10)
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

## Figure 3c
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


## Figure 3d

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

## Figure 3e

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
pdf('Fig3e.sankey.pdf', width = 8.27, height = 5)

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
