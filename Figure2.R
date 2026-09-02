## Figure 2a/2d

# Load pre-processed expression and metadata
RNAseq <- readRDS('RNAseq.final.rds')
ATACseq <- readRDS('ATACseq.final.rds')

inputfile <- 'gene.expression.tsv'
genefile <- 'GRCg7b.gene.alias'

# Read gene aliases and expression matrix
genealias <- read.delim(genefile, row.names = 1, header = F)
inputData <- read.delim(inputfile, row.names = 1) %>% as.matrix() 
# Transpose matrix for WGCNA: Rows = Samples, Columns = Genes
exprData <- t(inputData)

# 1. Data Pre-processing
# 1.1 Check for missing values and low-quality genes/samples
gsg <- goodSamplesGenes(exprData, verbose = 3)
if (!gsg$allOK) {
    exprData <- exprData[gsg$goodSamples, gsg$goodGenes]
}

# 1.2 Perform hierarchical clustering to detect potential outliers
sampleTree <- hclust(dist(exprData), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers")

# 1.3 Calculate Tissue Specificity Index (Tau) for each gene
tau_values <- apply(exprData, 2, calculate_tau)
hist(tau_values, breaks=50, col="skyblue", 
     xlab="Tau value", main="Distribution of Tissue Specificity (Tau)")
abline(v=0.8, col="red", lty=2)  # Reference threshold for high specificity

# 1.4 Stratified sampling: Select High Variable Genes (HVG) with tissue specificity
hvg <- names(head(sort(apply(exprData, 2, mad), decreasing=TRUE), n=10000)) 
tau_subset <- tau_values[hvg]
selected_genes <- names(tau_subset[tau_subset > 0.5]) 
exprData <- exprData[, selected_genes]

# 2. Network Construction
# 2.1 Soft-threshold selection for Scale-Free Topology
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(t(exprData), powerVector = powers, verbose = 5)

# Save soft power selection plots
pdf('WGCNA.gene.softPowers.pdf', width = 8.27, height = 8.27)
soft_power <- 8
plist <- plotSoftPowers(sft$fitIndices, soft_power = soft_power)
wrap_plots(plist, ncol=2)
dev.off()

# 2.2 Construct co-expression network and identify modules
net <- blockwiseModules(exprData, power = soft_power, 
                        TOMType = "unsigned", 
                        networkType = "unsigned",
                        minModuleSize = 15, 
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.075,
                        numericLabels = TRUE, 
                        maxBlockSize = 10000, 
                        saveTOMs = TRUE,
                        saveTOMFileBase = "geneTOM",
                        verbose = 3)

# 3. Module Analysis
# 3.1 Visualize gene modules using dendrograms and colors
library(ggthemes)
library(paletteer)
# Define a comprehensive custom color palette
cols <- unique(c(tableau_color_pal('Classic 20')(20), ...)) # (abbreviated for brevity)

moduleLabels <- net$colors
nx <- length(unique(moduleLabels))
mdcols <- setNames(cols[1:nx], sort(unique(moduleLabels)))
genecols <- setNames(mdcols[sprintf("%s",moduleLabels)], names(moduleLabels))
mergedColors <- cbind(Module=genecols[net$blockGenes[[1]]])

pdf('WGCNA.gene.dendroAndColors.pdf', width = 8.27, height = 4)
plotDendroAndColors(net$dendrograms[[1]], mergedColors, hang=0.03,
                    dendroLabels=FALSE, addGuide=TRUE, guideHang=0.05)
dev.off()

# 3.2 Calculate Module Membership (kME) and visualize using UMAP
MEs <- net$MEs
geneModuleMembership <- as.data.frame(cor(exprData, MEs, use = "p"))
colnames(geneModuleMembership) <- gsub("ME", "kME_", colnames(geneModuleMembership))

genemodules <- data.frame(gene=names(moduleLabels), row.names = names(moduleLabels)) %>%
    mutate(module = factor(moduleLabels, levels = unique(moduleLabels))) %>%
    mutate(color = genecols[gene]) %>%
    cbind(geneModuleMembership) %>%
    filter(module != '0') %>% droplevels()

# Calculate Topological Overlap Matrix (TOM) for network visualization
adjmat <- adjacency(exprData[, rownames(genemodules)], power = soft_power, type = "unsigned")
TOM <- TOMsimilarity(adjmat, TOMType = "unsigned")

pdf('WGCNA.full_geneNetwork.pdf', width = 12, height = 12)
plotModuleUMAP(genemodules, TOM = TOM)
dev.off()

## Figure 2b

############ Module-Trait associations

#' Elbow Method for K-Means Clustering
#' https://www.geeksforgeeks.org/machine-learning/elbow-method-for-optimal-value-of-k-in-kmeans/
#' 
#' @param data A numeric matrix or data frame of input data
#' @param max_k Maximum number of clusters to evaluate (default: 10)
#' @param nstart Number of random starts for K-Means (default: 25)
#' @param seed Random seed for reproducibility (default: 123)
#' @param plot Logical, whether to plot the elbow curve (default: TRUE)
#' @return A list containing:
#'   - wcss: Vector of WCSS values for each k
#'   - k_values: Vector of tested k values
#'   - optimal_k: Suggested optimal k based on elbow method
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' data <- matrix(c(rnorm(100, mean = 0), rnorm(100, mean = 5)), ncol = 2)
#' 
#' # Run elbow method
#' result <- elbow_method(data, max_k = 8)
#' print(paste("Suggested optimal k:", result$optimal_k))
elbow_method <- function(data, max_k = 10, nstart = 25, seed = 123, plot = TRUE) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  if (max_k < 2) {
    stop("max_k must be at least 2")
  }
  
  # Convert data to matrix if needed
  data <- as.matrix(data)
  
  # Remove any rows with missing values
  data <- na.omit(data)
  
  # Initialize WCSS vector
  wcss <- numeric(max_k)
  k_values <- 1:max_k
  
  set.seed(seed)
  
  # Calculate WCSS for each k
  for (k in k_values) {
    if (k == 1) {
      # For k=1, WCSS is just total variance
      wcss[k] <- sum(scale(data, scale = FALSE)^2)
    } else {
      # Run K-Means
      km <- kmeans(data, centers = k, nstart = nstart)
      wcss[k] <- km$tot.withinss
    }
  }
  
  # Find the elbow point (optimal k)
  # Using the method that finds the point with maximum curvature
  first_deriv <- diff(wcss)
  second_deriv <- diff(first_deriv)
  curvature <- abs(second_deriv) / (1 + first_deriv[-length(first_deriv)]^2)^(3/2)
  optimal_k <- which.max(curvature) + 1  # +1 because diff reduces length by 1
  
  # Plot if requested
  if (plot) {
    plot(k_values, wcss, 
         type = "b", 
         pch = 19, 
         frame = FALSE, 
         xlab = "Number of clusters (k)", 
         ylab = "Within-Cluster Sum of Squares (WCSS)",
         main = "Elbow Method for Optimal k")
    abline(v = optimal_k, lty = 2, col = "red")
    legend("topright", 
           legend = c("WCSS", paste("Suggested k =", optimal_k)),
           col = c("black", "red"), 
           lty = c(1, 2), 
           pch = c(19, NA))
  }
  
  # Return results
  list(wcss = wcss, 
       k_values = k_values, 
       optimal_k = optimal_k)
}


datTrait1 <- RNAseq %>% cor(method='s') %>% as.data.frame
module.trait.cor1 <- cor(net$MEs, datTrait1[rownames(net$MEs), ], use = "p") 
rownames(module.trait.cor1) <- gsub('ME', 'g', rownames(module.trait.cor1))
module.trait.cor1 <- module.trait.cor1[rownames(module.trait.cor1) !="g0", ]

datTrait2 <- ATACseq %>% cor(method='s') %>% as.data.frame
module.trait.cor2 <- cor(net2$MEs, datTrait2[rownames(net2$MEs), ], use = "p") 
rownames(module.trait.cor2) <- gsub('ME', 'p', rownames(module.trait.cor2))
module.trait.cor2 <- module.trait.cor2[rownames(module.trait.cor2) !="p0", ]

module.trait.cor <- rbind(module.trait.cor1, module.trait.cor2)
# module.trait.cor[module.trait.cor > 0.7] <- 0.7
# module.trait.cor[module.trait.cor < -0.7] <- -0.7
result <- elbow_method(module.trait.cor, max_k = 15)

annotation_mod <- data.frame(module=rownames(module.trait.cor), 
                             row.names = rownames(module.trait.cor)) %>%
  mutate(modType=gsub("[0-9]+","",module)) %>%
  mutate(modType=c(g="RNA",p="ATAC")[modType]) %>%
  mutate(GWAS=ifelse(module %in% c(sprintf("g%s", mods), 
                                   sprintf("p%s", mods2)), "Y", "N"))
modcols <- c(
  setNames(mdcols, sprintf("g%s", names(mdcols))),
  setNames(pmdcols, sprintf("p%s", names(pmdcols)))
)

dev.off()

library(pheatmap)
pdf('WGCNA.module-trait-associations.pdf', width = 8.27, height = 9)
annotation_cols <- list(module=modcols, GWAS=c(Y="#FC4E07", N="grey80"), 
                        modType=c(RNA="#00AFBB", ATAC="#E7B800"))

phtm <- pheatmap::pheatmap(module.trait.cor, cellwidth=10, fontsize_col=10,
                           cellheight=3, fontsize_row=3, border_color=NA, 
                           annotation_row = annotation_mod, 
                           annotation_colors = annotation_cols, 
                           cutree_rows=result$optimal_k)

mord <- phtm$tree_row$labels[phtm$tree_row$order] 
tord <- phtm$tree_col$labels[phtm$tree_col$order] 
annotation_class <- cutree(phtm$tree_row, result$optimal_k) %>% as.data.frame %>%
  rename(meta='.') %>% mutate(annotation_mod[rownames(.),])
annotation_stat <- annotation_class %>% select(modType, meta, GWAS) %>% table
x <- annotation_stat[,,'Y']
x <- colnames(x)[which(x[1,] > 0 & x[2,] > 0)]
# final_mod <- annotation_class %>% filter(meta %in% x & GWAS=="Y") %>% rownames()
final_mod <- annotation_class %>% filter(GWAS=="Y") %>% rownames()
gene_mod <- grep("^g", final_mod, value = T)
peak_mod <- grep("^p", final_mod, value = T)

annotation_cols$meta <- setNames(
  c(ggthemes::tableau_color_pal("Green-Orange-Teal")(12),
    ggthemes::tableau_color_pal("Red-Blue-Brown")(12)
  )[1:result$optimal_k], 1:result$optimal_k
)

pm <- ComplexHeatmap::pheatmap(module.trait.cor, cellwidth=10, fontsize_col=10,
                               cellheight=3, fontsize_row=3, border_color=NA, 
                               annotation_row = annotation_class, 
                               annotation_colors = annotation_cols, 
                               cutree_rows=result$optimal_k, 
                               labels_row = NULL)
ha <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=match(final_mod, rownames(module.trait.cor)), labels=final_mod, labels_gp=grid::gpar(fontsize = 8)))
print(pm + ha)

m1 <- signif_matrix
colnames(m1) <- sprintf("g%s", colnames(m1))
m2 <- signif_matrix_peak
colnames(m2) <- sprintf("p%s", colnames(m2))
m <- cbind(m1, m2) %>% as.data.frame()
m <- m[, mord]
m2 <- ((m != "") %>% t) + 0
mc <- mp <- cbind(annotation_class[colnames(m), ], m2[colnames(m), ])
mp[, rownames(m)] <- t(t(mp[, rownames(m)]) / colSums(mp[, rownames(m)]))

mlc <- mc %>% group_by(across(meta)) %>% 
  summarise(across(rownames(m), sum, na.rm = TRUE)) %>% 
  pivot_longer(cols = -meta, names_to = "trait", 
               values_to = "value") %>%
  mutate(meta = factor(meta, levels = mc$meta %>% unique)) %>%
  mutate(trait = factor(trait, levels = traits$variable))
mlp <- mp %>% group_by(across(meta)) %>% 
  summarise(across(rownames(m), sum, na.rm = TRUE)) %>% 
  pivot_longer(cols = -meta, names_to = "trait", 
               values_to = "value") %>%
  mutate(meta = factor(meta, levels = mp$meta %>% unique)) %>%
  mutate(trait = factor(trait, levels = traits$variable))

tcols <- setNames(traits$tcol, traits$variable)
p1 <- ggplot(mlc, aes(y = value, x = meta, fill = trait)) +
  geom_col(position = "stack") +
  labs(y = "Count", x = "Meta") +
  theme_minimal() +
  scale_fill_manual(values = tcols) +
  theme(legend.position = "bottom")
p2 <- ggplot(mlp, aes(y = value, x = meta, fill = trait)) +
  geom_col(position = "stack") +
  labs(y = "Proportion", x = "Meta") +
  theme_minimal() +
  scale_fill_manual(values = tcols) +
  theme(legend.position = "bottom")
print(p1)
print(p2)

pheatmap(module.trait.cor[gene_mod,], 
         border_color=NA, cellwidth=10, cellheight=10, 
         fontsize=10)
pheatmap(module.trait.cor[peak_mod, tord], 
         cluster_cols = F, fontsize=10, 
         border_color=NA, cellwidth=10, cellheight=10)
dev.off()

## Figure 2c

library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggsankey)

sankdat <- p2g %>% filter(peak %in% names(peaks) & gene %in% names(genes)) %>%
  mutate(pm = peaks[peak], gm = genes[gene]) %>% 
  select(gm, pm) 

colmap <- colorRampPalette(c('white', brewer.pal(n = 3, name = "Reds")))(100)
pdf('WGCNA.RNA-ATAC.pdf', width = 18, height = 8.27)
pheatmap(-log10(gps), 
         cellwidth = 10, cellheight = 10, color = colmap, 
         annotation_row = annotation_row, annotation_col = anno_col1, 
         annotation_colors = list(module = gmc, group = groupcol),
         display_numbers = gss, main = 'Genes', 
         border_color = 'white', cluster_rows = F)
pheatmap(-log10(pps), 
         cellwidth = 10, cellheight = 10, color = colmap, 
         annotation_row = annotation_row, annotation_col = anno_col2, 
         annotation_colors = list(module = pmc, group = groupcol),
         display_numbers = pss, main = 'OCRs', 
         border_color = 'white', cluster_rows = F)
pm1 <- pheatmap::pheatmap(-log10(pd1), 
                          cellwidth = 10, cellheight = 10, color = colmap, 
                          annotation_row = annotation_row, 
                          annotation_col = anno_col1 %>% filter(module %in% gms), 
                          annotation_colors = list(module = gmc[gms], group = groupcol),
                          display_numbers = gss[, gms], main = 'Genes', 
                          border_color = 'white', cluster_rows = F)
pm2 <- pheatmap::pheatmap(-log10(pd2), 
                          cellwidth = 10, cellheight = 10, color = colmap, 
                          annotation_row = annotation_row, 
                          annotation_col = anno_col2 %>% filter(module %in% pms), 
                          annotation_colors = list(module = pmc[pms], group = groupcol),
                          display_numbers = pss[, pms], main = 'OCRs', 
                          border_color = 'white', cluster_rows = F)

n1 <- pm1$tree_col$labels[pm1$tree_col$order]
n2 <- pm2$tree_col$labels[pm2$tree_col$order]
nodex <- c(sprintf("g%s", n1), sprintf("p%s", n2))
ncols <- c(setNames(gmc[n1], sprintf("g%s", names(gmc[n1]))), 
           setNames(pmc[n2], sprintf("p%s", names(pmc[n2]))))

pltd <- sankdat %>% table 
pltd <- pltd[n1, n2]
m1 <- pltd / rowSums(pltd)
m2 <- t(t(pltd) / colSums(pltd))
pheatmap(sqrt(m1), cellwidth = 15, cellheight = 15, 
         color = colmap, display_numbers = pltd, border_color = 'white', 
         cluster_rows = F, cluster_cols = F, main = 'by row')
pheatmap(sqrt(m2), cellwidth = 15, cellheight = 15, 
         color = colmap, display_numbers = pltd, border_color = 'white', 
         cluster_rows = F, cluster_cols = F, main="by column")

rmIdd <- which(m1 < 0.1 & m2 < 0.1, arr.ind=T) %>% as.data.frame() %>%
  mutate(gm=rownames(m1)[gm], pm=colnames(m1)[pm])

sankdat2 <- sankdat %>% dplyr::select(gm, pm) %>%
  filter(gm %in% gms) %>%
  filter(pm %in% pms)

sankdat3 <- sankdat %>% dplyr::select(gm, pm) %>%
  filter(gm %in% gsub("g","",gene_mod)) %>%
  filter(pm %in% gsub("p","",peak_mod))

p1 <- sankdat2 %>%
  mutate(id=sprintf("%s-%s", gm, pm)) %>%
  filter(!id %in% sprintf("%s-%s", rmIdd$gm, rmIdd$pm)) %>%
  mutate(gm = sprintf("g%s", gm)) %>%
  mutate(pm = sprintf("p%s", pm)) %>%
  make_long(pm, gm) %>% 
  mutate(label= node) %>%
  mutate(node = factor(node, levels=rev(nodex))) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node, 
             label = label,
             fill = node)) +
  geom_sankey(
    flow.alpha = 0.9,      # 桑基条带的不透明度
    # space = 50,          # 桑基节点间的距离
    # smooth = 6,          # 桑基条带的弯曲度
    width = 0.1, alpha = 1 # 桑基节点的宽度和不透明度
  ) +  
  geom_sankey_label(size = 3, color = 'white', fill = NA) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = ncols) + 
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle(sprintf("GWAS OCR2Gene links (n=%s)", nrow(sankdat2)))

p2 <- sankdat3 %>%
  mutate(id=sprintf("%s-%s", gm, pm)) %>%
  filter(!id %in% sprintf("%s-%s", rmIdd$gm, rmIdd$pm)) %>%
  mutate(gm = sprintf("g%s", gm)) %>%
  mutate(pm = sprintf("p%s", pm)) %>%
  make_long(pm, gm) %>% 
  mutate(label= node) %>%
  mutate(node = factor(node, levels=rev(nodex))) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node, 
             label = label,
             fill = node)) +
  geom_sankey(
    flow.alpha = 0.9,      # 桑基条带的不透明度
    # space = 50,          # 桑基节点间的距离
    # smooth = 6,          # 桑基条带的弯曲度
    width = 0.1, alpha = 1 # 桑基节点的宽度和不透明度
  ) +  
  geom_sankey_label(size = 3, color = 'white', fill = NA) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = ncols[c(gene_mod, peak_mod)]) + 
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle(sprintf("Final OCR2Gene links (n=%s)", nrow(sankdat3)))

patchwork::wrap_plots(list(p1=p1, p2=p2), ncol=2)
dev.off()

## Figure 2d

## See above

## Figure 2e/2f/2g

## motifs in p7 + p11
stat <- read.delim('peakmodule.motif.stat', header = F) %>%
  filter(V1 %in% c('p7','p11')) %>%
  mutate(V5 = sprintf("%s:%s", V2, V3))

mat <- stat %>% select(V1, V5, V4) %>% 
  mutate(V4 = -V4 / log(10), TF = mlist[V5,'V1']) %>%
  pivot_wider(names_from = V1, values_from = V4) %>%
  filter(p7 > 5 | p11 > 5) %>% 
  mutate(
    group = ifelse(p7 > p11, "asc", "desc")  # 分类标签
  ) %>%
  arrange(group,
          if_else(group == "asc", -p7, p11)
  ) %>%
  select(-group) %>% 
  mutate(TF = factor(TF, levels = unique(TF)))

final_list <- mlist[mat$V5, ] %>%
  filter(!duplicated(V2))
xid <- rownames(final_list)
rownames(final_list) <- final_list$V1

c2t <- final_list %>% select(V1, V2) %>% unique()
rownames(c2t) <- c2t$V2 

mat <- mat %>% filter(V5 %in% xid)

p <- ggplot(mat, aes(x = TF)) +
  geom_segment(aes(xend = TF, y = 0, yend = p7), color = pcols['7'], size = 0.5) +
  geom_point(aes(y = p7), color = pcols['7'], size = 3) +
  geom_segment(aes(xend = TF, y = 0, yend = -p11), color = pcols['11'], size = 0.5) +
  geom_point(aes(y = -p11), color = pcols['11'], size = 3) +
  coord_flip() +
  theme_bw()
print(p)
dev.off()

mfile <- 'known.motifs'
lines <- readLines(mfile)
motif_lines <- grep("^>", lines, value = TRUE)
motif_lines <- gsub(">","",gsub("\t.*","",motif_lines)) 
homers <- read_homer(mfile)
pcms <- lapply(1:length(homers), function(i){
  m <- homers[[i]]
  new("pcm", mat=m@motif, name=motif_lines[i])
})
pcms <- pcms[match(c2t$V2, motif_lines)]
names(pcms) <- c2t[sapply(pcms, function(m){m@name}),'V1']
nms <- names(pcms)
pcms <- lapply(names(pcms), function(n){
  m <- pcms[[n]]
  m@name <- n
  m
})
names(pcms) <- nms 

motifs <- lapply(1:length(homers), function(i){
  m <- homers[[i]]
  new("pfm", mat=m@motif, name=motif_lines[i])
})
motifs <- motifs[match(c2t$V2, motif_lines)]
names(motifs) <- c2t[sapply(motifs, function(m){m@name}),'V1']
nms <- names(motifs)
motifs <- lapply(names(motifs), function(n){
  m <- motifs[[n]]
  m@name <- n
  m
})
names(motifs) <- nms 

motifx <- pcms
pdf('Fig2e-g.motifs.pdf', height = 8.27, width = 8.27)
len <- length(motifx)
df <- data.frame(x=.5, y=(seq.int(len)-.5)/len, 
                 width=.75, height=1/(len+1))
df$motif <- motifx
library(ggplot2)
ggplot(df, aes(x=x, y=y, width=width, height=height, motif=motif)) +
  geom_motif(use.xy = TRUE) + theme_bw() + xlim(0, 1) + ylim(0, 1)

motifStack(motifx, layout = "tree")
motifStack(motifx, layout="phylog", f.phylog=.15, f.logo=0.25)

## cluster the motifs
hc <- clusterMotifs(motifx)
## convert the hclust to phylog object
library(ade4)
phylog <- hclust2phylog(hc)
## reorder the motifs by the order of hclust
leaves <- names(phylog$leaves)
motifx <- motifx[leaves]

## extract the motif signatures
motifSig <- motifSignature(motifx, phylog, cutoffPval = 0.0001, min.freq=1)

## get the signatures from object of motifSignature
sig <- signatures(motifSig)
## get the group color for each signature
gpCol <- sigColor(motifSig)

n <- length(motifx)

library(RColorBrewer)
color <- brewer.pal(12, "Set3")
## plot the logo stack with pile style.
motifPiles(phylog=phylog, pfms=motifx, pfms2=sig,
           col.tree=rep(color, each=5),
           col.leaves=rep(rev(color), each=5),
           col.pfms2=gpCol,
           r.anno=c(0.02, 0.03, 0.04),
           col.anno=list(sample(colors(), n),
                         sample(colors(), n),
                         sample(colors(), n)),
           motifScale="logarithmic",
           plotIndex=TRUE)

p1 <- mat %>% mutate(TF=factor(TF,levels=rev(leaves))) %>%
  ggplot(aes(x = TF)) +
  geom_segment(aes(xend = TF, y = 0, yend = p7), color = pcols['7'], size = 0.5) +
  geom_point(aes(y = p7), color = pcols['7'], size = 3) +
  geom_segment(aes(xend = TF, y = 0, yend = -p11), color = pcols['11'], size = 0.5) +
  geom_point(aes(y = -p11), color = pcols['11'], size = 3) +
  coord_flip() +
  theme_bw()

# library(patchwork)
library(cowplot)
print(plot_grid(p1, ncol = 3))

tissues <- c('Liver','Heart')
pltd <- RNAseq[final_list[leaves, 'V4'], tissues]
pltd <- apply(pltd, 1, function(x){x/mean(x)}) %>% t
colnames(pltd) <- tissues
rownames(pltd) <- gnames[rownames(pltd),1]
library(ComplexHeatmap)
p2 <- pheatmap(pltd, cluster_rows = F, cluster_cols = F, cellwidth = 10)
print(p2)
dev.off()

## Figure 2h

plist <- c(
  scan("peakmodule7.list", what = character()),
  scan("peakmodule11.list", what = character())
)

p2g <- read.delim('peak2genes.final.tsv') %>% 
  rename(target=gene) %>% select(peak, target)
fimo <- read.delim('final_fimo.tsv', header = F) %>%
  filter(V1 %in% final_list$V3 & V3 %in% plist) %>%
  rename(peak = V3, TF = V2) %>%
  select(peak, TF)


tf2target <- left_join(fimo, p2g, by = "peak", relationship='many-to-many') %>% 
  filter(moduleLabels[TF] == moduleLabels[target]) %>%
  select(TF, target) %>% unique()


expd <- RNAseq[unique(c(tf2target$TF, tf2target$target)), tissues] %>% na.omit()

pie_colors <- readRDS('tissue.color.rds')[colnames(expd)] %>% 
  as.character()
# pie_colors <- pcols[c(7,'11')]

links <- tf2target
colnames(links) <- c('from', 'to')
nodes <- unique(c(links$from, links$to))
net3.data <- links[, c('from', 'to')]
colnames(net3.data) <- c('source', 'target')
net3 <- graph.data.frame(net3.data, nodes, directed=T)

E(net3)$color <- adjustcolor('grey', alpha=.5) 
E(net3)$width <- 1 
E(net3)$arrow.size <- .4
E(net3)$arrow.width <- .4

V(net3)$label.cex <- 0.4
V(net3)$label.font <- 3
mx <- genemodules[nodes,'module'] %>% as.character()
mx[!mx %in% c('17', '28')] <- NA
V(net3)$frame.color <- ifelse(is.na(mx), 'grey80', gcols[mx])
# V(net3)$color <- ""
V(net3)$label <- gnames[nodes, 1]
nsize <- log1p(degree(net3)+1)*2
nsize[nsize > 5] <- 5
V(net3)$size <- nsize + 2

lab.dist <- rep(0.2, vcount(net3))
names(lab.dist) <- nodes
lab.dist[nodes %in% tf2target$TF] <- 0

pie.values <- lapply(nodes, function(x){
  if(x %in% rownames(expd)){
    v <- as.numeric(expd[x, ])
    v 
  }else{
    c(1,1)
  }
})

nshape <- rep("pie", vcount(net3)) ## 'circle'
names(nshape) <- nodes 

V(net3)$shape <- nshape

pdf("GRNs.pdf", width=8.27, heigh=8.27, pointsize=10)
set.seed(123456)
plot(net3, layout = layout_with_fr,
     vertex.label.degree=0, 
     vertex.label.dist=lab.dist, 
     vertex.pie=pie.values, 
     vertex.pie.border=NA, 
     vertex.pie.color=list(pie_colors), 
     main="GRN"
)
set.seed(123456)
plot(net3, layout = layout_with_kk,
     vertex.label.degree=0, 
     vertex.label.dist=lab.dist, 
     vertex.pie=pie.values, 
     vertex.pie.border=NA, 
     vertex.pie.color=list(pie_colors), 
     main="GRN"
)
set.seed(123456)
plot(net3, layout = layout_nicely,
     vertex.label.degree=0, 
     vertex.label.dist=lab.dist, 
     vertex.pie=pie.values, 
     vertex.pie.border=NA, 
     vertex.pie.color=list(pie_colors), 
     main="GRN"
)
dev.off()

## Figure 2i

go <- clusterProfiler::enrichGO(
  OrgDb = org.Mm.eg.db,
  gene = glist,
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  keyType = 'ENSEMBL',
  pAdjustMethod = 'fdr',
  ont = "BP"
)

## Figure 2j

path <- 'gga00071'
kgs <- kegg %>% filter(V1 == path & V3 != "") %>% pull(V3) %>% unique
pltd <- rowZscores(RNAseq[kgs,] %>% as.matrix(), limit=T)

library(pals)
library(ComplexHeatmap)
pdf("KEGG-gga00071.pdf", width=8.27, heigh=8.27, pointsize=10)
colmap <- rev(brewer.spectral(100))
colmap <- colorRampPalette(c("#5E4FA2","white","#9E0142"))(100)
p <- ComplexHeatmap::pheatmap(pltd, color = colmap, cellwidth=10, cellheight=10,
                              fontsize=10, labels_row=gnames[rownames(pltd),1])
print(p)
dev.off()

## Figure 2k

pid <- '00071'
library(pathview)
pathview(gene.data = gene_data, 
         pathway.id = pid,
         # cpd.data = cpd_data, 
         species = "gga", 
         limit= list(gene = 1, cpd = 1), 
         bins = list(gene = 50, cpd = 50), 
         low = list(gene = "#5E4FA2", cpd = "white"), 
         mid = list(gene = "white", cpd = "orange"), 
         high = list(gene = "#9E0142", cpd = "red"),
         kegg.dir = "./kegg/",
         out.suffix = "chicken", 
         kegg.native = T, 
         gene.idtype = 'kegg') 

