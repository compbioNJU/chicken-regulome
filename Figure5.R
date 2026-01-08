## Figure 5a / 5d / 5f

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
pdf('Fig5.WGCNA.gene.softPowers.pdf', width = 8.27, height = 8.27)
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

pdf('Fig5a.WGCNA.gene.dendroAndColors.pdf', width = 8.27, height = 4)
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

pdf('Fig5f.WGCNA.full_geneNetwork.pdf', width = 12, height = 12)
plotModuleUMAP(genemodules, TOM = TOM)
dev.off()

# 4. GWAS Enrichment Analysis
# Use Hypergeometric test (phyper) to check overlap between modules and GWAS traits
traits <- readRDS('trait.meta.rds')
bed <- read.table('gwas.assoc.sig.full.annotation.bed', sep="\t")
N <- nrow(RNAseq) # Total background genes
pout <- NULL

# Iterate through traits and modules to calculate enrichment p-values
for(trait in rownames(traits)){
    # (Enrichment logic: calculates overlap between GWAS genes and WGCNA modules)
    # ...
  g1 <- bed %>% filter(V5 == trait) %>% pull(V7) %>%
    strsplit(split = '/') %>% unlist() %>% 
    gsub(";.*", "", .) %>% unique 
  n1 <- length(g1)
  for(mod in levels(genemodules$module)){
    g2 <- genemodules %>% filter(module %in% mod) %>% pull(gene)
    n2 <- length(g2)
    g3 <- intersect(g1, g2)
    n3 <- length(g3)
    p <- phyper(q = n3 - 1, 
                m = n1, 
                n = N - n1, 
                k = n2, 
                lower.tail = FALSE) 
    pout <- data.frame(id=trait, trait=traits[trait, 'variable'], 
                       module = mod, pvalue = p) %>%
      rbind(pout)
  }
}

# 5. Visualization of Enrichment Result
# Create heatmap showing -log10(p-value) of module-trait enrichment
library(ComplexHeatmap)
pltd <- -log10(pvalue_matrix) %>% as.matrix()
signif_matrix <- ... # (Assign stars * based on p-value significance)

pdf('Fig5d.WGCNA.enrichment.pdf', width = 16, height = 8.27)
pheatmap(pltd %>% sqrt(), annotation_colors = list(module = mdcols, group = groupcol),
         display_numbers = signif_matrix, cluster_rows = F)
dev.off()

# 6. Final Network Plotting
# Visualize specific modules and label hub genes/GWAS-related genes
pdf('Fig5f.WGCNA.geneNetwork.pdf', width = 12, height = 12)
plotModuleUMAP(genemodules, TOM = TOM, label_hubs = 1, label_genes = label_genes)
dev.off()

## Figure 5b
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
pdf('Fig5b.WGCNA.module-trait-associations.pdf', width = 8.27, height = 9)
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

## Figure 5c
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggsankey)

sankdat <- p2g %>% filter(peak %in% names(peaks) & gene %in% names(genes)) %>%
  mutate(pm = peaks[peak], gm = genes[gene]) %>% 
  select(gm, pm) 

colmap <- colorRampPalette(c('white', brewer.pal(n = 3, name = "Reds")))(100)
pdf('Fig5c.WGCNA.RNA-ATAC.pdf', width = 18, height = 8.27)
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
