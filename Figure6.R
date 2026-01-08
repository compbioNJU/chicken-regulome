library(WGCNA)
library(patchwork)
library(ggplot2)
library(universalmotif)
library(motifStack)
library(igraph)
library(foreach)
library(enrichplot)
library(clusterProfiler)
library(paletteer)
library(dplyr)
library(tidyr)
library(tidyverse)
library(cowplot)

## Figure 6a / 6b / 6c
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
pdf('Fig6a-c.motifs.pdf', height = 8.27, width = 8.27)
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

## Figure 6d
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

pdf("Fig6d.GRNs.pdf", width=8.27, heigh=8.27, pointsize=10)
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

## Figure 6e
go <- clusterProfiler::enrichGO(
  OrgDb = org.Mm.eg.db,
  gene = glist,
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  keyType = 'ENSEMBL',
  pAdjustMethod = 'fdr',
  ont = "BP"
)

## Figure 6f
path <- 'gga00071'
kgs <- kegg %>% filter(V1 == path & V3 != "") %>% pull(V3) %>% unique
pltd <- rowZscores(RNAseq[kgs,] %>% as.matrix(), limit=T)

library(pals)
library(ComplexHeatmap)
pdf("Fig6f.KEGG-gga00071.pdf", width=8.27, heigh=8.27, pointsize=10)
colmap <- rev(brewer.spectral(100))
colmap <- colorRampPalette(c("#5E4FA2","white","#9E0142"))(100)
p <- ComplexHeatmap::pheatmap(pltd, color = colmap, cellwidth=10, cellheight=10,
                              fontsize=10, labels_row=gnames[rownames(pltd),1])
print(p)
dev.off()

## Figure 6g
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
