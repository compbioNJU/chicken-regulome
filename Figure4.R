## Figure 4a


## Figure 4b


## Figure 4c


## Figure 4d


## Figure 4e / 4k
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


## Figure 4f


## Figure 4g
# from WU browser screenshot

## Figure 4h

data <- fromJSON(sprintf("chickengtex/%s.json", gene))$result
names(data$y) <- data$name

tissues <- c("Ileum","Lung","Spleen","Jejunum","Cecum",
             "Adipose","Heart","Duodenum","Liver","Muscle",
             "Bursa","Thymus","Small_intestine","Brain")
expdat <- data$y[tissues]

o <- lapply(expdat, median) %>% as.data.frame %>% mutate(gene=gname, .before=Ileum)
pdat <- rbind(pdat, o)

o <- lapply(expdat, mean) %>% as.data.frame %>% mutate(gene=gname, .before=Ileum)
pdat2 <- rbind(pdat2, o)

df <- tibble(
    tissue = rep(names(expdat), lengths(expdat)),
    expression = unlist(expdat, use.names = FALSE)
) 
df_all <- rbind(df_all, df %>% mutate(gene = gname))
ord <- df %>% group_by(tissue) %>% summarize(m = mean(expression)) %>% 
    arrange(-m) %>% pull(tissue)
df <- df %>% mutate(tissue = factor(tissue, levels = ord))

p <- ggviolin(df, x = "tissue", y = "expression", 
              title = gname, fill = "#00AFBB", # palette = "jco",
              add = "boxplot", add.params = list(fill = "white")) + theme(
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)  # 垂直排列
              ) 
print(p)

## Figure 4i
raw <- fromJSON('APAF1.eQTLs.json')$result$data %>%
    filter(tissue %in% tissues & tssDistance >= 0 & 
               tssDistance < genelen)
data <- raw %>% filter(tissue == 'Spleen') %>% 
    select(rsId, tssDistance, alleleFreq, 
           pvalue, effectSize, alt, ref) 
SNP <- data$tssDistance

legends <- list(list(labels=c("ref", "alt"), 
                     fill=c("#d581b7","#83d3ad")))

sample2.gr$color <- ifelse(data$effectSize > 0, "#d60b6f","#26755d")
sample2.gr$score <- -log10(data$pvalue)

pie_data <- table(raw$tissue) %>% melt
colnames(pie_data) <- c('category', 'value')
pie_data <- pie_data %>% mutate(percentage=value * 100 / sum(value)) %>%
    mutate(label_pos = cumsum(percentage) - percentage / 2)
p <- ggplot(pie_data, aes(x = 2, y = percentage, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    geom_text(aes(y = label_pos, label = value),
              color = "white", size = 3, fontface = "bold") +
    coord_polar("y", start = 0) +
    xlim(0.5, 2.5) + theme_void() +
    theme(legend.position = "bottom")
print(p)

## Figure 4j
sample.gr <- sample2.gr <- GRanges("chr1", IRanges(SNP, width=1, names=data$rsId))
sample.gr$score <- NULL 
sample.gr$label <- NULL
sample.gr$node.label.col <- NULL
sample.gr$value1 <- data$alleleFreq*100
sample.gr$value2 <- 100 - sample.gr$value1 
sample.gr$color <- rep(list(c("#d581b7","#83d3ad")), length(SNP))
sample.gr$border <- "gray30"
lolliplot(sample.gr, features, type="pie")

lolliplot(list(A=sample.gr, B=sample2.gr), 
          list(x=features, y=features), 
          type=c("pie", "pin"), legend=legends)




## Figure 4l
# from WU browser screenshot

## Figure 4m


