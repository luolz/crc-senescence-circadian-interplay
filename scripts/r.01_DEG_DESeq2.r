rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists('./01_DEG_DESeq2')) {
  dir.create('./01_DEG_DESeq2')
}
setwd('./01_DEG_DESeq2')
library(tidyverse)
library(DESeq2)


count <- read.csv('../00_rawdata/01.count_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = F)
count <- round(count, digits = 0)
group <- read.csv('../00_rawdata/01.group_TCGA.COADREAD.csv', row.names = 1, check.names = F)
identical(colnames(count), rownames(group))
table(group$Type)
group$Type <- factor(group$Type, levels = c('Normal', 'Tumor'))

dds <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~Type)
dds <- dds[rowMeans(counts(dds)) > 1,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
save.image(file = 'DEGresult_TvsN.RData')
load('DEGresult_TvsN.RData')

res <- results(dds, contrast = c('Type', 'Tumor', 'Normal'))
head(res)

DEG <- as.data.frame(res)
DEG <- na.omit(DEG)

logFC_cutoff <- 2
DEG$change <- as.factor(ifelse(DEG$padj< 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                               ifelse(DEG$log2FoldChange > logFC_cutoff, 'UP', 'DOWN'), 'NOT'))
table(DEG$change)
# DOWN   NOT    UP 
# 1062 15357   994 
DEG_sig <- subset(DEG, DEG$change %in% c('UP', 'DOWN'))
dim(DEG_sig)
# 2056    7

DEG_write <- cbind(GeneSymbol = rownames(DEG), DEG)
DEG_write <- DEG_write[,-c(2,4,5)]
write.csv(DEG_write, file = '01.DEG_all.csv', quote = F, row.names = F)

DEG_sig_write <- cbind(GeneSymbol = rownames(DEG_sig), DEG_sig)
DEG_sig_write <- DEG_sig_write[,-c(2,4,5)]
write.csv(DEG_sig_write, file = '01.DEG_sig.csv', quote = F, row.names = F)


# Volcano plot ---------
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggrepel)

dat_rep <- rbind(
  head(DEG_sig[order(DEG_sig$log2FoldChange, decreasing = TRUE), ], 10),
  head(DEG_sig[order(DEG_sig$log2FoldChange, decreasing = FALSE), ], 10)
)

volcano_plot <- ggplot(data = DEG,
                       aes(x = log2FoldChange,
                           y = -log10(padj),
                           color = change)) +
  scale_color_manual(values = c('#66C2A5', 'darkgray', '#FFD700')) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0,1,5,10,20,50,100,200)) +
  geom_point(shape = 6,size = 1.5, alpha = 0.4, na.rm = TRUE) +
  theme_bw(base_size = 12, base_family = 'Times') +
  geom_vline(xintercept = c(-2, 2),
             lty = 4,
             col = '#00000050',
             lwd = 0.6) +
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = '#00000050',
             lwd = 0.6) +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.title = element_text(face = 'bold', color = 'black', size = 15),
        legend.text = element_text(face = 'bold', color = 'black', family = 'Times', size = 13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = 'bold', color = 'black', size = 15),
        axis.text.y = element_text(face = 'bold', color = 'black', size = 15),
        axis.title.x = element_text(face = 'bold', color = 'black', size = 15),
        axis.title.y = element_text(face = 'bold', color = 'black', size = 15),
        text = element_text(family = 'Times')) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, 'lines'),
    min.segment.length = 0,
    point.padding = unit(0.8, 'lines'),
    segment.color = 'black',
    show.legend = FALSE,
    family = 'Times') +
  labs(x = 'log2 (Fold Change)',
       y = '-log10(padj)')

pdf('01.volcano.pdf', w = 7, h = 5, family = 'Times')
volcano_plot
dev.off()
png('01.volcano.png', w = 7, h = 5, family = 'Times', units = 'in', res = 600)
volcano_plot
dev.off()

# Distribution as heatmap --------
library(ComplexHeatmap)
dat_rep1 <- rbind(
  head(DEG_sig[order(DEG_sig$log2FoldChange, decreasing = TRUE), ], 10),
  head(DEG_sig[order(DEG_sig$log2FoldChange, decreasing = FALSE), ], 10)
)
group1 <- group[order(group$Type, decreasing = TRUE),,drop=F]
x <- read.csv('../00_rawdata/01.fpkmlog2_TCGA.COADREAD_mRNA.csv', row.names = 1,check.names = F)
mat <- t(scale(t(x)))
mat <- pmin(pmax(mat, -1), 1)
mat <- mat[rownames(dat_rep1), rownames(group1)] %>% as.matrix()

pdf('01.heatmap.pdf', w = 8, h = 10, family = 'Times')
densityHeatmap(mat, title = 'Distribution as heatmap', ylab = ' ', height = unit(2.5, 'cm')) %v%
  HeatmapAnnotation(Group = group1$Type, col = list(Group = c('Tumor' = '#fed53b', 'Normal' = '#66C2A5'))) %v%
  Heatmap(mat,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = FALSE,
          show_row_names = TRUE,
          name = 'expression',
          height = unit(12, 'cm'),
          cluster_rows = FALSE,
          col = colorRampPalette(c('#66C2A5', 'white', '#fed53b'))(100))
dev.off()

png('01.heatmap.png', w = 8, h = 10, units = 'in', res = 600, family = 'Times')
densityHeatmap(mat, title = 'Distribution as heatmap', ylab = ' ', height = unit(2.5, 'cm')) %v%
  HeatmapAnnotation(Group = group1$Type, col = list(Group = c('Tumor' = '#fed53b', 'Normal' = '#66C2A5'))) %v%
  Heatmap(mat,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = FALSE,
          show_row_names = TRUE,
          name = 'expression',
          height = unit(12, 'cm'),
          cluster_rows = FALSE,
          col = colorRampPalette(c('#66C2A5', 'white', '#fed53b'))(100))
dev.off()
