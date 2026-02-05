rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists("03_enrichment")) {dir.create("03_enrichment")}
setwd("03_enrichment")

library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(writexl)

#富集
hub_gene_symbol <- read.csv("../02_venn/01.Intersection_gene.csv")
gene_symbol <- clusterProfiler::bitr(
  geneID = hub_gene_symbol$symbol,  #需要修改
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

enrich_go <- enrichGO(gene = gene_symbol$ENTREZID, 
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID", 
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = TRUE)
GO_ALL <- enrich_go@result
GO_ALL <- GO_ALL[GO_ALL$pvalue < 0.05,]
BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count,decreasing = T),]
CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count,decreasing = T),]
MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count,decreasing = T),]
paste0("得到",dim(GO_ALL)[[1]],"个结果，其中",dim(BP)[[1]],"个生物学过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
# "得到595个结果，其中552个生物学过程，12个细胞组分，31个分子功能"

library(writexl)
write.csv(GO_ALL,file = "00.go_ALL.csv")
write.csv(as.data.frame(BP), '00.GO_BP.csv')
write.csv(as.data.frame(CC), '00.GO_CC.csv')
write.csv(as.data.frame(MF), '00.GO_MF.csv')

kk <- enrichKEGG(gene = gene_symbol$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1, 
                 qvalueCutoff = 1
)

kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

kegg_result <- kk@result
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05,]
dim(hh)
# 12     9
hh <- hh[order(hh$Count,decreasing = T),]
hh <- na.omit(hh)    
write.csv(hh, '02.KEGG.csv')
dim(kk@result)

ekegg_all <- kk
ego_all <- enrich_go

pal <- c('#7bc4e2', '#acd372', '#fbb05b', '#ed6ca4')
pal <- c('#eaa052', '#b74147', '#90ad5b', '#23929c')
pal <- c('#c3e1e6', '#f3dfb7', '#dcc6dc', '#96c38e')

# 使用setReadable将GO富集结果中的基因ID转换为基因符号
ego_readable <- setReadable(ego_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
GO <- as.data.frame(ego_readable)

# 使用setReadable将KEGG富集结果中的基因ID转换为基因符号
ekegg_readable <- setReadable(ekegg_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG <- as.data.frame(ekegg_readable)

# 选择top5，根据pvalue排序，相同pvalue值取富集到基因最多的通路
use_pathway <- group_by(GO, ONTOLOGY) %>%
  top_n(5, wt = -pvalue) %>%
  # group_by(p.adjust) %>%
  # top_n(1, wt = Count) %>%
  rbind(
    top_n(KEGG, 5, -pvalue) %>%
      # group_by(p.adjust) %>%
      # top_n(1, wt = Count) %>%
      mutate(ONTOLOGY = 'KEGG')
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, pvalue) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

use_pathway <- use_pathway %>%
  mutate(geneID = sapply(str_split(geneID, "/"), function(x) paste(head(x, 20), collapse = "/")))
# 构造左侧标记数据
# 左侧分类标签和基因数量点图的宽度
width <- 0.5
# x 轴长度
xaxis_max <- max(-log10(use_pathway$pvalue)) + 1
# 左侧分类标签数据
# 先统计每个 ONTOLOGY 的数量
ontology_table <- table(use_pathway$ONTOLOGY)
ontology_df <- as.data.frame(ontology_table)
colnames(ontology_df) <- c("ONTOLOGY", "n")

# 计算 ymin 和 ymax
ontology_df$ymax <- cumsum(ontology_df$n)
ontology_df$ymin <- c(0, head(ontology_df$ymax, -1)) + 0.6
ontology_df$ymax <- ontology_df$ymax + 0.4

# 添加 xmin 和 xmax
ontology_df$xmin <- -3 * width
ontology_df$xmax <- -2 * width

# 最后得到 rect.data
rect.data <- ontology_df

# 绘制富集通路图
library(ggprism)
plot_enrichment <- function() {
  use_pathway_top_genes <- use_pathway %>%
    group_by(ONTOLOGY, Description) %>%
    slice_head(n = 10) %>%  # 每个通路最多选取10个基因
    ungroup()
  
  p <- use_pathway_top_genes %>%
    ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
    geom_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(
      aes(x = 0.05, label = Description),
      hjust = 0, size = 5
    ) +
    geom_text(
      aes(x = 0.1, label = geneID, colour = ONTOLOGY),
      hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic',
      show.legend = FALSE
    ) +
    # 基因数量
    geom_point(
      aes(x = -width, size = Count),
      shape = 21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name = 'Count', range = c(5, 12)) +
    # 分类标签 使用 geom_rect 替代 geom_round_rect
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes = FALSE,
      angle = 90,   # 旋转90度
      hjust = 0.5,  # 水平居中
      vjust = 0.5,  # 垂直居中
      size = 4
    ) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    labs(y = NULL) +
    scale_fill_manual(name = 'Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max, 2),
      expand = expansion(c(0, 0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  return(p)
}
enrichment_plot <- plot_enrichment()
print(enrichment_plot)

ggsave("03.GO_KEGG_dot.pdf", width = length(use_pathway$ONTOLOGY)/4*3, height = length(use_pathway$ONTOLOGY)/20*9, enrichment_plot,family = "Times")
png("03.GO_KEGG_dot.png", width = length(use_pathway$ONTOLOGY)/4*3, height = length(use_pathway$ONTOLOGY)/20*9,units = "in",res = 600,family = "Times")
print(enrichment_plot)
dev.off()
save.image("GO_KEGG.RData")





