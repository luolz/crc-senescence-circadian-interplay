rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists("09_GSEA")) {dir.create("09_GSEA")}
setwd("09_GSEA")

# library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(tidyverse)
library(GseaVis)
library(DESeq2)
#高低风险组------------------------------------
library(GSEABase)
library(GSVA)
library(msigdbr)
library(limma)

riskScore<- read.csv('../05_risk_train/03.risk.csv')
riskScore$group <- as.factor(riskScore$risk)
colnames(riskScore)[colnames(riskScore) == "id"] <- "sample"
risk <- riskScore[ , c("sample", "group")]
table(risk$group)
# high risk  low risk 
# 72       452 
risk$group <- as.factor(risk$group)
risk$group <- gsub("high risk", "High", risk$group)
risk$group <- gsub("low risk", "Low", risk$group)

dat <- read.csv(file = '../00_rawdata/01.count_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = FALSE)
dat <- dat[, colnames(dat) %in% risk$sample]

risk$group <-factor(risk$group)

# DEseq2 ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = dat, colData = risk, design = ~ group)
dds_norm <- vst(dds)
# # 提取标准化后的数据------
Controlized_counts <- assay(dds_norm)
# write.csv(Controlized_counts,file = "00.Controlized_counts.csv")
dds$group <- relevel(dds$group, ref ="Low")   #指定对照组
dds$group
dds <- DESeq(dds)#正式进行差异分析
res <- results(dds,contrast = c("group","High","Low")) 
res = res[order(res$padj),]
gene_df <- res %>%as.data.frame()
gene_df$rowname<- rownames(gene_df)
gene <- gene_df$rowname
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=gene_df$log2FoldChange,
                      SYMBOL = gene_df$rowname)
data_all <- merge(gene_df,gene,by="SYMBOL") %>% na.omit(.)
data_all_sort <- arrange(data_all, desc(logFC))


gmt <- read.gmt("/data/nas1/zhuxuying/common/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt") #读gmt文件

genelist <- data_all_sort$logFC # 把foldchange按照从大到小提取出来
names(genelist) <- data_all_sort$SYMBOL # 给上面提取的foldchange加上对应上ENTREZID
genelist

## GSEA分析
res_GSEA <- GSEA(genelist, TERM2GENE = gmt, pvalueCutoff = 1, eps = 0)
sortGESA <- data.frame(res_GSEA)
sortGESA <- sortGESA[order(sortGESA$pvalue),]#按照enrichment score从高到低排序
if(length(sortGESA$ID) > 5) {
  paths <- rownames(sortGESA[c(1: 5), ]) # 选取需要展示的通路ID
}else {
  paths <- sortGESA$ID
}

## 保存GSEA结果
write.csv(sortGESA, file = paste0("gsea_res.csv"))
sortGESA <- sortGESA[sortGESA$pvalue < 0.05, ]
dim(sortGESA) #79 11
write.csv(sortGESA, file = paste0("gsea_res_sig.csv"))

## 多通路展示
p <- gseaNb(object = res_GSEA,
            geneSetID = paths,
            subPlot = 2,
            termWidth = 45,
            # legend.position = c(0.72,0.8),
            addPval = T,
            rmHt = F,
            pvalX = 0.99,
            pvalY = 0.99, 
            newGsea = F,
            curveCol = c("#B1485B",
                         "#3288bdA9",
                         "#FF6F61",
                         "#FFCB65",
                         "#66BFBF"))

pdf(file = paste0("01.multi_paths_gsea.pdf"), w=13, h=9)
print(p)
dev.off()
png(file = paste0("01.multi_paths_gsea.png"), w=13, h=9, units = 'in', res = 600)
print(p)
dev.off()

##############单因素GSEA
library(magrittr)
library(stringr)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

gene <- readRDS('../05_risk_train/inputgene.rds')
hub_gene <- gene

exprSet <- dat

all.gene <- AnnotationDbi::select(org.Hs.eg.db, rownames(exprSet), "ENTREZID", "SYMBOL")
exprSet <- merge(exprSet, all.gene, by.x = "row.names", by.y  = "SYMBOL") %>% 
  dplyr::distinct(ENTREZID, .keep_all = T) %>% na.omit()
rownames(exprSet) <- NULL
exprSet <- exprSet %>% tibble::column_to_rownames(var = "ENTREZID") %>% dplyr::select(-Row.names)

hub.gene = AnnotationDbi::select(org.Hs.eg.db, hub_gene, "ENTREZID", "SYMBOL")
hub_gene <- hub.gene$ENTREZID
hub_dat <- exprSet[hub_gene, ] %>% as.data.frame()
hub_dat2 <- cbind(GeneID = rownames(hub_dat), hub_dat)
write.csv(hub_dat2, file = "02.Hub_gene.csv")

library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(ggrepel)
library(enrichplot)
library(aplot)

gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          text = element_text(face = "bold", family = "Times"))
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank(),
                              text = element_text(face = "bold", family = "Times")) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_brewer(palette = "Set1")
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

# 设置biomaRt连接
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gmt <- read.gmt("/data/nas1/zhuxuying/common/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt") #读gmt文件
gene_symbols <- unique(gmt$gene)
# 检索基因信息
genes <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id'),
               filters = 'hgnc_symbol',
               values = gene_symbols,
               mart = ensembl)
gmt <- gmt %>%
  left_join(genes, by = c("gene" = "hgnc_symbol"))

df.kegg = gmt[,c(1,3)]
colnames(df.kegg) <- c("gs_name","entrez_gene")

df.exp <- exprSet %>% t %>% as.data.frame() 

if (!dir.exists("03.GSEA")) {dir.create("03.GSEA")}
if (!dir.exists("03.GSEA_Res")) {dir.create("03.GSEA_Res")}
set.seed(783)
lapply(1:length(hub_gene), function(i){
  hub = hub_dat2$GeneID[i] %>% as.character()
  hub.exp = df.exp[[hub]]
  hub.cor = cor(df.exp, hub.exp, method = "spearman") %>% as.data.frame %>% na.omit
  hub.coreff = hub.cor[[1]]
  names(hub.coreff) = rownames(hub.cor)
  hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
  res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.05, seed = 1, pAdjustMethod = "BH", eps = 0)
  write.table(res.kegg, file = paste0("03.GSEA_Res/0", i, ".", hub.gene$SYMBOL[which(hub.gene$ENTREZID==hub)], ".res.xls"),
              sep = "\t", row.names = T, col.names = T, quote = F)
  top.kegg = res.kegg@result
  top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]$Description[1:5]
  p = gsea.plot(res.kegg, top.kegg, hub.gene$SYMBOL[which(hub.gene$ENTREZID==hub)])
  fn1 = paste0("03.GSEA/", sprintf("%02d",i), ".", hub.gene$SYMBOL[which(hub.gene$ENTREZID==hub)], ".png")
  fn2 = paste0("03.GSEA/", sprintf("%02d",i), ".", hub.gene$SYMBOL[which(hub.gene$ENTREZID==hub)], ".pdf")
  ggsave(fn1, p, width = 12, height = 6, units = "in", limitsize = 300)
  ggsave(fn2, p, width = 12, height = 6, units = "in", limitsize = 300)
  return(0)
})




