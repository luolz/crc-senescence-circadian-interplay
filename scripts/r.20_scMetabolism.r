rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./20_scMetabolism")){
  dir.create("./20_scMetabolism")
}
setwd("./20_scMetabolism")

load('../16_scRNA/04_SingleR/UMAP_singleR.Rdata')

# 加载R包
library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)

sc.metabolism.SeuratV5 <- function (obj, method = "VISION", imputation = F, ncores = 2, 
                                    metabolism.type = "KEGG") 
{
  countexp <- GetAssayData(obj, layer='counts')
  countexp <- data.frame(as.matrix(countexp))
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                       package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                           package = "scMetabolism")
  if (metabolism.type == "KEGG") {
    gmtFile <- signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if (metabolism.type == "REACTOME") {
    gmtFile <- signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  if (imputation == F) {
    countexp2 <- countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)
  }
  cat("Start quantify the metabolism activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  }
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                           nCores = ncores, plotStats = F)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  obj@assays$METABOLISM$score <- signature_exp
  obj
}

Idents(UMAP) <- "celltype"
countexp.Seurat <- sc.metabolism.SeuratV5(obj = UMAP,  #Seuratde单细胞object
                                          method = "AUCell", 
                                          imputation = F, 
                                          ncores = 2, 
                                          metabolism.type = "KEGG")

# 需要修改DimPlot.metabolism中的UMAP大小写
DimPlot.metabolismV5 <- function (obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, 
                                  size = 1) 
{
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  if (dimention.reduction.type == "umap") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
    umap.loc <- obj@reductions$umap@cell.embeddings
    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, 
    ]))
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = umap_1, 
                                                y = umap_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("UMAP 1") + ylab("UMAP 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  if (dimention.reduction.type == "tsne") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, 
    ]))
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, 
                                                y = tSNE_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("tSNE 1") + ylab("tSNE 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  plot
}

# check一下有哪些pathway
pathways <- countexp.Seurat@assays$METABOLISM$score
head(rownames(pathways))

write.csv(pathways, "01.pathways.csv")
input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:20]
write.csv(input.pathway, "02.input.pathway_TOP20.csv")

p <- DotPlot.metabolism(obj = countexp.Seurat,
                        pathway = input.pathway,
                        phenotype = "celltype", #这个参数需按需修改
                        norm = "y")
ggsave(p,file = "02.DotPlot.png",width = 6,height = 6)
ggsave(p,file = "02.DotPlot.pdf",width = 6,height = 6, family = "Times")

# # BoxPlot
# # 由于开发者默认吧obj定义为countexp.Seurat，所以还需要重新命名一下
# countexp.Seurat <- countexp.Seurat
# # input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:2]
# input.pathway<-c("Nitrogen metabolism",
#                  "Galactose metabolism")
# p <- BoxPlot.metabolism(obj = countexp.Seurat,
#                         pathway = input.pathway, 
#                         phenotype = "celltype", #这个参数需按需修改
#                         ncol = 1)
# ggsave(p,file = "03_BoxPlot.png",width = 13,height = 8)
# ggsave(p,file = "03_BoxPlot.pdf",width = 13,height = 8, family = "Times")
