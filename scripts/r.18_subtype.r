rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./18_subtype")){
  dir.create("./18_subtype")
}
setwd("./18_subtype")
load('../16_scRNA/04_SingleR/UMAP_singleR.Rdata')
library(Seurat)
table(UMAP$celltype)
# 关键细胞是fibroblasts
sce.car <- UMAP[,UMAP$celltype%in%c("B cells"),]
## PCA降维分析
## scale标准化——很多基因的reads数多达上千，但是有很多基因reads很少（ScaleData作用：是将全部基因reads数据压缩为0-1之间的数值）
sce.car.nor.sca <- ScaleData(sce.car)
sce.car.norm.pca <- RunPCA(sce.car.nor.sca, features = VariableFeatures(object = sce.car.nor.sca), npcs = 50)
## 首先找到最佳聚类数
## JackStraw函数的功能/作用：
## 随机排列数据的子集，并计算这些“随机”基因的预测PCA分数。
## 然后将“随机”基因的PCA分数与观察到的PCA分数进行比较，以确定统计显著性。
## 最终结果是每个基因与每个主成分相关的p值
## 简而言之: JackStraw()函数可以计算出每个主成分中各基因的P值，用于判断哪些主成分更具有统计学意义
sce.car.norm.pca <- JackStraw(sce.car.norm.pca, num.replicate = 100, dims = 50) # num.replicate : 要执行的复制采样数；dims : 计算重要性的PC数
## ScoreJackStraw()函数：用于量化主成分的显著性强度，富含低P值基因较多的主成分更有统计学意义
sce.car.norm.pca <- ScoreJackStraw(sce.car.norm.pca, dims = 1:50)
system.time(save(sce.car.norm.pca, file = "sce.car.norm.pca.Jack.Rdata"))

# load('./sce.car.norm.pca.Jack.Rdata')
## JackStrawPlot()函数：可视化比较每个主成分的 p 值分布和均匀分布（虚线）
plot_pca <- JackStrawPlot(sce.car.norm.pca, dims = 1:50)
plot_pca
pdf("01.pca_cluster.pdf", width = 12, height = 8,family = "Times")
plot_pca
dev.off()
png("01.pca_cluster.png", width = 12, height = 8,res = 600,units = 'in',family = "Times")
plot_pca
dev.off()

plot_elbow <- ElbowPlot(sce.car.norm.pca, ndims = 50)
plot_elbow
pdf("02.pca_sd.pdf", width = 6, height = 5,family = "Times")
plot_elbow
dev.off()
png("02.pca_sd.png", width = 6, height = 5,res = 600,units = 'in',family = "Times")
plot_elbow
dev.off()

# load('./sce.car.norm.pca.Jack.Rdata')
## FindNeighbors()函数作用：计算邻接距离
sce.car.norm.pca.clu <- FindNeighbors(sce.car.norm.pca, dims = 1:30)
## FindClusters()函数作用：以迭代方式将细胞分组在一起,resolution参数表示区分细胞群的分辨率，resolution越大，分的群也就越多。
sce.car.norm.pca.clu <- FindClusters(sce.car.norm.pca.clu, resolution = 0.1) # 0.4

#------------- non-linear dimensional reduction (UMAP) --------
UMAP<-RunTSNE(sce.car.norm.pca.clu,dims = 1:30) # tSNE

col = c("#E7298A","#66A61E","#E6AB02","#A6761D","#666666",
        "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33"
        ,"#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB"
        ,"#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",
        "#E41A1C","#377EB8","#984EA3","#FF7F00","#FFFF33","#A65628",
        "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
        "#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",
        "#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
        "#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1B9E77","#D95F02"
)

pdf("03.UMAP-30-res.pdf",height = 5,width = 6, family = "Times")
DimPlot(UMAP, reduction = "tsne",label = T,cols = col)
dev.off()
png("03.UMAP-30-res.png",height = 5,width = 6,units = 'in',res = 600, family = "Times")
DimPlot(UMAP, reduction = "tsne",label = T,cols = col)
dev.off()

pdf('04.UMAP.group.pdf',w=10,h=5,family = "Times")
DimPlot(UMAP,reduction = 'tsne',split.by  = "group",cols = col)
dev.off()
png('04.UMAP.group.png',w=10,h=5,units = 'in',res = 600, family = "Times")
DimPlot(UMAP,reduction = 'tsne',split.by  = "group",cols = col)
dev.off()

system.time(save(UMAP, file = "sce.car_singleR.Rdata"))  

load("sce.car_singleR.Rdata")
library(Seurat)
library(monocle)
library(tidyverse)
s_cmb<- UMAP #仅关键细胞数据
DefaultAssay(s_cmb) <- "RNA"
matrix = GetAssayData(s_cmb, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(matrix),gene_short_name=rownames(matrix))
rownames(feature_ann) = rownames(matrix)
fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = s_cmb@meta.data
pd = new("AnnotatedDataFrame", data =sample_ann)
cds = newCellDataSet(matrix, phenoData =pd, featureData=fd, expressionFamily=negbinomial.size())
rm(s_cmb)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
WGCNA::collectGarbage()

mycds2 <- reduceDimension(cds, max_components = 2,verbose = T,
                          method = 'log')
save(mycds2, file = "mycds_order.RData")

library(monocle)
# devtools::load_all('/data/nas2/software/monocle/')
# 拟时间轴轨迹构建和在拟时间内排列细胞
load( "mycds_order.RData")
mycds3 <- orderCells(mycds2) 
save(mycds3, file = "mycds_order3.RData")

load( "mycds_order3.RData")
library(patchwork)
pdf("05.monocle.pdf",width = 10,height = 8,onefile=FALSE,family='Times')
f1<-plot_cell_trajectory(mycds3, color_by = "Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  theme(text = element_text(size = 14))
f2<-plot_cell_trajectory(mycds3 , color_by = "group",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("#80B1D3","#FB8072")) +
  theme(text = element_text(size = 14))+
  facet_wrap(~group,nrow = 1)
f3<-plot_cell_trajectory(mycds3, color_by = "seurat_clusters",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("#FB8072","#80B1D3","#FDB462","#337AB7"
                                ,"#BC80BD","#00BFFF","#FFED6F","#1B9E77","#D95F02","#7570B3"
                                ,"#A3D8A3","#FF6F61","#56BCD3","#8E9BAE","#D09075","#FFD8B1"
  )) +
  theme(text = element_text(size = 14))

f4<-plot_cell_trajectory(mycds3, color_by = "State",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("#FB8072","#80B1D3","#FDB462","#FCCDE5","#FF69B4"
                                ,"#BC80BD","#00BFFF","#FFED6F","#1B9E77","#D95F02","#7570B3"
                                ,"#A3D8A3","#FF6F61","#56BCD3","#8E9BAE","#D09075","#FFD8B1"
                                ,"#337AB7")) +
  theme(text = element_text(size = 14))

f1+f3+f4+f2
dev.off()

png("05.monocle.png",width = 10,height = 8,units = 'in',res = 600, family = "Times")
f1<-plot_cell_trajectory(mycds3, color_by = "Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  theme(text = element_text(size = 14))
f2<-plot_cell_trajectory(mycds3 , color_by = "group",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("#80B1D3","#FB8072")) +
  theme(text = element_text(size = 14))+
  facet_wrap(~group,nrow = 1)
f3<-plot_cell_trajectory(mycds3, color_by = "seurat_clusters",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("#FB8072","#80B1D3","#FDB462","#337AB7"
                                ,"#BC80BD","#00BFFF","#FFED6F","#1B9E77","#D95F02","#7570B3"
                                ,"#A3D8A3","#FF6F61","#56BCD3","#8E9BAE","#D09075","#FFD8B1"
                                )) +
  theme(text = element_text(size = 14))

f4<-plot_cell_trajectory(mycds3, color_by = "State",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("#FB8072","#80B1D3","#FDB462","#FCCDE5","#FF69B4"
                                ,"#BC80BD","#00BFFF","#FFED6F","#1B9E77","#D95F02","#7570B3"
                                ,"#A3D8A3","#FF6F61","#56BCD3","#8E9BAE","#D09075","#FFD8B1"
                                ,"#337AB7")) +
  theme(text = element_text(size = 14))

f1+f3+f4+f2
dev.off()

load( "mycds_order3.RData")
gene <- readRDS('../05_risk_train/inputgene.rds')

library(tidyverse)
my_cds_subset = mycds3
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene,],# num_clusters = 2, # add_annotation_col = ac,
                                                 cluster_rows = F,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

pdf('06.pseudotime_hubgene.pdf',w=4,h=3,family='Times')
print(my_pseudotime_cluster)
dev.off()
png('06.pseudotime_hubgene.png',w=4,h=3,units = 'in',res = 300,family='Times')
print(my_pseudotime_cluster)
dev.off()
