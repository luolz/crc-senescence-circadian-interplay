# 0. 数据整理--------###########单细胞部分是112服务器跑的，考虑到版本等问题假设结果不一致请到112服务器找到同目录进行
rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists('./16_scRNA')){
  dir.create('./16_scRNA')
}
setwd('./16_scRNA/')

library(Seurat)
library(stringr)

temp1 <- fread('GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz',data.table = F)
temp1[1:4,1:4]
rownames(temp1)=temp1[,1]
temp1=temp1[,-1]
scRNA <- CreateSeuratObject(counts = temp1,
                            min.features = 200,
                            min.cells = 3, 
                            project = "scRNA_project",)

table(scRNA$orig.ident)
head(scRNA@meta.data)

scRNA$group<-ifelse(scRNA$orig.ident=='SMC01-N'|scRNA$orig.ident=='SMC02-N'|
                      scRNA$orig.ident=='SMC03-N'|scRNA$orig.ident=='SMC04-N'|
                      scRNA$orig.ident=='SMC05-N'|scRNA$orig.ident=='SMC06-N'|
                      scRNA$orig.ident=='SMC07-N'|scRNA$orig.ident=='SMC08-N'|
                      scRNA$orig.ident=='SMC09-N'|scRNA$orig.ident=='SMC10-N'
                    ,'Normal','Tumor')
table(scRNA$group)
# Normal  Tumor 
# 16404  47285
###线粒体基因表达占比(正则表达式，表示以MT-开头；scRNA[["percent.mt"]]这种写法会在meta.data矩阵加上一列)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  

###计算核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
head(scRNA@meta.data)
# orig.ident nCount_RNA nFeature_RNA group Patients_ID percent.mt percent.rb
# SMC01-T_AAACCTGCATACGCCG    SMC01-T      38052         4866 Tumor     SMC01-T  12.921791   23.12625
# SMC01-T_AAACCTGGTCGCATAT    SMC01-T      33750         5268 Tumor     SMC01-T   8.761481   23.15556
# SMC01-T_AAACCTGTCCCTTGCA    SMC01-T       7356         1714 Tumor     SMC01-T  19.711800   34.99184
# SMC01-T_AAACGGGAGGGAAACA    SMC01-T       3752         1229 Tumor     SMC01-T   9.541578   30.46375
# SMC01-T_AAACGGGGTATAGGTA    SMC01-T      23991         3914 Tumor     SMC01-T  17.314826   12.14205
# SMC01-T_AAAGATGAGGCCGAAT    SMC01-T      15662         3319 Tumor     SMC01-T   7.278764   25.60337

system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))

# 1. 质控------------------------------------
rm(list = ls())
setwd(('/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR/16_scRNA'))
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
setwd("./01_QC")
load('../01.scRNA_orig.Rdata')
length(colnames(scRNA))  #63689 细胞数量
length(rownames(scRNA))  #25655基因数量
# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
theme.set = theme(
  axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))
# 质控前小提琴图
library(patchwork)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       # group.by = "group",
                       cols = c("#FF8BA0","#A3D8A3","#FF6F61","#56BCD3","#8E9BAE",
                                "#D09075","#FFD8B1","#FFB6C1","#5CB85C","#D9534F",
                                "#FFCB65","#66BFBF","#FFD8B1","#FFB6C1","#FF8BA0",
                                "#87CEEB","#A3D8A3","#DDA0DD","#ADD8E6","#3288bdA9",
                                "#337AB7","#BA55D3","#FF6F61","#D53E4FA9","#FF6F61",
                                "#FFCB65","#66BFBF","#3D9970","#4EB5D6","#8778A0",
                                "#FF8BA0","#A3D8A3","#FF6F61","#56BCD3","#8E9BAE",
                                "#D09075","#FFD8B1","#FFB6C1","#5CB85C","#D9534F",
                                "#337AB7","#FFD700","#FFA500","#00FFFF","#87CEFA"),
                                         # type.by=group,
                       pt.size = 0,
                       features = plot.featrures[i])+ theme.set + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin
pdf("01.vlnplot_before_qc.pdf", width = 15, height = 6,family = "Times")
violin
dev.off()
png("01.vlnplot_before_qc.png", width = 15, height = 6,res = 600,units = 'in',family = "Times")
violin
dev.off()
dev.off()
## 质控后
## 设置质控标准
# minGene = quantile(scRNA$nFeature_RNA, .05) #486
# maxGene = quantile(scRNA$nFeature_RNA, .95) #4589
# maxUMI = quantile(scRNA$nCount_RNA, .95) #21935

pctmt=10
minGene <- 200
maxGene <- 2000
maxUMI <- 10000

scRNA <- subset(scRNA, subset =
                  nCount_RNA < maxUMI &
                  nFeature_RNA > minGene &
                  nFeature_RNA < maxGene &
                  percent.mt < pctmt)

length(colnames(scRNA))  #36625 细胞数量
length(rownames(scRNA))  #25655  基因数量

# 质控后小提琴图
# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       # group.by = "group",
                       cols = c("#FF8BA0","#A3D8A3","#FF6F61","#56BCD3","#8E9BAE",
                                "#D09075","#FFD8B1","#FFB6C1","#5CB85C","#D9534F",
                                "#FFCB65","#66BFBF","#FFD8B1","#FFB6C1","#FF8BA0",
                                "#87CEEB","#A3D8A3","#DDA0DD","#ADD8E6","#3288bdA9",
                                "#337AB7","#BA55D3","#FF6F61","#D53E4FA9","#FF6F61",
                                "#FFCB65","#66BFBF","#3D9970","#4EB5D6","#8778A0",
                                "#FF8BA0","#A3D8A3","#FF6F61","#56BCD3","#8E9BAE",
                                "#D09075","#FFD8B1","#FFB6C1","#5CB85C","#D9534F",
                                "#337AB7","#FFD700","#FFA500","#00FFFF","#87CEFA"),
                       pt.size = 0,
                       features = plot.featrures[i]) + theme.set+ NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)  
violin
pdf("02.vlnplot_after_qc.pdf", width = 15, height = 6,family = "Times")
violin
dev.off()
png("02.vlnplot_after_qc.png", width = 15, height = 6,res = 600,units = 'in',family = "Times")
violin
dev.off()
dev.off()

system.time(save(scRNA, file = "scRNA_qc.Rdata"))

# 2. Seurat v3的标准整合流程-----------------------------
## NormalizeData函数的作用：消除不同细胞测序深度的影响
rm(list = ls())
setwd(('/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR/16_scRNA'))
if (! dir.exists("./02_Integrate")){
  dir.create("./02_Integrate")
}
setwd("./02_Integrate")
load('../01_QC/scRNA_qc.Rdata')

library(Seurat)
library(tidyverse)

combined <- NormalizeData(scRNA)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

scRNA <- combined
scRNA <- JoinLayers(scRNA)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA), 10)

pdf("01.variable_feature.pdf",width = 8,height = 6,family = "Times")
plot1<-VariableFeaturePlot(scRNA)
plot2<-LabelPoints(plot = plot1,points = top10,repel = T)
print(plot2)
dev.off()

png("01.variable_feature.png",width = 8,height = 6,units = 'in',res = 600, family = "Times")
plot1<-VariableFeaturePlot(scRNA)
plot2<-LabelPoints(plot = plot1,points = top10,repel = T)
print(plot2)
dev.off()

system.time(save(scRNA , file = "scRNA_intergrated.Rdata"))

# 3. PCA UMAP降维聚类-------------------
rm(list = ls())
setwd(('/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR/16_scRNA'))
if (! dir.exists("./03_PCA")){
  dir.create("./03_PCA")
}
setwd("./03_PCA")
load('../02_Integrate/scRNA_intergrated.Rdata')

## PCA降维分析
## scale标准化——很多基因的reads数多达上千，但是有很多基因reads很少（ScaleData作用：是将全部基因reads数据压缩为0-1之间的数值）
scRNA.nor.sca <- ScaleData(scRNA)
scRNA.norm.pca <- RunPCA(scRNA.nor.sca, features = VariableFeatures(object = scRNA.nor.sca), npcs = 50)
## 首先找到最佳聚类数
## JackStraw函数的功能/作用：
## 随机排列数据的子集，并计算这些“随机”基因的预测PCA分数。
## 然后将“随机”基因的PCA分数与观察到的PCA分数进行比较，以确定统计显著性。
## 最终结果是每个基因与每个主成分相关的p值
## 简而言之: JackStraw()函数可以计算出每个主成分中各基因的P值，用于判断哪些主成分更具有统计学意义
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50) # num.replicate : 要执行的复制采样数；dims : 计算重要性的PC数
## ScoreJackStraw()函数：用于量化主成分的显著性强度，富含低P值基因较多的主成分更有统计学意义
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
system.time(save(scRNA.norm.pca, file = "scRNA.norm.pca.Jack.Rdata"))

# load('./scRNA.norm.pca.Jack.Rdata')
## JackStrawPlot()函数：可视化比较每个主成分的 p 值分布和均匀分布（虚线）
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot_pca
pdf("01.pca_cluster.pdf", width = 12, height = 8,family = "Times")
plot_pca
dev.off()
png("01.pca_cluster.png", width = 12, height = 8,res = 600,units = 'in',family = "Times")
plot_pca
dev.off()

plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
pdf("02.pca_sd.pdf", width = 6, height = 5,family = "Times")
plot_elbow
dev.off()
png("02.pca_sd.png", width = 6, height = 5,res = 600,units = 'in',family = "Times")
plot_elbow
dev.off()

load('./scRNA.norm.pca.Jack.Rdata')
## FindNeighbors()函数作用：计算邻接距离
scRNA.norm.pca.clu <- FindNeighbors(scRNA.norm.pca, dims = 1:30)
## FindClusters()函数作用：以迭代方式将细胞分组在一起,resolution参数表示区分细胞群的分辨率，resolution越大，分的群也就越多。
scRNA.norm.pca.clu <- FindClusters(scRNA.norm.pca.clu, resolution = 0.1) # 0.4

##non-linear dimensional reduction (UMAP) --------
UMAP<-RunTSNE(scRNA.norm.pca.clu,dims = 1:30) # UMAP

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

pdf("03.UMAP-res.pdf",height = 5,width = 6, family = "Times")
DimPlot(UMAP, reduction = "tsne",label = T,cols = col)
dev.off()
png("03.UMAP-res.png",height = 5,width = 6,units = 'in',res = 600, family = "Times")
DimPlot(UMAP, reduction = "tsne",label = T,cols = col)
dev.off()

pdf('04.UMAP_group.pdf',w=10,h=5,family = "Times")
DimPlot(UMAP,reduction = 'tsne',split.by  = "group",cols = col)
dev.off()
png('04.UMAP_group.png',w=10,h=5,units = 'in',res = 600, family = "Times")
DimPlot(UMAP,reduction = 'tsne',split.by  = "group",cols = col)
dev.off()
system.time(save(UMAP, file = "UMAP.Rdata")) 

# 4. 细胞类型注释------
rm(list = ls())
setwd(('/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR/16_scRNA'))
if (! dir.exists("./04_SingleR")){
  dir.create("./04_SingleR")
}
setwd("./04_SingleR")
load('../03_PCA/UMAP.Rdata')

library(ggplot2)
library(SingleR)
library(Seurat)
##参考文献：PMID: 40565588
gene <- c(
  'MS4A1','CD79A','CD79B', ## T cells
  'CD3D','CD3E', ## B cells
  'IGHG1','FCRL5','MZB1', ## Plasma cells
  'KRT18','CLDN4','KRT19','EPCAM', ## Epithelial cells
  'RAMP2','PLVAP','VWF','PECAM1', ## Endothelial cells
  'C1QB','C1QA','CD14','CD68', ## Macrophages
  'FCN1','S100A8', ## Monocytes
  'PLD4','TCL1A','LILRA4','IRF4','PTCRA', ## pDC
  'KIT','CPA3','TPSAB1', ## Mast cells
  'CDH19','PLP1', ## Schwann cells
  # 'RGS5','NOTCH3', ## Pericytes
  'FBLN1','COL1A2','COL1A1' ## Fibroblasts
)

theme.set = theme(
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12,angle=45 ,hjust=1,face = "bold", family = "Times"),
  axis.text.y = element_text(size = 12,  face = "bold", family = "Times"),
  legend.text = element_text(size = 14, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

pdf(file="01.DotPlot.pdf",width=8,height=5,family='Times')
DotPlot(UMAP,features = gene)+theme.set+ 
  scale_color_gradientn(colours = c('white',"#A3D8A3","#3D9970"))
dev.off()

png(file="01.DotPlot.png",width=8,height=5,family='Times',units='in',res=600)
DotPlot(UMAP,features = gene)+theme.set+ 
  scale_color_gradientn(colours = c('white',"#A3D8A3","#3D9970"))
dev.off()

# gene <- c(
#   'MS4A1','CD79A','CD79B', ## T cells
#   'CD3D','CD3E', ## B cells
#   'IGHG1','FCRL5','MZB1', ## Plasma cells
#   'KRT18','CLDN4','KRT19','EPCAM', ## Epithelial cells
#   'RAMP2','PLVAP','VWF','PECAM1', ## Endothelial cells
#   'C1QB','C1QA','CD14','CD68', ## Macrophages
#   'FCN1','S100A8', ## Monocytes
#   'PLD4','TCL1A','LILRA4','IRF4','PTCRA', ## pDC
#   'KIT','CPA3','TPSAB1', ## Mast cells
#   'CDH19','PLP1', ## Schwann cells
#   'FBLN1','COL1A2','COL1A1' ## Fibroblasts
# )

new.cluster.ids = c(
  "0" = "B cells",  
  "1" = "B cells",  
  "2" = "Macrophages/Monocytes",  
  "3" = "T cells", 
  "4" = "Epithelial cells",   
  "5" = "Fibroblasts",   
  "6" = "Plasma cells",  
  "7" = "Endothelial cells", 
  "8" = "pDC", 
  '9' = "B cells", 
  '10' = "Schwann cells",
  '11' = "Mast cells"
  ) 

UMAP <- RenameIdents(UMAP, new.cluster.ids)
UMAP$celltype = Idents(UMAP)
table(UMAP$celltype)

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

pdf("02.UMAP-res.pdf",height = 5,width = 7, family = "Times")
DimPlot(UMAP, reduction = "tsne",label = T,cols = col)
dev.off()
png("02.UMAP-res.png",height = 5,width = 7,units = 'in',res = 600, family = "Times")
DimPlot(UMAP, reduction = "tsne",label = T,cols = col)
dev.off()

pdf('03.UMAP-group.pdf',w=10,h=5,family = "Times")
DimPlot(UMAP,reduction = 'tsne',split.by  = "group",cols = col)
dev.off()
png('03.UMAP-group.png',w=10,h=5,units = 'in',res = 600, family = "Times")
DimPlot(UMAP,reduction = 'tsne',split.by  = "group",cols = col)
dev.off()

library(scCustomize)
library(stringr)
genes_to_check=str_to_upper(unique(gene))
genes_to_check
th=theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(UMAP, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_all_markers
data<-p_all_markers$data
data<-p_all_markers$data
colnames(data)
colnames(data)<-c("AverageExpression_unscaled",
                  "Precent Expressed",
                  "Features",
                  "celltype",
                  "Average Expression")
unique(data$`Precent Expressed`)
table(data$celltype)
####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(0,8))+theme_bw()+
  coord_flip()+
  scale_fill_gradient(low = "white", high = "#3D9970")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        # legend.margin=margin(t= 0.5, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=12,angle = 40,
                                    vjust = 1, hjust=1,family = "Times"),#x轴
        axis.text.y  = element_text(color="black",size=12,family = "Times"),#y轴
        legend.text = element_text(size =10,color="black",family = "Times"),#图例
        legend.title = element_text(size =10,color="black",family = "Times"),#图例
        axis.title.y=element_text(vjust=1,size=16,family = "Times")
  )+labs(x=" ",y = "")
p
ggsave(file="04.dotplot_celltype.pdf",height= 6,width = 12)
ggsave(file="04.dotplot_celltype.png",height= 6,width = 12)

system.time(save(UMAP, file = "UMAP_singleR.Rdata"))
