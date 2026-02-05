rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./21_irGSEA")){
  dir.create("./21_irGSEA")
}
setwd("./21_irGSEA")

load('../16_scRNA/04_SingleR/UMAP_singleR.Rdata')
# install packages from Bioconductor
# install packages from CRAN

library(irGSEA)
library(Seurat)
# library(SeuratData)
library(RcppML)
library(AUCell)
library(doMC)
library(scCustomize)

sc_dataset <- UMAP
# Check
UMAP_celltype <- DimPlot(sc_dataset, reduction ="tsne",
                         group.by="celltype",label = T);UMAP_celltype
Idents(sc_dataset) <- sc_dataset$celltype
scCustomize::DimPlot_scCustom(sc_dataset, figure_plot = TRUE)


sc_dataset <- SeuratObject::UpdateSeuratObject(sc_dataset)
sc_dataset2 <- CreateSeuratObject(counts = CreateAssay5Object(GetAssayData(sc_dataset,
                                                                           assay = "RNA", 
                                                                           slot="counts")),
                                  meta.data = sc_dataset[[]])
sc_dataset2 <- NormalizeData(sc_dataset2)
sc_dataset2 <- irGSEA.score(object = sc_dataset2, assay = "RNA",
                            slot = "data", seeds = 123, 
                            #ncores = 1,
                            min.cells = 3, min.feature = 0,
                            custom = F, geneset = NULL, msigdb = T,
                            species = "Homo sapiens", 
                            category = "H",  
                            subcategory = NULL, 
                            geneid = "symbol",
                            method = c("AUCell","UCell","singscore",
                                       "ssgsea", "JASMINE", "viper"),
                            aucell.MaxRank = NULL, 
                            ucell.MaxRank = NULL,
                            kcdf = 'Gaussian')


result.dge <- irGSEA.integrate(object = sc_dataset2,
                               group.by = "celltype",
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))

# 查看RRA识别的在多种打分方法中都普遍认可的差异基因集
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)

pathway <- as.data.frame(geneset.show)
write.csv(pathway,"01.pathway.csv")

#####热图
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                               method = "RRA", 
                               top = 10,
                               show.geneset = NULL)

heatmap.plot


#####气泡图      
bubble.plot <- irGSEA.bubble(object = result.dge,
                             method = "RRA",
                             top = 10,
                             show.geneset = geneset.show)
bubble.plot

pdf("02.irGSEA.pdf",w = 12,h = 10,family = "Times")

print(bubble.plot)
dev.off()

png("02.irGSEA.png",w = 12,h = 10,units = "in",res = 600,family = "Times")

print(bubble.plot)
dev.off()





