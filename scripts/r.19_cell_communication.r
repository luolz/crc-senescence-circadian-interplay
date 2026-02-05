rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./19_cell_communication")){
  dir.create("./19_cell_communication")
}
setwd("./19_cell_communication")

##Normal---------------------
if (! dir.exists("./Normal")){dir.create("./Normal")}
setwd("./Normal")

load('/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR/16_scRNA/04_SingleR/UMAP_singleR.Rdata')
library(Seurat)
library(dplyr)
library(CellChat)
library(ggalluvial)
scRNA <- UMAP
table(scRNA$celltype)
data.input <- scRNA@assays$RNA$counts # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$group) 

cell.use = rownames(meta)[meta$group == "Normal"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)

## CellChat的输入需要的matrix和meta已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将"celltype"设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量
groupSize

## 基于配体-受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用人

library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  # 展示互作记录
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use # 将提取出来的数据库内容加载到cellchat对象中

cellchat <- subsetData(cellchat, features = NULL) # 如果有感兴趣的基因可以填写到features这里
cellchat <- identifyOverExpressedGenes(cellchat) # 首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat) # 寻找高表达的通路
cellchat <- smoothData(cellchat, adj = PPI.human) # 将基因映射为蛋白

## 计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat) # 有个参数为raw.use = T，意味着上一步基因映射为蛋白这一步就没用了，计算胞间通讯概率时会用原始数据，如果为F，就会用映射的蛋白做计算
df.net <- subsetCommunication(cellchat) # 将细胞通讯预测结果以数据框的形式取出

write.csv(df.net, 'df_net.csv')

## 信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

# ### 细胞间相互作用--------------------------
pdf('01.Number_of_interactions.pdf',w=6,h=5,family='Times')
netVisual_heatmap(cellchat)
dev.off()
png('01.Number_of_interactions.png',w=6,h=5,family='Times',units = "in",res = 600,)
netVisual_heatmap(cellchat)
dev.off()

pdf('02.Net_circle_number.pdf',w=7,h=7,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('02.Net_circle_number.png',w=7,h=7,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('02.Net_circle_weight.pdf',w=7,h=7,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('02.Net_circle_weight.png',w=7,h=7,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

pdf("03.Net_circle_number_hubcell.pdf",width = 20,height = 20,family='Times')
mat <- cellchat@net$count

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png("03.Net_circle_number_hubcell.png",width = 20,height = 20,family='Times',units='in',res=600)
mat <- cellchat@net$count

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf("03.Net_circle_weight_hubcell.pdf",width = 20,height = 20,family='Times')
mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png("03.Net_circle_weight_hubcell.png",width = 20,height = 20,family='Times',units='in',res=600)
mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

save(cellchat, file = "cellchat.RData")

### 配体-受体的配对-----------------------------------
library(ggplot2)
levels(cellchat@idents) # 5个

pdf('04.Single_circle（B cells）.pdf', family='Times', width = 5, height = 5)
p <- netVisual_bubble(cellchat, sources.use = 1,targets.use = c(2:9), remove.isolate = FALSE)
p 
dev.off()

png('04.Single_circle（B cells）.png', family='Times', width = 5, height = 5, units = 'in', res = 600)
p <- netVisual_bubble(cellchat, sources.use = 1,targets.use = c(2:9), remove.isolate = FALSE)
p 
dev.off()

##Disease---------------------
rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./19_cell_communication")){
  dir.create("./19_cell_communication")
}
setwd("./19_cell_communication")

# timing----

##IMN---------------------
if (! dir.exists("./Tumor")){dir.create("./Tumor")}
setwd("./Tumor")

load('/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR/16_scRNA/04_SingleR/UMAP_singleR.Rdata')

scRNA <- UMAP
table(scRNA$celltype)
data.input <- scRNA@assays$RNA$counts # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$group)
cell.use = rownames(meta)[meta$group == "Tumor"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)

## CellChat的输入需要的matrix和meta已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将"celltype"设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量
groupSize

## 基于配体-受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用人

library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  # 展示互作记录
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use # 将提取出来的数据库内容加载到cellchat对象中

cellchat <- subsetData(cellchat, features = NULL) # 如果有感兴趣的基因可以填写到features这里
cellchat <- identifyOverExpressedGenes(cellchat) # 首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat) # 寻找高表达的通路
cellchat <- smoothData(cellchat, adj = PPI.human) # 将基因映射为蛋白

## 计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat) # 有个参数为raw.use = T，意味着上一步基因映射为蛋白这一步就没用了，计算胞间通讯概率时会用原始数据，如果为F，就会用映射的蛋白做计算
df.net <- subsetCommunication(cellchat) # 将细胞通讯预测结果以数据框的形式取出

write.csv(df.net, 'df.net.csv')

## 信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

pdf('01.Number_of_interactions.pdf',w=6,h=5,family='Times')
netVisual_heatmap(cellchat)
dev.off()
png('01.Number_of_interactions.png',w=6,h=5,family='Times',units = "in",res = 600,)
netVisual_heatmap(cellchat)
dev.off()

pdf('02.Net_circle_number.pdf',w=7,h=7,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('02.Net_circle_number.png',w=7,h=7,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('02.Net_circle_weight.pdf',w=7,h=7,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('02.Net_circle_weight.png',w=7,h=7,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

pdf("03.Net_circle_number_hubcell.pdf",width = 20,height = 20,family='Times')
mat <- cellchat@net$count

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png("03.Net_circle_number_hubcell.png",width = 20,height = 20,family='Times',units='in',res=600)
mat <- cellchat@net$count

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf("03.Net_circle_weight_hubcell.pdf",width = 20,height = 20,family='Times')
mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png("03.Net_circle_weight_hubcell.png",width = 20,height = 20,family='Times',units='in',res=600)
mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)        #设置画板，两行三列 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

save(cellchat, file = "cellchat.RData")
### 配体-受体的配对-----------------------------------
library(ggplot2)
levels(cellchat@idents) 

pdf('04.Single_circle（B cells）.pdf', family='Times', width = 5, height = 5)
p <- netVisual_bubble(cellchat, sources.use = 1,targets.use = c(2:10), remove.isolate = FALSE)
p
dev.off()

png('04.Single_circle（B cells）.png', family='Times', width = 5, height = 5, units = 'in', res = 600)
p <- netVisual_bubble(cellchat, sources.use = 1,targets.use = c(2:10), remove.isolate = FALSE)
p 
dev.off()
