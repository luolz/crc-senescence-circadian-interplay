rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists("10_ssGSEA")) {dir.create("10_ssGSEA")}
setwd("10_ssGSEA")

library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(immunedeconv)
library(GSVA)
library(DESeq2)
expr <- read.csv(file = '../00_rawdata/01.fpkmlog2_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = FALSE)

riskScore <- read.csv('../05_risk_train/03.risk.csv')
risk<-riskScore[,c(1,9)]
colnames(risk) <- c("sample","group")

fpkm <-expr[,risk$sample]
train_group<- risk

gene_set <- read.table("/data/nas1/zhuxuying/common/mmc3.txt",header = T,sep ="\t")
# expr[,1:802] <- as.data.frame(lapply(expr[,1:802],as.numeric))
dat.final <- as.matrix(fpkm)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.csv(ssgsea_score,
          file = "01.ssgsea_result_cell.csv",
          quote = F)
# 热图----
colnames(train_group) <- c('sample', 'group')
train_group<-train_group[order(train_group$group),]

ssgsea_score<-read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1)
ssgsea_score<-ssgsea_score[,train_group$sample]
annotation_col<-as.data.frame(train_group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)

library(pheatmap)
color.key<-c("#66C2A5", "white","#FFD700")
ann_colors<-list(Group = c('low risk'="#66C2A5",'high risk'="#FFD700"))
p <- pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  # border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
p
png(filename = "01.heatmap.png", height = 5, width = 8, units = 'in',res = 600,family='Times')
p
dev.off()
pdf(file = "01.heatmap.pdf", height = 5,width = 8,family='Times')
p
dev.off()
# 差异免疫细胞鉴定----
colnames(train_group) <- c('sample', 'group')
group_case<-train_group[train_group$group=='high risk',]
group_case<-as.character(group_case$sample)
group_control<-train_group[train_group$group=='low risk',]
group_control<-as.character(group_control$sample)

tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1) 
tiics_result <- tiics_result[, train_group$sample] %>% as.matrix()
#pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
pvalue = padj <- matrix(0, nrow(tiics_result), 1)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_case],
                                       tiics_result[i, group_control])$p.value
  # log2FoldChange[i, 1] = mean(tiics_result[i, group_case]) - 
  #   mean(tiics_result[i, group_control])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(#log2FoldChange, 
  pvalue, 
  padj,
  row.names = rownames(tiics_result))
# control <- signif(apply(tiics_result[rownames(rTable), group_control], #signif 4位小数，mean平均值
#                     1,
#                     mean), 4)
# case <- signif(apply(tiics_result[rownames(rTable), group_case], 
#                      1, 
#                      mean), 4)
# rTable <- data.frame(control, 
#                      case,
#                      rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$pvalue<0.05),]
dim(diff_Table)
##7  4
write.csv(rTable,
          file = "02.cell_all_wilcox_test.csv",
          quote = F,
          row.names = F)
write.csv(diff_Table,
          file = "02.cell_diff_wilcox_test.csv",
          quote = F,
          row.names = F)
# 箱线图----
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

cell.data <- data.frame(Immune_Cell=rownames(tiics_result), 
                        tiics_result, 
                        pvalue=rTable$pvalue)
plot.cell <- cell.data[which(cell.data$pvalue<0.05),]
diff_tiics <- rownames(plot.cell)
violin_dat <- gather(plot.cell, key=Group, value=score, -c("Immune_Cell","pvalue"))
group_control<-gsub("-",".",group_control)
violin_dat$Group <- ifelse(violin_dat$Group%in% group_control,
                           "Low", "High") 
violin_dat$Group <- factor(violin_dat$Group, levels = c("Low", "High"))
violin_dat <- violin_dat[,-2]
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  # geom_violin(trim=T,color=alpha('black',alpha = 0.5)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  # stat_boxplot(geom="errorbar",
  #              width=0.1,
  #              position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = 21,
               outlier.fill = "black",
               outlier.size = 0.5,
  )+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  # geom_point(aes(fill = Group),
  #            size = 0.05,
  #            position = position_dodge(0.9))+
  scale_fill_manual(values= c("#66C2A5","#FFD700"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     method = "wilcox.test", #没有直接用差异分析的结果，但检验方法是一样的
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=60,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        text=element_text(family = 'Times'),
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs
ggsave('03.ssgsea_Box.pdf',boxplot_diff_TIICs,w=6,h=5)
ggsave('03.ssgsea_Box.png',boxplot_diff_TIICs,w=6,h=5)

# 细胞相关性----
tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1)
#tiics_result <- tiics_result[, group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
#tiics_result <- tiics_result[rownames(expr),]

diff <- read.csv('./02.cell_diff_wilcox_test.csv')
tiics_result <- tiics_result[,diff$immune_cell]

library(ggcorrplot)
cor_res <- round(cor(tiics_result, method = 'spearman'), 3)

p.mat <- ggcorrplot::cor_pmat(tiics_result)

write.csv(cor_res, '04.cor_res.csv')
write.csv(p.mat, '04.cor_res-p.csv')
library(corrplot)
library(RColorBrewer)
library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(ggExtra)
col1 = colorRampPalette(colors =c("#66C2A5","gray","#FFD700"),space="Lab") # space参数选择使用RGB或者CIE Lab颜色空间

pdf(file = "04.cor_heatmap.pdf",width = 6,height = 6,family="Times")
corrplot(corr = cor_res, 
         p.mat = p.mat,
         method = "pie",
         type = "upper", 
         tl.pos = "lt", tl.cex = 1, tl.col = "black", tl.srt = 60, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1, 
         col = col1(20)) # 用*作为显著性标签，sig.level参数表示0.05用*表示。0.01用**。0.001用***。pch.cex = 0.8,用于设置显著性标签的字符大小。
corrplot(corr = cor_res, 
         type="lower", 
         add=TRUE,
         method="number", 
         tl.pos = "n", 
         cl.pos = "n",
         diag=FALSE,
         number.cex = 1,
         col = col1(20))
dev.off()

png(file = "04.cor_heatmap.png",width =6,height =6,family="Times",units = 'in', res = 600)
corrplot(corr = cor_res, 
         p.mat = p.mat,
         method = "pie",
         type = "upper", 
         tl.pos = "lt", tl.cex = 1, tl.col = "black", tl.srt = 60, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1, 
         col = col1(20)) # 用*作为显著性标签，sig.level参数表示0.05用*表示。0.01用**。0.001用***。pch.cex = 0.8,用于设置显著性标签的字符大小。
corrplot(corr = cor_res, 
         type="lower", 
         add=TRUE,
         method="number", 
         tl.pos = "n", 
         cl.pos = "n",
         diag=FALSE,
         number.cex = 1,
         col = col1(20))
dev.off()

#  生物标志物与关键免疫细胞相关性----
expr1<-fpkm
model_gene <- readRDS('../05_risk_train/inputgene.rds')
expr1 <- expr1[model_gene,] %>% t %>% as.data.frame()

tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1) 
tiics_result <- tiics_result[, train_group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[rownames(expr1),]

diff <- read.csv('./02.cell_diff_wilcox_test.csv')
tiics_result <- tiics_result[,diff$immune_cell]

identical(rownames(expr1), rownames(tiics_result))

cor_r <- cor(expr1,tiics_result,method = "spearman") 
cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(expr1)))

cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"05.correlation_cor.csv")

#棒棒糖图--------------
cor_dat.p0.05 <- cor_dat[cor_dat$Pvalue<0.05&abs(cor_dat$Correlation)>0.3,]
correlation.dat <- cor_dat.p0.05[order(cor_dat.p0.05$Correlation),]

for (i in 1:4) {
  gene.name <- model_gene[i]
  correlation <- correlation.dat[correlation.dat$gene==gene.name,]
  p0<-ggdotchart(correlation, x = "cell", y = "Correlation",
                 dot.size ='Correlation',
                 # filled.contour(p.value,col=Lab.palette(20)),
                 color ='Pvalue',
                 # label='sig',
                 #  font.label = list( size = 9),
                 #ylab="mgp = c(3, 2, 0)",
                 # font.label = list(size=10,face="bold",color='black',position="dodge",
                 #                   vjust=0.5),
                 #  y.label=0.7,
                 # palette = palette,
                 # 按照cyl填充颜色
                 # palette = palette, # 修改颜色
                 sorting = "descending",
                 add = "segments",                             # 添加棒子
                 rotate = TRUE,
                 ggtheme = theme_pubr(),                        # 改变主题
                 ylab="Correlation",
                 xlab='',
                 # dot.size = 6 ,
                 title=gene.name
  )+scale_colour_gradient( high = "#FFD700",low = "#66C2A5")
  p10<-p0+theme(legend.position = "right",
                panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
    theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
          axis.text.x =element_text(size=12,family = "Times", face = "bold"),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=12,family = "Times", face = "bold"),
          plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
          legend.text = element_text(size = 16, family = "Times"),
          legend.title = element_text(size = 18, family = "Times",face = "bold"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  p10
  ggsave(file=paste0('05.correlation_',gene.name,'.png'), height = 6, width = 8, p10)
  ggsave(file=paste0('05.correlation_',gene.name,'.pdf'), height = 6, width = 8, p10) 
}
