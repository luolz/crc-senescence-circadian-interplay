rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists('./13_maftools')){
  dir.create('./13_maftools')
}
setwd('./13_maftools/')

library(TCGAmutations)
library(tidyverse)

maf1 <- TCGAmutations::tcga_load(study = "COAD")
maf2 <- TCGAmutations::tcga_load(study = "READ")

all.maf <- maftools::merge_mafs(list(maf1, maf2))
maf <- all.maf
## 计算tmb值
x = tmb(maf = maf)

riskScore <- read.csv('../05_risk_train/03.risk.csv')
risk<-riskScore[,c(1,9)]
colnames(risk) <- c("sample","risk")
High.sample <- risk$sample[which(risk$risk == "high risk")]
Low.sample <- risk$sample[which(risk$risk == "low risk")]

maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
maf_sample<- merge(maf_sample, risk, by = 'sample')

sample <- subset(maf_sample)$barcode
maf <- subsetMaf(maf, tsb = sample)

# 自定义变异类型的颜色
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Set3')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample, ]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
pdf(file = "01.oncoplot.high.pdf",family = "Times", height = 6, width = 8)
oncoplot(maf = maf.high,colors = vc_cols, top = 20)
dev.off()

png(file = "01.oncoplot.high.png", family = "Times", height = 6, width = 8, units = "in", res = 600)
oncoplot(maf = maf.high,colors = vc_cols,top = 20)
dev.off()

low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
pdf(file = "02.oncoplot.low.pdf",family = "Times", height = 6, width = 8)
oncoplot(maf = maf.low,colors = vc_cols, top = 20)
dev.off()

png(file = "02.oncoplot.low.png", family = "Times", height = 6, width = 8, units = "in", res = 600)
oncoplot(maf = maf.low,colors = vc_cols, top = 20)
dev.off()

