rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists('./00_rawdata')) {
  dir.create('./00_rawdata')
}
setwd('./00_rawdata')

disease <- 'COADREAD'

###TCGA-------------
library(readr)
library(readxl)
library(tidyverse)
library(tibble)

#coad---------------不用跑这块，太慢了
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
length(query$results[[1]]$file_id)
GDCdownload(query,directory = "GDCdata")
expquery <- GDCprepare(query,directory = "GDCdata",summarizedExperiment = T)
save(expquery,file = "COAD.gdc_2025.5.rda")

##########可以从这里跑
load('/data/nas2/TCGA/COAD/COAD.gdc_2025.5.rda')
rowdata <- rowData(expquery) %>% as.data.frame()
table(rowdata$gene_type)

probe2symbol <- rowdata %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  dplyr::select(ID = gene_id, symbol = gene_name)
sample <- data.frame(sample = expquery$sample,
                     barcode = expquery$barcode)
duplicated_samples1 <- sample %>%
  group_by(sample) %>%
  filter(n() == 1)
duplicated_samples2 <- sample %>%
  group_by(sample) %>%
  filter(n() > 1)
if(nrow(duplicated_samples2)>1){
  duplicated_samples2 <- duplicated_samples2%>%
    arrange(sample,
            desc(str_split(duplicated_samples2$barcode, "-", simplify = TRUE)[,5]),
            desc(str_split(duplicated_samples2$barcode, "-", simplify = TRUE)[,6]))
  duplicated_samples2 <- duplicated_samples2[!duplicated(duplicated_samples2$sample),]
}
barcode <- c(duplicated_samples1$barcode,duplicated_samples2$barcode)


## 删除重复样本
count <- setNames(as.data.frame(assay(expquery, "unstranded")), expquery$barcode)
count <- count[,barcode]
colnames(count) <- substr(colnames(count), 1, 16)

#count <- count[,!duplicated(colnames(count))]
range(count)
dat_count <- count
dat_count$ID <- rownames(dat_count) %>% as.character()

##id转换
dat_count <- dat_count %>%
  inner_join(probe2symbol, by = 'ID') %>%
  dplyr::select(ID, symbol, everything()) %>% 
  mutate(rowMean = rowMeans(.[grep('TCGA', names(.))])) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(symbol, .keep_all = T) %>% 
  dplyr::select(-rowMean)  

probeid <- dat_count$ID
rownames(dat_count) <- dat_count$symbol
dat_count[,c('ID', 'symbol')] <- NULL
dim(dat_count) 
#[1] 19938   514


## 样本分组
phenotype <- colData(expquery) %>% as.data.frame()
group <- data.frame(sample = phenotype$sample,
                    sample_type = phenotype$sample_type)
group <- group[!duplicated(group),]
group <- group[substr(group$sample, 16, 16) == "A", ]
table(group$sample_type)
group <- group %>%
  subset(sample_type %in% c('Primary Tumor', 'Solid Tissue Normal')) %>% 
  transform(Type = ifelse(sample_type == 'Solid Tissue Normal', 'Normal', 'Tumor')) %>%
  with(data.frame(row.names = sample, Type))
table(group$Type)
#Normal  Tumor 
#41    455 
sample <- Reduce(intersect, list(rownames(group), colnames(dat_count)))
group.write.code <- group[sample,,drop=F] %>% dplyr::arrange(Type)
table(group.write.code$Type)
# Normal  Tumor 
# 41    455 
#write.csv(group.write.code, file = paste0('../00_rawdata/01.group_TCGA.',disease,'.csv'), quote=F) 

count.write.code <- dat_count[rownames(group.write.code)]
identical(colnames(count.write.code), rownames(group.write.code))
dim(count.write.code)
#19938   496
#write.csv(count.write.code,file = paste0('../00_rawdata/01.count_TCGA.',disease,'_mRNA.csv'), quote = F, row.names = T)


# ## 保留有生存数据的样本
survival <- subset(phenotype,select =c(sample,vital_status,days_to_death,days_to_last_follow_up))
survival <- survival[!duplicated(survival),]
table(survival$vital_status)
survival <- survival[survival$vital_status != 'Not Reported',]
survival <- survival %>%
  transmute(sample,
            OS = as.integer(vital_status == "Dead"),
            OS.time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)) %>%
  filter(!is.na(OS.time) & OS.time > 0)
survival <- data.frame(row.names = survival$sample,
                       OS = survival$OS,
                       OS.time = survival$OS.time)

tumor <- rownames(group.write.code[group.write.code$Type == 'Tumor', , drop = FALSE])
survival.sample <- Reduce(intersect, list(tumor, rownames(survival)))
survival.write.code <- survival[survival.sample,]
dim(survival.write.code)
#387      2
#write.csv(survival.write.code, file = paste0('../00_rawdata/01.survival_TCGA.',disease,'.csv'), quote=F) 






##fpkm
fpkm <- setNames(as.data.frame(assay(expquery, "fpkm_unstrand")), expquery$barcode)
fpkm <- fpkm[probeid,barcode]
colnames(fpkm) <- substr(colnames(fpkm), 1, 16)
#fpkm <- fpkm[,!duplicated(colnames(fpkm))]
range(fpkm)
dat_fpkm <- log2(fpkm + 1)
range(dat_fpkm)
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm <- dat_fpkm %>%
  inner_join(probe2symbol, by = 'ID')
rownames(dat_fpkm) <- dat_fpkm$symbol
dat_fpkm[,c('ID', 'symbol')] <- NULL
dim(dat_fpkm) 
## 19938   514

fpkm.write.coad <- dat_fpkm[rownames(count.write.code), colnames(count.write.code)]
identical(colnames(fpkm.write.coad), rownames(group.write.code))
identical(colnames(fpkm.write.coad), colnames(count.write.code))
identical(rownames(fpkm.write.coad), rownames(count.write.code))
dim(fpkm.write.coad)
#  19938   496
#write.csv(fpkm.write.coad, file = paste0('../00_rawdata/01.fpkmlog2_TCGA.',disease,'_mRNA.csv'), row.names = T, quote = F)

#read------------
# library(dplyr)
# library(TCGAbiolinks)
# library(SummarizedExperiment)
# 
# query <- GDCquery(
#   project = "TCGA-READ", 
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )
# length(query$results[[1]]$file_id)
# GDCdownload(query,directory = "GDCdata")
# expquery <- GDCprepare(query,directory = "GDCdata",summarizedExperiment = T)
# save(expquery,file = "READ.gdc_2025.5.rda") 

rm(list = setdiff(ls(), c("count.write.code", "fpkm.write.coad",'group.write.code','survival.write.code','disease')))
load('/data/nas2/TCGA/READ/READ.gdc_2025.5.rda')
rowdata <- rowData(expquery) %>% as.data.frame()
table(rowdata$gene_type)

probe2symbol <- rowdata %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  dplyr::select(ID = gene_id, symbol = gene_name)

sample <- data.frame(sample = expquery$sample,
                     barcode = expquery$barcode)
duplicated_samples1 <- sample %>%
  group_by(sample) %>%
  filter(n() == 1)
duplicated_samples2 <- sample %>%
  group_by(sample) %>%
  filter(n() > 1)
if(nrow(duplicated_samples2)>1){
  duplicated_samples2 <- duplicated_samples2%>%
    arrange(sample,
            desc(str_split(duplicated_samples2$barcode, "-", simplify = TRUE)[,5]),
            desc(str_split(duplicated_samples2$barcode, "-", simplify = TRUE)[,6]))
  duplicated_samples2 <- duplicated_samples2[!duplicated(duplicated_samples2$sample),]
}
barcode <- c(duplicated_samples1$barcode,duplicated_samples2$barcode)


count <- setNames(as.data.frame(assay(expquery, "unstranded")), expquery$barcode)
count <- count[,barcode]
colnames(count) <- substr(colnames(count), 1, 16)
range(count)

##id转换
dat_count <- count
dat_count$ID <- rownames(dat_count) %>% as.character()
dat_count <- dat_count %>%
  inner_join(probe2symbol, by = 'ID') %>% 
  dplyr::select(ID, symbol, everything()) %>% 
  mutate(rowMean = rowMeans(.[grep('TCGA', names(.))])) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(symbol, .keep_all = T) %>% 
  dplyr::select(-rowMean)  

probeid <- dat_count$ID
rownames(dat_count) <- dat_count$symbol
dat_count[,c('ID', 'symbol')] <- NULL
dim(dat_count) 
##19938   177


## 样本分组
phenotype <- colData(expquery) %>% as.data.frame()
group <- data.frame(sample = phenotype$sample,
                    sample_type = phenotype$sample_type)
group <- group[substr(group$sample, 16, 16) == "A", ]
table(group$sample_type)
group <- group %>%
  subset(sample_type %in% c('Primary Tumor', 'Solid Tissue Normal')) %>% 
  transform(Type = ifelse(sample_type == 'Solid Tissue Normal', 'Normal', 'Tumor')) %>%
  with(data.frame(row.names = sample, Type))
table(group$Type)
#Normal  Tumor 
#10    163
sample <- Reduce(intersect, list(rownames(group), colnames(dat_count)))
group.write.read <- group[sample,,drop=F] %>% dplyr::arrange(Type)
table(group.write.read$Type)
# Normal  Tumor 
# 10    163
#write.csv(group.write.read, file = paste0('../00_rawdata/01.group_TCGA.',disease,'.csv'), quote=F) 
#rownames(group) <- substr(rownames(group), 1, 16)
count.write.read <- dat_count[rownames(group.write.read)]
identical(colnames(count.write.read), rownames(group.write.read))
dim(count.write.read)
#19938   173
#write.csv(count.write.read,file = paste0('../00_rawdata/01.count_TCGA.',disease,'_mRNA.csv'), quote = F, row.names = T)


# ## 保留有生存数据的样本
survival <- subset(phenotype,select =c(sample,vital_status,days_to_death,days_to_last_follow_up))
table(survival$vital_status)
survival <- survival[survival$vital_status != 'Not Reported',]
survival <- survival %>%
  transmute(sample,
            OS = as.integer(vital_status == "Dead"),
            OS.time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)) %>%
  filter(!is.na(OS.time) & OS.time > 0)
survival <- data.frame(row.names = survival$sample,
                       OS = survival$OS,
                       OS.time = survival$OS.time)

tumor <- rownames(group[group$Type == 'Tumor', , drop = FALSE])
survival.sample <- Reduce(intersect, list(tumor, rownames(survival)))
survival.write.read <- survival[survival.sample,]
dim(survival.write.read)
#137         2
#write.csv(survival.write.read, file = paste0('../00_rawdata/01.survival_TCGA.',disease,'.csv'), quote=F) 






##fpkm
fpkm <- setNames(as.data.frame(assay(expquery, "fpkm_unstrand")), expquery$barcode)
fpkm <- fpkm[probeid,barcode]
colnames(fpkm) <- substr(colnames(fpkm), 1, 16)
#fpkm <- fpkm[,!duplicated(colnames(fpkm))]
range(fpkm)
dat_fpkm <- log2(fpkm + 1)
range(dat_fpkm)
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm <- dat_fpkm %>%
  inner_join(probe2symbol, by = 'ID')
rownames(dat_fpkm) <- dat_fpkm$symbol
dat_fpkm[,c('ID', 'symbol')] <- NULL
dim(dat_fpkm) 
## 19938   177

fpkm.write.read <- dat_fpkm[rownames(count.write.read), colnames(count.write.read)]
identical(colnames(fpkm.write.read), rownames(group.write.read))
identical(colnames(fpkm.write.read), colnames(count.write.read))
identical(rownames(fpkm.write.read), rownames(count.write.read))
dim(fpkm.write.read)
#  19938   173
#write.csv(fpkm.write.read, file = paste0('../00_rawdata/01.fpkmlog2_TCGA.',disease,'_mRNA.csv'), row.names = T, quote = F)

fpkm.write.all <- merge(fpkm.write.coad, fpkm.write.read, by=0) %>% column_to_rownames('Row.names')
count.write.all <- merge(count.write.code, count.write.read, by=0) %>% column_to_rownames('Row.names')
group.write.all <- rbind(group.write.code, group.write.read)
survival.write.all <- rbind(survival.write.code, survival.write.read)
identical(colnames(fpkm.write.all), rownames(group.write.all))
identical(colnames(fpkm.write.all), colnames(count.write.all))
identical(rownames(fpkm.write.all), rownames(count.write.all))


write.csv(fpkm.write.all, file = paste0('../00_rawdata/01.fpkmlog2_TCGA.',disease,'_mRNA.csv'), row.names = T, quote = F)
write.csv(survival.write.all, file = paste0('../00_rawdata/01.survival_TCGA.',disease,'.csv'), quote=F) 
write.csv(count.write.all,file = paste0('../00_rawdata/01.count_TCGA.',disease,'_mRNA.csv'), quote = F, row.names = T)
write.csv(group.write.all, file = paste0('../00_rawdata/01.group_TCGA.',disease,'.csv'), quote=F) 
###TCGA





#4 ###GSE12945---------
library(GEOquery)
library(Biobase)
gset<-getGEO('GSE12945',
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gset
gpl<-getGEO("GPL96",destdir = '/data/nas1/liky/database/GPL/')

gpl<-Table(gpl)
colnames(gpl)


colnames.sym <- intersect(colnames(gpl), c("Symbol", "Gene Symbol","gene_assignment","SYMBOL"))
colnames.sym
if (colnames.sym %in% c("Symbol", "Gene Symbol","SYMBOL")) {
  select_col <- 1
} else if ("gene_assignment" %in% colnames.sym) {
  select_col <- 2
} else {
  select_col <- NA 
}

gpl$Symbol <- gpl[[colnames.sym]]
library(tidyverse)
probe2symobl <- gpl %>%
  select(ID, Symbol) %>%
  mutate(Symbol = trimws(str_split(Symbol, "///", simplify = TRUE)[, select_col])) %>%
  mutate(Symbol = trimws(str_split(Symbol, "//", simplify = TRUE)[, select_col])) %>%
  mutate(Symbol = trimws(str_split(Symbol, "/", simplify = TRUE)[, select_col]))

probe2symobl=probe2symobl[probe2symobl$Symbol!='',]
colnames(probe2symobl) <- c('ID', 'Symbol')
probe2symobl[1:20,]

dat_expr<-expr
dat_expr$ID<-rownames(dat_expr)
dat_expr$ID<-as.character(dat_expr$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat_expr<-dat_expr %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select(Symbol,everything())%>%     
  mutate(rowMean=rowMeans(.[,-1]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(Symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
range(dat_expr, na.rm = TRUE)
range(dat_expr)

#dat_expr_log2 <- log2(dat_expr + 1)
dat_expr_log2 <- dat_expr
range(dat_expr_log2, na.rm = TRUE)

pd<-pData(gset[[1]])
colnames(pd)
table(pd$source_name_ch1)
#pd <- pd[pd$source_name_ch1 != 'Frozen tissue of non tumoral colorectal mucosa',]
survival_va <- tibble(sample = rownames(pd), OS=pd$`Death:ch1`, OS.time=pd$`OverallSurvival_months:ch1`)
survival_va <- na.omit(survival_va)
table(survival_va$OS)
#survival_va$OS <- ifelse(survival_va$OS == 'death', 1, 0)

#survival_va$OS <- as.numeric(gsub("os.event: ", "", survival_va$OS))
#survival_va$OS.time <- as.numeric(gsub("os.delay \\(months\\): ", "", survival_va$OS.time))
survival_va$OS <- as.numeric(survival_va$OS)
survival_va$OS.time <- as.numeric(survival_va$OS.time)
survival_va <- survival_va %>% na.omit()

survival_va$OS.time <- as.integer(survival_va$OS.time*30)
survival_va <- survival_va[survival_va$OS.time>0,]
survival_va <- column_to_rownames(survival_va, var = "sample")

common_samples <- intersect(colnames(dat_expr_log2), rownames(survival_va))
dat_expr_log2 <- dat_expr_log2[, common_samples]
survival_va <- survival_va[common_samples, ]


write.csv(dat_expr_log2, file = '02.exprlog2_GSE12945.csv', quote = F, row.names = T)
write.csv(survival_va, file = '02.survival_GSE12945.csv', quote = F, row.names = T)
###GSE12945





