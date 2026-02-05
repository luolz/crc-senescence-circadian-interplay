rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists("11_checkpoint")){dir.create("11_checkpoint")}
setwd("11_checkpoint")

library(tibble)
library(magrittr)
library(corrplot)

datExp <- read.csv(file = '../00_rawdata/01.fpkmlog2_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = FALSE)
riskScore <- read.csv('../05_risk_train/03.risk.csv')
risk<-riskScore[,c(1,9)]

colnames(risk) <- c('sample','risk')
risk <- risk[order(risk$risk,decreasing = T),]
datExp <- datExp[,risk$sample]

stimuGenes <- read.csv('/data/nas1/zhuxuying/common/check_point.csv',header = T)
stimuGenes <- stimuGenes$gene
# stimuGenes <- read.table("./immue_MHC.txt",header = T)$gene %>% gsub("-","_",.)#21
stimuGenes_exp <- subset(datExp,rownames(datExp) %in% stimuGenes) %>% na.omit()#
data <- stimuGenes_exp %>% t %>% as.data.frame()%>% rownames_to_column(.,var = "sample")
dat1 <- merge(risk,data,by = "sample")
library(reshape2)
dat2 <- melt(dat1,id = c("sample","risk")) 
dat2 <- dat2[order(dat2$risk,decreasing = T),]
colnames(dat2) <- c('sampele','group','gene','expr')

library(rstatix)
stat_res <- dat2 %>%
  group_by(gene) %>%
  wilcox_test(expr ~ group) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_res

write.csv(stat_res,"01.Wilcoxon_res.csv")
diff_res<-stat_res[which(stat_res$p<0.05),]
dim(diff_res)#11
write.csv(diff_res,"01.Wilcoxon_diffres.csv")
# 差异柱状 --------------------------------------------------------------------
dat3 <- dat2[dat2$gene %in% diff_res$gene, ]

# dat3$group<- factor(dat3$group)
dat3$group<- factor(dat3$group,levels = c("high risk","low risk"))
colnames(dat3)
library(RColorBrewer)
color3 <- brewer.pal(4,"Set3") 
dat3<-transform(dat3,pos=as.numeric(group),
                pos_adj=ifelse(group == "low risk",-0.22,0.22))
# library(gghalves)
library(ggplot2)
p <-  ggplot(dat3,
             aes(x=gene,y=expr,
                 fill=group,
                 outlier.shape = NA,
                 bxp.errorbar = T)) +
  geom_boxplot(width=0.8,
               alpha=1,
               position = position_dodge(0.9))+
  scale_fill_manual(values = c("#FFD700","#66C2A5")) +
  annotate(geom = "text", x = diff_res$gene, y = 9, size = 3, family = "Times",
           label =as.character(diff_res$p.signif)) +
  theme_bw(base_family = "Times")+
  theme(legend.position = "top")+
  theme(axis.title.x =element_text(size=12,family = "Times", face = "bold"),
        axis.text.x =element_text(angle=45,size=12,hjust = 1,family = "Times", face = "bold"),
        axis.title.y =element_text(size=15,family = "Times", face = "bold"),
        axis.text.y=element_text(size=10,family = "Times", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.title=element_text(size=15) , legend.text=element_text(size=15))+
  labs(x="",y="Expression")
p
ggsave(width=8,height=6,'02.Train_Expr.pdf')
ggsave(width=8,height=6,'02.Train_Expr.png')

