rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./15_TIDE")){
  dir.create("./15_TIDE")
}
setwd("./15_TIDE")

#将文档上传到TIDE官网可得到结果文件
library(magrittr)
library(dplyr)
library(plyr)
expr <- read.csv(file = '../00_rawdata/01.fpkmlog2_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = FALSE)

riskScore <- read.csv('../05_risk_train/03.risk.csv')
riskScore<-riskScore[,c(1,9)]
colnames(riskScore) <- c("sample","group")

fpkm <-expr[,riskScore$sample]
group<- riskScore

# 数据准备
tide_dat <- fpkm
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

tide_dat <- apply(tide_dat,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
# TIDE要求的数据处理
rownmean <- apply(fpkm,1,mean)
tide_dat2 <- sweep(fpkm, 1, rownmean)
# colnames(tide_dat2) <- gsub("-","",colnames(tide_dat2))
dim(tide_dat2)
write.table(tide_dat2, file ="01.tide_dat.xls",sep = "\t",quote = F,row.names = T)

#结果文件——————————————————————————————————————————————————————————————————————————————————
riskScore <- read.csv('../05_risk_train/03.risk.csv')
risk<-riskScore[,c(1,8:9)]
colnames(risk) <- c("id","riskScore","risk")
tide_result <- read.csv("tide.csv",header = T)
colnames(tide_result)

#TIDE
tide_result2 <- subset(tide_result, select = c("Patient", "TIDE"))
tide_plot_dat <- data.frame(Patient=risk$id,
                            riskScore=risk$riskScore,
                            group=risk$risk)
plot.dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
colnames(plot.dat)
##散点图线性拟合
p1<-ggplot(plot.dat,aes(x= riskScore,y=TIDE))+geom_point(aes( color = group), size=4)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='red')+
  stat_cor(data=plot.dat, method = "spearman")+
  scale_color_manual(values = c("#FFD700", "#66C2A5"))+
  scale_fill_manual(values = c("#FFD700", "#66C2A5")) +theme_bw()+
  theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
        axis.text.x =element_text(size=10,family = "Times", face = "bold"),
        axis.title.y =element_text(size=15,family = "Times", face = "bold"),
        axis.text.y=element_text(size=10,family = "Times", face = "bold"),
        legend.title=element_text(size=10,family = "Times", face = "bold") , 
        legend.text=element_text(size=10,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="riskScore",y="TIDE")+
  theme(legend.position = "bottom")  ##修改图例位置
##箱线图显著性检验
my_comparisons = list( c("high risk", "low risk"))
library(ggsci)
p2<-ggplot(plot.dat, aes(x=group, y=TIDE,fill=group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",cex=6)+
  scale_fill_npg()+
  scale_fill_manual(values=c("#FFD700", "#66C2A5")) + 
  labs(x = "", y = "", title = "") + 
  #theme_bw() + 
  # geom_text(aes(label = B, vjust = 1.1, hjust = -0.5, angle = 45), show_guide = FALSE) + 
  theme(panel.grid =element_blank()) +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())+#去除背景
  guides(fill='none')  ##去除图例

##图片结合
library(patchwork)  
##ggarrange(p1,p2,ncol = 2,nrow =1,widths = c(2,0.3),heights = c(1,1))  ##y轴没对应
P1<-p1+p2+plot_layout(widths = c(2, 0.3))
P1
ggsave('02.TIDE.boxplot.pdf',width=6,height=6)
ggsave('02.TIDE.boxplot.png',width=6,height=6)

