rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists("08_clinical")){dir.create("08_clinical")}
setwd("08_clinical")

library(readr)
library(tidyverse)
library(rstatix)
phenotype<-read.table('../07_prog_model/01.clinical.csv',header = T,sep = ",")

train_phenotype<-phenotype
colnames(train_phenotype)
colnames(train_phenotype) <- c('sample','Gender','Age','Stage','T_stage','M_stage','OS','OS.time','riskScore','risk')
rownames(train_phenotype) <- train_phenotype$sample
train_phenotype <- train_phenotype[,-1]
write.csv(train_phenotype,'00.train_phenotype_survival_use.csv')

# 2 临床信息风险评分的差异分析----------
# row.names(risk_data) <- risk_data$id
risk_data <- train_phenotype[,c('OS','OS.time','riskScore')]
risk_data <- risk_data[row.names(train_phenotype),]
train_phenotype <- train_phenotype[,-c(6:9)]
identical(row.names(train_phenotype),row.names(risk_data))

gene<-"riskScore"
TT<-colnames(train_phenotype)

#循环
for (i in 3) {
  gene_data<-risk_data[,3] %>% as.data.frame()
  colnames(gene_data)<-"gene"
  gene_data$Sample<-rownames(risk_data)
  for(j in TT){
    phe<-train_phenotype[,j] %>% as.data.frame()
    colnames(phe)<-"phe"
    phe$Sample<-rownames(train_phenotype)
    phe<-phe %>% dplyr::filter(phe!="")
    data<-merge(gene_data,phe,by="Sample")
    stat.test <- data%>%
      wilcox_test(gene ~ phe)%>%
      adjust_pvalue(method = 'fdr')
    stat.test$p <- ifelse(stat.test$p < 0.001, "***",
                          ifelse(stat.test$p < 0.01, "**",
                                 ifelse(stat.test$p < 0.05, "*", 'ns')))
    write.csv(stat.test, paste0("Clinical_riskScore","_",j,".csv"))
    stat.test <- stat.test %>%
      add_xy_position(x='phe', dodge = 1, step.increase = 0.5)
    nu<-data$phe %>% unique() %>% as.data.frame()
    
    if (nrow(nu) < 3) {
      p1 <- ggplot(data, aes(phe, gene))+
        geom_boxplot(aes(fill = phe))+
        #geom_jitter()+
        scale_fill_manual(values =c( "#66C2A5","#FFD700"))+
        stat_pvalue_manual(stat.test,
                           label = "p",
                           tip.length = 0.01,
                           size = 4,
                           hide.ns = F) +
        labs(title=paste0("Clinical character-",j), x=j, y = paste0("The expression of riskScore"),size=20)+
        theme_bw()+
        theme(plot.title = element_text(family = "Times",hjust =0.5,colour="black",face="bold",size=15),
              axis.text.x=element_text(family = "Times",angle=0,hjust=0.5,colour="black",face="bold",size=12),
              axis.text.y=element_text(family = "Times",hjust=0.5,colour="black",face="bold",size=12),
              axis.title.x=element_text(family = "Times",size=16,face="bold"),
              axis.title.y=element_text(family = "Times",size=16,face="bold"),
              legend.text=element_text(family = "Times",face = "bold", hjust = 0.5,colour="black", size=12),
              legend.title = element_text(family = "Times",face = "bold", size = 12),
              legend.position = "right",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      p1
      ggsave(filename = paste0("Clinical_riskScore","_",j,'_bar.pdf') , p1, w=6, h=8)
      ggsave(filename = paste0("Clinical_riskScore","_",j,'_bar.png'), p1, w=6, h=8, dpi = 600)
    }
    if (nrow(nu) > 2){
      p2 <- ggplot(data, aes(phe, gene))+
        geom_boxplot(aes(fill = phe))+
        #geom_jitter()+
        scale_fill_manual(values =c("#66C2A5","#FFD700","#8FB4DC","#AC99D2"))+
        stat_pvalue_manual(stat.test,
                           label = "p",
                           tip.length = 0.01,
                           size = 4,
                           hide.ns = F) +
        labs(title=paste0("Clinical character-",j), x=j, y = paste0("The expression of riskScore"),size=20)+
        theme_bw()+
        theme(plot.title = element_text(family = "Times",hjust =0.5,colour="black",face="bold",size=15),
              axis.text.x=element_text(family = "Times",angle=0,hjust=0.5,colour="black",face="bold",size=12),
              axis.text.y=element_text(family = "Times",hjust=0.5,colour="black",face="bold",size=12),
              axis.title.x=element_text(family = "Times",size=16,face="bold"),
              axis.title.y=element_text(family = "Times",size=16,face="bold"),
              legend.text=element_text(family = "Times",face = "bold", hjust = 0.5,colour="black", size=12),
              legend.title = element_text(family = "Times",face = "bold", size = 12),
              legend.position = "right",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      p2
      ggsave(filename = paste0("Clinical_riskScore","_",j,'_bar.pdf') , p2, w=6, h=8)
      ggsave(filename = paste0("Clinical_riskScore","_",j,'_bar.png'), p2, w=6, h=8, dpi = 600)}
  }
  print(i)
}
