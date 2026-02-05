rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists('./06_risk_test_GSE12945')) {
  dir.create('./06_risk_test_GSE12945')
}
setwd('./06_risk_test_GSE12945')
library(readxl)
library(readr)
library(tidyverse)
library(glmnet)

survival <- read.csv(file = '../00_rawdata/02.survival_GSE12945.csv', row.names = 1)
test_rowdata <- read.csv(file = '../00_rawdata/02.exprlog2_GSE12945.csv', row.names = 1, check.names = F)
range(test_rowdata)
#test_rowdata <- log2(test_rowdata + 1)
test_dat <- t(test_rowdata) %>% as.data.frame()

inputgene <- readRDS('../05_risk_train/inputgene.rds')
train_center <- readRDS('../05_risk_train/train_center.rds')
train_scale <- readRDS('../05_risk_train/train_scale.rds')
test_scaled <- test_dat[,inputgene] %>% na.omit()
test_scaled <- scale(test_scaled, center = train_center, scale = train_scale)
test_scaled <- scale(test_scaled)
test_input <- merge(survival,test_scaled,by=0) %>% column_to_rownames(var = 'Row.names')


#predict----
rsf_model <- readRDS('../05_risk_train/rsf_model.rds')
risk_out <- cbind(id=rownames(test_input),
                  test_input,
                  riskScore_out = as.numeric(predict(rsf_model,newdata = test_input[,inputgene])$predicted))
risk_out <- risk_out[risk_out$OS.time>0,]

#risk_out$risk_out <- ifelse(risk_out$riskScore_out>median(risk_out$riskScore_out),1,0)
dat_score <- risk_out[,c('OS','OS.time','riskScore_out')]
res.cut <- surv_cutpoint(dat_score, time = 'OS.time', event = 'OS', minprop = 0.1,
                         variables = 'riskScore_out')
cutpoint <- res.cut$cutpoint$cutpoint
cutpoint
# 32.40562
risk_out$risk_out <- ifelse(risk_out$riskScore_out > cutpoint, 1, 0)

risk2 <- risk_out
risk2$risk_out <- ifelse(risk2$risk_out == 1, 'high risk', 'low risk')
table(risk2$risk_out)
# high risk  low risk 
# 17        45 
write.csv(risk2, file = '01.risk_out.csv', quote = FALSE, row.names = FALSE)

library(ggplot2)
library(ggthemes)
library(Ipaper)
median(risk_out$riskScore_out)
riskScore <- risk_out$riskScore

risk_dis <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                             y=riskScore_out,
                             color = factor(risk_out,
                                            levels = c(1, 0),
                                            labels = c('High Risk', 'Low Risk')))) +
  geom_point() +
  scale_color_manual(values = c('#FFD700', '#66C2A5')) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,100,200,300,400,500,600,700,800,900,1000,1100)],
                   labels = c(1,100,200,300,400,500,600,700,800,900,1000,1100),
                   expand = c(0.02, 0)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==0),]) + 0.5, lty = 2) +
  geom_hline(yintercept = cutpoint, lty = 2) +
  labs(x = 'Patients (increasing risk score)', y = 'Risk Score', title = 'GSE12945') +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(1, 0),
        legend.background = element_rect(color = 'black', size = 0.3),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

surv_stat <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                              y = OS.time / 365,
                              color = factor(OS,
                                             levels = c(0, 1),
                                             labels = c('Alive', 'Dead')))) +
  geom_point() +
  scale_color_manual(values = c('#FFD700', '#66C2A5')) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,100,200,300,400,500,600,700,800,900,1000,1100)],
                   labels = c(1,100,200,300,400,500,600,700,800,900,1000,1100),
                   expand = c(0.02, 0)) +
  ylim(c(0, 10)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==0),]) + 0.5, lty = 2) +
  labs(x = 'Patients (increasing risk score)', y = 'Survival time (Years)', title = '') +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0, 1),
        legend.background = element_rect(color = 'black', size = 0.3),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)
ggsave(filename = '02.Test_risk_and_survival_Distribution.pdf', p, width = 6, height = 7)
ggsave(filename = '02.Test_risk_and_survival_Distribution.png', p, width = 6, height = 7, units = 'in', dpi = 600)


## KM------
kmfit <- survfit(Surv(OS.time, OS) ~ risk_out, data = risk_out)

train_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = FALSE, 
                                      legend.labs = c('Low risk', 'High risk'), 
                                      legend.title = 'Risk score', 
                                      title = 'GSE12945', 
                                      font.main = c(15, 'bold'), 
                                      risk.table = TRUE, 
                                      risk.table.col = 'strata', 
                                      linetype = 'strata', 
                                      ggtheme = theme_bw(), 
                                      palette = c('#66C2A5','#FFD700'))

train_survival_median

pdf(file = '03.Test_KM.pdf', width = 5, height = 5, family = 'Times')
print(train_survival_median)
dev.off()

png(file = '03.Test_KM.png', width = 5, height = 5, units = 'in', res = 600, family = 'Times')
print(train_survival_median)
dev.off()


## ROC------
multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$OS.time,
                        status = risk_score_table$OS,
                        marker = risk_score_table$riskScore,
                        predict.time = single_time,
                        method = 'KM')
    data.frame('True_positive' = for_ROC$TP, 'False_positive' = for_ROC$FP,
               'Cut_values' = for_ROC$cut.values, 'Time' = rep(single_time, length(for_ROC$TP)),
               'AUC' = rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365,730,1095),
                           risk_score_table = risk_out)
table(for_multi_ROC$AUC)
a <- unique(data.frame(Time = for_multi_ROC$Time, AUC = for_multi_ROC$AUC))
a
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365), 5][1], 2)
auc_y2 <- round(for_multi_ROC[which(for_multi_ROC$Time==730), 5][1], 2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095), 5][1], 2)
ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                  y=True_positive,
                                  label=Cut_values,
                                  color=Time)) +
  scale_color_manual(breaks = c('365','730','1095'),
                      labels = c('1 years', '2 years', '3 years'),
                      values = c('#66C2A5','#FFD700',"#8DA0CB")) +
  geom_roc(labels = FALSE, stat = 'identity') +
  style_roc() +
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = 'GSE12945') +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
        text=element_text(family='Times')) +
  annotate('text', x=0.75, y=c(0.25, 0.15, 0.05), label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)), paste('AUC of 2 years =', format(auc_y2,nsmall=2)), paste('AUC of 3 years =', format(auc_y3,nsmall=2))))
ROC
pdf('04.Test_ROC.pdf', w=5,h=4,family='Times')
print(ROC)
dev.off()
png('04.Test_ROC.png',w=5,h=4,units='in',res=600,family='Times')
print(ROC)
dev.off()
dev.off()
