rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists('./05_risk_train')) {
  dir.create('./05_risk_train')
}
setwd('./05_risk_train')
library(readxl)
library(readr)
library(tidyverse)

group <- read.csv(file = '../00_rawdata/01.group_TCGA.COADREAD.csv', row.names = 1)
Tumor.sample <- group[group$Type == 'Tumor', , drop = FALSE]

gene <- read.csv('../02_venn/01.Intersection_gene.csv')
gene <- unique(gene$symbol)

dat_fpkm <- read.csv(file = '../00_rawdata/01.fpkmlog2_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = FALSE)
train_dat <- dat_fpkm[rownames(dat_fpkm) %in% gene, rownames(Tumor.sample)]
train_dat <- t(train_dat) %>% as.data.frame()

survival <- read.csv(file = '../00_rawdata/01.survival_TCGA.COADREAD.csv', row.names = 1)
train_dat <- merge(survival, train_dat, by = 0)
train_dat <- column_to_rownames(train_dat, var = 'Row.names') %>% na.omit()

train_data <- train_dat[train_dat[['OS.time']] > 0, ]

##Univariate Cox regression-----------
colnames(train_data) <- gsub('-', '_', colnames(train_data))
train_data <- train_data[, colSums(train_data != 0) > 0]
gene <- setdiff(colnames(train_data), c('OS.time', 'OS'))
library(survival)
library(survminer)
res.all <- data.frame()
for (i in 1:length(gene)) {
  gene_list <- gene[i]
  cox_data <- as.formula(paste0('Surv(OS.time, OS)~', gene_list))
  cox <- coxph(cox_data, data = train_data)
  coxSummary <- summary(cox)
  outTab <- data.frame(
    id = gene_list,
    HR = coxSummary$conf.int[,'exp(coef)'],
    HR.95L = coxSummary$conf.int[,'lower .95'],
    HR.95H = coxSummary$conf.int[,'upper .95'],
    pvalue = coxSummary$coefficients[,'Pr(>|z|)']
  )
  cox_zph <- cox.zph(cox)
  cox_table <- cox_zph$table[-nrow(cox_zph$table), ]
  outTab$PH.p <- cox_table[3]
  res.all <- rbind(res.all, outTab)
}
res.all <- res.all[order(as.numeric(res.all$HR), decreasing = TRUE), ]
write.csv(res.all, file = '01.uniCox_all.csv', row.names = FALSE, quote = FALSE)
res_results_0.05 <- res.all %>%
  dplyr::filter(as.numeric(pvalue) < 0.05,
                as.numeric(PH.p) > 0.05,
                as.numeric(HR) != 1) %>% na.omit()
write.csv(res_results_0.05,
          file = '01.uniCox_filtered.csv',
          quote = FALSE,
          row.names = FALSE)


res_results_plot <- res_results_0.05
#res_results_plot <- head(res_results_plot,10)
res_results_plot$HR <- as.numeric(res_results_plot$HR)
res_results_plot$HR.95L <- as.numeric(res_results_plot$HR.95L)
res_results_plot$HR.95H <- as.numeric(res_results_plot$HR.95H)
res_results_plot$pvalue <- as.numeric(res_results_plot$pvalue)
hz <- paste(round(res_results_plot$HR, 3),
            '(', round(res_results_plot$HR.95L, 3),
            '-', round(res_results_plot$HR.95H, 3), ')', sep = '')
tabletext <- cbind(c(NA, 'Gene', res_results_plot$id),
                   c(NA, 'P value', ifelse(res_results_plot$pvalue < 0.001,
                                         '< 0.001',
                                         round(res_results_plot$pvalue, 4))),
                   c(NA, 'Hazard Ratio(95% CI)', hz))

library(forestplot)
p <- forestplot(labeltext = tabletext,
                graph.pos = 2,
                is.summary = c(FALSE, TRUE, rep(FALSE, nrow(tabletext) - 2)),
                col = fpColors(box = '#FFD700', lines = '#66C2A5', zero = 'gray50'),
                mean = c(NA, NA, res_results_plot$HR),
                lower = c(NA, NA, res_results_plot$HR.95L),
                upper = c(NA, NA, res_results_plot$HR.95H),
                boxsize = 0.15, lwd.ci = 2,
                ci.vertices.height = 0.05, ci.vertices = TRUE,
                zero = 1, lwd.zero = 2,
                colgap = unit(5, 'mm'),
                xticks = c(0, 1, 2, 3, 4),
                lwd.xaxis = 2,
                lineheight = unit(0.7, 'cm'),
                graphwidth = unit(0.5, 'npc'),
                cex = 1.2, fn.ci_norm = fpDrawNormalCI,
                hrzl_lines = list('3' = gpar(col = 'black', lty = 1, lwd = 2)),
                txt_gp = fpTxtGp(label = gpar(cex = 1.2, fontface = 'bold', fontfamily = 'Times'),
                                 summary = gpar(cex = 1.4, fontfamily = 'Times'),
                                 ticks = gpar(cex = 1.2, fontface = 'bold', fontfamily = 'Times'),
                                 xlab = gpar(cex = 1.3, fontface = 'bold', fontfamily = 'Times'),
                                 title = gpar(cex = 1.5, fontface = 'bold', fontfamily = 'Times')),
                xlab = 'OS Hazard Ratio',
                grid = TRUE,
                clip = c(0, 4))

pdf(file = '01.univariate_cox_forest.pdf', height = 5, width = 11, onefile = FALSE)
print(p)
dev.off()

png(file = '01.univariate_cox_forest.png', height = 5, width = 11, units = 'in', res = 600, family = 'Times')
print(p)
dev.off()
dev.off()


###RSF---------------
inputgene <- res_results_0.05$id
set.seed(2)
train_scaled <- scale(train_data[,inputgene])
train_center <- attr(train_scaled, 'scaled:center')
train_scale  <- attr(train_scaled, 'scaled:scale')
saveRDS(train_center,'train_center.rds')
saveRDS(train_scale,'train_scale.rds')

train_input <- cbind(train_data[,c('OS','OS.time')],train_scaled)
library(randomForestSRC)
rsf_model <- rfsrc(Surv(OS.time, OS) ~ ., data = train_input,
                   ntree = 1000, mtry = 4)
imp <- vimp(rsf_model)$importance
imp
saveRDS(rsf_model,'rsf_model.rds')
saveRDS(inputgene,'inputgene.rds')

risk <- cbind(id = rownames(train_data),
              train_data[,c('OS','OS.time',inputgene)],
              riskScore = as.numeric(predict(rsf_model,newdata = train_input[,inputgene])$predicted))


#risk$risk=as.vector(ifelse(riskScore>median(riskScore),'1','0'))
dat_score <- risk[,c('OS','OS.time','riskScore')]
res.cut <- surv_cutpoint(dat_score, time = 'OS.time', event = 'OS', minprop = 0.1,
                         variables = 'riskScore')
cutpoint <- res.cut$cutpoint$cutpoint
cutpoint
# 29.6248
risk$risk <- ifelse(risk$riskScore > cutpoint, 1, 0)

risk2 <- risk
risk2$risk <- ifelse(risk2$risk == 1, 'high risk', 'low risk')
table(risk2$risk)
# high risk  low risk 
# 168       356 
write.csv(risk2, file = '03.risk.csv', quote = FALSE, row.names = FALSE)

library(ggplot2)
library(ggthemes)
library(Ipaper)
riskScore <- risk$riskScore
median(riskScore)

risk_dis <- ggplot(risk, aes(x = reorder(id, riskScore),
                             y = riskScore,
                             color = factor(risk,
                                            levels = c(1, 0),
                                            labels = c('High Risk', 'Low Risk')))) +
  geom_point() +
  scale_color_manual(values = c('#FFD700', '#66C2A5')) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900,1000,1100)],
                   labels = c(1,100,200,300,400,500,600,700,800,900,1000,1100),
                   expand = c(0.02, 0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == 0),]) + 0.5, lty = 2) +
  geom_hline(yintercept = cutpoint, lty = 2) +
  labs(x = 'Patients (increasing risk score)', y = 'Risk Score', title = 'TCGA-COADREAD') +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(1, 0),
        legend.background = element_rect(color = 'black', size = 0.3),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

surv_stat <- ggplot(risk, aes(x = reorder(id, riskScore),
                              y = OS.time / 365,
                              color = factor(OS,
                                             levels = c(0, 1),
                                             labels = c('Alive', 'Dead')))) +
  geom_point() +
  scale_color_manual(values = c('#FFD700', '#66C2A5')) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900,1000,1100)],
                   labels = c(1,100,200,300,400,500,600,700,800,900,1000,1100),
                   expand = c(0.02, 0)) +
  ylim(c(0, 15)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == 0),]) + 0.5, lty = 2) +
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
ggsave(filename = '03.Train_risk_and_survival_Distribution.pdf', p, width = 6, height = 7)
ggsave(filename = '03.Train_risk_and_survival_Distribution.png', p, width = 6, height = 7, units = 'in', dpi = 600)


## KM------
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = risk)

train_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = FALSE, 
                                      legend.labs = c('Low risk', 'High risk'), 
                                      legend.title = 'Risk score', 
                                      title = 'TCGA-COADREAD', 
                                      font.main = c(15, 'bold'), 
                                      risk.table = TRUE, 
                                      risk.table.col = 'strata', 
                                      linetype = 'strata', 
                                      ggtheme = theme_bw(), 
                                      palette = c('#66C2A5','#FFD700' ))

train_survival_median

pdf(file = '04.Train_KM.pdf', width = 5, height = 5, family = 'Times')
print(train_survival_median)
dev.off()

png(file = '04.Train_KM.png', width = 5, height = 5, units = 'in', res = 600, family = 'Times')
print(train_survival_median)
dev.off()
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
                           risk_score_table = risk)
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
  labs(title = 'TCGA-COADREAD') +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
        text=element_text(family='Times')) +
  annotate('text', x=0.75, y=c(0.25, 0.15, 0.05), label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)), paste('AUC of 2 years =', format(auc_y2,nsmall=2)), paste('AUC of 3 years =', format(auc_y3,nsmall=2))))
ROC
pdf('05.Train_ROC.pdf', w=5,h=4,family='Times')
print(ROC)
dev.off()
png('05.Train_ROC.png',w=5,h=4,units='in',res=600,family='Times')
print(ROC)
dev.off()

#单因素PH
library(readxl)
cox <- res_results_0.05
diff_expr_clinical <- train_data
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))
univ_models <- lapply(univ_formulas,
                      function(x) coxph(x, data = diff_expr_clinical))
if(!dir.exists("00.PH")){dir.create("00.PH")}
setwd("00.PH")
ph_results_list <- list()
for (gene in cox$id) {
  test.ph <- cox.zph(univ_models[[gene]])
  
  # 保存PH检验图
  plot_zph <- ggcoxzph(test.ph)
  ggsave(paste0(gene, "-ggcoxzph.pdf"), arrangeGrob(grobs = plot_zph), height = 5, width = 5)
  
  # 保存PH检验结果为CSV
  ph_result_df <- test.ph$table
  ph_result_df <- as.data.frame(ph_result_df)
  ph_result_df$gene <- gene # 使用循环变量gene作为基因名称
  ph_result <- ph_result_df[1,c(4,3)]
  # 将数据框写入CSV文件
  ph_results_list[[gene]] <- ph_result
}

ph_results_df <- do.call(rbind, ph_results_list)
write.csv(ph_results_df, file = "01.ph_results.csv", row.names = FALSE)







