rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists("07_prog_model")){dir.create("07_prog_model")}
setwd("07_prog_model")

#phe
load('/data/nas2/TCGA/COAD/COAD.gdc_2025.5.rda')
phenotype1 <- colData(expquery) %>% as.data.frame()
phenotype1 <- data.frame(sample = phenotype1$sample,
                         gender = phenotype1$gender,
                         age = phenotype1$age_at_index,
                         stage = phenotype1$ajcc_pathologic_stage,
                         T_stage = phenotype1$ajcc_pathologic_t,
                         # N_stage = phenotype1$ajcc_pathologic_n,
                         M_stage = phenotype1$ajcc_pathologic_m)
phenotype1 <- phenotype1[substr(phenotype1$sample, 16, 16) == "A", ]
colnames(phenotype1)
phenotype1 <- na.omit(phenotype1)

load('/data/nas2/TCGA/READ/READ.gdc_2025.5.rda')
phenotype2 <- colData(expquery) %>% as.data.frame()
phenotype2 <- data.frame(sample = phenotype2$sample,
                         gender = phenotype2$gender,
                         age = phenotype2$age_at_index,
                         stage = phenotype2$ajcc_pathologic_stage,
                         T_stage = phenotype2$ajcc_pathologic_t,
                         # N_stage = phenotype2$ajcc_pathologic_n,
                         M_stage = phenotype2$ajcc_pathologic_m)
phenotype2 <- phenotype2[substr(phenotype2$sample, 16, 16) == "A", ]
colnames(phenotype2)
phenotype2 <- na.omit(phenotype2)

colnames(phenotype1) == colnames(phenotype2)

cli_dat <- rbind(phenotype1,phenotype2)
cli_dat <- as.data.frame(cli_dat)
cli_dat[cli_dat == ""] <- NA
cli_dat <- na.omit(cli_dat)
cli_dat <- distinct(cli_dat)
rownames(cli_dat)<-cli_dat$sample

riskScore <- read.csv('../05_risk_train/03.risk.csv')
risk<-riskScore[,c(1:3,8:9)]
colnames(risk) <- gsub("id","sample",colnames(risk))

clinical <- merge(cli_dat,risk,by = "sample")
rownames(clinical)<-clinical$sample

table(clinical$stage)
table(clinical$T_stage)
# table(clinical$N_stage)
table(clinical$gender)
table(clinical$M_stage)
#标准化stage
clinical$stage <- gsub("not reported", NA, clinical$stage)
clinical$stage <- gsub("Stage I[A-D]*$", "Stage I", clinical$stage)
clinical$stage <- gsub("Stage II[A-D]*$", "Stage II", clinical$stage)
clinical$stage <- gsub("Stage III[A-D]*$", "StageIII/IV", clinical$stage)
clinical$stage <- gsub("Stage IV[A-D]*$", "StageIII/IV", clinical$stage)####Stage IV在多因素中na影响构建列线图，在此去除
table(clinical$stage)
#标准化T
clinical$T_stage <- gsub("T1[a-d]*$", "T1", clinical$T_stage)
clinical$T_stage <- gsub("T2[a-d]*$", "T2", clinical$T_stage)
clinical$T_stage <- gsub("T3[a-d]*$", "T3", clinical$T_stage)
clinical$T_stage <- gsub("T4[a-d]*$", "T4", clinical$T_stage)
clinical$T_stage <- gsub("Tis", NA, clinical$T_stage)
table(clinical$T_stage)
# #标准化M
clinical$M_stage <- gsub("M1[a-d]*$", "M1", clinical$M_stage)
clinical$M_stage <- gsub("MX", NA, clinical$M_stage)
table(clinical$M_stage)
# #标准化N
# clinical$N_stage <- gsub("N0.*$", "N0", clinical$N_stage)
# clinical$N_stage <- gsub("N1[a-d]*$", "N1/2", clinical$N_stage)
# clinical$N_stage <- gsub("N2(a|b)?$", "N1/2", clinical$N_stage)
# clinical$N_stage <- gsub("NX", NA, clinical$N_stage)
# table(clinical$N_stage)
clinical$age <- ifelse(clinical$age >= 60, "Age >= 60", "Age < 60")

clinical<-na.omit(clinical)#405
write.csv(clinical, file = '01.clinical.csv', quote = FALSE, row.names = FALSE)

## 1 单因素Cox----------
train_risk_clinical <- clinical[,-10]
dim(train_risk_clinical)
colnames(train_risk_clinical)
colnames(train_risk_clinical) <-c("sample","Gender","Age","Stage",
                                  "T_stage","M_stage","OS","OS.time","riskScore")


train_risk_clinical$Age <- factor(train_risk_clinical$Age)
train_risk_clinical$Gender <- factor(train_risk_clinical$Gender)
train_risk_clinical$Stage <- factor(train_risk_clinical$Stage)
train_risk_clinical$T_stage <- factor(train_risk_clinical$T_stage)
train_risk_clinical$M_stage <- factor(train_risk_clinical$M_stage)

library(survival)
res.risk <- coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk <- c(res.risk$conf.int[-2], res.risk$coefficients[5])

res.Age <- coxph(Surv(time = OS.time, event = OS) ~ Age, data = train_risk_clinical) %>% summary
res.Age <- c(res.Age$conf.int[-2], res.Age$coefficients[5])

res.Gender <- coxph(Surv(time = OS.time, event = OS) ~ Gender, data = train_risk_clinical) %>% summary
res.Gender <- c(res.Gender$conf.int[-2], res.Gender$coefficients[5])

res.Stage <- coxph(Surv(time = OS.time, event = OS) ~ Stage, data = train_risk_clinical) %>% summary
res.Stage <- cbind(res.Stage$conf.int[,-2], res.Stage$coefficients[,5])

res.T_stage <- coxph(Surv(time = OS.time, event = OS) ~ T_stage, data = train_risk_clinical) %>% summary
res.T_stage <- cbind(res.T_stage$conf.int[,-2], res.T_stage$coefficients[,5])

res.M_stage <- coxph(Surv(time = OS.time, event = OS) ~ M_stage, data = train_risk_clinical) %>% summary
res.M_stage <- c(res.M_stage$conf.int[-2], res.M_stage$coefficients[5])

### ph假定检验------------
cox_results <- list()
ph_test_results <- list()
colnames_sum <- colnames(train_risk_clinical)
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
uniSigExp <- train_risk_clinical[,c("OS","OS.time",covariates)]

for (i in colnames(uniSigExp[,4:(ncol(uniSigExp))])) {
  cox <- coxph(Surv(OS.time, OS) ~ uniSigExp[,i], data = uniSigExp)
  cox_results[[i]] <- cox
  ph_test <- cox.zph(cox)
  ph_test_results[[i]] <- ph_test
}
ph_test_results

ph_test_data <- data.frame(
  Variable = character(),
  Time = numeric(),
  Schoenfeld_Residuals = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in names(ph_test_results)) {
  temp_data <- data.frame(
    Variable = i,
    Time = ph_test_results[[i]]$time,
    Schoenfeld_Residuals = ph_test_results[[i]]$y[, 1],  
    p_value = ph_test_results[[i]]$table[1, 3]
  )
  ph_test_data <- rbind(ph_test_data, temp_data)
}

p_value_summary <- unique(ph_test_data[, c("Variable", "p_value")])

pdf("01.PH_Schoenfeld.pdf",width = 4,height = 4)
ggplot(p_value_summary, aes(x = reorder(Variable, -p_value), y = p_value)) +
  geom_segment(aes(xend =Variable, yend = 0), color = "#66C2A5") +  # 绘制线段
  geom_point(size = 3, color = "#FFD700") +  # 绘制点
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +  # 添加阈值线
  annotate(
    "text",
    x = Inf,  # 将标注放在图的最右侧
    y = 0.05,  # 标注的位置与阈值线对齐
    label = "p > 0.05",  # 标注文本
    vjust = -1.2,  # 调整垂直位置
    hjust = 1,  # 调整水平位置
    color = "black",  # 文本颜色
    size = 3  # 文本大小
  ) +
  labs(
    title = "Proportional Hazards Assumption Test",
    x = "phenotype",
    y = "p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,color = "black"))  # 调整 x 轴标签角度
dev.off()

png("01.PH_Schoenfeld.png",width = 1600,height = 1200,res=300)
ggplot(p_value_summary, aes(x = reorder(Variable, -p_value), y = p_value)) +
  geom_segment(aes(xend =Variable, yend = 0), color = "#66C2A5") +  # 绘制线段
  geom_point(size = 3, color = "#FFD700") +  # 绘制点
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +  # 添加阈值线
  annotate(
    "text",
    x = Inf,  # 将标注放在图的最右侧
    y = 0.05,  # 标注的位置与阈值线对齐
    label = "p > 0.05",  # 标注文本
    vjust = -1.2,  # 调整垂直位置
    hjust = 1,  # 调整水平位置
    color = "black",  # 文本颜色
    size = 3  # 文本大小
  ) +
  labs(
    title = "Proportional Hazards Assumption Test",
    x = "phenotype",
    y = "p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,color = "black"))  # 调整 x 轴标签角度
dev.off()

p_value_summary_filter <- p_value_summary[p_value_summary$p_value > 0.05,]
table(p_value_summary_filter$Variable)

### cox绘图------------
# res.ref <- c(1,1,1,NA)
res <- rbind(res.risk, res.Age, res.Gender, res.Stage, res.T_stage, res.M_stage) %>% as.data.frame()
rownames(res)

# 添加列名
colnames(res) <- c("coefficients", "conf.int.lower","conf.int.upper",  "p.value")

# 添加行名
table(train_risk_clinical$T_stage)
rownames(res) <- c("riskScore","Age(>=60 VS <60)","Gender(male VS female)", 
                   "Stage(II VS I)", "Stage(III/IV VS I)",
                   "T_stage(T2 VS T1)", "T_stage(T3 VS T1)","T_stage(T4 VS T1)",
                   "M_stage(M1 VS M0)")

# 将结果转换为数据框（如果还不是）
res <- as.data.frame(res)
# res <- res[,c(1:4)] #突然多了几列，删除一下
# 添加指标列
res$Indicators <- rownames(res)

colnames(res) <- c("hr","low","up","pv","Indicator")
res$p <- signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] <- NA
res$Indicator <- factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.csv(res2, file = "01.Univariate_cox.csv")
res2 <- subset(res2, select = -c(Indicator))

library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.05, "< 0.05", round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))

# nrow(tabletext) + 1
library(forestplot)
{p <- forestplot(labeltext = tabletext, 
                 graph.pos=4,  #为Pvalue箱线图所在的位置
                 is.summary = c(TRUE,FALSE,FALSE,FALSE,
                                FALSE,FALSE,FALSE,
                                FALSE,FALSE,FALSE,
                                FALSE,FALSE,
                                rep(FALSE, 10)), #该处用于调整分割线（加粗的那个）
                 col=fpColors(box="#FFD700", lines="#66C2A5", zero = "gray50"),
                 mean=c(NA,res2$HR),
                 lower=c(NA,res2$HR.95L), #95%置信区间下限
                 upper=c(NA,res2$HR.95H), #95%置信区间上限
                 boxsize=0.15,lwd.ci=3,   #箱子大小，线的宽度
                 ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                 zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
                 colgap=unit(3,"mm"),    #列间隙
                 xticks = c(0,2,4,6,8,10), #横坐标刻度
                 lwd.xaxis=2,            #X轴线宽
                 lineheight = unit(1.2,"cm"), #固定行高
                 graphwidth = unit(.6,"npc"), #图在表中的宽度比例
                 cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                 hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
                 # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                 #                # "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                 #                 "14" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
                 # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
                 #fpTxtGp函数中的cex参数设置各个组件的大小
                 txt_gp=fpTxtGp(label=gpar(cex=1),
                                ticks=gpar(cex=0.8, fontface = "bold"),
                                xlab=gpar(cex = 1, fontface = "bold"),
                                title=gpar(cex = 1.25, fontface = "bold")),
                 xlab="Hazard Ratio",
                 grid = T,
                 # title = "Univariate",
                 clip = c(0,9)) # 垂直于x轴的网格线，对应每个刻度
}
p
pdf(file = '01.Uncoxforest.pdf', family = "Times", height = 10, width = 12, onefile = F)
print(p)
dev.off()

png(file = '01.Uncoxforest.png', family = "Times", height = 10, width = 12, units = 'in',res = 600)
print(p)
dev.off()
dev.off()
## 2 多因素cox-------------------
### ph假定检验------------
library(survminer)
cox_more <- coxph(Surv(time = OS.time, event = OS) ~ riskScore + Age + Stage + T_stage + M_stage, data = train_risk_clinical)
cox_zph <- cox.zph(cox_more)
summary(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]      #PH假定检验，p>0.05,都过了

cox_results <- list()
ph_test_results <- list()
uniSigExp <- train_risk_clinical[,c("OS","OS.time","Stage","T_stage","M_stage","riskScore","Age")]

for (i in colnames(uniSigExp[,3:(ncol(uniSigExp))])) {
  cox <- coxph(Surv(OS.time, OS) ~ riskScore + Age + Stage+ T_stage + M_stage, data = uniSigExp)
  cox_results[[i]] <- cox
  ph_test <- cox.zph(cox)
  ph_test_results[[i]] <- ph_test
}
ph_test_results

ph_test_data <- data.frame(
  Variable = character(),
  Time = numeric(),
  Schoenfeld_Residuals = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in names(ph_test_results)) {
  temp_data <- data.frame(
    Variable = i,
    Time = ph_test_results[[i]]$time,
    Schoenfeld_Residuals = ph_test_results[[i]]$y[, 1],
    p_value = ph_test_results[[i]]$table[i,3]
  )
  ph_test_data <- rbind(ph_test_data, temp_data)
}

p_value_summary <- unique(ph_test_data[, c("Variable", "p_value")])

pdf("02.PH_Schoenfeld.pdf",width = 4,height = 4)
ggplot(p_value_summary, aes(x = reorder(Variable, -p_value), y = p_value)) +
  geom_segment(aes(xend =Variable, yend = 0), color = "#66C2A5") +  # 绘制线段
  geom_point(size = 3, color = "#FFD700") +  # 绘制点
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +  # 添加阈值线
  annotate(
    "text",
    x = Inf,  # 将标注放在图的最右侧
    y = 0.05,  # 标注的位置与阈值线对齐
    label = "p > 0.05",  # 标注文本
    vjust = -1.2,  # 调整垂直位置
    hjust = 1,  # 调整水平位置
    color = "black",  # 文本颜色
    size = 3  # 文本大小
  ) +
  labs(
    title = "Proportional Hazards Assumption Test",
    x = "phenotype",
    y = "p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,color = "black"))  # 调整 x 轴标签角度
dev.off()

png("02.PH_Schoenfeld.png",width = 1200,height = 1200,res=300)
ggplot(p_value_summary, aes(x = reorder(Variable, -p_value), y = p_value)) +
  geom_segment(aes(xend =Variable, yend = 0), color = "#66C2A5") +  # 绘制线段
  geom_point(size = 3, color = "#FFD700") +  # 绘制点
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +  # 添加阈值线
  annotate(
    "text",
    x = Inf,  # 将标注放在图的最右侧
    y = 0.05,  # 标注的位置与阈值线对齐
    label = "p > 0.05",  # 标注文本
    vjust = -1.2,  # 调整垂直位置
    hjust = 1,  # 调整水平位置
    color = "black",  # 文本颜色
    size = 3  # 文本大小
  ) +
  labs(
    title = "Proportional Hazards Assumption Test",
    x = "phenotype",
    y = "p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,color = "black"))  # 调整 x 轴标签角度
dev.off()

### 多因素绘图------------
res.mul <- coxph(Surv(time = OS.time, event = OS) ~ riskScore + Age + Stage+ T_stage + M_stage, data = train_risk_clinical) %>% summary

res.mul <- cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame() #
rownames(res.mul)

res.mul$Indicators <- c("riskScore","Age(>=60 VS <60)", 
                        "Stage(II VS I)", "Stage(III/IV VS I)",
                        "T_stage(T2 VS T1)", "T_stage(T3 VS T1)","T_stage(T4 VS T1)",
                        "M_stage(M1 VS M0)")
# 添加行名
rownames(res.mul) <- res.mul$Indicators
# 将结果转换为数据框
res.mul <- as.data.frame(res.mul)

# 添加列名
#colnames(res.mul) <- c("conf.int.lower", "conf.int.upper", "coefficients", "p.value")
colnames(res.mul) <- c("hr","low","up","pv","Indicator")

res.mul$p <- signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] <- NA
res.mul$Indicator <- factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value = res.mul$pv,
                        HR = res.mul$hr,
                        HR.95L = res.mul$low,
                        HR.95H = res.mul$up,
                        Indicator = res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res, file = "02.Multivariate_cox.csv", row.names = T)

library(tidyr)
hz <- paste(round(multi_res$HR, 3), "(", round(multi_res$HR.95L, 3), "-", round(multi_res$HR.95H, 3), ")", sep = "")
hz
tabletext <- cbind(c(NA, rownames(multi_res)),
                   c("P value", ifelse(multi_res$p.value<0.05, "< 0.05", round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)", hz))
tabletext
{p <- forestplot(labeltext=tabletext, 
                 graph.pos=4,  #为Pvalue箱线图所在的位置
                 is.summary = c(TRUE,FALSE,FALSE,FALSE,
                                FALSE,FALSE,FALSE,
                                FALSE,FALSE,FALSE),
                 col=fpColors(box="#FFD700", lines="#66C2A5", zero = "gray50"),
                 mean=c(NA,multi_res$HR),
                 lower=c(NA,multi_res$HR.95L), #95%置信区间下限
                 upper=c(NA,multi_res$HR.95H), #95%置信区间上限
                 boxsize=0.15,lwd.ci=3,   #箱子大小，线的宽度
                 ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                 zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
                 colgap=unit(5,"mm"),    #列间隙
                 xticks = c(0,2,4,6), #横坐标刻度
                 lwd.xaxis=2,            #X轴线宽
                 lineheight = unit(1.2,"cm"), #固定行高
                 graphwidth = unit(.6,"npc"), #图在表中的宽度比例
                 cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                 hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
                 # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                 #                 # "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                 #                 "13" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
                 # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
                 #fpTxtGp函数中的cex参数设置各个组件的大小
                 txt_gp=fpTxtGp(label=gpar(cex=1),
                                ticks=gpar(cex=0.8, fontface = "bold"),
                                xlab=gpar(cex = 1, fontface = "bold"),
                                title=gpar(cex = 1.25, fontface = "bold")),
                 xlab="Hazard Ratio",
                 grid = T,
                 # title = "Multivariate",
                 clip = c(0,6)) # 垂直于x轴的网格线，对应每个刻度
}
p
pdf(file = '02.Mulcoxforest.pdf',family = "Times", height = 6, width = 12, onefile = F)
print(p)
dev.off()

png(file = '02.Mulcoxforest.png',family = "Times", height = 6, width = 12, units = 'in',res = 600)
print(p)
dev.off()
dev.off()
## 3 列线图-----------------------------
library(survival)
library(regplot)
library(rms)

train_phenotype3 <- train_risk_clinical
train_phenotype3$id <- row.names(train_phenotype3)
train_phenotype3 <- train_phenotype3[,c("sample","OS.time","OS","riskScore","Age","M_stage")]
rownames(train_phenotype3) <- train_phenotype3$sample # 将第1列改变为行名
train_phenotype3 <- train_phenotype3[, -1] # 删除原来的第1列（现在已变为行名）

#绘制列线图
res.cox=coxph(Surv(OS.time, OS) ~ . , data = train_phenotype3)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=train_phenotype3[12,],
             rank="sd",
             failtime = c(365,730,1095),
             prfail = F)
pdf("03.Nomogram_line_points_.pdf", family = "Times", width = 10, height = 6)
print(nom1)
dev.off()

nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(train_phenotype3, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
outTab <- outTab[-1, ]
write.csv(outTab, file="03.nomoRisk.csv", row.names=T)

## 4 绘制校准曲线-----------------------------
d <- datadist(rt)
options(datadist='d')

pdf(file="04.Calibration.pdf", family = "Times", w=5, h=5)
f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=365)
cal <- rms::calibrate(f, cmethod="KM", method="boot", u=365, m=(nrow(rt)/3), B=20)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=2, col="#66C2A5", sub=F)

f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=730)
cal <- rms::calibrate(f, cmethod="KM", method="boot", u=730, m=(nrow(rt)/3), B=20)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=2, col="#FFD700", sub=F, add=T)

f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1095)
cal <- rms::calibrate(f, cmethod="KM", method="boot", u=1095, m=(nrow(rt)/3), B=20)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=2, col="#8DA0CB", sub=F, add=T)
legend('bottomright', c('1-year', '2-years', '3-years'),
       col=c("#66C2A5","#FFD700","#8DA0CB"), lwd=2, bty = 'n')
dev.off()

png(file="04.Calibration.png", family = "Times", w=1500, h=1500, units = "px", res = 300)
f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=365)
cal <- rms::calibrate(f, cmethod="KM", method="boot", u=365, m=(nrow(rt)/3), B=20)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=2, col="#66C2A5", sub=F)

f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=730)
cal <- rms::calibrate(f, cmethod="KM", method="boot", u=730, m=(nrow(rt)/3), B=20)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=2, col="#FFD700", sub=F, add=T)

f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1095)
cal <- rms::calibrate(f, cmethod="KM", method="boot", u=1095, m=(nrow(rt)/3), B=20)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=2, col="#8DA0CB", sub=F, add=T)
legend('bottomright', c('1-year', '2-years', '3-years'),
       col=c("#66C2A5","#FFD700","#8DA0CB"), lwd=2, bty = 'n')
dev.off()


## 5 ROC曲线-----------------------------
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
                           risk_score_table = train_phenotype3)
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
                     labels = c('1 year', '2 years', '3 years'),
                     values = c('#66C2A5','#FFD700',"#8DA0CB")) +
  geom_roc(labels = FALSE, stat = 'identity') +
  style_roc() +
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = '') +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
        text=element_text(family='Times')) +
  annotate('text', x=0.75, y=c(0.25, 0.15, 0.05), label = c(paste('AUC of 1 year =', format(auc_y1,nsmall=2)), paste('AUC of 2 years =', format(auc_y2,nsmall=2)), paste('AUC of 3 years =', format(auc_y3,nsmall=2))))
ROC

pdf('05.ROC.pdf', w=5,h=4,family='Times')
print(ROC)
dev.off()

png('05.ROC.png',w=5,h=4,units='in',res=600,family='Times')
print(ROC)
dev.off()

