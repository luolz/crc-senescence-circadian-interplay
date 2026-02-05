rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists("12_Drug_sensitivity")){dir.create("12_Drug_sensitivity")}
setwd("12_Drug_sensitivity")

library(tidyverse)
library(lance)
library(oncoPredict)
library(ggplot2)

GDSC2_expr <- readRDS('/data/nas1/zhuxuying/common/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_res <- readRDS('/data/nas1/zhuxuying/common/GDSC2_Res.rds')

train_data <- read.csv(file = '../00_rawdata/01.fpkmlog2_TCGA.COADREAD_mRNA.csv', row.names = 1, check.names = FALSE)

group <- read.csv('../05_risk_train/03.risk.csv')
group <- subset(group, select = c(id, risk))
colnames(group) <- c("sample","group")
# group$risk <- ifelse(group$risk == 0, "High_risk", "Low_risk")
rownames(group) <- group$sample
#group$group <- as.factor(group$group)
sub_group <- group

train_data <- train_data[,group$sample]
identical(colnames(train_data), group$sample)

calcPhenotype(trainingExprData = GDSC2_expr, 
              trainingPtype = GDSC2_res, 
              testExprData = as.matrix(train_data), 
              batchCorrect = 'eb', 
              powerTransformPhenotype = TRUE, 
              removeLowVaryingGenes = 0.2, 
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' 
)

## 1.1 药物IC50高低风险组wilcox检验-----------------------------------------------------
library(data.table)
library(dplyr)
IC50<-fread('calcPhenotype_Output/DrugPredictions.csv')
IC50 <- IC50 %>% column_to_rownames("V1")
library(stringr)
colnames(IC50) <- str_replace(colnames(IC50), "_.*", "")
IC50$sample <- rownames(IC50)

sub_group$group <- as.factor(sub_group$group)

dat.IC50 <- merge(IC50, sub_group, by = "sample")
dat.IC50_2 <- dat.IC50 %>% 
  pivot_longer(
    cols = -c("sample", "group"),
    names_to = "drug",
    values_to = "Score"
  )
colnames(dat.IC50_2)
# dat.IC50_2 <- separate(dat.IC50_2, drug, into = c('drug', 'id'), sep = '_')
# dat.IC50_2 <- dat.IC50_2[,-4]

library(rstatix)
## 差异分析
stat_res <- dat.IC50_2 %>% 
  group_by(drug) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p") #190

DE.res <- stat_res[which(stat_res$p < 0.05), ] # 73个药物显著

## 导出数据
write.csv(stat_res, file = '01.stat.IC50.csv')
write.csv(DE.res, file = '01.DE.IC50.csv')
#library(xlsx)

### 数据处理-----
## 计算每组的评分总和
drug_sensitivity_summary <- dat.IC50_2 %>%
  group_by(drug, group) %>%
  summarise(total_score = sum(Score, na.rm = TRUE)) %>%
  pivot_wider(names_from = group, values_from = total_score) %>%
  rename_with(~ c("drug", "high_risk_sum", "low_risk_sum"), 
              .cols = c(drug, `high risk`, `low risk`))

library(dplyr)
stat_res_df <- as.data.frame(stat_res)

print(colnames(stat_res_df))
print(colnames(drug_sensitivity_summary))

final_result <- stat_res_df %>%
  dplyr::left_join(drug_sensitivity_summary, by = "drug") %>%
  dplyr::mutate(
    sensitivity_direction = dplyr::case_when(
      p > 0.05 ~ "No sensitive",
      high_risk_sum > low_risk_sum ~ "Sensitive to low risk",
      TRUE ~ "Sensitive to high risk"
    )
  ) %>%
  dplyr::select(drug, p, p.adj, p.signif, high_risk_sum, low_risk_sum, sensitivity_direction)
table(final_result$sensitivity_direction)
# No sensitive Sensitive to high risk  Sensitive to low risk 
#   117                     68                      5 
write.csv(final_result, file = "03.Drug_Sensitivity_Direction_with_pvalue.csv", row.names = FALSE)

### 绘图-----
stat.test <- final_result
red <- data.frame(drug = stat.test$drug[stat.test$sensitivity_direction == "Sensitive to high risk"])
# add x and y
red$x <- c(rep(1, 25), rep(2, 25), rep(3, 18))
red$y <- c(seq(11, 35, 1), seq(11, 35, 1),seq(11, 28, 1))
###35=10+25, 16=10+6
# organize blue drug
blue <- data.frame(drug = stat.test$drug[stat.test$sensitivity_direction == "Sensitive to low risk"])
# add x and y
blue$x <- c(rep(1, 5))
blue$y <- c(seq(11, 15, 1))

# organize gray drug
gray <- data.frame(drug = stat.test$drug[stat.test$sensitivity_direction == "No sensitive"])
# add x and y
gray$x <- c(rep(1, 25), rep(2, 25),rep(3, 25), rep(4, 25), rep(5, 17))
gray$y <- c(seq(11, 35, 1), seq(11, 35, 1),seq(11, 35, 1), seq(11, 35, 1),seq(11, 27, 1))

###---- IC50 fig ----####
# p0 data 
plot.p0 <- stat.test
num <- table(plot.p0$sensitivity_direction) %>% data.frame()

# organzie p0 data
p0.data <- data.frame(Parties = factor(c("Sensitive to high risk", "No sensitive", "Sensitive to low risk"),
                                       levels = c("Sensitive to high risk", "No sensitive", "Sensitive to low risk")),
                      seats = c(num$Freq[2], num$Freq[1], num$Freq[3]),   #需要与table的数量一致
                      colors = c("#FFD700", "grey", "#66C2A5"),
                      stringsAsFactors = FALSE)

# p0 fig
library(ggpol)
p0 <- ggplot(p0.data) + 
  geom_parliament(aes(seats = seats, fill = Parties), color = "black") + 
  scale_fill_manual(values = p0.data$colors, labels = p0.data$Parties) +
  coord_fixed() + 
  theme_void()+
  ggtitle("Medicinal Sensity")+
  theme(plot.title = element_text(color = 'black', size = 24, face = "bold", hjust = 0.5),
        text = element_text(family = 'Times'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(color = 'black', size = 12, face = "bold"),
        legend.title = element_text(color = 'black', size = 18, face = "bold"),
        legend.position = "top")
p0

# read red data file
# data.p1 <- read.csv('04.plot_red.csv')
data.p1 <- red
table(data.p1$x)
range(data.p1$y)

# p1 fig add red table
p1 <- ggplot(data.p1, aes(x = x, y = y)) +
  geom_point(color = 'white') +
  geom_text(label = data.p1$drug, check_overlap = T,
            vjust = 0.5, size = 4, hjust = 0, nudge_x = 0, alpha = 0.9, col = '#FFD700')+
  theme_void()+
  xlim(1, 4) +######3=2+1
  ylim(11, 35) +
  theme(text = element_text(angle = 0, hjust = 0, family = 'Arial'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "#FFD700", fill = NA, size = 1))
p1

# read gray data file
data.p2 <- gray
table(data.p2$x)
range(data.p2$y)

# p2 fig add gray table
p2 <- ggplot(data.p2, aes(x = x, y = y)) +
  geom_point(color = 'white') +
  geom_text(label = data.p2$drug, check_overlap = T,
            vjust = 0.5, size = 4, hjust = 0, nudge_x = 0, alpha = 0.9, col = 'grey') +
  theme_void() +
  xlim(1, 8) +
  ylim(11, 35) +
  theme(text = element_text(angle = 0, hjust = 0, family = 'Arial'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, size = 1))
p2

# read gray data file
data.p3 <- blue
table(data.p3$x)
range(data.p3$y)

# p3 fig add blue table
p3 <- ggplot(data.p3, aes(x = x, y = y)) +
  geom_point(color = 'white') +
  geom_text(label = data.p3$drug, check_overlap = T,
            vjust = 0.5, size = 4, hjust = 0, nudge_x = 0, alpha = 0.9, col='#66C2A5')+
  theme_void()+
  xlim(1, 3) +
  ylim(11, 20) +
  theme(text = element_text(angle = 0, hjust = 0, family='Arial'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "#66C2A5", fill = NA, size = 1))
p3

library(ggpubr)
# arrange fig 
p <- ggarrange(p3, p2, # 要排版的图形
               ncol = 2, nrow = 1, # 行数和列数为2
               widths = c(2, 8)) # 宽度为2：1
p
pp <- ggarrange(p, p1, # 要排版的图形
                ncol = 2.5, nrow = 1, # 行数和列数为2
                widths = c(2, 1)) # 宽度为2：1
pp
ppp <- ggarrange(p0, pp, # 要排版的图形s
                 ncol = 1, nrow = 2, # 行数和列数为2
                 widths = c(2, 1)) # 宽度为2：1
ppp

# save ic50 fig
ggsave(filename = "04.IC50.pdf", ppp, width = 16, height = 9, family = 'Times')
ggsave(filename = "04.IC50.png", ppp, width = 16, height = 9, units = "in", dpi = 600, bg = "white")


# ## 2 差异分析-----
# ## 根据p值排序导出
# library(reshape)
# write.csv(stat_res,'05.stat_res.csv')
# dif_stat_res<-stat_res[stat_res$p<0.05,]
# dif_stat_res<-dif_stat_res[order(dif_stat_res$p,decreasing = F),]
# write.csv(dif_stat_res,'05.dif_stat_res.csv')
# 
# ## 小提琴图
# colnames(DE.res)
# DE.res_1 <- DE.res
# DE.res_1 <- DE.res_1[order(DE.res_1$p), ]
# DE.res_1 <- DE.res_1[1 : 10, ] #数字表示后面展示的药物
# 
# # violin.cibersort1 <- dat.IC50_2
# violin.cibersort1 <- dat.IC50_2[dat.IC50_2$drug %in% DE.res_1$drug, ]
# # violin.cibersort1 <- separate(violin.cibersort1, drug, into = c('drug', 'id'), sep = '_')
# class(violin.cibersort1$group) # 检查分组是否为因子，如果不是要转换成因子
# 
# library(ggpubr)
# p1 <- ggplot(violin.cibersort1, aes(x = drug, y = Score, fill = group)) +
#   geom_violin(trim=F, color="black", aes(fill = group)) + #绘制小提琴图, “color”设置小提琴图的轮廓线的颜色(不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
#   #"trim"如果为TRUE(默认值), 则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
#   stat_boxplot(geom = "errorbar",
#                width = 0.1,
#                position = position_dodge(0.5)) +
#   geom_boxplot(aes(x = drug, y = Score, fill = group),
#                width = 0.3,
#                position = position_dodge(0.5),
#                outlier.shape = NA,
#                outlier.colour = NA)+ #绘制箱线图，此处width=0.1控制箱线图的宽窄
#   scale_fill_manual(values = c("#337AB7", "#87CEEB"), name = "Group")+
#   labs(title = "IC50", x = "", y = "IC50", size = 20) +
#   stat_compare_means(data = violin.cibersort1,
#                      mapping = aes(group = group),
#                      label = "p.signif",
#                      method = 'wilcox.test',
#                      paired = F) +
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
#         axis.text.x = element_text(angle = 45, hjust=1, colour = "black", face = "bold", size = 12),
#         axis.text.y = element_text(hjust = 0.5, colour ="black", face="bold", size=12),
#         axis.title.x = element_text(size = 20, face = "bold"),
#         axis.title.y = element_text(size = 16, face = "bold"),
#         legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
#         legend.title = element_text(face = "bold", size = 12),
#         legend.position = "top",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# p1
# ggsave(filename = '05.IC50_violin_top10.pdf', family = "Times", p1, w=8, h=6)
# ggsave(filename = '05.IC50_violin_top10.png', p1, w=8, h=6, dpi = 600)