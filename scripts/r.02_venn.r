rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (!dir.exists('./02_venn')) {
  dir.create('./02_venn')
}
setwd('./02_venn')

gene1 <- read.csv('../01_DEG_DESeq2/01.DEG_sig.csv', row.names = 1)
gene2 <- read.csv('../00_rawdata/ACR.csv') %>%
  unique() %>% .[[1]] %>% trimws()

inter <- data.frame(symbol = intersect(rownames(gene1), gene2))
write.csv(inter, file = '01.Intersection_gene.csv', row.names = FALSE, quote = FALSE)

library(ggvenn)
mydata <- list('DEGs' = rownames(gene1), 'ACR' = gene2)
pdf('01.venn.pdf', w = 6, h = 5, family = 'Times')
ggvenn(mydata, c('DEGs', 'ACR'),
       fill_color = c('#5CB85C', '#FFA500'),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 1,
       text_size = 4,
       stroke_color = '#7e530e',
       stroke_linetype = 'solid',
       set_name_color = c('#5CB85C', '#FFA500'),
       text_color = 'black')
dev.off()

png('01.venn.png', w = 6, h = 5, units = 'in', res = 600, family = 'Times')
ggvenn(mydata, c('DEGs', 'ACR'),
       fill_color = c('#5CB85C', '#FFA500'),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 1,
       text_size = 4,
       stroke_color = '#7e530e',
       stroke_linetype = 'solid',
       set_name_color = c('#5CB85C', '#FFA500'),
       text_color = 'black')
dev.off()
