rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if (! dir.exists("./17_Hubcell")){
  dir.create("./17_Hubcell")
}
setwd("./17_Hubcell")
load('../16_scRNA/04_SingleR/UMAP_singleR.Rdata')

######01.比例图########
table(UMAP$group)
##tumor
scRNA_CC<-subset(UMAP,group=='Tumor')
phe=scRNA_CC@meta.data
table(phe$celltype)
cell_type_CC <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_CC)<-c('Celltype','Tumor')

##Healthy 
scRNA_Normal<-subset(UMAP,group=='Normal')
phe=scRNA_Normal@meta.data
table(phe$celltype)
cell_type_Normal <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_Normal)<-c('Celltype','Normal')

cell_type<-merge(cell_type_Normal,cell_type_CC,by='Celltype')

i<-1
while(i<10){ 
  cell_type$All[i]<-sum(cell_type$CC[i],cell_type$Normal[i])
  i<-i+1
}

write.csv(cell_type,'01.cell_type.csv',row.names = F)

##换个形式
all.count <- data.frame(celltype=UMAP$celltype, sample=UMAP$group)
table(all.count$celltype)

library(reshape2)
plot.celltype = table(all.count$celltype, all.count$sample) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype = cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)[2:5] = c("Count_Normal","Count_Tumor","Proportion_Normal","Proportion_Tumor")
plot.celltype <- plot.celltype[order(plot.celltype$Count_Normal,decreasing = T),]
write.csv(plot.celltype, "CellCount.csv", row.names = F,quote = F)

dat.celltype = plot.celltype[c(1,4:5)]

colnames(dat.celltype)[2:3] = c("Tumor","Normal")
dat.celltype[2:3] = 100 * dat.celltype[2:3]

dat.plot = melt(dat.celltype, id.vars = "Cell", variable.name = "Sample", value.name = "value")
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
library(IOBR)

color2 <- c('#ABD0F1','#E56F5E','#92B4C8','#F19685','#F6C957','#FFB77F','#FBE8D5','#EEA599',"#cab2d6","#377EB8")
dat.plot$label = round(dat.plot$value, 2) %>% paste0(.,"%")
a <- ggplot(data = dat.plot, mapping = aes(x = value, y = Cell)) + 
  geom_text(aes(x = value + 1, label = label), hjust = 0) +
  geom_bar(aes(fill = Cell), stat = "identity") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values = palettes(category = "random",22,show_col=F)) + 
  # scale_fill_manual(values = col) + 
  xlim(c(0,70)) +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 16) + 
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"), legend.position = "none",
        panel.grid = element_blank())
a
pdf("01.CellProportion.pdf", width = 9, height = 4,family = "Times")
a
dev.off()
png("01.CellProportion.png", width = 9, height = 4,res = 600,units = 'in',family = "Times")
a
dev.off()
dev.off()


#3. 生物标志物在各细胞类型的表达情况------
gene <- readRDS('../05_risk_train/inputgene.rds')##
vln.df=as.data.frame(UMAP[["RNA"]]$data[gene,])
vln.df$symbol=rownames(vln.df)
vln.df=melt(vln.df,id="symbol")
colnames(vln.df)[c(2,3)]=c("sample","exp")
head(vln.df)
anno=UMAP@meta.data[,c("group","celltype")]
anno$sample <- rownames(anno)
vln.df=inner_join(vln.df,anno,by="sample")
vln.df$group <- factor(vln.df$group,levels = c('Tumor','Normal'))

my_comparisons <- list(c('Tumor','Normal'))
table(vln.df$celltype)
vln.dat <- vln.df
celltype <- as.data.frame(table(vln.dat$celltype))[1]
#celltype <- celltype[-10,] %>% as.data.frame()
colnames(celltype) <- 'Var1'
plots <- list()
for (i in c(1:nrow(celltype))) {
  plot.dat <- vln.dat[which(vln.dat$celltype==celltype$Var1[i]),]
  position=c(max(plot.dat$exp)+0.5,max(plot.dat$exp)+1,max(plot.dat$exp+1.5))
  
  p <- ggplot(plot.dat,aes(x=group,y=exp,fill=group))+
    geom_violin(trim = F,color='black')+
    stat_boxplot(geom = 'errorbar',
                 width=0.1,
                 position = position_dodge(0.9))+
    geom_boxplot(width=0.4,
                 position = position_dodge(0.9),
                 outlier.shape = NA,fill='white')+
    scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
    labs(title=celltype$Var1[i], x="", y = "",size=20) +
    geom_signif(comparisons = my_comparisons,
                test = t.test,
                map_signif_level = T,
                y_position = position,textsize = 2.5)+
    ggtitle(plot.dat$celltype[1])+
    ylim(x=c(0,max(position+0.3))) +
    theme_bw()+
    theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=10),
          axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10),
          axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
          axis.title.x=element_text(size=16,face="bold"),
          axis.title.y=element_text(size=16,face="bold"),
          legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
          legend.title = element_text(face = "bold", size = 12),
          legend.position = "top",
          panel.grid = element_line(color = 'gray',size = 0.2),
          panel.border = element_rect(colour = 'black',fill = NA,size = 1))+
    facet_wrap(~symbol,ncol = 11,scales = 'free')
  p
  plots[[i]] <- p
}
p
library(gridExtra)
pdf(file = '02.gene.solo.pdf',w=16,h=10,family = "Times")
grid.arrange(grobs = plots, ncol = 3)
dev.off()

png(file = '02.gene.solo.png',w=16,h=10,units = 'in',res = 300,family = "Times")
grid.arrange(grobs = plots, ncol = 3)
dev.off()

#关键细胞选作B cells
