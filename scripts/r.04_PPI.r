rm(list = ls())
setwd("/data/nas1/zhuxuying/00.BRRU/CRC/15.PX-ACR")
if(!dir.exists("04_PPI")){dir.create("04_PPI")}
setwd("./04_PPI")

# 加载R包
library(circlize)
library(grid)
library(graphics)
library(ComplexHeatmap)
library(dplyr)
library(readr)
select=dplyr::select

# 下载ppi网络信息 0.15
# options(timeout = 9999999)
# get.string(gene$symbol, score = 400, species = 9606)

string_data <- read_tsv('string_interactions_short.tsv') %>%###0.15
  dplyr::select(c("preferredName_A"="#node1", "preferredName_B" = "node2", "score" = "combined_score")) %>%
  filter(score>=0.15)

p_fun <- function(string_data,select_dregee=1){
  
  # 步骤1：初始化过滤过程
  filtered_data <- string_data  # 先使用原始数据
  has_changes <- TRUE  # 设置初始条件以进入循环
  
  # 步骤2：迭代过滤，直到不再有变化
  while(has_changes) {
    # 获取当前数据中的所有节点
    all_nodes <- c(filtered_data$preferredName_A, filtered_data$preferredName_B)
    
    # 计算每个节点的度数（连接数）
    node_degrees <- as.data.frame(table(all_nodes))
    colnames(node_degrees) <- c("node", "degree")
    
    # 找出度数≥2的节点
    nodes_to_keep <- node_degrees$node[node_degrees$degree >= select_dregee]
    
    # 过滤数据，只保留度数≥2的节点之间的连接
    new_filtered_data <- filtered_data %>%
      filter(preferredName_A %in% nodes_to_keep & 
               preferredName_B %in% nodes_to_keep)
    
    # 检查是否有变化
    if(nrow(new_filtered_data) == nrow(filtered_data)) {
      has_changes <- FALSE  # 没有变化，退出循环
    } else {
      filtered_data <- new_filtered_data  # 更新数据，继续迭代
    }
  }
  
  # 步骤3：检查结果
  print(paste("原始数据行数:", nrow(string_data)))
  print(paste("过滤后数据行数:", nrow(filtered_data)))
  
  # 最终确认所有保留的节点都有度数≥2
  final_nodes <- c(filtered_data$preferredName_A, filtered_data$preferredName_B)
  final_degrees <- as.data.frame(table(final_nodes))
  colnames(final_degrees) <- c("node", "degree")
  
  # 检查是否有度数<2的节点
  low_degree_nodes <- final_degrees$node[final_degrees$degree < select_dregee]
  if(length(low_degree_nodes) > 0) {
    print(paste0("警告：仍有度数<",select_dregee,"的节点:"))
    print(low_degree_nodes)
  } else {
    print(paste0("成功：所有节点的度数都≥",select_dregee))
  }
  
  print(paste0('包含节点：', length(unique(c(filtered_data$preferredName_A,filtered_data$preferredName_B)))))
  print(paste0('包含边：', nrow(filtered_data)))
  
  library(ComplexHeatmap)
  library(circlize)
  # 准备数据
  chord_data <- data.frame(
    from = filtered_data$preferredName_A,
    to = filtered_data$preferredName_B,
    value = filtered_data$score
  )
  
  # 获取唯一的基因名称
  all_genes <- unique(c(chord_data$from, chord_data$to))
  
  # 计算每个基因的连接度（与其他节点的连接数量）
  gene_connectivity <- table(c(chord_data$from, chord_data$to))
  gene_connectivity <- gene_connectivity[all_genes]  # 确保顺序一致
  
  # 按照连接度从高到低排序基因
  sorted_genes <- names(sort(gene_connectivity, decreasing = TRUE))
  
  # 创建连接度的颜色渐变 - 使用ComplexHeatmap的colorRamp2函数
  connectivity_range <- range(gene_connectivity)
  connectivity_col_fun <- colorRamp2(
    breaks = c(connectivity_range[1], mean(connectivity_range), connectivity_range[2]),
    colors = c("blue", "green", "red")
  )
  sector_colors <- connectivity_col_fun(as.numeric(gene_connectivity))
  names(sector_colors) <- all_genes
  
  # 创建score的颜色渐变函数
  score_range <- range(chord_data$value)
  score_col_fun <- colorRamp2(
    breaks = c(score_range[1], mean(score_range), score_range[2]),
    colors = c("lightblue", "purple", "darkred")
  )
  
  par(mar = c(1, 1, 3, 1))  # 增加顶部边距
  
  # 清除任何现有的circos布局
  circos.clear()
  
  # 设置基于扇区数量的较小间隔
  n_sectors <- length(all_genes)
  gap_degree <- min(2, 360/(n_sectors * 3))  # 自适应间隔大小
  
  # 设置圆形布局参数 - 从12点钟方向开始
  circos.par(
    gap.degree = gap_degree,
    start.degree = 90,  # 90度是12点钟方向
    # 调整画布大小，留出更多空间给标签
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1.2, 1.2)  # 上方留出更多空间
  )
  
  # 创建弦图，使用order参数按连接度排序
  chordDiagram(
    chord_data,
    grid.col = sector_colors,  # 根据连接度的颜色
    col = score_col_fun,  # 连线根据score着色
    transparency = 0.5,
    link.lwd = 1.5,
    link.lty = 1,
    link.sort = TRUE,
    link.decreasing = TRUE,
    annotationTrack = "grid",  # 只显示网格线，不显示标签
    order = sorted_genes  # 按照连接度排序基因
  )
  
  # 添加更外围的标签
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(
        x = mean(c(get.cell.meta.data("xlim"))),
        y = get.cell.meta.data("ylim")[1] + 1.2,  # 将标签放置得更靠外
        labels = sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        # cex = 0.3
        cex = 0.7
      )
    },
    bg.border = NA
  )
  title("", line = -2)
  # 添加标题 - 放置得更高
  title(
    "Protein-Protein Interaction Network",
    line = 0.3,  # 增加此值会将标题放置得更高
    cex.main = 1.2,  # 稍微增大标题字体
    font.main = 2    # 粗体标题
  )
  
  # 使用ComplexHeatmap的Legend函数创建图例
  # 为连接度创建图例
  connectivity_legend <- Legend(
    title = "Degree",
    col_fun = connectivity_col_fun,
    at = round(seq(connectivity_range[1], connectivity_range[2], length = 5)),
    labels = round(seq(connectivity_range[1], connectivity_range[2], length = 5)),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  
  # 为score值创建图例
  score_legend <- Legend(
    title = "Interaction Score",
    col_fun = score_col_fun,
    at = round(seq(score_range[1], score_range[2], length = 5), 3),
    labels = round(seq(score_range[1], score_range[2], length = 5), 3),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  
  # 将两个图例打包在一起
  legend_list <- packLegend(connectivity_legend, score_legend, direction = "vertical")
  
  # 绘制图例
  pushViewport(viewport(x = 0.85, y = 0.15, width = 0.3, height = 0.3))
  draw(legend_list)
  popViewport()
  
  # 重置circlize参数
  circos.clear()
  
}

# [1] "原始数据行数: 27"
# [1] "过滤后数据行数: 27"
# [1] "成功：所有节点的度数都≥1"
# [1] "包含节点：10"
# [1] "包含边：27"
library(Cairo)
CairoPDF("01.PPI.pdf", width=8, height=8)
p_fun(string_data)
dev.off()

CairoPNG("01.PPI.png", width=8, height=8, units="in", res=300)
p_fun(string_data)
dev.off()



