# 加载必要的库
library(ggplot2)
library(dplyr)
library(stringr)

# 读取数据
data <- read.table("selcount.ONT.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("Chromosome", "Start", "End", "Reads")

# 提取染色体名称中的数字部分，并转换为数值
data <- data %>%
  mutate(Chromosome_Number = as.numeric(str_extract(Chromosome, "\\d+")))

# 将区间转换为更友好的格式（例如 "0-1Mb"）
data <- data %>%
  mutate(Interval = paste0(Start / 1e6, "-", End / 1e6, "Mb"))

# 按染色体数字部分排序
data <- data %>%
  arrange(Chromosome_Number)

# 将 Chromosome 列转换为因子，并按照数字顺序设置水平
data$Chromosome <- factor(data$Chromosome, levels = unique(data$Chromosome))

# 创建一个新的列，用于控制横坐标标签的显示
data <- data %>%
  group_by(Chromosome) %>%
  mutate(Interval_Index = row_number()) %>%
  ungroup()

# 自定义横坐标标签：每隔 10 个区间显示一个标签
custom_labels <- function(breaks) {
  labels <- ifelse(breaks %% 10 == 1, data$Interval[breaks], "")
  return(labels)
}

# 按染色体分面绘制柱状图
ggplot(data, aes(x = Interval_Index, y = Reads, fill = Chromosome)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +  # 每行显示 4 个染色体
  scale_x_continuous(
    breaks = seq(1, max(data$Interval_Index)),  # 设置所有区间为断点
                 labels = custom_labels  # 自定义标签显示
    ) +
  labs(
        title = "Reads Distribution per 1 Mb Interval",
        x = "Interval (Mb)",
        y = "Read Count"
      ) +
  theme_minimal() +
  theme(
        axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        strip.text = element_text(size = 10, face = "bold")  # 调整分面标题字体
        )
    
    # 保存图像
ggsave("reads_distribution_plot.png", width = 14, height = 8, dpi = 300)
print("Plot saved to reads_distribution_plot.png")
    
# 載入必要的套件
library(ggplot2)
library(readr)

setwd("D:/R")
# 讀取數據
data <- read_tsv("chromosome_bin_counts2.tsv")

# 繪製條形圖
ggplot(data, aes(x = bin, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_wrap(~ chromosome) +
  labs(x = "Bin", y = "Count", title = "Chromosome Bin Counts") +
  theme_minimal()
