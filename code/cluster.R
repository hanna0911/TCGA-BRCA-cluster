# 设置工作目录
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # 设置当前R语言源文件目录作为工作目录, reference: https://blog.csdn.net/deephill/article/details/104818274

# 读取、预处理数据
GeneMatrix <- read.table(file = "./GeneMatrix.txt", header = TRUE, sep = "\t", row.names = 1) # 读取GeneMatrix.txt数据, reference: https://www.cnblogs.com/caiyishuai/p/9729554.html
GeneMatrix <- t(scale(GeneMatrix)) # 数据标准化，矩阵转置

# 层次聚类
distance <- dist(GeneMatrix, method = "euclidean") # 计算距离
fit <- hclust(distance, method = "average") # 聚类，reference: https://blog.csdn.net/weixin_39692037/article/details/110289511
pdf("./cluster.pdf", width = 25.0, height = 10.0) # 保存图片，reference: https://blog.csdn.net/weixin_41929524/article/details/100152389
plot(fit, hang = -1, cex = 0.3) # 作图
dev.off()

# 热点图
pdf("./heatmap.pdf", width = 15.0, height = 15.0)
heatmap(as.matrix(distance)) # 聚类并作图，reference: https://www.cnblogs.com/payton/p/4240824.html
dev.off()

# PCA
# reference: https://www.jianshu.com/p/c3b9653ed605
pca <- prcomp(GeneMatrix, center = FALSE, scale = FALSE) # prcomp函数横行是样本，矩阵转置
pdf("./pca_screeplot.pdf", width = 20.0, height = 12.0)
screeplot(pca, type = "barplot", npcs = length(pca$sdev), main = "Scree Plot", xlab = "Principal Component") # 可视化
dev.off()

# 决定PCA主成分数量
pdf("./pca_cumulative.pdf", width = 20.0, height = 12.0)
plot(summary(pca)$importance[3:3, 1:389], ylim = c(0.0, 1.0), xlab = "Principal Component", ylab = "Cumulative Proportion", type = "l", las = 1) # 选取Cumulative Proportion作图
abline(v = 172, lwd = 1, lty = 3, col = "blue")
abline(h = 0.95, lwd = 1, lty = 3, col = "blue")
axis(side = 2, at = c(0.95), las = 2, labels = c("0.95"))
axis(side = 1, at = c(172), labels = c("172"))
dev.off()

pca.data <- data.frame(pca$x[,1:172]) # 选取前172个主成分

pca_distance <- dist(pca.data, method = "euclidean") # 计算距离
pca_fit <- hclust(pca_distance, method = "average") # 聚类
pdf("./pca_cluster.pdf", width = 25.0, height = 10.0) # 保存图片
plot(pca_fit, hang = -1, cex = 0.3) # 作图
dev.off()

# 热点图
pdf("./pca_heatmap.pdf", width = 15.0, height = 15.0)
heatmap(as.matrix(pca_distance)) # 聚类并作图
dev.off()

# 与clinical_data ER_Status_nature2012分类方式进行对比
# 导入clinical_data数据，并保留ER_Status_nature2012列
clinical_data <- read.table(file = "./clinical_data.txt", header = TRUE, sep = "\t", row.names = 1)
ER_Status_nature2012 <- clinical_data[7:7]

# 处理ER_Status_nature2021的行名（病例名）
for (i in 1:771)
row.names(ER_Status_nature2012)[i] <- chartr("-", ".", row.names(ER_Status_nature2012)[i])

# 选取GeneMatrix与clinical_data病例的交集
ER_Status_nature2012$patient <- row.names(ER_Status_nature2012)
ER_Status_nature2012 <- subset(ER_Status_nature2012, patient %in% row.names(GeneMatrix))
ER_Status_nature2012 <- ER_Status_nature2012[1:1]

# 将ER_Status_nature2012指标由chr类型转为factor类型再转为numeric类型
library(dplyr)
ER_Status_nature2012 <- mutate_at(ER_Status_nature2012, vars(ER_Status_nature2012), as.factor) # reference: https://www.codenong.com/9251326/
ER_Status_nature2012$ER_Status_nature2012 <- as.numeric(ER_Status_nature2012$ER_Status_nature2012)
ER_Status_nature2012 <- scale(ER_Status_nature2012)

# 层次聚类
ER_distance <- dist(ER_Status_nature2012, method = "euclidean") # 计算距离
ER_fit <- hclust(ER_distance, method = "average") # 聚类
pdf("./ER_cluster.pdf", width = 25.0, height = 10.0) # 保存图片
plot(ER_fit, hang = -1, cex = 0.3) # 作图
dev.off()

# 热点图
pdf("./ER_heatmap.pdf", width = 15.0, height = 15.0)
heatmap(as.matrix(ER_distance)) # 聚类并作图
dev.off()
