#增强子课题分析差异表达基因
library('DESeq2')
library('RColorBrewer')
library(tximport)
library('corrplot')
library(stringr)
library("magrittr")
library(tidyverse)
library(hrbrthemes)
library(ggtext)
library(gridExtra)
library(ggplot2)
library(lattice)
library(ggpubr)
library(patchwork)
library(scales)
library(tidyr)
#install.packages("ggrepel")
#install.packages("matrixStats")
library(matrixStats)
library(ggrepel)
require(pals)
library(rdist)
library(rpart)
library(ggthemes)
library(ComplexHeatmap)
library(ggdendro)
#利用DESeq2构建表达量矩阵
library('DESeq2')
library('RColorBrewer')
library(ggnewscale)
library(clusterProfiler)
library(enrichplot)
library(org.Twheat.eg.db)
library(simplifyEnrichment)
library(VennDiagram)
filePath <- "D:/高粱/RNA-seq/RNA-seq/"
setwd("D:/高粱/RNA-seq/")
s2c <- read.table("design_matrix.txt", header = TRUE, sep='\t',stringsAsFactors=FALSE)
sampleNames <- s2c[,1]
countData.list <- sapply(sampleNames, function(x) read.table(file=paste0(filePath, x, "/abundance.tsv"), header=T, sep="\t"), simplify=F)
countData.df <- do.call("cbind", countData.list)
colsToKeep <- c(1,grep("est_count", names(countData.df)))
ct <- countData.df[,colsToKeep]
#构建count矩阵
names(ct) <- c("transcript_id", sampleNames)
ct[,2:10] <- round(ct[,2:10])
write.csv(ct,"count_matrix.csv",row.names=F)
#构建tpm矩阵
colsToKeep <- c(1,grep("tpm", names(countData.df)))
ct <- countData.df[,colsToKeep]
names(ct) <- c("transcript_id", sampleNames)
write.csv(ct,"TPM_matrix.csv",row.names=F)

#质控
#绘制相关性热图
library('corrplot')
tpm <- read.csv('../TPM_matrix.csv', header = T, row.names = 1)
corr <- cor(tpm, method = 'spearman') 
#corr <- cor(tpm, method = 'pearson')
pdf("spearman_correlation.pdf")
corrplot(corr, type = 'upper', tl.col = 'black', order = 'hclust', tl.srt = 45, addCoef.col = 'white')	#并不能直接计算两两样本之间的相关系数，而需要通过R内置函数cor计算，corrplot真正做的是把cor计算得到的样本相关系数矩阵映射成不同的颜色和样式
# 映射成热图
# type='upper'：只显示右上角相关系数矩阵
# tl.col='black'：字体颜色黑色
# order='hclust'：使用层次聚类算法
# tl.srt = 45：x轴标签倾斜45度
# addCoef.col='white'：添加相关系数数值，颜色白色
dev.off()

#绘制热图
options(stringsAsFactors = F)
a = read.csv("TPM_matrix.csv",header=T,row.name=1)
dim(a)
dat = a
library(stringr)
ac = data.frame(group=str_split(colnames(dat),'_',simplify = T)[,1])
rownames(ac) = colnames(dat)
M=cor(log(dat+1))
pdf("heatmap.pdf")
pheatmap::pheatmap(M, annotation_col = ac, display_numbers=F)
dev.off()

#合并转录本表达量为基因表达量
setwd(("D:/高粱/RNA-seq"))
#设置kallisto生成的路径
base_dir <-("D:/高粱/RNA-seq/RNA-seq")
#获取所有的simple_id,具体可以看一些file.path命令的使用
sample_id <- dir(file.path(base_dir))
sample_id
## [1] "WT_1" "WT_2" "WT_3" Mutation_1"
## [5] "Mutation_2" Muatation_3"
#获取结果文件的路径
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs
##                                                  WT_1 
## "/home/617/sleuth_analysis/kallisto_qaunt_output/WT_1" 
##                                                  WT_2
## "/home/617/sleuth_analysis/kallisto_qaunt_output/WT_2" 
##                                                  WT_3 
## "/home/617/sleuth_analysis/kallisto_qaunt_output/WT_3" 
##                                                  Mutation_1 
## "/home/617/sleuth_analysis/kallisto_qaunt_output/Mutation_1" 
##                                                  Mutation_2
## "/home/617/sleuth_analysis/kallisto_qaunt_output/Mutation_2" 
##                                                  Mutation_3
## "/home/617/sleuth_analysis/kallisto_qaunt_output/Mutation_3"

#读取实验设计表
s2c <- read.table("design_matrix.txt", header = TRUE, sep='\t',stringsAsFactors=FALSE)
s2c
##           sample condition reps
## 1 WT_1        WT    1
## 2 WT_2        WT    2
## 3 WT_3        WT    3
## 4 Mutation_1        Mutation    1
## 5 Mutation_2        Mutation    2
## 6 Mutation_3        Mutation    3
#与路径合并
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
##           sample condition reps
## 1 WT_1        WT    1
## 2 WT_2        WT    2
## 3 WT_3        WT    3
## 4 Mutation_1        Mutation    1
## 5 Mutation_2        Mutation    2
## 6 Mutation_3        Mutation    3
##                                                                 path
## 1 /home/617/sleuth_analysis/kallisto_qaunt_output/WT_1
## 2 /home/617/sleuth_analysis/kallisto_qaunt_output/WT_2
## 3 /home/617/sleuth_analysis/kallisto_qaunt_output/WT_3
## 4 /home/617/sleuth_analysis/kallisto_qaunt_output/Mutation_1
## 5 /home/617/sleuth_analysis/kallisto_qaunt_output/IMutation_2
## 6 /home/617/sleuth_analysis/kallisto_qaunt_output/Mutation_3

#获取转录本ID信息
#方法三：自己制作相应的txt文件第一列为target_id，第二列为对应的gene（可以用ID或者是Symbol）
t2g <- read.table("t2g.txt", header = TRUE, stringsAsFactors=FALSE)

library(sleuth)
#so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, extra_bootstrap_summary = TRUE)

#如果想要将一个基因的不同转录本合并，可以按照以下模式导入
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, num_cores = 1,extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE,gene_mode=T,aggregation_column = 'gene')
#so <- sleuth_prep(s2c, ~ condition, num_cores = 1,target_mapping = t2g,aggregation_column = 'gene', extra_bootstrap_summary=TRUE,read_bootstrap_tpm=TRUE, gene_mode = TRUE)

# HDF5. File accessibility. Unable to open file.

#summarizing bootstraps
#.Error in process_bootstrap(i, samp_name, kal_path, num_transcripts, est_counts_sf[[i]],  : 
#  File h5 has no bootstraps.Please generate bootstraps using "kallisto quant -b".File D:/ZLH/RNA-seq/C0-2_result/abundance.h5 has no bootstraps.Please generate bootstraps using "kallisto quant -b".
#此外: Warning message:
#In sleuth_prep(s2c, ~condition, target_mapping = t2g, num_cores = 1,  :
#  829 target_ids are missing annotations for the aggregation_column: gene.
#These target_ids will be dropped from the gene-level analysis.
#If you did not expect this, check your 'target_mapping' table for missing values.

#两两差异比较
library(stringr)
a <- list(c("V0","V26N0"),c("V0","V26N6"), c("V26N0","V26N6"))

for (xxx in a) {
  s2b <- dplyr::filter(s2c, condition == xxx[1] | condition == xxx[2])
  s2o <- sleuth_prep(s2b, ~ condition, 
                     num_cores = 1,target_mapping = t2g,
                     aggregation_column = 'gene', 
                     extra_bootstrap_summary=TRUE,
                     read_bootstrap_tpm=TRUE,
                     gene_mode = TRUE)
  s2o <- sleuth_fit(s2o)
  s2o <- sleuth_fit(s2o, formula = ~ 1, fit_name = "reduced")
  s2o_lrt <- sleuth_lrt(s2o, "reduced", "full")
  models(s2o_lrt)
  #LRT检验结果
  sleuth_table <- sleuth_results(s2o_lrt, 'reduced:full', 'lrt', show_all = FALSE)
  table(sleuth_table[,"qval"] < 0.05)
  sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
  head(sleuth_significant)
  write.csv(sleuth_table,str_c(xxx[1], "vs", xxx[2],"_sleuth_gene_level.csv"),row.names=TRUE,quote=TRUE)
  write.csv(sleuth_significant,str_c(xxx[1], "vs", xxx[2],"_sleuth_significant_gene_level.csv"),row.names=TRUE,quote=TRUE)
  #导入基因表达量TPM
  sleuth_matrix <- sleuth_to_matrix(s2o_lrt, 'obs_norm', 'tpm')
  head(sleuth_matrix)
  write.csv(sleuth_matrix, str_c(xxx[1],"vs", xxx[2],"_sleuth_tpm_norm_gene_level.csv"),row.names=TRUE,quote=TRUE)
}

#获得差异表达倍数
#a <- list(c("V0","V26N0"),c("V0","V26N6"), c("V26N0","V26N6"))
a <- list(c("V0","V26N6"))
for (xxx in a) {
  s2b <- dplyr::filter(s2c, condition == xxx[1] | condition == xxx[2])
  so <- sleuth_prep(s2b, 
                    ~ condition, 
                    target_mapping = t2g, 
                    read_bootstrap_tpm = TRUE,
                    aggregation_column = 'gene',
                    gene_mode = TRUE,
                    extra_bootstrap_summary = TRUE,
                    transformation_function = function(x) log2(x + 0.5)) %>%
    sleuth_fit()
  models(so)
  oe <- sleuth_wt(so, 
                  which_beta = 'conditionV26N6')
  sleuth_results_oe <- sleuth_results(oe, 
                                      test = 'conditionV26N6', 
                                      show_all = TRUE)
  write.csv(sleuth_results_oe, str_c(xxx[1],"vs", xxx[2],"_sleuth_log2FC.csv"),row.names=TRUE,quote=TRUE)
}



###绘制火山图
#绘制火山图
setwd("D:/高粱/RNA-seq")
#leafCvsrootC
dat<-read.csv("CL6vsCR6_sleuth_significant_gene_level.csv",header=T)
head(dat)
#根据阈值FC>=1和P<0.05筛选显著差异表达基因，为后面火山图的点表上颜色
loc_up <- intersect(which(dat$qval < 0.05),which(dat$b >=1))
loc_down <- intersect(which(dat$qval < 0.05),which(dat$b <= (-1)))
logFC_cutoff <- 1 
dat$sig = as.factor(ifelse(dat$qval < 0.05 & abs(dat$b) > logFC_cutoff,
                           ifelse(dat$b > logFC_cutoff ,'UP','DOWN'),'NOT'))#添加sig这栏up,down,not的行数据
this_tile <- paste0('CL6 vs CR6 Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(dat[which(dat$sig == 'UP'),]),
                    '\nThe number of down gene is ',nrow(dat[which(dat$sig == 'DOWN'),])) #设置标题，体现阈值、上下调基因数量
#基本散点图
p0 <- ggplot(data=dat, aes(x=-log10(qval), y=b, color=sig)) +
  geom_point(alpha=0.4, size=1.5) +
  theme_bw(base_size=15)+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "up")+
  xlab("-log10(qval)") + ylab("log2FC") +
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('#1170aa','#a2a2a2','#d62728')) +
  geom_vline(xintercept=-log(0.05,10) ,linetype=4) +
  geom_hline(yintercept=c(-1,1) ,linetype=4 )


#leaf
dat<-read.csv("CL6vsPL6_sleuth_significant_gene_level.csv",header=T)
head(dat)
#根据阈值FC>=1和P<0.05筛选显著差异表达基因，为后面火山图的点表上颜色
loc_up <- intersect(which(dat$qval < 0.05),which(dat$b >=1))
loc_down <- intersect(which(dat$qval < 0.05),which(dat$b <= (-1)))
logFC_cutoff <- 1 
dat$sig = as.factor(ifelse(dat$qval < 0.05 & abs(dat$b) > logFC_cutoff,
                           ifelse(dat$b > logFC_cutoff ,'UP','DOWN'),'NOT'))#添加sig这栏up,down,not的行数据
this_tile <- paste0('CL6 vs PL6 Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(dat[which(dat$sig == 'UP'),]),
                    '\nThe number of down gene is ',nrow(dat[which(dat$sig == 'DOWN'),])) #设置标题，体现阈值、上下调基因数量
#基本散点图
p1 <- ggplot(data=dat, aes(x=-log10(qval), y=b, color=sig)) +
  geom_point(alpha=0.4, size=1.5) +
  theme_bw(base_size=15)+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "up")+
  xlab("-log10(qval)") + ylab("log2FC") +
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('#1170aa','#a2a2a2','#d62728')) +
  geom_vline(xintercept=-log(0.05,10) ,linetype=4) +
  geom_hline(yintercept=c(-1,1) ,linetype=4 )

#root
dat<-read.csv("CR6vsPR6_sleuth_significant_gene_level.csv",header=T)
head(dat)
#根据阈值FC>=1和P<0.05筛选显著差异表达基因，为后面火山图的点表上颜色
loc_up <- intersect(which(dat$qval < 0.05),which(dat$b >=1))
loc_down <- intersect(which(dat$qval < 0.05),which(dat$b <= (-1)))
logFC_cutoff <- 1 
dat$sig = as.factor(ifelse(dat$qval < 0.05 & abs(dat$b) > logFC_cutoff,
                           ifelse(dat$b > logFC_cutoff ,'UP','DOWN'),'NOT'))#添加sig这栏up,down,not的行数据
this_tile <- paste0('CR6 vs PR6 Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(dat[which(dat$sig == 'UP'),]),
                    '\nThe number of down gene is ',nrow(dat[which(dat$sig == 'DOWN'),])) #设置标题，体现阈值、上下调基因数量
#基本散点图
p2 <- ggplot(data=dat, aes(x=-log10(qval), y=b, color=sig)) +
  geom_point(alpha=0.4, size=1.5) +
  theme_bw(base_size=15)+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "up")+
  xlab("-log10(qval)") + ylab("log2FC") +
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('#1170aa','#a2a2a2','#d62728')) +
  geom_vline(xintercept=-log(0.05,10) ,linetype=4) +
  geom_hline(yintercept=c(-1,1) ,linetype=4 )
#组图
grid.arrange(p0,p1,p2,
             nrow=2,ncol=2)     %>%  ggsave("差异表达基因火山图.pdf",.,width=210,height=210, units="mm")

###根和叶差异表达基因韦恩图
leaf <- read.csv("CL6vsPL6_sleuth_significant_gene_level.csv",header=T)
leaf_up <- leaf[which((leaf$qval < 0.05) & (leaf$b >= 1)),]
leaf_down <- leaf[which(leaf$qval < 0.05) & (leaf$b <= (-1)),]
root <- read.csv("CR6vsPR6_sleuth_significant_gene_level.csv",header=T)
root_up <- root[which((root$qval < 0.05) & (root$b >= 1)),]
root_down <- root[which(root$qval < 0.05) & (root$b <= (-1)),]
write.csv(leaf_up,"leaf_up.csv")
write.csv(leaf_down,"leaf_down.csv")
write.csv(root_up,"root_up.csv")
write.csv(root_down,"root_down.csv")

#出图
venn_list1 <- list(leaf_up = leaf_up$target_id, leaf_down = leaf_down$target_id, root_up = root_up$target_id, root_down = root_down$target_id)
#V0vsV26N0_up  "#A6CEE3"
#V0vsV26N0_down  "#1F78B4"
#V26N0vsV26N6_up  "#B2DF8A"
#V26N0vsV26N6_down  "#33A02C"
#V0vsV26N6  "#FDBF6F"
#V0vsV26N6_down  "#FF7F00"
#V0vsV26N0_V26N0vsV26N6
venn.plot1<-venn.diagram(venn_list1,  filename = NULL,
                         fill = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"),  alpha = 0.50,
                         cat.col = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"),  cat.cex = 1,  cat.fontfamily = 'serif',
                         col = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"), cex = 1, fontfamily = 'serif')
pdf(file="leafCvsT_rootCvsT.pdf")
grid.draw(venn.plot1)
dev.off()
#提取交集元素
inter <- get.venn.partitions(venn_list1)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'leafCvsT_rootCvsT_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#可视化组蛋白修饰peak数量
setwd("D:/高粱/ChIP")
peak <- read.csv("统计.csv")
head(peak)
number <- unique(peak[1:63,c(2,3,4,9)])
number$peak_number <- as.numeric(number$peak_number)
number$type <- gsub("CL","",number$sample)
number$type <- gsub("CR","",number$type)
number$type <- gsub("PL","",number$type)
number$type <- gsub("PR","",number$type)
number$group <- paste0(number$组织,"_",number$对照或处理)
#绘制柱形图
p1 <- ggplot(number,mapping = aes(x=type,y=peak_number,fill=group))+
  geom_bar(stat='identity',position='dodge',show.legend = TRUE) +
  labs(y = 'Peak number') +
  scale_fill_manual(values=c("#6388b4", "#ffae34", "#ef6f6a", "#8cc2ca", "#55ad89", "#c3bc3f", "#bb7693", "#baa094", "#a9b5ae", "#767676")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 1.2),
        axis.ticks.x=element_blank(),
        legend.position="top") +
  ggtitle("Interaction_number") +
  guides(fill = guide_legend(title = 'Type')) +
  theme_bw()
#小提琴图展示FRiP值
number <- unique(peak[1:63,c(2,3,4,8)])
number$FRiP <- as.numeric(number$FRiP)*100
number$type <- gsub("CL","",number$sample)
number$type <- gsub("CR","",number$type)
number$type <- gsub("PL","",number$type)
number$type <- gsub("PR","",number$type)
number$group <- paste0(number$组织,"_",number$对照或处理)
p2 <- ggviolin(number, x = "type", y = "FRiP", fill = "type",
         bxp.errorbar=T, 
         outlier.shape = NA,
         #add = c("mean_sd"), error.plot = "crossbar") +
         add = c("boxplot")) +
  #scale_fill_manual(values=c("#aec7e8","c5b0d5","#ffbb78","#98df8a","#ff9896","#c49c94","#f7b6d2","lightgrey")) +
  #scale_fill_manual(values=c("#B2DF8A","#FDBF6F","#CAB2D6","#FB9A99","#1F78B4","#33A02C","#A6CEE3","#E31A1C","#FF7F00","#6A3D9A","#FFFF99","#B15928")) +
  scale_fill_manual(values=c(brewer.pal(8,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))) + #设置填充色
  ggtitle("FRiP") + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 1)) +
  xlab(NULL) + ylab("value") +
  theme (legend.position = 'none')
#小提琴图展示比对率值
peak <- read.csv("统计.csv")
number <- unique(peak[1:63,c(2,3,4,7)])
number$alignment <- as.numeric(number$alignment)*100
number$type <- gsub("CL","",number$sample)
number$type <- gsub("CR","",number$type)
number$type <- gsub("PL","",number$type)
number$type <- gsub("PR","",number$type)
number$group <- paste0(number$组织,"_",number$对照或处理)
p3 <- ggviolin(number, x = "type", y = "alignment", fill = "type",
               bxp.errorbar=T, 
               outlier.shape = NA,
               #add = c("mean_sd"), error.plot = "crossbar") +
               add = c("boxplot")) +
  #scale_fill_manual(values=c("#aec7e8","c5b0d5","#ffbb78","#98df8a","#ff9896","#c49c94","#f7b6d2","lightgrey")) +
  #scale_fill_manual(values=c("#B2DF8A","#FDBF6F","#CAB2D6","#FB9A99","#1F78B4","#33A02C","#A6CEE3","#E31A1C","#FF7F00","#6A3D9A","#FFFF99","#B15928")) +
  scale_fill_manual(values=c(brewer.pal(8,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))) + #设置填充色
  ggtitle("alignment") + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 1)) +
  xlab(NULL) + ylab("value") +
  theme (legend.position = 'none')

#组图
grid.arrange(p1,p0,p2,p3,
             nrow=2,ncol=2)     %>%  ggsave("质控值.pdf",.,width=210,height=210, units="mm")


#####根据基因表达量对基因进行分类
setwd("D:/高粱/RNA-seq/final_RNAseq")
tpm <- read.csv("all_tpm_expression.csv")
head(tpm)
row.names(tpm) <- tpm[,1]
CL <- tpm[,2:4]
CR <- tpm[,5:7]
PL <- tpm[,8:10]
PR <- tpm[,11:13]
CL$average <- rowMeans(CL)
CR$average <- rowMeans(CR)
PL$average <- rowMeans(PL)
PR$average <- rowMeans(PR)
CR$average <- rowMeans(CR)
PL$average <- rowMeans(PL)
PR$average <- rowMeans(PR)
#删除tpm小于0.1的基因
CL <- CL[which(rowMaxs(as.matrix(CL[,1:3])) > 0.1),]
CR <- CR[which(rowMaxs(as.matrix(CR[,1:3])) > 0.1),]
PL <- PL[which(rowMaxs(as.matrix(PL[,1:3])) > 0.1),]
PR <- PR[which(rowMaxs(as.matrix(PR[,1:3])) > 0.1),]
head(CL)
head(CR)
head(PL)
head(PR)
##将基因分成下25%、50%、75%和上25%
#定义函数，设置函数，四分位分组
Quantile<-function(x){
  ifelse(x>quantile(x,.75),"Q1",ifelse(x>quantile(x,.5),"Q2",ifelse(x>quantile(x,.25),"Q3","Q4")))
}
CL$cluster<-apply(CL, 2, Quantile)
CR$cluster<-apply(CR, 2, Quantile)
PL$cluster<-apply(PL, 2, Quantile)
PR$cluster<-apply(PR, 2, Quantile)
#加上基因的位置
gene <- read.table("geneR1.bed")
row.names(gene) <- gene[,4]
head(gene)
CL_pos <- merge(CL,gene,by="row.names",all.x=F)
CR_pos <- merge(CR,gene,by="row.names",all.x=F)
PL_pos <- merge(PL,gene,by="row.names",all.x=F)
PR_pos <- merge(PR,gene,by="row.names",all.x=F)
write.csv(CL_pos,"CL.csv")
write.csv(CR_pos,"CR.csv")
write.csv(PL_pos,"PL.csv")
write.csv(PR_pos,"PR.csv")


#####分析H3K27me3,A可视化peak宽度分布，B可视化落在基因区的数量，C圈图可视化全基因组H3K27me3的分布
setwd("D:/高粱/ChIP/peaks")
###合并重复的peaks
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
cat CRK27me3_rep1.peaks.bed CRK27me3_rep2.peaks.bed | sort -k1,1 -k2n,2 | bedtools merge -i - -c 4 -o collapse >CRK27me3_merged.bed
cat CLK27me3_rep1.peaks.bed CLK27me3_rep2.peaks.bed | sort -k1,1 -k2n,2 | bedtools merge -i - -c 4 -o collapse >CLK27me3_merged.bed
###合并相邻1kb的peak
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
bedtools merge -i CRK27me3_merged.bed -d 1000 >CRK27me3_mergedR1.bed
bedtools merge -i CLK27me3_merged.bed -d 1000 >CLK27me3_mergedR1.bed

##绘制长度分布直方图
#准备数据
cl <- read.table("CLK27me3_merged.bed")
cl$length <- cl$V3-cl$V2
a1 <- cl[which(cl$length >= 50000),]
a2 <- cl[which(cl$length < 50000 & cl$length >= 20000),]
a3 <- cl[which(cl$length < 20000 & cl$length >= 10000),]
a4 <- cl[which(cl$length < 10000 & cl$length >= 5000),]
a5 <- cl[which(cl$length < 5000 & cl$length >= 2000),]
a6 <- cl[which(cl$length < 2000 & cl$length >= 1000),]
a7 <- cl[which(cl$length < 1000 & cl$length >= 0),]
a1$group <- rep(">=50",nrow(a1))
a2$group <- rep("20-50",nrow(a2))
a3$group <- rep("10-20",nrow(a3))
a4$group <- rep("5-10",nrow(a4))
a5$group <- rep("2-5",nrow(a5))
a6$group <- rep("1-2",nrow(a6))
a7$group <- rep("0-1",nrow(a7))
cl <- rbind(a1,a2,a3,a4,a5,a6,a7)
cl$type <- rep("CL",nrow(cl))
head(cl)
CR <- read.table("CRK27me3_merged.bed")
CR$length <- CR$V3-CR$V2
a1 <- CR[which(CR$length >= 50000),]
a2 <- CR[which(CR$length < 50000 & CR$length >= 20000),]
a3 <- CR[which(CR$length < 20000 & CR$length >= 10000),]
a4 <- CR[which(CR$length < 10000 & CR$length >= 5000),]
a5 <- CR[which(CR$length < 5000 & CR$length >= 2000),]
a6 <- CR[which(CR$length < 2000 & CR$length >= 1000),]
a7 <- CR[which(CR$length < 1000 & CR$length >= 0),]
a1$group <- rep(">=50",nrow(a1))
a2$group <- rep("20-50",nrow(a2))
a3$group <- rep("10-20",nrow(a3))
a4$group <- rep("5-10",nrow(a4))
a5$group <- rep("2-5",nrow(a5))
a6$group <- rep("1-2",nrow(a6))
a7$group <- rep("0-1",nrow(a7))
CR <- rbind(a1,a2,a3,a4,a5,a6,a7)
CR$type <- rep("CR",nrow(CR))
head(CR)
len1 <- rbind(cl,CR)
len1$cluster <- rep("MG0",nrow(len1))
CL2 <- read.table("CLK27me3_mergedR1.bed")
CL2$length <- CL2$V3-CL2$V2
a1 <- CL2[which(CL2$length >= 50000),]
a2 <- CL2[which(CL2$length < 50000 & CL2$length >= 20000),]
a3 <- CL2[which(CL2$length < 20000 & CL2$length >= 10000),]
a4 <- CL2[which(CL2$length < 10000 & CL2$length >= 5000),]
a5 <- CL2[which(CL2$length < 5000 & CL2$length >= 2000),]
a6 <- CL2[which(CL2$length < 2000 & CL2$length >= 1000),]
a7 <- CL2[which(CL2$length < 1000 & CL2$length >= 0),]
a1$group <- rep(">=50",nrow(a1))
a2$group <- rep("20-50",nrow(a2))
a3$group <- rep("10-20",nrow(a3))
a4$group <- rep("5-10",nrow(a4))
a5$group <- rep("2-5",nrow(a5))
a6$group <- rep("1-2",nrow(a6))
a7$group <- rep("0-1",nrow(a7))
CL2 <- rbind(a1,a2,a3,a4,a5,a6,a7)
CL2$type <- rep("CL",nrow(CL2))
head(CL2)
CR2 <- read.table("CRK27me3_mergedR1.bed")
CR2$length <- CR2$V3-CR2$V2
a1 <- CR2[which(CR2$length >= 50000),]
a2 <- CR2[which(CR2$length < 50000 & CR2$length >= 20000),]
a3 <- CR2[which(CR2$length < 20000 & CR2$length >= 10000),]
a4 <- CR2[which(CR2$length < 10000 & CR2$length >= 5000),]
a5 <- CR2[which(CR2$length < 5000 & CR2$length >= 2000),]
a6 <- CR2[which(CR2$length < 2000 & CR2$length >= 1000),]
a7 <- CR2[which(CR2$length < 1000 & CR2$length >= 0),]
a1$group <- rep(">=50",nrow(a1))
a2$group <- rep("20-50",nrow(a2))
a3$group <- rep("10-20",nrow(a3))
a4$group <- rep("5-10",nrow(a4))
a5$group <- rep("2-5",nrow(a5))
a6$group <- rep("1-2",nrow(a6))
a7$group <- rep("0-1",nrow(a7))
CR2 <- rbind(a1,a2,a3,a4,a5,a6,a7)
CR2$type <- rep("CR",nrow(CR2))
head(CR2)
len2 <- rbind(CL2,CR2)
len2$cluster <- rep("MG1000",nrow(len2))
head(len2)
len1 <- len1[,-4]
len <- rbind(len1,len2)
len$length <- as.numeric(len$length)
head(len)
#统计数量
total_number <- aggregate(len$V1,list(len$group,len$type,len$cluster),length)
total_number$Group.1 <- factor(total_number$Group.1, levels=c("0-1","1-2","2-5","5-10",
                                                              "10-20","20-50",">=50"), ordered=TRUE)
head(total_number)
#绘图
pdf("2AR2.pdf")
ggplot(total_number,mapping = aes(x=Group.1,y=x,fill=Group.3))+
  geom_bar(stat='identity',position='dodge',show.legend = TRUE) +
  labs(x = 'Length (kb) of H3K27me3 regions', y = 'Number of regions') +
  scale_fill_manual(values=brewer.pal(8,"Dark2")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1.2),
        axis.ticks.x=element_blank(),
        legend.position="top") +
  ggtitle("H3K27me3") +
  guides(fill = guide_legend(title = 'Type')) +
  geom_text(aes(label = x), size = 3, vjust =-0.2, position = position_dodge(0.9)) +
  theme_bw() +
  facet_grid(. ~Group.2)
dev.off()
#没用下面这个图
#hist(log10(cl$length),main="CL",col="dark red",xlab="Length (bp) of H3K27me3 regions (log10)")
ggplot(data = len, mapping = aes(x = length/1000,fill = cluster)) +
  scale_fill_manual(values = brewer.pal(9,"Set1")) +
  #scale_fill_manual(values = c("#FF7F00","#6A3D9A")) +
  geom_histogram( bins = 60, position="dodge") + facet_grid(. ~ type) +
  geom_vline(aes(xintercept=median( length/1000, na.rm=T)),
             color="red", linetype="dashed", size=0.75) +
  theme_bw() +
  xlab("Length (kb) of H3K27me3 regions") + ylab("Number of regions") +
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  #coord_cartesian(xlim = c(0,60))
  xlim(c(0,60))
#标出中位数
cl <- read.table("CLK27me3_merged.bed")
cl$length <- cl$V3-cl$V2
CR <- read.table("CRK27me3_merged.bed")
CR$length <- CR$V3-CR$V2
CL2 <- read.table("CLK27me3_mergedR1.bed")
CL2$length <- CL2$V3-CL2$V2
CR2 <- read.table("CRK27me3_mergedR1.bed")
CR2$length <- CR2$V3-CR2$V2
median(cl$length)  #1357
median(CR$length)  #1892.5
median(CL2$length)  #1497
median(CR2$length)  #2243

###### B 看落在基因上k27me3的比列
#####求CL CR落在基因上peak的比列
####MG0
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
sort -k1,1 -k2,2n CLK27me3_merged.bed | bedtools closest -D ref -t all -mdb all -a - -b \
/public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed  >CLK27me3_peak_gene.txt
sort -k1,1 -k2,2n CRK27me3_merged.bed | bedtools closest -D ref -t all -mdb all -a - -b \
/public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed  >CRK27me3_peak_gene.txt
####读入数据
####CR
cr27 <- read.table("CRK27me3_peak_gene.txt")
cr27$type <- rep("CR",nrow(cr27))
###不能去重复
#cr27$fz<-abs(cr27[,11])
#cr27<-cr27[order(cr27[,1],cr27[,4],cr27[,13]),]
#cr27<-cr27[!duplicated(cr27[,4]),]
#cr27<-cr27[,-13]
#cr27<-cr27[order(cr27[,1],cr27[,2],cr27[,11]),]
####CL
cl27 <- read.table("CLK27me3_peak_gene.txt")
cl27$type <- rep("CL",nrow(cl27))
###不能去重复
#cl27$fz<-abs(cl27[,11])
#cl27<-cl27[order(cl27[,1],cl27[,4],cl27[,13]),]
#cl27<-cl27[!duplicated(cl27[,4]),]
#cl27<-cl27[,-13]
#cl27<-cl27[order(cl27[,1],cl27[,2],cl27[,11]),]
#合并
k27 <-rbind(cr27,cl27)
head(k27)
ge <- k27[which(k27$V11 == 0),]
ne <- k27[which(k27$V11 != 0),]
ne$V11 <- rep("no_genes",nrow(ne))
###统计同一个修饰同时标记多少个基因
ge_number <- aggregate(ge$V1,list(ge$V4,ge$V11,ge$type),length)
head(ge_number)
##统计含含不同基因的peak的数量
ge_number <- aggregate(ge_number$ Group.1,list(ge_number$x,ge_number$Group.2,ge_number$Group.3),length)
head(ge_number)
###统计没有修饰的基因数,要去重
head(ne)
ne<-ne[!duplicated(ne[,4]),]
ne_number <- aggregate(ne$V1,list(ne$V11,ne$type),length)
head(ne_number)
##合并
ge_number <- ge_number[,-2]
colnames(ge_number) <- colnames(ne_number)
to_number <- rbind(ge_number,ne_number)
to_number$Group.1 <- gsub("no_genes",0,to_number$Group.1)
to_number$Group.1 <- as.numeric(to_number$Group.1)
head(to_number)
#分为0,1,2,3,4,5-10,10-20,>20
low <- to_number[which(to_number$Group.1 <5),]
a1 <- to_number[which(to_number$Group.1 >= 5 & to_number$Group.1 < 10),]
a1$Group.1 <- rep("5-10",nrow(a1))
a1 <- aggregate(a1$x,list(a1$Group.1,a1$Group.2),sum)
a2 <- to_number[which(to_number$Group.1 >= 10 & to_number$Group.1 <20),]
a2$Group.1 <- rep("10-20",nrow(a2))
a2 <- aggregate(a2$x,list(a2$Group.1,a2$Group.2),sum)
a3 <- to_number[which(to_number$Group.1 >= 20),]
a3$Group.1 <- rep(">=20",nrow(a3))
a3 <- aggregate(a3$x,list(a3$Group.1,a3$Group.2),sum)
to1 <- rbind(low,a1,a2,a3)
head(to1)
to1$type <- rep("MG0",nrow(to1))


####MG1000
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
bedtools merge -i CRK27me3_merged.bed -d 1000 -c 4 -o collapse >CRK27me3_mergedR1.bed
bedtools merge -i CLK27me3_merged.bed -d 1000 -c 4 -o collapse >CLK27me3_mergedR1.bed
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
sort -k1,1 -k2,2n CLK27me3_mergedR1.bed | bedtools closest -D ref -t all -mdb all -a - -b \
/public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed  >CLK27me3_peak_geneR1.txt
sort -k1,1 -k2,2n CRK27me3_mergedR1.bed | bedtools closest -D ref -t all -mdb all -a - -b \
/public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed  >CRK27me3_peak_geneR1.txt
####读入数据
####CR
cr27 <- read.table("CRK27me3_peak_geneR1.txt")
cr27$type <- rep("CR",nrow(cr27))
###不能去重复
#cr27$fz<-abs(cr27[,11])
#cr27<-cr27[order(cr27[,1],cr27[,4],cr27[,13]),]
#cr27<-cr27[!duplicated(cr27[,4]),]
#cr27<-cr27[,-13]
#cr27<-cr27[order(cr27[,1],cr27[,2],cr27[,11]),]
####CL
cl27 <- read.table("CLK27me3_peak_geneR1.txt")
cl27$type <- rep("CL",nrow(cl27))
###不能去重复
#cl27$fz<-abs(cl27[,11])
#cl27<-cl27[order(cl27[,1],cl27[,4],cl27[,13]),]
#cl27<-cl27[!duplicated(cl27[,4]),]
#cl27<-cl27[,-13]
#cl27<-cl27[order(cl27[,1],cl27[,2],cl27[,11]),]
#合并
k27 <-rbind(cr27,cl27)
head(k27)
ge <- k27[which(k27$V11 == 0),]
ne <- k27[which(k27$V11 != 0),]
ne$V11 <- rep("no_genes",nrow(ne))
###统计同一个修饰同时标记多少个基因
ge_number <- aggregate(ge$V1,list(ge$V4,ge$V11,ge$type),length)
head(ge_number)
##统计含含不同基因的peak的数量
ge_number <- aggregate(ge_number$ Group.1,list(ge_number$x,ge_number$Group.2,ge_number$Group.3),length)
head(ge_number)
###统计没有修饰的基因数,要去重
head(ne)
ne<-ne[!duplicated(ne[,4]),]
ne_number <- aggregate(ne$V1,list(ne$V11,ne$type),length)
head(ne_number)
##合并
ge_number <- ge_number[,-2]
colnames(ge_number) <- colnames(ne_number)
to_number <- rbind(ge_number,ne_number)
to_number$Group.1 <- gsub("no_genes",0,to_number$Group.1)
to_number$Group.1 <- as.numeric(to_number$Group.1)
head(to_number)
#分为0,1,2,3,4,5-10,10-20,>20
low <- to_number[which(to_number$Group.1 <5),]
a1 <- to_number[which(to_number$Group.1 >= 5 & to_number$Group.1 < 10),]
a1$Group.1 <- rep("5-10",nrow(a1))
a1 <- aggregate(a1$x,list(a1$Group.1,a1$Group.2),sum)
a2 <- to_number[which(to_number$Group.1 >= 10 & to_number$Group.1 <20),]
a2$Group.1 <- rep("10-20",nrow(a2))
a2 <- aggregate(a2$x,list(a2$Group.1,a2$Group.2),sum)
a3 <- to_number[which(to_number$Group.1 >= 20),]
a3$Group.1 <- rep(">=20",nrow(a3))
a3 <- aggregate(a3$x,list(a3$Group.1,a3$Group.2),sum)
to2 <- rbind(low,a1,a2,a3)
head(to2)
to2$type <- rep("MG1000",nrow(to2))
#合并
to <- rbind(to1,to2)
#画图
to$Group.1 <- factor(to$Group.1, levels=c("0","1","2","3","4","5-10","10-20",">=20"), ordered=TRUE)
pdf("2c.pdf")
ggplot(to,mapping = aes(x=Group.1,y=x,fill=type))+
  geom_bar(stat='identity',position='dodge',show.legend = TRUE) +
  labs(x = 'Number of genes in H3K27me3 region', y = 'Number of regions') +
  #scale_fill_manual(values=brewer.pal(8,"Paired")) +
  scale_fill_manual(values=c("#1F78B4","#A6CEE3")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1.2),
        axis.ticks.x=element_blank(),
        legend.position="top") +
  ggtitle("H3K27me3") +
  guides(fill = guide_legend(title = 'Type')) +
  geom_text(aes(label = x), size = 3, vjust =-0.2, position = position_dodge(0.9)) +
  theme_bw() +
  facet_grid(. ~Group.2)
dev.off()

########Circos,不做，换种方式展示
#绘制H3K27me3的circos图
##Ciros
library(circlize)
library(stringr)
library(grid)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
#circos.clear()
setwd("D:/高粱/ChIP/CIRCOS")
seq_stat<-read.csv("3chromosome_length.csv",header=T,row.names=1)
seq_stat$seq_start <- as.numeric(seq_stat$seq_start)
seq_stat$seq_end <- as.numeric(seq_stat$seq_end)
circle_size = unit(1, "snpc")
circos.par(gap.degree = 2)
circos.genomicInitialize(seq_stat, plotType = 'axis')
circos.track(
  ylim = c(0, 0.5), track.height = 0.08, bg.border = NA, bg.col =rep(c("#CAB2D6", "#B2DF8A", "#FDBF6F"), 7),
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    seq_ID = CELL_META$sector.index
    circos.text(mean(xlim), mean(ylim), seq_ID, cex = 0.7, col = 'black', facing = 'inside', niceFacing = FALSE)
  } )

#基因密度 200000一个滑窗
#计算基因密度，将参考基因组拆成10000个滑窗
bedops --chop 100000 --stagger 100000 -x  /public/home/chaohe/sorghum/RNA-seq/db/Sorghum_bicolor.genome_table.txt > gene_100kb.bed
bedmap --echo --count --delim '\t'  gene_100kb.bed  /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed > count_gene_100kb.txt
density<-read.table("count_gene_100kb.txt",header=T)
color_assign <- colorRamp2(breaks = c(1, 10, 21), 
                           col = c('#FEE8C8', '#FC8D59', '#D7301F'))
circos.genomicTrackPlotRegion(
  density, track.height = 0.12, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )
####对H3K27me3的peak进行定量
#ChIP量化
#featurecount计算count值-安装subread
#conda install -c bioconda subread
#https://jaist.dl.sourceforge.net/project/subread/subread-2.0.2/subread-2.0.2-Linux-x86_64.tar.gz
#tar -zxvf subread-2.0.2-Linux-x86_64.tar.gz
#cd subread-2.0.2-Linux-x86_64
#cd bin
#ls
##产生saf文件
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
cut -f 1 b.txt | uniq | while read i; do
awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' "$i".peaks.bed >"$i".saf
done
#featurecount计算count值
cut -f 1,2 b.txt | uniq | while read i j; do
/public/home/chaohe/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -P -B -C -T 4 \
-a "$i".saf \
-F SAF \
-o "$i"_counts_subread.txt \
/public/home/chaohe/sorghum/chip/align/"$j".final.bam
done
#narrow
cd /public/home/chaohe/sorghum/chip/macs2/narrow-spp-0.00001
cut -f 1 a.txt | uniq | while read i; do
awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' "$i".peaks.bed >"$i".saf
done
cut -f 1,2 a.txt | uniq | while read i j; do
/public/home/chaohe/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -P -B -C -T 4 \
-a "$i".saf \
-F SAF \
-o "$i"_counts_subread.txt \
/public/home/chaohe/sorghum/chip/align/"$j".final.bam
done
#TPM量化每个peak
#定义公式
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
#计算TPM值
name=c("CLK27me3_rep1","CLK27me3_rep2","CRK27me3_rep1","CRK27me3_rep2")
for (i in name){
  pea<-read.table(paste0(i,"_counts_subread.txt"),skip=1)
  head(pea)
  colnames(pea)<-pea[1,]
  pea<-pea[-1,]
  count<-as.data.frame(as.numeric(pea[,7]))
  #Length<-as.data.frame(as.numeric(pea[,6]))
  Length<-as.data.frame(as.numeric(pea$Length))
  tpms <- apply(count, 2, function(x) tpm(x, Length))
  peak_tpms<-cbind(pea,tpms)
  colnames(peak_tpms)<-c("Geneid","Chr","Start","End","Strand","Length",
                         "counts","TPM")
  peak_tpms <- peak_tpms[,c(2,3,4,8)]
  colnames(peak_tpms)<-c("seq_ID","seq_start","seq_end","value1")
  write.csv(peak_tpms,paste0(i,"peak_tpms.csv"),row.names=F)
}
###CRK27me31
CR1 <- read.csv("CRK27me3_rep2peak_tpms.csv")
CR1$value1 <- log10(CR1$value1+1)
head(CR1)
#genescore$value1 <- log2(genescore$value1+0.01)
circos.genomicTrack(
  CL1,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#1B9E77',border = '#1B9E77')
  } )
#画图
setwd("/public/home/chaohe/sorghum/chip/align/RPKM")
seq_stat<-read.csv("3chromosome_length.csv",header=T,row.names=1)
seq_stat$seq_start <- as.numeric(seq_stat$seq_start)
seq_stat$seq_end <- as.numeric(seq_stat$seq_end)
circle_size = unit(1, "snpc")
circos.par(gap.degree = 2)
circos.genomicInitialize(seq_stat, plotType = 'axis')
circos.track(
  ylim = c(0, 0.5), track.height = 0.08, bg.border = NA, bg.col =rep(c("#CAB2D6", "#B2DF8A", "#FDBF6F"), 7),
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    seq_ID = CELL_META$sector.index
    circos.text(mean(xlim), mean(ylim), seq_ID, cex = 0.7, col = 'black', facing = 'inside', niceFacing = FALSE)
  } )
density<-read.table("count_gene_100kb.txt",header=T)
color_assign <- colorRamp2(breaks = c(1, 10, 21), 
                           col = c('#FEE8C8', '#FC8D59', '#D7301F'))
circos.genomicTrackPlotRegion(
  density, track.height = 0.12, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )
cr1 <- read.csv("CRK27me3_rep1.bedGraph")
colnames(cr1) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cr1,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#D95F02',border = '#D95F02')
  } )
cr2 <- read.csv("CRK27me3_rep1.bedGraph")
colnames(cr2) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cr2,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#7570B3',border = '#7570B3')
  } )
cl1 <- read.csv("CLK27me3_rep1.bedGraph")
colnames(cl1) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cl1,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#66A61E',border = '#E78AC3')
  } )
cl2 <- read.csv("CLK27me3_rep1.bedGraph")
colnames(cl2) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cl2,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#E7298A',border = '#FFD92F')
  } )
bsub  -J circos -n 8 -o %J.3.out -e %J.3.err  -q smp -R "rusage[mem=300GB]" "Rscript RR.r"

###利用Sushi做
#BiocManager::install("Sushi")
library('Sushi')
Sushi_data = data(package = 'Sushi')
data(list = Sushi_data$results[,3])
###用bw文件做，bigwig转bedgraph
##H3K27me3
cd /public/home/chaohe/sorghum/chip/align/RPKM
bigWigToBedGraph  CRK27me3_rep1_rpkm.bw  CRK27me3_rep1.bedGraph 
bigWigToBedGraph  CRK27me3_rep2_rpkm.bw  CRK27me3_rep2.bedGraph 
bigWigToBedGraph  CLK27me3_rep1_rpkm.bw  CLK27me3_rep1.bedGraph 
bigWigToBedGraph  CLK27me3_rep2_rpkm.bw  CLK27me3_rep2.bedGraph 
sed -i "/super/d" CRK27me3_rep1.bedGraph
sed -i "/super/d" CRK27me3_rep2.bedGraph
sed -i "/super/d" CLK27me3_rep1.bedGraph 
sed -i "/super/d" CLK27me3_rep2.bedGraph 

#合并2kb
bedtools merge -i CRK27me3_rep1.bedGraph -d 2000 -c 4 -o sum >CRK27me3_rep1.bed
CRK27me3_rep1.bedGraph <- read.table("CRK27me3_rep1.bedGraph")
plotBedgraph(CRK27me3_rep1.bedGraph,transparency=.50,overlay=TRUE,rescaleoverlay=TRUE)

###利用pyGenomeTracks画
conda activate py3
#conda install -c bioconda pygenometracks
module load pyGenomeTracks/3.5
cd /public/home/chaohe/sorghum/chip/align/RPKM
make_tracks_file --trackFiles C*K27me3_rep1_rpkm.bw C*K27me3_rep2_rpkm.bw /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed -o H3K27me3.ini
bsub  -J 242-1 -n 8 -o %J.3.out -e %J.3.err  -q smp -R "rusage[mem=80GB]" "pyGenomeTracks --tracks H3K27me3R1.ini --BED region.bed -o image.pdf"
pyGenomeTracks --tracks H3K27me3R1.ini --region 2:1-12,874,906 -o 2.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 5:373,594-8,960,825 -o 5-1.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 5:58,788,817-71,854,669 -o 5-2.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 8:53,472,765-61,966,666 -o 8.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 2:11,031,103-11,192,785 -o 2-2.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 5:4,688,193-5,143,220 -o 2-5-1.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 5:68,664,052-68,937,714 -o 2-5-2.pdf 
pyGenomeTracks --tracks H3K27me3R1.ini --region 8:61,454,778-61,683,380 -o 2-8.pdf 

#Figure 3B
module load pyGenomeTracks/3.5
cd /public/home/chaohe/sorghum/chip/align/RPKM
make_tracks_file --trackFiles CLK27me3_rep1_rpkm.bw CLK27me3_rep2_rpkm.bw CLK36me3_rep1_rpkm.bw CLK36me3_rep2_rpkm.bw /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed -o K27me3-K26me3.ini
pyGenomeTracks --tracks K27me3-K26me3.ini --region 2:10,488,076-12,686,713 -o 3-2.pdf 
pyGenomeTracks --tracks K27me3-K26me3.ini --region 5:3,065,051-5,836,587 -o 3-5-1.pdf 
pyGenomeTracks --tracks K27me3-K26me3.ini --region 5:66,879,551-69,451,087 -o 3-5-2.pdf 
pyGenomeTracks --tracks K27me3-K26me3.ini --region 8:57,682,802-59,900,710 -o 3-8.pdf 
make_tracks_file --trackFiles CLKH2AZ_rep1_rpkm.bw CLKH2AZ_rep2_rpkm.bw CLK36me3_rep1_rpkm.bw CLK36me3_rep2_rpkm.bw /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed -o KH2AZ-K26me3.ini
pyGenomeTracks --tracks KH2AZ-K26me3.ini --region 4:65,227,990-65,393,504 -o 3-4.pdf 
pyGenomeTracks --tracks KH2AZ-K26me3.ini --region 6:53,765,402-54,060,844 -o 3-6.pdf 
pyGenomeTracks --tracks KH2AZ-K26me3.ini --region 9:58,032,216-58,318,687 -o 3-9.pdf 

#Figure S3B
module load pyGenomeTracks/3.5
cd /public/home/chaohe/sorghum/chip/align/RPKM
#make_tracks_file --trackFiles CRK27me3_rep1_rpkm.bw CRK27me3_rep2_rpkm.bw CRK36me3_rep1_rpkm.bw CRK36me3_rep2_rpkm.bw /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed -o CRK27me3-K26me3.ini
pyGenomeTracks --tracks CRK27me3-K26me3.ini --region 2:10,488,076-12,686,713 -o CR3-2.pdf 
pyGenomeTracks --tracks CRK27me3-K26me3.ini --region 5:3,065,051-5,836,587 -o 3-CR5-1.pdf 
pyGenomeTracks --tracks CRK27me3-K26me3.ini --region 5:66,879,551-69,451,087 -o CR3-5-2.pdf 
pyGenomeTracks --tracks CRK27me3-K26me3.ini --region 8:57,682,802-59,900,710 -o CR3-8.pdf 
#make_tracks_file --trackFiles CRKH2AZ_rep1_rpkm.bw CRKH2AZ_rep2_rpkm.bw CRK36me3_rep1_rpkm.bw CRK36me3_rep2_rpkm.bw /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed -o CRKH2AZ-K26me3.ini
pyGenomeTracks --tracks CRKH2AZ-K26me3.ini --region 4:65,227,990-65,393,504 -o CR3-4.pdf 
pyGenomeTracks --tracks CRKH2AZ-K26me3.ini --region 6:53,765,402-54,060,844 -o CR3-6.pdf 
pyGenomeTracks --tracks CRKH2AZ-K26me3.ini --region 9:58,032,216-58,318,687 -o CR3-9.pdf 

#20230224补充，可视化FigS2中用到的水稻、大麦和玉米的全基因组H3K27me3的分布
module load pyGenomeTracks/3.5
cd /public/home/chaohe/sorghum/other_species
#rice
make_tracks_file --trackFiles Nip_young_leaf_H3K27me3_Rep1.rpkm.bw Nip_young_leaf_H3K27me3_Rep2.rpkm.bw rice.bed -o rice.ini
pyGenomeTracks --tracks rice.ini --BED rice.bed -o rice.pdf
#barley
make_tracks_file --trackFiles ERX654615.final.rpkm.bw ERX654616.final.rpkm.bw ERX654617.final.rpkm.bw barley.bed -o barley.ini
pyGenomeTracks --tracks barley.ini --BED barley.bed -o barley.pdf 
#maize
make_tracks_file --trackFiles SRX1895365.final.rpkm.bw SRX1895366.final.rpkm.bw SRX1895367.final.rpkm.bw -o maize.ini
pyGenomeTracks --tracks maize.ini --BED maize.bed -o maize.pdf 
#####更新，可视化FigS2r1中水稻和拟南芥的全基因组H3K27me3的分布，用peak来做，后面有更新
cd /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001
#rice
perl /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001/writeChromInfoBed.pl /public/home/zhruan/ATAC-rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
mv chromInfo.bed  rice_genome_table.txt
#bedops --chop 100000 --stagger 100000 -x  rice_genome_table.txt > rice_gene_100kb.bed
#bedmap --echo --count --delim '\t'  rice_gene_100kb.bed  SRX7426637.broadPeak > count_rice_rep1_1ookb.txt
#bedmap --echo --count --delim '\t'  rice_gene_100kb.bed  SRX7426638.broadPeak > count_rice_rep2_1ookb.txt
#ara
perl writeChromInfoBed.pl /public/home/zhruan/ATAC-aradopisis2/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
mv chromInfo.bed  ara_genome_table.txt
#bedops --chop 100000 --stagger 100000 -x  ara_genome_table.txt > ara_gene_100kb.bed
#bedmap --echo --count --delim '\t'  ara_gene_100kb.bed  SRX7734769.broadPeak > count_ara_rep1_100kb.txt
#bedmap --echo --count --delim '\t'  ara_gene_100kb.bed  SRX7734770.broadPeak > count_ara_rep2_100kb.txt
#rice
####水稻
cd /public/home/chaohe/sorghum/other_species
###提取落在peak的reads
#bw转bedgraph
module load BEDTools/2.27
bigWigToBedGraph  Nip_young_leaf_H3K27me3_Rep1.rpkm.bw  rice_rep1.bedGraph 
bigWigToBedGraph  Nip_young_leaf_H3K27me3_Rep2.rpkm.bw  rice_rep2.bedGraph 
awk '{print $1"\t"$3}' /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001/rice_genome_table.txt >rice_genome_table.txt
bedtools  intersect  -a rice_rep1.bedGraph -b SRX7426637.broadPeak | sort -k1,1 -k2,2n > rice_rep1_reads.bed
bedtools  intersect  -a rice_rep2.bedGraph -b SRX7426638.broadPeak | sort -k1,1 -k2,2n > rice_rep2_reads.bed
bedGraphToBigWig rice_rep1_reads.bed rice_genome_table.txt rice_rep1_reads.bw
bedGraphToBigWig rice_rep2_reads.bed rice_genome_table.txt rice_rep2_reads.bw
####拟南芥
bigWigToBedGraph  SRX7734769.final.rpkm.bw  ara_rep1.bedGraph 
bigWigToBedGraph  SRX7734770.final.rpkm.bw  ara_rep2.bedGraph 
awk '{print $1"\t"$3}' /public/home/chaohe/sorghum/chip/macs2/broad_rep-spp-0.00001/ara_genome_table.txt >ara_genome_table.txt
bedtools  intersect  -a ara_rep1.bedGraph -b SRX7734769.broadPeak | sort -k1,1 -k2,2n > ara_rep1_reads.bed
bedtools  intersect  -a ara_rep2.bedGraph -b SRX7734770.broadPeak | sort -k1,1 -k2,2n > ara_rep2_reads.bed
perl -p -i -e 's/Chr//g' ara_rep1_reads.bed
perl -p -i -e 's/Chr//g' ara_rep2_reads.bed
bedGraphToBigWig ara_rep1_reads.bed ara_genome_table.txt ara_rep1_reads.bw
bedGraphToBigWig ara_rep2_reads.bed ara_genome_table.txt ara_rep2_reads.bw
#绘图
conda activate py36
module load pyGenomeTracks/3.5
awk '{print $1"\t1\t"$2}' rice_genome_table.txt >rice_gene.bed
make_tracks_file --trackFiles rice_rep1_reads.bw rice_rep2_reads.bw rice_gene.bed -o rice.ini
pyGenomeTracks --tracks rice.ini --BED rice_gene.bed -o riceR1.pdf
#ara
awk '{print $1"\t1\t"$2}' ara_genome_table.txt >ara_gene.bed
make_tracks_file --trackFiles ara_rep1_reads.bw ara_rep2_reads.bw ara_genome_table.txt -o ara.ini
pyGenomeTracks --tracks ara.ini --BED ara_genome_table.txt -o ara.pdf
#####################更新，只展示长达大于10 kb的信号
#水稻 提取大于10kb的peaks
#conda install -c bioconda ucsc-bigwigtobedgraph
module load OpenSSL/1.0.2h-foss-2016b
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
awk '{if(($3-$2) > 10000) print $0}' SRX7426637.broadPeak | bedtools  intersect  -a rice_rep1.bedGraph -b - | sort -k1,1 -k2,2n > riceR1_rep1_reads.bed
awk '{if(($3-$2) > 10000) print $0}' SRX7426638.broadPeak | bedtools  intersect  -a rice_rep1.bedGraph -b - | sort -k1,1 -k2,2n > riceR1_rep2_reads.bed
bedGraphToBigWig riceR1_rep1_reads.bed rice_genome_table.txt riceR1_rep1_reads.bw
bedGraphToBigWig riceR1_rep2_reads.bed rice_genome_table.txt riceR1_rep2_reads.bw
#拟南芥 提取大于10kb的peaks
awk '{if(($3-$2) > 10000) print $0}' SRX7734769.broadPeak | bedtools  intersect  -a ara_rep1.bedGraph -b - | sort -k1,1 -k2,2n > araR1_rep1_reads.bed
awk '{if(($3-$2) > 10000) print $0}' SRX7734770.broadPeak | bedtools  intersect  -a ara_rep2.bedGraph -b - | sort -k1,1 -k2,2n > araR1_rep2_reads.bed
perl -p -i -e 's/Chr//g' araR1_rep1_reads.bed
perl -p -i -e 's/Chr//g' araR1_rep2_reads.bed
bedGraphToBigWig araR1_rep1_reads.bed ara_genome_table.txt araR1_rep1_reads.bw
bedGraphToBigWig araR1_rep2_reads.bed ara_genome_table.txt araR1_rep2_reads.bw
#绘图
conda activate py36
module load pyGenomeTracks/3.5
pyGenomeTracks --tracks riceR1.ini --BED rice_gene.bed -o rice_10kb.pdf
#ara
pyGenomeTracks --tracks araR1.ini --BED ara_gene.bed -o ara_10kb.pdf

###计算水稻这套数据的重复性
multiBigwigSummary bins -b Nip_young_leaf_H3K27me3_Rep1.rpkm.bw Nip_young_leaf_H3K27me3_Rep2.rpkm.bw -o results.npz 
plotCorrelation \
-in results.npz \
--corMethod pearson --skipZeros \
--plotTitle "Person Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o heatmap_PersonCorr_readCounts.png   \
--outFileCorMatrix SpearmanCorr_readCounts.tab
###重复性不好，重新找数据,这套数据重复性可以
multiBigwigSummary bins -b SRX7426661.final.rpkm.bw SRX7426662.final.rpkm.bw SRX7426671.final.rpkm.bw SRX7426672.final.rpkm.bw -o results.npz 
plotCorrelation \
-in results.npz \
--corMethod pearson --skipZeros \
--plotTitle "Person Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o heatmap_PersonCorr_readCounts.png   \
--outFileCorMatrix SpearmanCorr_readCounts.tab
##################水稻选择root这套数据SRX7426661  SRX7426662来自SRP238435重新分析
cd /public/home/chaohe/sorghum/other_species
conda install -c bioconda ucsc-bigwigtobedgraph
module load BEDTools/2.27
bigWigToBedGraph  SRX7426661.final.rpkm.bw  rice_rep1.bedGraph 
bigWigToBedGraph  SRX7426662.final.rpkm.bw  rice_rep2.bedGraph 
awk '{if(($3-$2) > 10000) print $0}' SRX7426661.broadPeak | bedtools  intersect  -a rice_rep1.bedGraph -b - | sort -k1,1 -k2,2n > riceR1_rep1_reads.bed
awk '{if(($3-$2) > 10000) print $0}' SRX7426662.broadPeak | bedtools  intersect  -a rice_rep1.bedGraph -b - | sort -k1,1 -k2,2n > riceR1_rep2_reads.bed
module load OpenSSL/1.0.2h-foss-2016b
bedGraphToBigWig riceR1_rep1_reads.bed rice_genome_table.txt riceR1_rep1_reads.bw
bedGraphToBigWig riceR1_rep2_reads.bed rice_genome_table.txt riceR1_rep2_reads.bw
#绘图
conda activate py36
module load pyGenomeTracks/3.5
pyGenomeTracks --tracks riceR1.ini --BED rice_gene.bed -o rice_10kb.pdf
####


#Fgiure 3c,表达量热图
setwd("D:/高粱/ChIP/文章图片/Figure 3/表达量热图")
rna <- read.csv("D:/高粱/RNA-seq/final_RNAseq/all_tpm_expression.csv",row.names=1)
head(rna)
cl <- rna[,1:3]
cl$average <- rowMeans(cl)
head(cl)
#kh2az
h2az <- read.table("Heatmap1sortedRegions_CL_CLKH2AZ_TSS.bed")
h2az$order <- row.names(h2az)
rownames(h2az) <- h2az[,4]
head(h2az)
h2az_rna <- merge(h2az,cl,by="row.names",all.x=F)
h2az_map <- h2az_rna[,c(1,14,15,18)]
h2az_map$order <- as.numeric(h2az_map$order)
h2az_map <- h2az_map[order(h2az_map$order),]
head(h2az_map)
row.names(h2az_map) <- h2az_map[,1]
h2az_number <- aggregate(h2az_map$Row.names,list(h2az_map$V13),length)
head(h2az_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",9121),rep("cluster2",9547),rep("cluster3",13804)))
row.names(annotation_row) <- rownames(h2az_map)
h2az_map<- as.matrix(h2az_map[,-1:-3])
h2az_map[,1] <- log2(h2az_map[,1]+1)
pdf("h2az.pdf",width=2,height=3)
p1 <- pheatmap(h2az_map, name = "expression", show_rownames = F,#不显示行名
              show_colnames = F,#不显示列名
              scale = "column",
              col = colorRampPalette(c("blue","white","red"))(100),
              #annotation_col = annotation_col, #列注释信息
              annotation_row = annotation_row,#行注释信息
              row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
              #column_split = annotation_col$group,#列截断
              annotation_names_row = F,#不显示行注释信息
              annotation_names_col = F ,#不显示列注释信息
              column_title = "CL",#不显示列标题
              cluster_rows = FALSE, 
              row_title = "KH2AZ")#不显示行标题
dev.off()
#K9ac
K9ac <- read.table("Heatmap1sortedRegions_CL_CLK9ac_TSS.bed")
K9ac$order <- row.names(K9ac)
rownames(K9ac) <- K9ac[,4]
head(K9ac)
K9ac_rna <- merge(K9ac,cl,by="row.names",all.x=F)
K9ac_map <- K9ac_rna[,c(1,14,15,18)]
K9ac_map$order <- as.numeric(K9ac_map$order)
K9ac_map <- K9ac_map[order(K9ac_map$order),]
head(K9ac_map)
row.names(K9ac_map) <- K9ac_map[,1]
K9ac_number <- aggregate(K9ac_map$Row.names,list(K9ac_map$V13),length)
head(K9ac_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",5594),rep("cluster2",8163),rep("cluster3",18715)))
row.names(annotation_row) <- rownames(K9ac_map)
K9ac_map<- as.matrix(K9ac_map[,-1:-3])
K9ac_map[,1] <- log2(K9ac_map[,1]+1)
pdf("k9ac.pdf",width=2,height=3)
p2 <- pheatmap(K9ac_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CL",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "k9ac")#不显示行标题
dev.off()

#K27ac
K27ac <- read.table("Heatmap1sortedRegions_CL_CLK27ac_TSS.bed")
K27ac$order <- row.names(K27ac)
rownames(K27ac) <- K27ac[,4]
head(K27ac)
K27ac_rna <- merge(K27ac,cl,by="row.names",all.x=F)
K27ac_map <- K27ac_rna[,c(1,14,15,18)]
K27ac_map$order <- as.numeric(K27ac_map$order)
K27ac_map <- K27ac_map[order(K27ac_map$order),]
head(K27ac_map)
row.names(K27ac_map) <- K27ac_map[,1]
K27ac_number <- aggregate(K27ac_map$Row.names,list(K27ac_map$V13),length)
head(K27ac_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",4921),rep("cluster2",7967),rep("cluster3",19584)))
row.names(annotation_row) <- rownames(K27ac_map)
K27ac_map<- as.matrix(K27ac_map[,-1:-3])
K27ac_map[,1] <- log2(K27ac_map[,1]+1)
pdf("K27ac.pdf",width=2,height=3)
p3 <- pheatmap(K27ac_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CL",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K27ac")#不显示行标题
dev.off()

#K43
K43 <- read.table("Heatmap1sortedRegions_CL_CLK43_TSS.bed")
K43$order <- row.names(K43)
rownames(K43) <- K43[,4]
head(K43)
K43_rna <- merge(K43,cl,by="row.names",all.x=F)
K43_map <- K43_rna[,c(1,14,15,18)]
K43_map$order <- as.numeric(K43_map$order)
K43_map <- K43_map[order(K43_map$order),]
head(K43_map)
row.names(K43_map) <- K43_map[,1]
K43_number <- aggregate(K43_map$Row.names,list(K43_map$V13),length)
head(K43_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",4119),rep("cluster2",10215),rep("cluster3",18138)))
row.names(annotation_row) <- rownames(K43_map)
K43_map<- as.matrix(K43_map[,-1:-3])
K43_map[,1] <- log2(K43_map[,1]+1)
pdf("K43.pdf",width=2,height=3)
p4 <- pheatmap(K43_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CL",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K43")#不显示行标题
dev.off()

#K4me2
K4me2 <- read.table("Heatmap1sortedRegions_CL_CLK4me2_TSS.bed")
K4me2$order <- row.names(K4me2)
rownames(K4me2) <- K4me2[,4]
head(K4me2)
K4me2_rna <- merge(K4me2,cl,by="row.names",all.x=F)
K4me2_map <- K4me2_rna[,c(1,14,15,18)]
K4me2_map$order <- as.numeric(K4me2_map$order)
K4me2_map <- K4me2_map[order(K4me2_map$order),]
head(K4me2_map)
row.names(K4me2_map) <- K4me2_map[,1]
K4me2_number <- aggregate(K4me2_map$Row.names,list(K4me2_map$V13),length)
head(K4me2_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",7393),rep("cluster2",13064),rep("cluster3",12015)))
row.names(annotation_row) <- rownames(K4me2_map)
K4me2_map<- as.matrix(K4me2_map[,-1:-3])
K4me2_map[,1] <- log2(K4me2_map[,1]+1)
pdf("K4me2.pdf",width=2,height=3)
p5 <- pheatmap(K4me2_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CL",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K4me2")#不显示行标题
dev.off()

#K36me3
K36me3 <- read.table("Heatmap1sortedRegions_CL_CLK36me3_TSS.bed")
K36me3$order <- row.names(K36me3)
rownames(K36me3) <- K36me3[,4]
head(K36me3)
K36me3_rna <- merge(K36me3,cl,by="row.names",all.x=F)
K36me3_map <- K36me3_rna[,c(1,14,15,18)]
K36me3_map$order <- as.numeric(K36me3_map$order)
K36me3_map <- K36me3_map[order(K36me3_map$order),]
head(K36me3_map)
row.names(K36me3_map) <- K36me3_map[,1]
K36me3_number <- aggregate(K36me3_map$Row.names,list(K36me3_map$V13),length)
head(K36me3_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",6362),rep("cluster2",6783),rep("cluster3",19327)))
row.names(annotation_row) <- rownames(K36me3_map)
K36me3_map<- as.matrix(K36me3_map[,-1:-3])
K36me3_map[,1] <- log2(K36me3_map[,1]+1)
pdf("K36me3.pdf",width=2,height=3)
p6 <- pheatmap(K36me3_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CL",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K36me3")#不显示行标题
dev.off()

#K27me3
K27me3 <- read.table("Heatmap1sortedRegions_CL_CLK27me3_TSS.bed")
K27me3$order <- row.names(K27me3)
rownames(K27me3) <- K27me3[,4]
head(K27me3)
K27me3_rna <- merge(K27me3,cl,by="row.names",all.x=F)
K27me3_map <- K27me3_rna[,c(1,14,15,18)]
K27me3_map$order <- as.numeric(K27me3_map$order)
K27me3_map <- K27me3_map[order(K27me3_map$order),]
head(K27me3_map)
row.names(K27me3_map) <- K27me3_map[,1]
K27me3_number <- aggregate(K27me3_map$Row.names,list(K27me3_map$V13),length)
head(K27me3_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",2053),rep("cluster2",3658),rep("cluster3",26761)))
row.names(annotation_row) <- rownames(K27me3_map)
K27me3_map<- as.matrix(K27me3_map[,-1:-3])
K27me3_map[,1] <- log2(K27me3_map[,1]+1)
pdf("K27me3.pdf",width=2,height=3)
p7 <- pheatmap(K27me3_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CL",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K27me3")#不显示行标题
dev.off()

#####################################################################根
setwd("D:/高粱/ChIP/文章图片/Figure 3/根表达量热图")
rna <- read.csv("D:/高粱/RNA-seq/final_RNAseq/all_tpm_expression.csv",row.names=1)
head(rna)
CR <- rna[,4:6]
CR$average <- rowMeans(CR)
head(CR)
#kh2az
h2az <- read.table("Heatmap1sortedRegions_CR_CRKH2AZ_TSS.bed")
h2az$order <- row.names(h2az)
rownames(h2az) <- h2az[,4]
head(h2az)
h2az_rna <- merge(h2az,CR,by="row.names",all.x=F)
h2az_map <- h2az_rna[,c(1,14,15,18)]
h2az_map$order <- as.numeric(h2az_map$order)
h2az_map <- h2az_map[order(h2az_map$order),]
head(h2az_map)
row.names(h2az_map) <- h2az_map[,1]
h2az_number <- aggregate(h2az_map$Row.names,list(h2az_map$V13),length)
head(h2az_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",10727),rep("cluster2",11918),rep("cluster3",9804)))
row.names(annotation_row) <- rownames(h2az_map)
h2az_map<- as.matrix(h2az_map[,-1:-3])
h2az_map[,1] <- log2(h2az_map[,1]+1)
pdf("h2azCR.pdf",width=2,height=3)
p1 <- pheatmap(h2az_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "KH2AZ")#不显示行标题
dev.off()
#K9ac
K9ac <- read.table("Heatmap1sortedRegions_CR_CRK9ac_TSS.bed")
K9ac$order <- row.names(K9ac)
rownames(K9ac) <- K9ac[,4]
head(K9ac)
K9ac_rna <- merge(K9ac,CR,by="row.names",all.x=F)
K9ac_map <- K9ac_rna[,c(1,14,15,18)]
K9ac_map$order <- as.numeric(K9ac_map$order)
K9ac_map <- K9ac_map[order(K9ac_map$order),]
head(K9ac_map)
row.names(K9ac_map) <- K9ac_map[,1]
K9ac_number <- aggregate(K9ac_map$Row.names,list(K9ac_map$V13),length)
head(K9ac_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",4883),rep("cluster2",8521),rep("cluster3",19045)))
row.names(annotation_row) <- rownames(K9ac_map)
K9ac_map<- as.matrix(K9ac_map[,-1:-3])
K9ac_map[,1] <- log2(K9ac_map[,1]+1)
pdf("k9acCR.pdf",width=2,height=3)
p2 <- pheatmap(K9ac_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "k9ac")#不显示行标题
dev.off()

#K27ac
K27ac <- read.table("Heatmap1sortedRegions_CR_CRK27ac_TSS.bed")
K27ac$order <- row.names(K27ac)
rownames(K27ac) <- K27ac[,4]
head(K27ac)
K27ac_rna <- merge(K27ac,CR,by="row.names",all.x=F)
K27ac_map <- K27ac_rna[,c(1,14,15,18)]
K27ac_map$order <- as.numeric(K27ac_map$order)
K27ac_map <- K27ac_map[order(K27ac_map$order),]
head(K27ac_map)
row.names(K27ac_map) <- K27ac_map[,1]
K27ac_number <- aggregate(K27ac_map$Row.names,list(K27ac_map$V13),length)
head(K27ac_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",2828),rep("cluster2",7052),rep("cluster3",22569)))
row.names(annotation_row) <- rownames(K27ac_map)
K27ac_map<- as.matrix(K27ac_map[,-1:-3])
K27ac_map[,1] <- log2(K27ac_map[,1]+1)
pdf("K27acCR.pdf",width=2,height=3)
p3 <- pheatmap(K27ac_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K27ac")#不显示行标题
dev.off()

#K43
K43 <- read.table("Heatmap1sortedRegions_CR_CRK43_TSS.bed")
K43$order <- row.names(K43)
rownames(K43) <- K43[,4]
head(K43)
K43_rna <- merge(K43,CR,by="row.names",all.x=F)
K43_map <- K43_rna[,c(1,14,15,18)]
K43_map$order <- as.numeric(K43_map$order)
K43_map <- K43_map[order(K43_map$order),]
head(K43_map)
row.names(K43_map) <- K43_map[,1]
K43_number <- aggregate(K43_map$Row.names,list(K43_map$V13),length)
head(K43_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",4027),rep("cluster2",9556),rep("cluster3",18866)))
row.names(annotation_row) <- rownames(K43_map)
K43_map<- as.matrix(K43_map[,-1:-3])
K43_map[,1] <- log2(K43_map[,1]+1)
pdf("K43CR.pdf",width=2,height=3)
p4 <- pheatmap(K43_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K43")#不显示行标题
dev.off()

#K4me2
K4me2 <- read.table("Heatmap1sortedRegions_CR_CRK4me2_TSS.bed")
K4me2$order <- row.names(K4me2)
rownames(K4me2) <- K4me2[,4]
head(K4me2)
K4me2_rna <- merge(K4me2,CR,by="row.names",all.x=F)
K4me2_map <- K4me2_rna[,c(1,14,15,18)]
K4me2_map$order <- as.numeric(K4me2_map$order)
K4me2_map <- K4me2_map[order(K4me2_map$order),]
head(K4me2_map)
row.names(K4me2_map) <- K4me2_map[,1]
K4me2_number <- aggregate(K4me2_map$Row.names,list(K4me2_map$V13),length)
head(K4me2_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",13004),rep("cluster2",9731),rep("cluster3",9714)))
row.names(annotation_row) <- rownames(K4me2_map)
K4me2_map<- as.matrix(K4me2_map[,-1:-3])
K4me2_map[,1] <- log2(K4me2_map[,1]+1)
pdf("K4me2CR.pdf",width=2,height=3)
p5 <- pheatmap(K4me2_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K4me2")#不显示行标题
dev.off()

#K36me3
K36me3 <- read.table("Heatmap1sortedRegions_CR_CRK36me3_TSS.bed")
K36me3$order <- row.names(K36me3)
rownames(K36me3) <- K36me3[,4]
head(K36me3)
K36me3_rna <- merge(K36me3,CR,by="row.names",all.x=F)
K36me3_map <- K36me3_rna[,c(1,14,15,18)]
K36me3_map$order <- as.numeric(K36me3_map$order)
K36me3_map <- K36me3_map[order(K36me3_map$order),]
head(K36me3_map)
row.names(K36me3_map) <- K36me3_map[,1]
K36me3_number <- aggregate(K36me3_map$Row.names,list(K36me3_map$V13),length)
head(K36me3_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",6548),rep("cluster2",6652),rep("cluster3",19249)))
row.names(annotation_row) <- rownames(K36me3_map)
K36me3_map<- as.matrix(K36me3_map[,-1:-3])
K36me3_map[,1] <- log2(K36me3_map[,1]+1)
pdf("K36me3CR.pdf",width=2,height=3)
p6 <- pheatmap(K36me3_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K36me3")#不显示行标题
dev.off()

#K27me3
K27me3 <- read.table("Heatmap1sortedRegions_CR_CRK27me3_TSS.bed")
K27me3$order <- row.names(K27me3)
rownames(K27me3) <- K27me3[,4]
head(K27me3)
K27me3_rna <- merge(K27me3,CR,by="row.names",all.x=F)
K27me3_map <- K27me3_rna[,c(1,14,15,18)]
K27me3_map$order <- as.numeric(K27me3_map$order)
K27me3_map <- K27me3_map[order(K27me3_map$order),]
head(K27me3_map)
row.names(K27me3_map) <- K27me3_map[,1]
K27me3_number <- aggregate(K27me3_map$Row.names,list(K27me3_map$V13),length)
head(K27me3_number)
annotation_row = data.frame(
  cluster = c(rep("cluster1",1913),rep("cluster2",3563),rep("cluster3",26973)))
row.names(annotation_row) <- rownames(K27me3_map)
K27me3_map<- as.matrix(K27me3_map[,-1:-3])
K27me3_map[,1] <- log2(K27me3_map[,1]+1)
pdf("K27me3_CR.pdf",width=2,height=3)
p7 <- pheatmap(K27me3_map, name = "expression", show_rownames = F,#不显示行名
               show_colnames = F,#不显示列名
               scale = "column",
               col = colorRampPalette(c("blue","white","red"))(100),
               #annotation_col = annotation_col, #列注释信息
               annotation_row = annotation_row,#行注释信息
               row_split = annotation_row$cluster,#行截断（按照pathway，不像之前随机）
               #column_split = annotation_col$group,#列截断
               annotation_names_row = F,#不显示行注释信息
               annotation_names_col = F ,#不显示列注释信息
               column_title = "CR",#不显示列标题
               cluster_rows = FALSE, 
               row_title = "K27me3")#不显示行标题
dev.off()

######绘制质控图
setwd("D:/高粱/ChIP")
qc <- read.csv("质控统计统计R1.csv",row.names = 1)
head(qc)
#绘制柱形图
#col <- colorRampPalette(brewer.pal(8,'Accent'))(16)
qc <- qc[,c(1,9,10)]
qc<- qc %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer( cols =  c("FRiP":"SPOT"),
                names_to = 'stage',
                values_to = 'expr')
col1<- brewer.pal(8,'Dark2')
col2 <- brewer.pal(8,"Accent")
pdf("signal2noise.pdf")
ggboxplot(qc, x = "stage", y = "expr", fill = "stage",
          bxp.errorbar=T,# 
          #add = c("mean_sd"), error.plot = "crossbar") +
          add = c("dotplot"),
          add.params = list( size = 0.4,alpha = 0.9),
          legend="none",
                outlier.shape = NA) +
  #coord_cartesian(ylim = ylim1) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2")) +
  scale_fill_manual(values=c("#c2e59c","#be93c5")) +
  theme(plot.title = element_text(hjust = 0.5,size=10),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  xlab(NULL) + ylab("Signal to noise rate")
dev.off()

###TTS-TES差异显著性分析
setwd("D:/高粱/ChIP/文章图片/Figure 4/TSS-TES")
name <- c("CRK4me2_CLK4me2_down","CRK4me2_CLK4me2_up","CRK9ac_CLK9ac_down",
          "CRK9ac_CLK9ac_up","CRK27ac_CLK27ac_down","CRK27ac_CLK27ac_up","CRK27me3_CLK27me3_down",
          "CRK27me3_CLK27me3_up","CRK36me3_CLK36me3_down","CRK36me3_CLK36me3_up",
          "CRK43_CLK43_down","CRK43_CLK43_up","CRKH2AZ_CLKH2AZ_down","CRKH2AZ_CLKH2AZ_up")
plot_list <- list()
j=1
for (i in name) {
  x <- read.csv(paste0(i,".tab"))
  x <- t(x[,c(-2,-78:-152)])
  colnames(x) <- c("bins",strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2])
  x <- x[-1,]
  x <- as.data.frame(x)
  head(x)
  x<- x %>% 
    rownames_to_column(var = 'sample') %>% 
    pivot_longer( cols =  c(strsplit(i,"_")[[1]][1]:strsplit(i,"_")[[1]][2]),
                  names_to = 'stage',
                  values_to = 'expr')
  x$expr <- as.numeric(x$expr)
  head(x)
  my_comparisons <- list(c(strsplit(i,"_")[[1]][1], strsplit(i,"_")[[1]][2]))
  p1 <- ggboxplot(x, x = "stage", y = "expr", fill = "stage",
                  bxp.errorbar=T, 
                  outlier.shape = NA) + 
    #add = c("mean_sd"), error.plot = "crossbar") +
    #add = c("boxplot")) +
    #scale_fill_manual(values=c("#aec7e8","c5b0d5","#ffbb78","#98df8a","#ff9896","#c49c94","#f7b6d2","lightgrey")) +
    #scale_fill_manual(values=c("#B2DF8A","#FDBF6F","#CAB2D6","#FB9A99","#1F78B4","#33A02C","#A6CEE3","#E31A1C","#FF7F00","#6A3D9A","#FFFF99","#B15928")) +
    scale_fill_manual(values=c(brewer.pal(8,"Paired"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))) + #设置填充色
    ggtitle(i) + 
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text.x=element_text()) +
    xlab(NULL) + ylab("value") +
    theme (legend.position = 'none') +
    stat_compare_means(comparisons = my_comparisons)
  plot_list[[j]] <- p1
  j = j+1
}

grid.arrange(plot_list[[1]],plot_list[[2]], 
             plot_list[[3]],plot_list[[4]], 
             plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],
             plot_list[[9]], plot_list[[10]],
             plot_list[[11]],plot_list[[12]],
             plot_list[[13]],plot_list[[14]],
             nrow=5,ncol=3)     %>%  ggsave("TTS-TES_DEGs_significant.pdf",.,width=210,height=550, units="mm")



#################分析差异组蛋白修饰基因，以基因为单位计算count只后构建DEseq2矩阵，随后分析差异表达基因
##构建每种修饰的count矩阵，同时，分析根和叶中的差异修饰基因
setwd("D:/高粱/ChIP/定量")
#构建tpm矩阵
a <- list.files(path="D:/高粱/ChIP/定量", pattern="_counts_subread.txt", all.files=F, full.names=F)
dir = paste("./",a,sep="")     
n = length(dir)  
merge.data = read.table(file = dir[1],header=TRUE,dec = ".")
merge.data$count <- merge.data[,7]
merge.data <- merge.data[,-7]
merge.data$type <- rep(dir[1],nrow(merge.data))
for (i in 2:n){
  new.data = read.table(file = dir[i], header=TRUE, dec = ".",skip=1)
  new.data$count <- new.data[,7]
  new.data <- new.data[,-7]
  new.data$type <- rep(dir[i],nrow(new.data))
  merge.data = rbind(merge.data,new.data)
}
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
#计算TPM值
total_tpm <- data.frame()
for (i in 1:n)  {
  pea<-merge.data[which(merge.data$type == dir[i]),]
  head(pea)
  count <- as.data.frame(pea$count)
  Length <- as.data.frame(pea$Length)
  tpmm <- apply(count, 2, function(x) tpm(x, Length))
  peak_tpms <-cbind(pea,tpmm)
  peak_tpms$tpm <- peak_tpms[,9]
  peak_tpms <- peak_tpms[,-9]
  total_tpm <- rbind(total_tpm,peak_tpms)
}
#长变宽
total_tpmR1 <- total_tpm[,-7]
toal <- tidyr::spread(data=total_tpmR1,key=type,value=tpm)
head(toal)
write.csv(toal,"chip_tpm_matrix.csv")
###构建count矩阵
head(merge.data)
toal_count <- tidyr::spread(data=merge.data,key=type,value=count)
head(toal_count)
write.csv(toal_count,"chip_count_matrix.csv")

###差异peak分析，用Deseq2的方法做
library(DESeq2)
library(dplyr)
## 导入TPM数据矩阵
#countdata <- read.csv("H3K4me3_count.csv", row.names = 1)
countdata <- read.csv("chip_count_matrix.csv", row.names = 1)
row.names(countdata) <- countdata[,1]
countdata <- countdata[,-1:-6]
head(countdata)
tpmdata <- read.csv("chip_TPM_matrix.csv", row.names = 1)
row.names(tpmdata) <- tpmdata[,1]
tpmdata <- tpmdata[,-1:-6]
## 过滤在所有重复样本中小于1的基因
countdata = countdata[rowMax(as.matrix(countdata)) > 1,]
#导入样本注释信息
coldata  <- read.csv("col_data.csv",row.names = 1)
coldata$sample <- coldata$condition
head(coldata)
#检查数据Counts文件与coldata数据是否匹配
all(rownames(coldata) %in% colnames(countdata))  
all(rownames(coldata) == colnames(countdata))
#PCA分析
rsem.in <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
rsem.de <- DESeq(rsem.in)
dds<-rsem.in
rsem.rlog <- rlog(rsem.de)
pdf("Deseq2_PCR1.pdf")
pcaData <- plotPCA(rsem.rlog, intgroup=c("condition"), returnData=TRUE)
#plot(pcaData[,1:2],pch=19,col=pal_npg("nrc")(12))
#text(pcaData[,1],pcaData[,2]+0.2,row.names(pcaData),cex=0.8)
ggplot(data = pcaData, aes(x = PC1, y = PC2, label = condition, shape = condition)) +
  geom_point(aes(color = condition), size = 3) +  #根据样本坐标绘制二维散点图
  geom_text_repel(aes(label = condition,color=condition),
                  #color = "gray20",
                  #data = subset(pca_sample, samples %in% pointsToLabel),
                  force = 10) +
  scale_shape_manual(values = c(1:60)) +
  #scale_color_d3("category20") +
  #scale_linetype_manual(values = c('twodash', 'longdash', 'dashed')) 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) 
dev.off()
##差异分析
# 制作差异矩阵
##CLK27ac和CRK27ac
a1 <- countdata[,c(1,2,15,16)]
c1 <- tpmdata[,c(1,2,15,16)]
head(a1)
b1 <- coldata[c(1,2,15,16),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"CLK27ac_CRk27ac_diff_allR1.csv")  # 导出所有的差异文件

##CLK9ac和CRK9ac
a2 <- countdata[,c(11,12,25,26)]
c2 <- tpmdata[,c(11,12,25,26)]
head(a2)
b2 <- coldata[c(11,12,25,26),]
head(b2)
dds <-  DESeqDataSetFromMatrix(countData = a2,colData = b2,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c2,by="row.names",all.x=F)
head(diff_tpm)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
write.csv(diff_tpm,"CLK9acvsCRK9ac_diff_allR1.csv")  # 导出所有的差异文件

##CLK27me3vsCRK27me3
a1 <- countdata[,c(3,4,17,18)]
c1 <- tpmdata[,c(3,4,17,18)]
head(a1)
b1 <- coldata[c(3,4,17,18),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"CLK27me3vsCRK27me3_diff_allR1.csv")  # 导出所有的差异文件

##CLK36me3vsCRK36me3
a1 <- countdata[,c(5,6,19,20)]
c1 <- tpmdata[,c(5,6,19,20)]
head(a1)
b1 <- coldata[c(5,6,19,20),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"CLK36me3vsCRK36me3_diff_allR1.csv")  # 导出所有的差异文件

##CLK4me3vsCRK4me3
a1 <- countdata[,c(7,8,21,22)]
c1 <- tpmdata[,c(7,8,21,22)]
head(a1)
b1 <- coldata[c(7,8,21,22),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"CLK4me3vsCRK4me3_diff_allR1.csv")  # 导出所有的差异文件

##CLK4me2vsCRK4me2
a1 <- countdata[,c(9,10,23,24)]
c1 <- tpmdata[,c(9,10,23,24)]
head(a1)
b1 <- coldata[c(9,10,23,24),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"CLK4me2vsCRK4me2_diff_allR1.csv")  # 导出所有的差异文件

##CLKH1AZvsCRKH1AZ
a1 <- countdata[,c(13,14,27,28)]
c1 <- tpmdata[,c(13,14,27,28)]
head(a1)
b1 <- coldata[c(13,14,27,28),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"CLKH1AZvsCRKH1AZ_diff_allR1.csv")  # 导出所有的差异文件

###########################PEG处理根与对照的差异peak分析
setwd("D:/高粱/ChIP/定量/PRvsCR")
##PRK27ac和CRK27ac
a1 <- countdata[,c(43,44,15,16)]
c1 <- tpmdata[,c(43,44,15,16)]
head(a1)
b1 <- coldata[c(43,44,15,16),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PRK27ac_CRk27ac_diff_allR1.csv")  # 导出所有的差异文件

##PRK9ac和CRK9ac
a2 <- countdata[,c(53,54,25,26)]
c2 <- tpmdata[,c(53,54,25,26)]
head(a2)
b2 <- coldata[c(53,54,25,26),]
head(b2)
dds <-  DESeqDataSetFromMatrix(countData = a2,colData = b2,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c2,by="row.names",all.x=F)
head(diff_tpm)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
write.csv(diff_tpm,"PRK9acvsCRK9ac_diff_allR1.csv")  # 导出所有的差异文件

##PRK27me3vsCRK27me3
a1 <- countdata[,c(45,46,17,18)]
c1 <- tpmdata[,c(45,46,17,18)]
head(a1)
b1 <- coldata[c(45,46,17,18),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1000,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
diff[which(diff$padj < 0.05),]
head(diff)
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PRK27me3vsCRK27me3_diff_allR1.csv")  # 导出所有的差异文件


fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- data.frame(apply(fpkm,2,fpkmToTpm), check.names = FALSE)
head(tpms)
colSums(tpms)


##PRK36me3vsCRK36me3
a1 <- countdata[,c(47,48,19,20)]
c1 <- tpmdata[,c(47,48,19,20)]
head(a1)
b1 <- coldata[c(47,48,19,20),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PRK36me3vsCRK36me3_diff_allR1.csv")  # 导出所有的差异文件

##PRK4me3vsCRK4me3
a1 <- countdata[,c(49,50,21,22)]
c1 <- tpmdata[,c(49,50,21,22)]
head(a1)
b1 <- coldata[c(49,50,21,22),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PRK4me3vsCRK4me3_diff_allR1.csv")  # 导出所有的差异文件

##PRK4me2vsCRK4me2
a1 <- countdata[,c(51,52,23,24)]
c1 <- tpmdata[,c(51,52,23,24)]
head(a1)
b1 <- coldata[c(51,52,23,24),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PRK4me2vsCRK4me2_diff_allR1.csv")  # 导出所有的差异文件

##PRKH1AZvsCRKH1AZ
a1 <- countdata[,c(55,56,27,28)]
c1 <- tpmdata[,c(55,56,27,28)]
head(a1)
b1 <- coldata[c(55,56,27,28),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PPKH1AZvsCRKH1AZ_diff_allR1.csv")  # 导出所有的差异文件

###########################PEG处理叶与对照的差异peak分析
##PLK27ac和CLK27ac
a1 <- countdata[,c(29,30,1,2)]
c1 <- tpmdata[,c(29,30,1,2)]
head(a1)
b1 <- coldata[c(29,30,1,2),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PLK27ac_CLk27ac_diff_allR1.csv")  # 导出所有的差异文件

##PLK9ac和CLK9ac
a2 <- countdata[,c(39,40,11,12)]
c2 <- tpmdata[,c(39,40,11,12)]
head(a2)
b2 <- coldata[c(39,40,11,12),]
head(b2)
dds <-  DESeqDataSetFromMatrix(countData = a2,colData = b2,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c2,by="row.names",all.x=F)
head(diff_tpm)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
write.csv(diff_tpm,"PLK9acvsCLK9ac_diff_allR1.csv")  # 导出所有的差异文件

##PLK27me3vsCLK27me3
a1 <- countdata[,c(31,32,3,4)]
c1 <- tpmdata[,c(31,32,3,4)]
head(a1)
b1 <- coldata[c(31,32,3,4),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=TRUE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PLK27me3vsCLK27me3_diff_allR1.csv")  # 导出所有的差异文件

##PLK36me3vsCLK36me3
a1 <- countdata[,c(33,34,5,6)]
c1 <- tpmdata[,c(33,34,5,6)]
head(a1)
b1 <- coldata[c(33,34,5,6),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PLK36me3vsCLK36me3_diff_allR1.csv")  # 导出所有的差异文件

##PLK4me3vsCLK4me3
a1 <- countdata[,c(35,36,7,8)]
c1 <- tpmdata[,c(35,36,7,8)]
head(a1)
b1 <- coldata[c(35,36,7,8),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PLK4me3vsCLK4me3_diff_allR1.csv")  # 导出所有的差异文件

##PLK4me2vsCLK4me2
a1 <- countdata[,c(37,38,9,10)]
c1 <- tpmdata[,c(37,38,9,10)]
head(a1)
b1 <- coldata[c(37,38,9,10),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PLK4me2vsCLK4me2_diff_allR1.csv")  # 导出所有的差异文件

##PLKH1AZvsCLKH1AZ
a1 <- countdata[,c(41,42,13,14)]
c1 <- tpmdata[,c(41,42,13,14)]
head(a1)
b1 <- coldata[c(41,42,13,14),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=FALSE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
write.csv(diff_tpm,"PLKH1AZvsCLKH1AZ_diff_allR1.csv")  # 导出所有的差异文件



#####差异表达基因和差异修饰基因的关系,弃用
#导入差异表达基因
library(data.table)
library(viridis)
deg <- read.csv("D:/高粱/RNA-seq/final_RNAseq/CLvsCR_DEGs.csv",row.names=1)
degs <- read.csv("D:/高粱/RNA-seq/final_RNAseq/CLvsCR.csv",row.names=1) ##计算所有基因的差异倍数，是差异表达基因就用Deseq2的结果，如果CL和CR都小于0.5则倍数为0，如果分母为0，则就去log2分子，如果分子为0，则分子设为1，其余就取log2(CR/CL)
degs$Row.names <- row.names(degs)
#deg <- deg[which(deg$padj < 0.05),]
deg <- deg[,c(-1,-3,-4,-5,-7:-14)]
colnames(deg) <- c("RNA_log2FC","RNA_padj")
head(deg)
##CLvsCR
name <- c("K27ac","K27me3","K36me3","K4me3","K4me2","K9ac","KH1AZ")
m = 1
n=15
j=1
plot_list <- list()
for (i in name) {
  his <- read.csv(paste0("CL",i,"vsCR",i,"_diff_all.csv"),row.names=1)
  his <- his[which(his$padj < 0.05),]
  row.names(his) <- his[,1]
  his <- his[,c(-1,-2,-4,-5,-6,-8:-11)]
  colnames(his) <- c("histone_log2FC","histone_padj")
  his <- his[which(his$histone_log2FC >1 | his$histone_log2FC < -1),]
  head(his)
  #融合
  deg_his1 <- merge(deg,his,"row.names",all.x=T)
  deg_his2 <- merge(his,deg,"row.names",all.x=T)
  deg_his2 <- deg_his2[,c(1,4,5,2,3)]
  deg_his <- rbind(deg_his1,deg_his2)
  deg_his <- unique(deg_his)
  deg_his[is.na(deg_his)] <- 0
  #deg_his <- na.omit(deg_his)
  up <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1) & deg_his$histone_log2FC > 1),]
  down <- deg_his[which(deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC < -1),]
  up$type <- rep("Up",nrow(up))
  down$type <- rep("Down",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(deg_his$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(deg_his,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ###将no分为3类，no1是rna倍数未知，no2是hisone倍数未知，no3是都已知
  no3 <- no[which(no$RNA_log2FC != 0 & no$histone_log2FC != 0),]
  no1 <- no[which(no$RNA_log2FC == 0),]
  no1 <- merge(no1,degs,by="Row.names",all.x=F)
  head(no1)
  no1$RNA_log2FC <- no1$log2FoldChange
  no1$RNA_padj <- no1$padj
  no1 <- no1[,c(-7:-10)]
  no2 <- no[which(no$histone_log2FC == 0),]
  ot <- read.csv("chip_tpm_matrix.csv",row.names=1)
  row.names(ot) <- ot[,1]
  ot <- ot[,c(-1:-6)]
  to <- ot[,c(m,m+1,n,n+1)]
  m=m+2
  n=n+2
  to$Row.names <- row.names(to)
  to$CL <- (to[,1] + to[,2])/2
  to$CR <- (to[,3] + to[,4])/2
  to1 <- to[which(to$CL < 0.5 & to$CR < 0.5),]
  to1$log2FoldChange <- 0
  to1$padj <- 1
  to2 <- to[which(to$CL == 0 & to$CR > 0.5),]
  to2$log2FoldChange <- log2(to2$CR/(to2$CL+1))
  to2$padj <- 1
  to3 <- to[which(to$CL >0.5 & to$CR == 0),]
  to3$log2FoldChange <- log2((to3$CR+1)/to3$CL)
  to3$padj <- 1
  t1 <- rbind(to1,to2,to3)
  t1 <- t1[,-1:-4]
  x<- list(t1$Row.names)
  y <- list(to$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  t2 <- merge(to,mm,by="Row.names",is.all=FALSE)
  t2$log2FoldChange <- log2(t2$CR/t2$CL)
  t2$padj <- 1
  t2 <- t2[,-2:-5]
  row.names(t2) <- t2[,1]
  head(t2)
  tm <- rbind(t1,t2)
  no2 <- merge(no2,tm,by="Row.names",all.x=F)
  no2$histone_log2FC <- no2$log2FoldChange
  no2$histone_padj <- no2$padj
  no2 <- no2[,c(-7:-10)]
  tn <- unique(rbind(up,down,no1,no2,no3))
  up <-  tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC > 1),]
  down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC < -1),]
  up$type <- rep("Up",nrow(up))
  down$type <- rep("Down",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(tn$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(tn,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  this_tile <- paste0('The number of up gene is ',nrow(up),',down gene is ',nrow(down),",R value is ",cor.test(ttt$RNA_log2FC,ttt$histone_log2FC)$estimate[[1]])
  #绘图
  head(ttt)
  p1 <- ggplot(ttt, aes(RNA_log2FC, histone_log2FC)) + 
    geom_point(aes(colour = factor(type))) +
    #geom_point() +
    scale_colour_manual(values=c("#FF7F00","lightgrey","#FF7F00")) +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(size = 12,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          #axis.ticks = element_line(size = .8),
          #去除图例标题
          legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    theme_bw() + 
    labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
    geom_hline(aes(yintercept=0), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed")
  plot_list[[j]] <- p1
  j=j+1
}

##CLvsCR，负相关的两外着色，R 大于1 H 小于-1，弃用
name <- c("K27ac","K27me3","K36me3","K4me3","K4me2","K9ac","KH1AZ")
m = 1
n=15
j=1
qlot_list <- list()
for (i in name) {
  his <- read.csv(paste0("CL",i,"vsCR",i,"_diff_all.csv"),row.names=1)
  his <- his[which(his$padj < 0.05),]
  row.names(his) <- his[,1]
  his <- his[,c(-1,-2,-4,-5,-6,-8:-11)]
  colnames(his) <- c("histone_log2FC","histone_padj")
  his <- his[which(his$histone_log2FC >1 | his$histone_log2FC < -1),]
  head(his)
  #融合
  deg_his1 <- merge(deg,his,"row.names",all.x=T)
  deg_his2 <- merge(his,deg,"row.names",all.x=T)
  deg_his2 <- deg_his2[,c(1,4,5,2,3)]
  deg_his <- rbind(deg_his1,deg_his2)
  deg_his <- unique(deg_his)
  deg_his[is.na(deg_his)] <- 0
  #deg_his <- na.omit(deg_his)
  up <- deg_his[which(((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1) & deg_his$histone_log2FC > 1) | ((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1) & deg_his$histone_log2FC < -1)),]
  down <- deg_his[which(((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1) & deg_his$histone_log2FC < -1) | ((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1) & deg_his$histone_log2FC > 1)),]
  up$type <- rep("Up",nrow(up))
  down$type <- rep("Down",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(deg_his$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(deg_his,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ###将no分为3类，no1是rna倍数未知，no2是hisone倍数未知，no3是都已知
  no3 <- no[which(no$RNA_log2FC != 0 & no$histone_log2FC != 0),]
  no1 <- no[which(no$RNA_log2FC == 0),]
  no1 <- merge(no1,degs,by="Row.names",all.x=F)
  head(no1)
  no1$RNA_log2FC <- no1$log2FoldChange
  no1$RNA_padj <- no1$padj
  no1 <- no1[,c(-7:-10)]
  no2 <- no[which(no$histone_log2FC == 0),]
  ot <- read.csv("chip_tpm_matrix.csv",row.names=1)
  row.names(ot) <- ot[,1]
  ot <- ot[,c(-1:-6)]
  to <- ot[,c(m,m+1,n,n+1)]
  m=m+2
  n=n+2
  to$Row.names <- row.names(to)
  to$CL <- (to[,1] + to[,2])/2
  to$CR <- (to[,3] + to[,4])/2
  to1 <- to[which(to$CL < 0.5 & to$CR < 0.5),]
  to1$log2FoldChange <- 0
  to1$padj <- 1
  to2 <- to[which(to$CL == 0 & to$CR > 0.5),]
  to2$log2FoldChange <- log2(to2$CR/(to2$CL+1))
  to2$padj <- 1
  to3 <- to[which(to$CL >0.5 & to$CR == 0),]
  to3$log2FoldChange <- log2((to3$CR+1)/to3$CL)
  to3$padj <- 1
  t1 <- rbind(to1,to2,to3)
  t1 <- t1[,-1:-4]
  x<- list(t1$Row.names)
  y <- list(to$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  t2 <- merge(to,mm,by="Row.names",is.all=FALSE)
  t2$log2FoldChange <- log2(t2$CR/t2$CL)
  t2$padj <- 1
  t2 <- t2[,-2:-5]
  row.names(t2) <- t2[,1]
  head(t2)
  tm <- rbind(t1,t2)
  no2 <- merge(no2,tm,by="Row.names",all.x=F)
  no2$histone_log2FC <- no2$log2FoldChange
  no2$histone_padj <- no2$padj
  no2 <- no2[,c(-7:-10)]
  tn <- unique(rbind(up,down,no1,no2,no3))
  #up <-  tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC > 1),]
  #down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC < -1),]
  #up$type <- rep("Up",nrow(up))
  #down$type <- rep("Down",nrow(down))
  #up <- tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC < -1),]
  #down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC > 1),]
  up <- tn[which((tn$RNA_log2FC > 1 & tn$histone_log2FC > 1) | (tn$RNA_log2FC < -1 & tn$histone_log2FC < -1)),]
  down <- tn[which((tn$RNA_log2FC > 1 & tn$histone_log2FC < -1) | (tn$RNA_log2FC < -1 & tn$histone_log2FC > 1)),]
  up$type <- rep("positive",nrow(up))
  down$type <- rep("negative ",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(tn$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(tn,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  this_tile <- paste0('positive: ',nrow(up),',negative: ',nrow(down),",R value is ",round(cor.test(ttt$RNA_log2FC,ttt$histone_log2FC)$estimate[[1]],digits=3))
  #绘图
  head(ttt)
  p1 <- ggplot(ttt, aes(RNA_log2FC, histone_log2FC)) + 
    geom_point(aes(colour = factor(type))) +
    #geom_point() +
    scale_colour_manual(values=c("#4E5BF4","lightgrey","#E7298A")) +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(size = 12,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          #axis.ticks = element_line(size = .8),
          #去除图例标题
          legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    theme_bw() + 
    labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
    geom_hline(aes(yintercept=0), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed")
  qlot_list[[j]] <- p1
  j=j+1
}

#组图
grid.arrange(qlot_list[[1]],qlot_list[[2]], 
             qlot_list[[3]],qlot_list[[4]], 
             qlot_list[[5]],qlot_list[[6]],
             qlot_list[[7]],
             nrow=4,ncol=2)     %>%  ggsave("DEGs_DMG_ScatterplotR1.pdf",.,width=210,height=297, units="mm")


#####差异表达基因和差异修饰基因的关系,最终用这个
#导入差异表达基因
library(data.table)
library(viridis)
library(gridExtra)
library(ggplot2)
library(lattice)
library(ggpubr)
library(patchwork)
library(scales)
library(tidyr)
####图片和作图用的数据（手动更改了数值，变为0的）放在了DEG_DMG_correlationR1，没有更改数值的原始数据放在了DEG_DMG_correlationR2
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation")
deg <- read.csv("all_tpm_expression.csv",row.names=1)
#degs <- read.csv("D:/高粱/RNA-seq/final_RNAseq/CLvsCR.csv",row.names=1) ##计算所有基因的差异倍数，是差异表达基因就用Deseq2的结果，如果CL和CR都小于0.5则倍数为0，如果分母为0，则就去log2分子，如果分子为0，则分子设为1，其余就取log2(CR/CL)
#degs$Row.names <- row.names(degs)
#deg <- deg[which(deg$padj < 0.05),]
#deg <- deg[which(deg$CL != "Unexpress" & deg$CR != "Unexpress"),]
deg <- deg[,c(7:8,9,10)]
colnames(deg) <- c("RNA_log2FC","RNA_padj","RNA_CL","RNA_CR")
head(deg)
##CLvsCR
name <- c("K27ac","K27me3","K36me3","K4me3","K4me2","K9ac","KH1AZ")
m = 1
n=15
j=1
plot_list <- list()
for (i in name) {
  his <- read.csv(paste0("CL",i,"vsCR",i,"_diff_allR2.csv"),row.names=1)
  #his <- his[which(his$CL != "Unexpress" & his$CR != "Unexpress"),]
  row.names(his) <- his[,1]
  his <- his[,c(3,7,12,13)]
  colnames(his) <- c("histone_log2FC","histone_padj","CL","CR")
  head(his)
  #融合
  deg_his1 <- merge(deg,his,"row.names",all.x=T)
  deg_his2 <- merge(his,deg,"row.names",all.x=T)
  deg_his2 <- deg_his2[,c(1,6,7,8,9,2,3,4,5)]
  deg_his <- rbind(deg_his1,deg_his2)
  deg_his <- unique(deg_his)
  deg_his[is.na(deg_his)] <- 0
  #deg_his <- deg_his[which((deg_his$CL != "Unexpress" & deg_his$CR != "Unexpress") | (deg_his$RNA_CL != "Unexpress" & deg_his$RNA_CL != "Unexpress")),]
  row.names(deg_his) <- deg_his[,1]
  #deg_his <- na.omit(deg_his)
  up <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1 & deg_his$histone_log2FC > 0.75 & deg_his$RNA_padj < 0.05) | (deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC < -0.75 & deg_his$RNA_padj < 0.05)),]
  up <- up[which(up$CL == "Express" | up$CR == "Express"),]
  up <- up[which((up$RNA_CL == "Express" | up$RNA_CR == "Express")),]
  down <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1 & deg_his$histone_log2FC < -0.75 & deg_his$RNA_padj < 0.05) | (deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC > 0.75 & deg_his$RNA_padj < 0.05)),]
  down <- down[which((down$CL == "Express" | down$CR == "Express")),]
  down <- down[which((down$RNA_CL == "Express" | down$RNA_CR == "Express")),]
  up$type <- rep("Positive",nrow(up))
  down$type <- rep("Negative",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(deg_his$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(deg_his,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ###将RNA_log2FC绝对值大于1和histone_log2FC绝对值大于0.75的设为0
  no1 <- no[which(no$RNA_log2FC > -1 & no$RNA_log2FC < 1 & no$histone_log2FC > -0.75 & no$histone_log2FC < 0.75),]
  no1$RNA_log2FC <- 0
  no1$histone_log2FC <- 0
  ###将unexcpected但显著差异且差异倍数大于阈值的设为0，没用显著性这个条件来筛
  no2 <- no[which((no$RNA_log2FC < -1 & no$histone_log2FC > 0.75) | (no$RNA_log2FC < -1 & no$histone_log2FC < -0.75) | (no$RNA_log2FC > 1 & no$histone_log2FC > 0.75) | (no$RNA_log2FC > 1 & no$histone_log2FC < -0.75)),]
  no2$RNA_log2FC <- 0
  no2$histone_log2FC <- 0 
  ###提取剩余的no
  x<- list(c(no1$Row.names,no2$Row.names))
  y <- list(no$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no3 <- merge(no,mm,by="Row.names",is.all=FALSE)
  ###重新合并no
  no <- no3
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  this_tile <- paste0('The number of up gene is ',nrow(up),',down gene is ',nrow(down),",R value is ",cor.test(ttt$RNA_log2FC,ttt$histone_log2FC)$estimate[[1]])
  #t.test(ttt$RNA_log2FC,ttt$histone_log2FC)
  #绘图
  head(ttt)
  p1 <- ggplot(ttt, aes(RNA_log2FC, histone_log2FC)) + 
    geom_point(aes(colour = factor(type))) +
    #geom_point() +
    scale_colour_manual(values=c("#4E5BF4","lightgrey","#E7298A")) +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(size = 6,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          #axis.ticks = element_line(size = .8),
          #去除图例标题
          legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    #theme_classic() +
    theme_bw() +
    theme(text = element_text(family = ,face='bold',size=6),
          axis.text = element_text(size = 6,face = 'bold'),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
    geom_hline(aes(yintercept=0.75), colour="#BB0000", linetype="dashed") +
    geom_hline(aes(yintercept=-0.75), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=-1), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed")
  plot_list[[j]] <- p1
  j=j+1
  #保存数据
  no <- rbind(no1,no2,no3)
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no)
  write.csv(ttt,paste0("CL",i,"vsCR",i,"_DEG_DMG_clusterd.csv"))
}
####保存数据
grid.arrange(plot_list[[1]],plot_list[[2]], 
             plot_list[[3]],plot_list[[4]], 
             plot_list[[5]],plot_list[[6]],
             plot_list[[7]],
             nrow=4,ncol=2)     %>%  ggsave("DEGs_DMG_correlationR1.pdf",.,width=210,height=297, units="mm")

#####GO富集分析，对激活类markerpositive 上下调基因分别做GO，抑制类marker negtive 上下调基因分别做GO；https://biit.cs.ut.ee/gprofiler/gost
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR1")
nn <- c("H2AZ_neg_down","H2AZ_neg_up","k4me2_neg_down","k4me2_neg_up",
        "k4me3_pos_down","k4me3_pos_up","k9ac_pos_down","k9ac_pos_up",
        "k27ac_pos_down","k27ac_pos_up","k36_pos_down","k36_pos_up",
        "k273_neg_down","k273_neg_up")
nn <- c("k4me2_neg_down","k4me2_neg_up",
        "k4me3_pos_down","k4me3_pos_up","k9ac_pos_down","k9ac_pos_up",
        "k27ac_pos_down","k27ac_pos_up","k36_pos_down","k36_pos_up",
        "k273_neg_down","k273_neg_up")
c <- data.frame()
for (i in nn) {
  a <- read.csv(paste0(i,".csv"),header=T)
  a <- a[,c(1,2,3,4,5,9)]
  head(a)
  name <- c("GO:BP","GO:CC","GO:MF")
  plot_list <- list()
  for (j in name) {
    b <- a[which(a$source == j),]
    if (nrow(b) != 0) {
      b <- b[order(b$intersection_size),]
      mat = GO_similarity(b$term_id)
      term_id <- row.names(mat)
      row.names(b) <- b$term_id
      b <- b[term_id,]
      b$term_name <- factor(b$term_name, levels=b$term_name, ordered=TRUE)
      b$term_id <- factor(b$term_id, levels=b$term_id, ordered=TRUE)
      b$type <- rep(paste(i,j,sep="_"),nrow(b))
      c <- rbind(c,b)
    }
  }
}
#write.csv(c,"GO_total_reduplicated.csv")
c <- read.csv("GO_total_reduplicated.csv",row.names=1)
dim(c)
unique(c$type)
head(c)
c <- unique(c)
c<- c %>% group_by(type) %>% top_n(n = 10, wt = -log10(adjusted_p_value))
c$term_name <- gsub("acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen","",c$term_name)
c1 <- c[which(c$source == "GO:BP"),]
c1 <- c1[order(c1$type),]
c1$term_name <- factor(c1$term_name, levels=unique(c1$term_name), ordered=TRUE)
#c1$term_id <- factor(c1$term_id, levels=unique(c1$term_id), ordered=TRUE)
p1 <- ggplot(c1,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#f80759",high="#A8C0FF")+
  ggtitle("BP") +
  labs(x=NULL,y=NULL)

c2 <- c[which(c$source == "GO:MF"),]
c2 <- c2[order(c2$type),]
c2$term_name <- factor(c2$term_name, levels=unique(c2$term_name), ordered=TRUE)
#c2$term_id <- factor(c2$term_id, levels=unique(c2$term_id), ordered=TRUE)
p2 <- ggplot(c2,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#f80759",high="#A8C0FF")+
  ggtitle("MF") +
  labs(x=NULL,y=NULL)

c3 <- c[which(c$source == "GO:CC"),]
c3 <- c3[order(c3$type),]
c3$term_name <- factor(c3$term_name, levels=unique(c3$term_name), ordered=TRUE)
#c3$term_id <- factor(c3$term_id, levels=unique(c3$term_id), ordered=TRUE)
p3 <- ggplot(c3,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#f80759",high="#A8C0FF")+
  ggtitle("CC") +
  labs(x=NULL,y=NULL)
#组图
grid.arrange(p1,p2,p3,
             nrow=3,ncol=1)     %>%  ggsave("GO_total_top.pdf",.,width=300,height=1800, limitsize = FALSE,units="mm")

#柱形图，不用
p1 <- ggplot(b,mapping = aes(x=intersection_size,y=term_id,fill=adjusted_p_value))+
  geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  #geom_bar(stat='identity',position = position_dodge(),show.legend = TRUE) +
  #scale_fill_manual(values=brewer.pal(10,"Paired")) +
  #scale_fill_manual(values=c("#CAB2D6","#6A3D9A")) +
  #scale_colour_gradient(low = "#659999", high = "#f4791f") +
  scale_fill_gradient2(low="#3a1c71",mid="#d76d77",high="#ffaf7b") +
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  ggtitle(i) +
  xlab(NULL) + ylab(j) +
  geom_text(aes(label = term_name), size = 3, hjust = 0, vjust = 0, position = "stack")
plot_list[[m]] <- p1
m+1

#####KEGG富集分析，对激活类markerpositive 上下调基因分别做GO，抑制类marker negtive 上下调基因分别做GO；https://biit.cs.ut.ee/gprofiler/gost
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR1")
nn <- c("H2AZ_neg_down","H2AZ_neg_up","k4me2_neg_down","k4me2_neg_up",
        "k4me3_pos_down","k4me3_pos_up","k9ac_pos_down","k9ac_pos_up",
        "k27ac_pos_down","k27ac_pos_up","k36_pos_down","k36_pos_up",
        "k273_neg_down","k273_neg_up")
b <- data.frame()
for (i in nn) {
  a <- read.csv(paste0(i,".csv"),header=T)
  #a <- a[which(a$source=="KEGG"),c(1,2,3,4,5,9)]
  a <- a[which(a$source=="KEGG"),]
  a$type <- rep(i,nrow(a))
  b <- rbind(b,a)
}
write.csv(b,"total_KEGG.csv",row.names = F)

##绘制气泡图
library(tidyverse)
#pdf("KEGG_all.pdf",family="MT")
b <- b[order(b$term_name,b$type),]
b$term_name <- factor(b$term_name,levels=unique(b$term_name),ordered=TRUE)
ggplot(b,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  labs(x=NULL,y=NULL)
dev.off()
write.csv(data,"KEGG_total.csv")
write.csv(top,"KEGG_top10.csv")

##绘制Spot Matrix
ke <- b[,c(2,5,6,7)]
ke <- ke[which(ke$term_name != "Metabolic pathways"),]
#ke_wide <- spread(data=ke,key=type,value=intersection_size)
#ke_wide[is.na(ke_wide)] <- 0
#head(ke_wide)
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(legend.position = "top"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  #scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))

colors <- pals::parula(10)[c(1,4,7,9)]
pdf("KEGG_all.pdf",family = "ArialMT",width=12,height = 10)
ggplot(ke, aes(type, term_name)) + spot.theme + ggtitle("KEGG") +
  geom_point(colour = "grey20", aes(size = 1)) +
  #geom_point(colour = "white", aes(size = 0.8)) +
  geom_point( aes(colour = adjusted_p_value,size = 0.81*intersection_size)) +
  scale_color_gradient(low=colors[[4]],high=colors[[3]]) +
  labs(x=NULL,y=NULL)
dev.off()

ggplot(ke, aes(type, term_name)) + spot.theme + ggtitle("KEGG") +
  geom_point(colour = "gray20",    aes(size = 1)) +
  geom_point(colour = colors[[3]], aes(size = intersection_size))



#####计算CR CL差异表达基因的metaplot的显著性
setwd("D:/高粱/宇宁/Figure4R1-A")
tobeCopy <- list.files(".", pattern="*.tab")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}

##################PRvsCR###############
setwd("D:/高粱/ChIP/定量/PRvsCR/data/PRvsCR")
deg <- read.csv("../all_tpm_expression_pr_vs_cr.csv",row.names=1)
deg <- deg[,c(7:8,9,10)]
colnames(deg) <- c("RNA_log2FC","RNA_padj","RNA_PR","RNA_CR")
head(deg)
##CRvsCR
name <- c("K27ac","K27me3","K36me3","K4me3","K4me2","K9ac","KH1AZ")
m = 1
n=15
j=1
plot_list <- list()
for (i in name) {
  his <- read.csv(paste0("PR",i,"vsCR",i,"_diff_allR1.csv"),row.names=1)
  #his <- his[which(his$CR != "Unexpress" & his$CR != "Unexpress"),]
  row.names(his) <- his[,1]
  his <- his[,c(3,7,12,13)]
  colnames(his) <- c("histone_log2FC","histone_padj","PR","CR")
  head(his)
  #融合
  deg_his1 <- merge(deg,his,"row.names",all.x=TRUE)
  deg_his2 <- merge(his,deg,"row.names",all.x=TRUE)
  deg_his2 <- deg_his2[,c(1,6,7,8,9,2,3,4,5)]
  deg_his <- rbind(deg_his1,deg_his2)
  deg_his <- unique(deg_his)
  deg_his[is.na(deg_his)] <- 0
  #deg_his <- deg_his[which((deg_his$PR != "Unexpress" & deg_his$CR != "Unexpress") | (deg_his$RNA_PR != "Unexpress" & deg_his$RNA_PR != "Unexpress")),]
  row.names(deg_his) <- deg_his[,1]
  #deg_his <- na.omit(deg_his)
  up <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1 & deg_his$histone_log2FC > 0.75 & deg_his$RNA_padj < 0.05) | (deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC < -0.75 & deg_his$RNA_padj < 0.05)),]
  up <- up[which(up$PR == "Express" | up$CR == "Express"),]
  up <- up[which((up$RNA_PR == "Express" | up$RNA_CR == "Express")),]
  down <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1 & deg_his$histone_log2FC < -0.75 & deg_his$RNA_padj < 0.05) | (deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC > 0.75 & deg_his$RNA_padj < 0.05)),]
  down <- down[which((down$PR == "Express" | down$CR == "Express")),]
  down <- down[which((down$RNA_PR == "Express" | down$RNA_CR == "Express")),]
  up$type <- rep("Positive",nrow(up))
  down$type <- rep("Negative",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(deg_his$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(deg_his,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ###将RNA_log2FC绝对值大于1和histone_log2FC绝对值大于0.75的设为0
  no1 <- no[which(no$RNA_log2FC > -1 & no$RNA_log2FC < 1 & no$histone_log2FC > -0.75 & no$histone_log2FC < 0.75),]
  no1$RNA_log2FC <- 0
  no1$histone_log2FC <- 0
  ###将unexcpected但显著差异且差异倍数大于阈值的设为0，没用显著性这个条件来筛
  no2 <- no[which((no$RNA_log2FC < -1 & no$histone_log2FC > 0.75) | (no$RNA_log2FC < -1 & no$histone_log2FC < -0.75) | (no$RNA_log2FC > 1 & no$histone_log2FC > 0.75) | (no$RNA_log2FC > 1 & no$histone_log2FC < -0.75)),]
  no2$RNA_log2FC <- 0
  no2$histone_log2FC <- 0 
  ###提取剩余的no
  x<- list(c(no1$Row.names,no2$Row.names))
  y <- list(no$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no3 <- merge(no,mm,by="Row.names",is.all=FALSE)
  ###重新合并no
  no <- no3
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  this_tile <- paste0('The number of up gene is ',nrow(up),',down gene is ',nrow(down),",R value is ",cor.test(ttt$RNA_log2FC,ttt$histone_log2FC)$estimate[[1]])
  t.test(ttt$RNA_log2FC,ttt$histone_log2FC)
  i
  #绘图
  head(ttt)
  p1 <- ggplot(ttt, aes(RNA_log2FC, histone_log2FC)) + 
    geom_point(aes(colour = factor(type))) +
    #geom_point() +
    scale_colour_manual(values=c("#4E5BF4","lightgrey","#E7298A")) +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(size = 6,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          #axis.ticks = element_line(size = .8),
          #去除图例标题
          legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    #theme_PRassic() +
    theme_bw() +
    theme(text = element_text(family = ,face='bold',size=6),
          axis.text = element_text(size = 6,face = 'bold'),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
    geom_hline(aes(yintercept=0.75), colour="#BB0000", linetype="dashed") +
    geom_hline(aes(yintercept=-0.75), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=-1), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed")
  plot_list[[j]] <- p1
  j=j+1
  #保存数据
  no <- rbind(no1,no2,no3)
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no)
  write.csv(ttt,paste0("PR",i,"vsCR",i,"_DEG_DMG_clusterd.csv"))
}
####保存数据
grid.arrange(plot_list[[1]],plot_list[[2]], 
             plot_list[[3]],plot_list[[4]], 
             plot_list[[5]],plot_list[[6]],
             plot_list[[7]],
             nrow=4,ncol=2)     %>%  ggsave("DEGs_DMG_correlationR1.pdf",.,width=210,height=297, units="mm")


#####/PRvsCR，GO富集分析，对激活类markerpositive 上下调基因分别做GO，
#抑制类marker negtive 上下调基因分别做GO；https://biit.cs.ut.ee/gprofiler/gost
setwd("D:/高粱/ChIP/定量/PRvsCR/data/PRvsCR")
nn <- c("H2AZ_neg_down","H2AZ_neg_up","k4me2_neg_down","k4me2_neg_up",
        "k43_pos_down","k43_pos_up","k9ac_pos_down","k9ac_pos_up",
        "k27ac_pos_down","k27ac_pos_up","k36me3_pos_up")
c <- data.frame()
for (i in nn) {
  a <- read.csv(paste0(i,".csv"),header=T)
  a <- a[,c(1,2,3,4,5,9)]
  head(a)
  name <- c("GO:BP","GO:CC","GO:MF")
  plot_list <- list()
  for (j in name) {
    b <- a[which(a$source == j),]
    if (nrow(b) > 1) {
      b <- b[order(b$intersection_size),]
      mat = GO_similarity(b$term_id)
      term_id <- row.names(mat)
      row.names(b) <- b$term_id
      b <- b[term_id,]
      b$term_name <- factor(b$term_name, levels=b$term_name, ordered=TRUE)
      b$term_id <- factor(b$term_id, levels=b$term_id, ordered=TRUE)
      b$type <- rep(paste(i,j,sep="_"),nrow(b))
      c <- rbind(c,b)
    }
  }
}
write.csv(c,"GO_total_reduplicated.csv")
c <- read.csv("GO_total_reduplicated.csv",row.names=1)
dim(c)
unique(c$type)
head(c)
c <- unique(c)
c$term_name <- gsub("acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen","",c$term_name)
c1 <- c[which(c$source == "GO:BP"),]
c1$term_name <- factor(c1$term_name, levels=unique(c1$term_name), ordered=TRUE)
c1$term_id <- factor(c1$term_id, levels=unique(c1$term_id), ordered=TRUE)
p1 <- ggplot(c1,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  ggtitle("BP") +
  labs(x=NULL,y=NULL)

c2 <- c[which(c$source == "GO:MF"),]
c2$term_name <- factor(c2$term_name, levels=unique(c2$term_name), ordered=TRUE)
c2$term_id <- factor(c2$term_id, levels=unique(c2$term_id), ordered=TRUE)
p2 <- ggplot(c2,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  ggtitle("MF") +
  labs(x=NULL,y=NULL)

c3 <- c[which(c$source == "GO:CC"),]
c3$term_name <- factor(c3$term_name, levels=unique(c3$term_name), ordered=TRUE)
c3$term_id <- factor(c3$term_id, levels=unique(c3$term_id), ordered=TRUE)
p3 <- ggplot(c3,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  ggtitle("CC") +
  labs(x=NULL,y=NULL)
#组图
grid.arrange(p1,p2,p3,
             nrow=3,ncol=1)     %>%  ggsave("GO_PRvsCR_total.pdf",.,width=300,height=1000, limitsize = FALSE,units="mm")
######KEGG
setwd("D:/高粱/ChIP/定量/PRvsCR/data/PRvsCR")
nn <- c("H2AZ_neg_down","H2AZ_neg_up","k4me2_neg_down","k4me2_neg_up",
        "k43_pos_down","k43_pos_up","k9ac_pos_down","k9ac_pos_up",
        "k27ac_pos_down","k27ac_pos_up","k36me3_pos_up")
b <- data.frame()
for (i in nn) {
  a <- read.csv(paste0(i,".csv"),header=T)
  #a <- a[which(a$source=="KEGG"),c(1,2,3,4,5,9)]
  a <- a[which(a$source=="KEGG"),]
  a$type <- rep(i,nrow(a))
  b <- rbind(b,a)
}
write.csv(b,"total_KEGG.csv",row.names = F)
##绘制气泡图
library(tidyverse)
b <- b[order(b$term_name,b$type),]
b$term_name <- factor(b$term_name,levels=unique(b$term_name),ordered=TRUE)
##绘制Spot Matrix
ke <- b[,c(2,3,5,9,12)]
#ke <- ke[order(ke$term_name,ke$type),]
ke <- ke[order(ke$type),]
ke$term_name <- factor(ke$term_name,levels=unique(ke$term_name),ordered=TRUE)

#ke <- ke[which(ke$term_name != "Metabolic pathways"),]
#ke_wide <- spread(data=ke,key=type,value=intersection_size)
#ke_wide[is.na(ke_wide)] <- 0
#head(ke_wide)
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(legend.position = "top"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  #scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))

colors <- pals::parula(10)[c(1,4,7,9)]
pdf("KEGG_PRvsCR_all.pdf",family = "ArialMT",width=12,height = 10)
ggplot(ke, aes(type, term_name)) + spot.theme + ggtitle("KEGG") +
  geom_point(colour = "grey20", aes(size = 1)) +
  #geom_point(colour = "white", aes(size = 0.8)) +
  geom_point( aes(colour = adjusted_p_value,size = 0.81*intersection_size)) +
  scale_color_gradient(low=colors[[1]],high=colors[[3]]) +
  labs(x=NULL,y=NULL)
dev.off()



















#####计算PR CR差异表达基因的metaplot的显著性
setwd("D:/高粱/ChIP/定量/PRvsCR/data/metaplot")
tobeCopy <- list.files(".", pattern="*.tab")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}


setwd("D:/高粱/宇宁/PEG处理前后的修饰变化/PLvsCL")
deg <- read.csv("all_tpm_expression_pl_vs_cl.csv",row.names=1)
deg <- deg[,c(7:8,9,10)]
colnames(deg) <- c("RNA_log2FC","RNA_padj","RNA_pl","RNA_cl")
head(deg)
##CRvsCR
name <- c("K27ac","K27me3","K36me3","K4me3","K4me2","K9ac","KH1AZ")
m = 1
n=15
j=1
plot_list <- list()
for (i in name) {
  his <- read.csv(paste0("PL",i,"vsCL",i,"_diff_allR1.csv"),row.names=1)
  #his <- his[which(his$CR != "Unexpress" & his$CR != "Unexpress"),]
  row.names(his) <- his[,1]
  his <- his[,c(3,7,12,13)]
  colnames(his) <- c("histone_log2FC","histone_padj","PL","CL")
  head(his)
  #融合
  deg_his1 <- merge(deg,his,"row.names",all.x=TRUE)
  deg_his2 <- merge(his,deg,"row.names",all.x=TRUE)
  deg_his2 <- deg_his2[,c(1,6,7,8,9,2,3,4,5)]
  deg_his <- rbind(deg_his1,deg_his2)
  deg_his <- unique(deg_his)
  deg_his[is.na(deg_his)] <- 0
  #deg_his <- deg_his[which((deg_his$PR != "Unexpress" & deg_his$CR != "Unexpress") | (deg_his$RNA_PR != "Unexpress" & deg_his$RNA_PR != "Unexpress")),]
  row.names(deg_his) <- deg_his[,1]
  #deg_his <- na.omit(deg_his)
  up <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1 & deg_his$histone_log2FC > 0.75 & deg_his$RNA_padj < 0.05) | (deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC < -0.75 & deg_his$RNA_padj < 0.05)),]
  up <- up[which(up$PL == "Express" | up$CL == "Express"),]
  up <- up[which((up$RNA_pl == "Express" | up$RNA_cl == "Express")),]
  down <- deg_his[which((deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC > 1 & deg_his$histone_log2FC < -0.75 & deg_his$RNA_padj < 0.05) | (deg_his$histone_padj < 0.05 & deg_his$RNA_log2FC < -1 & deg_his$histone_log2FC > 0.75 & deg_his$RNA_padj < 0.05)),]
  down <- down[which((down$PL == "Express" | down$CL == "Express")),]
  down <- down[which((down$RNA_pl == "Express" | down$RNA_cl == "Express")),]
  up$type <- rep("Positive",nrow(up))
  down$type <- rep("Negative",nrow(down))
  x <- list(c(up$Row.names,down$Row.names))
  y<- list(deg_his$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no <- merge(deg_his,mm,by="Row.names",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  ###将RNA_log2FC绝对值大于1和histone_log2FC绝对值大于0.75的设为0
  no1 <- no[which(no$RNA_log2FC > -1 & no$RNA_log2FC < 1 & no$histone_log2FC > -0.75 & no$histone_log2FC < 0.75),]
  no1$RNA_log2FC <- 0
  no1$histone_log2FC <- 0
  ###将unexcpected但显著差异且差异倍数大于阈值的设为0，没用显著性这个条件来筛
  no2 <- no[which((no$RNA_log2FC < -1 & no$histone_log2FC > 0.75) | (no$RNA_log2FC < -1 & no$histone_log2FC < -0.75) | (no$RNA_log2FC > 1 & no$histone_log2FC > 0.75) | (no$RNA_log2FC > 1 & no$histone_log2FC < -0.75)),]
  no2$RNA_log2FC <- 0
  no2$histone_log2FC <- 0 
  ###提取剩余的no
  x<- list(c(no1$Row.names,no2$Row.names))
  y <- list(no$Row.names)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "Row.names"
  no3 <- merge(no,mm,by="Row.names",is.all=FALSE)
  ###重新合并no
  no <- no3
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  this_tile <- paste0('The number of up gene is ',nrow(up),',down gene is ',nrow(down),",R value is ",cor.test(ttt$RNA_log2FC,ttt$histone_log2FC)$estimate[[1]])
  t.test(ttt$RNA_log2FC,ttt$histone_log2FC)
  i
  #绘图
  head(ttt)
  p1 <- ggplot(ttt, aes(RNA_log2FC, histone_log2FC)) + 
    geom_point(aes(colour = factor(type))) +
    #geom_point() +
    scale_colour_manual(values=c("#4E5BF4","lightgrey","#E7298A")) +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(size = 6,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          #axis.ticks = element_line(size = .8),
          #去除图例标题
          legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    #theme_PRassic() +
    theme_bw() +
    theme(text = element_text(family = ,face='bold',size=6),
          axis.text = element_text(size = 6,face = 'bold'),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
    geom_hline(aes(yintercept=0.75), colour="#BB0000", linetype="dashed") +
    geom_hline(aes(yintercept=-0.75), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=-1), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed")
  plot_list[[j]] <- p1
  j=j+1
  #保存数据
  no <- rbind(no1,no2,no3)
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no)
  write.csv(ttt,paste0("PR",i,"vsCR",i,"_DEG_DMG_clusterd.csv"))
}

########/PLvsCL，GO富集分析，对激活类markerpositive 上下调基因分别做GO，
#抑制类marker negtive 上下调基因分别做GO；https://biit.cs.ut.ee/gprofiler/gost
setwd("D:/高粱/宇宁/PEG处理前后的修饰变化/PLvsCL")
nn <- c("H2AZ_neg_down","H2AZ_neg_up","k4me2_neg_down","k4me2_neg_up",
        "k4me3_pos_up","k9ac_pos_down","k9ac_pos_up",
        "k27ac_pos_down","k27ac_pos_up")
c <- data.frame()
for (i in nn) {
  a <- read.csv(paste0(i,".csv"),header=T)
  a <- a[,c(1,2,3,4,5,9)]
  head(a)
  name <- c("GO:BP","GO:CC","GO:MF")
  plot_list <- list()
  for (j in name) {
    b <- a[which(a$source == j),]
    if (nrow(b) > 2) {
      b <- b[order(b$intersection_size),]
      #mat = GO_similarity(b$term_id)
      #term_id <- row.names(mat)
      #row.names(b) <- b$term_id
      #b <- b[term_id,]
      b$term_name <- factor(b$term_name, levels=b$term_name, ordered=TRUE)
      b$term_id <- factor(b$term_id, levels=b$term_id, ordered=TRUE)
      b$type <- rep(paste(i,j,sep="_"),nrow(b))
      c <- rbind(c,b)
    }
  }
}
write.csv(c,"GO_total_reduplicated.csv")
c <- read.csv("GO_total_reduplicated.csv",row.names=1)
dim(c)
unique(c$type)
head(c)
c <- unique(c)
c$term_name <- gsub("acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen","",c$term_name)
c1 <- c[which(c$source == "GO:BP"),]
c1$term_name <- factor(c1$term_name, levels=unique(c1$term_name), ordered=TRUE)
c1$term_id <- factor(c1$term_id, levels=unique(c1$term_id), ordered=TRUE)
p1 <- ggplot(c1,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  ggtitle("BP") +
  labs(x=NULL,y=NULL)

c2 <- c[which(c$source == "GO:MF"),]
c2$term_name <- factor(c2$term_name, levels=unique(c2$term_name), ordered=TRUE)
c2$term_id <- factor(c2$term_id, levels=unique(c2$term_id), ordered=TRUE)
p2 <- ggplot(c2,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  ggtitle("MF") +
  labs(x=NULL,y=NULL)

c3 <- c[which(c$source == "GO:CC"),]
c3$term_name <- factor(c3$term_name, levels=unique(c3$term_name), ordered=TRUE)
c3$term_id <- factor(c3$term_id, levels=unique(c3$term_id), ordered=TRUE)
p3 <- ggplot(c3,aes(x=type,y=term_name))+
  geom_point(aes(size=`intersection_size`,
                 color=`adjusted_p_value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  ggtitle("CC") +
  labs(x=NULL,y=NULL)
#组图
grid.arrange(p1,p2,p3,
             nrow=3,ncol=1)     %>%  ggsave("GO_PLvsCL_total.pdf",.,width=300,height=1000, limitsize = FALSE,units="mm")
######KEGG
b <- data.frame()
for (i in nn) {
  a <- read.csv(paste0(i,".csv"),header=T)
  #a <- a[which(a$source=="KEGG"),c(1,2,3,4,5,9)]
  a <- a[which(a$source=="KEGG"),]
  a$type <- rep(i,nrow(a))
  b <- rbind(b,a)
}
write.csv(b,"total_KEGG.csv",row.names = F)
##绘制气泡图
library(tidyverse)
b <- b[order(b$term_name,b$type),]
b$term_name <- factor(b$term_name,levels=unique(b$term_name),ordered=TRUE)
##绘制Spot Matrix
ke <- b[,c(2,3,5,9,12)]
#ke <- ke[order(ke$term_name,ke$type),]
ke <- ke[order(ke$type),]
ke$term_name <- factor(ke$term_name,levels=unique(ke$term_name),ordered=TRUE)
#ke <- ke[which(ke$term_name != "Metabolic pathways"),]
#ke_wide <- spread(data=ke,key=type,value=intersection_size)
#ke_wide[is.na(ke_wide)] <- 0
#head(ke_wide)
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(legend.position = "top"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  #scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))
colors <- pals::parula(10)[c(1,4,7,9)]
pdf("KEGG_PLvsCL_all.pdf",family = "ArialMT",width=12,height = 10)
ggplot(ke, aes(type, term_name)) + spot.theme + ggtitle("KEGG") +
  geom_point(colour = "grey20", aes(size = 1)) +
  #geom_point(colour = "white", aes(size = 0.8)) +
  geom_point( aes(colour = adjusted_p_value,size = 0.81*intersection_size)) +
  scale_color_gradient(low=colors[[1]],high=colors[[3]]) +
  labs(x=NULL,y=NULL)
dev.off()

###gprofiler2
#install.packages("gprofiler2")
library(gprofiler2)
gostres1 <- gost(query = c("KEGG:00196"), organism = "sbicolor",multi_query = FALSE)
write.csv(gostres1$meta$genes_metadata$query$query_1$mapping$`KEGG:00196`,"KEGG00196.csv")

gostres2 <- gost(query = c("KEGG:00195"), organism = "sbicolor",multi_query = FALSE)
head(gostres2$result, 3)
write.csv(gostres2$meta$genes_metadata$query$query_1$mapping$`KEGG:00195`,"KEGG00195.csv")


#####计算Pl Cl差异表达基因的metaplot的显著性
setwd("D:/高粱/宇宁/PEG处理前后的修饰变化/PLvsCL")
tobeCopy <- list.files(".", pattern="*.tab")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}

#####核小体部分，计算Pl Cl差异表达基因的metaplot的显著性
setwd("D:/高粱/核小体/figure4_tab/C")
tobeCopy <- list.files(".", pattern="plvscl*")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}

tobeCopy <- list.files(".", pattern="prvscr*")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}

setwd("D:/高粱/核小体/figure4_tab/A")
tobeCopy <- list.files(".", pattern="*tab")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}

setwd("D:/高粱/核小体/figure4_tab/B")
tobeCopy <- list.files(".", pattern="*tab")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}

####H3K27me3的metaplot的pvalue
setwd("D:/高粱/H3K27me3Metaplot")
tobeCopy <- list.files(".", pattern="*tab")
for (i in tobeCopy) {
  te <- read.csv(i,sep="\t",header=T,row.names=1)
  te <- te[,-1]
  te <- t(te)
  te <- te[,-1]
  te <- te[!is.na(te[,1]),]
  print (t.test(te[,1],te[,2]))
  print (i)
}


########根叶中各种组蛋白标签富集变化两两重叠的关联图
install.packages("UpSetR")
all.genes <- unique(c(HSCR_5c3.DEG, HSCR_10c2.DEG, HSCR_20c7.DEG, HSCR_23c9.DEG, HSCR_1c11.DEG, HSCR_17c8.DEG))
length(all.genes)

DEG.UpSetR.df <- data.frame(Name=all.genes, `HSCR#5`=as.integer(all.genes %in% HSCR_5c3.DEG),
                            `HSCR#10`=as.integer(all.genes %in% HSCR_10c2.DEG),
                            `HSCR#20`=as.integer(all.genes %in% HSCR_20c7.DEG),
                            `HSCR#23`=as.integer(all.genes %in% HSCR_23c9.DEG),
                            `HSCR#1`=as.integer(all.genes %in% HSCR_1c11.DEG),
                            `HSCR#17`=as.integer(all.genes %in% HSCR_17c8.DEG)
)


#####分析基因长度、外显子长度、内含子长度分布
setwd("D:/高粱/ChIP/normalized_BW/genomeFile")
exon <- read.table("sorghum.exon.gtf")
exon$exon <- exon$V3-exon$V2
head(exon)
p1 <- ggplot(exon, mapping = aes(exon)) +
  geom_density(size=1.5,key_glyph = draw_key_path,color="#A6CEE3") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_colour_manual(values=brewer.pal(9,"Paired")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(size=1,colour = "black"),
        axis.text.x = element_text(color="black", 
                                   size=12),
        axis.text.y = element_text(color="black", 
                                   size=12)) +
  annotation_logticks()
intron <- read.table("sorghum.intron.bed")
intron$intron <- intron$V3-intron$V2
head(intron)
p2 <- ggplot(intron, mapping = aes(intron)) +
  geom_density(size=1.5,key_glyph = draw_key_path,color="#B2DF8A") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_colour_manual(values=brewer.pal(9,"Paired")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(size=1,colour = "black"),
        axis.text.x = element_text(color="black", 
                                   size=12),
        axis.text.y = element_text(color="black", 
                                   size=12)) +
  annotation_logticks()

gene <- read.table("geneR1.bed")
gene$gene <- gene$V3-gene$V2
head(gene)
p3 <- ggplot(gene, mapping = aes(gene)) +
  geom_density(size=1.5,key_glyph = draw_key_path,color="#FB9A99") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_colour_manual(values=brewer.pal(9,"Paired")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(size=1,colour = "black"),
        axis.text.x = element_text(color="black", 
                                   size=12),
        axis.text.y = element_text(color="black", 
                                   size=12)) +
  annotation_logticks()
####组图
grid.arrange(p1,p2,p3,
             nrow=3,ncol=1)     %>%  ggsave("exon_intron_gene_length.pdf",.,width=210,height=297, units="mm")


####根和叶组蛋白修饰变化的基因顺式作用元件分析（上升和下降分开做）
gene <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/差异修饰的基因的顺式调控元件motif富集分析/Sorghum_promoter.csv",header=F)
row.names(gene) <- gene$V4
head(gene)
##CLvsCR
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation")
name <- c("K27ac","K27me3","K36me3","K4me3","K4me2","K9ac","KH1AZ")
i=name[[1]]
for (i in name) {
  his <- read.csv(paste0("CL",i,"vsCR",i,"_diff_allR2.csv"),row.names=1)
  #his <- his[which(his$CL != "Unexpress" & his$CR != "Unexpress"),]
  row.names(his) <- his[,1]
  his <- his[,c(3,7,12,13)]
  colnames(his) <- c("histone_log2FC","histone_padj","CL","CR")
  head(his)
  up <- his[which((his$histone_padj < 0.05 & his$histone_log2FC > 0.75)),]
  up <- up[which(up$CL == "Express" | up$CR == "Express"),]
  down <- his[which((his$histone_padj < 0.05 & his$histone_log2FC < -0.75)),]
  down <- down[which((down$CL == "Express" | down$CR == "Express")),]
  ##上游2kb的序列
  up <- merge(up,gene,by="row.names",all=F)
  down <- merge(down,gene,by="row.names",all=F)
  up <- up[,c(6,7,8)]
  down <- down[,c(6,7,8)]
  write.csv(up,paste0("差异修饰的基因的顺式调控元件motif富集分析/CL",i,"vsCR",i,"_up.bed"))
  write.csv(down,paste0("差异修饰的基因的顺式调控元件motif富集分析/CL",i,"vsCR",i,"_down.bed"))
}
##提取序列
cd /public/home/chaohe/sorghum/chip/align/diff
perl -p -i -e 's/^M//g' *bed
perl -p -i -e 's/,/\t/g' *bed
perl -p -i -e 's/\"//g' *bed
for i in *bed
do
lsns -t net | awk 'NR == 1 {next} {print $2"\t"$3"\t"$4}' $i > R1_"$i"
done
module load BEDTools/2.27
for i in R1*
do
bedtools getfasta -fi /public/home/chaohe/db/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -bed $i -fo "$i".fa
done

#####干旱诱导上升的有几个修饰也做一下k9ac，k27ac，k4me3的motif富集分析
####PLvsCL
setwd("D:/高粱/ChIP/定量/PLvsCL")
name <- c("K27ac","K4me3","K9ac")
####标记express和Unexpress，两个重复都大于10则为express
for (i in name) {
  his <- read.csv(paste0("PL",i,"vsCL",i,"_diff_allR1.csv"),row.names=1)
  his1 <- his[which(his[,8] > 10 & his[,9] > 10),]
  his2 <- his[which((his[,8] <= 10 & his[,9] >= 10) | (his[,8] >= 10 & his[,9] <= 10) | (his[,8] <= 10 & his[,9] <= 10)),]
  his1$PL <- rep("Express",nrow(his1))
  his2$PL <- rep("Unexpress",nrow(his2))
  his <- rbind(his1,his2)
  his1 <- his[which(his[,10] > 10 & his[,11] > 10),]
  his2 <- his[which((his[,10] <= 10 & his[,11] >= 10) | (his[,10] >= 10 & his[,11] <= 10) | (his[,10] <= 10 & his[,11] <= 10)),]
  his1$CL <- rep("Express",nrow(his1))
  his2$CL <- rep("Unexpress",nrow(his2))
  his <- rbind(his1,his2)
  write.csv(his,paste0("PL",i,"vsCL",i,"_diff_allR2.csv"))
}
####调取位置
name <- c("K27ac","K4me3","K9ac")
for (i in name) {
  his <- read.csv(paste0("PL",i,"vsCL",i,"_diff_allR2.csv"),row.names=1)
  #his <- his[which(his$CL != "Unexpress" & his$CR != "Unexpress"),]
  row.names(his) <- his[,1]
  his <- his[,c(3,7,12,13)]
  colnames(his) <- c("histone_log2FC","histone_padj","PL","CL")
  head(his)
  up <- his[which((his$histone_padj < 0.05 & his$histone_log2FC > 0.75)),]
  up <- up[which(up$PL == "Express" | up$CL == "Express"),]
  down <- his[which((his$histone_padj < 0.05 & his$histone_log2FC < -0.75)),]
  down <- down[which((down$PL == "Express" | down$CL == "Express")),]
  ##上游2kb的序列
  up <- merge(up,gene,by="row.names",all=F)
  down <- merge(down,gene,by="row.names",all=F)
  up <- up[,c(6,7,8)]
  down <- down[,c(6,7,8)]
  write.csv(up,paste0("差异修饰的基因的顺式调控元件motif富集分析/PL",i,"vsCL",i,"_up.bed"))
  write.csv(down,paste0("差异修饰的基因的顺式调控元件motif富集分析/PL",i,"vsCL",i,"_down.bed"))
}
####PRvsCR
setwd("D:/高粱/ChIP/定量/PRvsCR")
name <- c("K27ac","K4me3","K9ac")
####标记express和Unexpress，两个重复都大于10则为express
for (i in name) {
  his <- read.csv(paste0("PR",i,"vsCR",i,"_diff_allR1.csv"),row.names=1)
  his1 <- his[which(his[,8] > 10 & his[,9] > 10),]
  his2 <- his[which((his[,8] <= 10 & his[,9] >= 10) | (his[,8] >= 10 & his[,9] <= 10) | (his[,8] <= 10 & his[,9] <= 10)),]
  his1$PL <- rep("Express",nrow(his1))
  his2$PL <- rep("Unexpress",nrow(his2))
  his <- rbind(his1,his2)
  his1 <- his[which(his[,10] > 10 & his[,11] > 10),]
  his2 <- his[which((his[,10] <= 10 & his[,11] >= 10) | (his[,10] >= 10 & his[,11] <= 10) | (his[,10] <= 10 & his[,11] <= 10)),]
  his1$CL <- rep("Express",nrow(his1))
  his2$CL <- rep("Unexpress",nrow(his2))
  his <- rbind(his1,his2)
  write.csv(his,paste0("PR",i,"vsCR",i,"_diff_allR2.csv"))
}
####调取位置
name <- c("K27ac","K4me3","K9ac")
for (i in name) {
  #his <- read.csv(paste0("D:/高粱/ChIP/定量/PRvsCR/data/PRvsCR/","PR",i,"vsCR",i,"_diff_allR1.csv"),row.names=1)
  his <- read.csv(paste0("PR",i,"vsCR",i,"_diff_allR2.csv"),row.names=1)
  row.names(his) <- his[,1]
  his <- his[,c(3,7,12,13)]
  colnames(his) <- c("histone_log2FC","histone_padj","PR","CR")
  head(his)
  up <- his[which((his$histone_padj < 0.05 & his$histone_log2FC > 0.75)),]
  up <- up[which(up$PR == "Express" | up$CR == "Express"),]
  down <- his[which((his$histone_padj < 0.05 & his$histone_log2FC < -0.75)),]
  down <- down[which((down$PR == "Express" | down$CR == "Express")),]
  ##上游2kb的序列
  up <- merge(up,gene,by="row.names",all=F)
  down <- merge(down,gene,by="row.names",all=F)
  up <- up[,c(6,7,8)]
  down <- down[,c(6,7,8)]
  write.csv(up,paste0("差异修饰的基因的顺式调控元件motif富集分析/PR",i,"vsCR",i,"_up.bed"))
  write.csv(down,paste0("差异修饰的基因的顺式调控元件motif富集分析/PR",i,"vsCR",i,"_down.bed"))
}

##提取序列
cd /public/home/chaohe/sorghum/chip/align/diff
perl -p -i -e 's/^M//g' *bed
perl -p -i -e 's/,/\t/g' *bed
perl -p -i -e 's/\"//g' *bed
for i in PR*bed PL*bed
do
lsns -t net | awk 'NR == 1 {next} {print $2"\t"$3"\t"$4}' $i > R1_"$i"
done
module load BEDTools/2.27
for i in R1_PR*bed R1_PL*bed
  do
bedtools getfasta -fi /public/home/chaohe/db/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -bed $i -fo "$i".fa
done

######绘制motif富集分析图
setwd("D:/高粱/ChIP/定量/PRvsCR/差异修饰的基因的顺式调控元件motif富集分析")
name=c("PRK4me3vsCRK4me3_up","PRK4me3vsCRK4me3_down",
       "PRK9acvsCRK9ac_up","PRK9acvsCRK9ac_down",
       "PRK27acvsCRK27ac_up","PRK27acvsCRK27ac_down")
b <- data.frame()
for (i in name) {
  a <- read.table(paste0("R1_",i,".tsv"),header=T)
  a$type <- rep(i,nrow(a))
  b <- rbind(b,a)
}
head(b)
setwd("D:/高粱/ChIP/定量/PLvsCL/差异修饰的基因的顺式调控元件motif富集分析")
name=c("PLK4me3vsCLK4me3_up","PLK4me3vsCLK4me3_down",
       "PLK9acvsCLK9ac_up","PLK9acvsCLK9ac_down",
       "PLK27acvsCLK27ac_up","PLK27acvsCLK27ac_down")
c <- data.frame()
for (i in name) {
  a <- read.table(paste0("R1_",i,".tsv"),header=T)
  a$type <- rep(i,nrow(a))
  c <- rbind(c,a)
}
head(c)
setwd("D:/高粱/ChIP/定量/PRvsCR/差异修饰的基因的顺式调控元件motif富集分析")
total <- rbind(b,c)
write.csv(total,"total_motif.csv",row.names=F)
####绘图
ge <- read.csv("total_motif.csv",header=T) 
head(ge)
topp<- ge %>% group_by(type) %>% top_n(n = 20, wt = X.TP)
top<- topp %>% group_by(type) %>% top_n(n = 5, wt = -log10(E.value))
#top <- ge[which(ge$rank < 6),]
top <- top[order(top$type,-log10(top$E.value)),]
top$symbol <- paste(top$Family_Name,top$consensus,top$DBID,sep="_")
top$symbol <- factor(top$symbol, levels=unique(top$symbol),ordered=TRUE)
head(top)
top[top == 0] <- 1e-250
top$E.value <- -log(top$E.value,10)
top <- top[which(top$E.value > 0),]
#pdf("top5_motif_enriched_stressed.pdf",width=12)
ggplot(top,mapping = aes(x=E.value,y=symbol,fill=consensus))+
  #geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  geom_bar(stat='identity',position = position_dodge(),show.legend = TRUE) +
  #labs(y = 'Number of peak region') +
  #scale_fill_manual(values=brewer.pal(10,"Paired")) +
  scale_fill_igv()+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),legend.position = "none") + 
  guides(fill = guide_legend(title = 'Type')) +
  xlab("-log10(E_value)") +
  geom_text(aes(label = X.TP), size = 3, hjust = 0.5, vjust = 1) +
  ylab(NULL) + 
  #coord_cartesian(xlim = c(0,800)) +
  facet_grid(~type,scales = "free_x")
#dev.off()
pdf("top5_motif_enriched_stressed_气泡图.pdf")
ggplot(top,aes(x=type,y=symbol))+
  geom_point(aes(size=`X.TP`,
                 color=`E.value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#ffede3",high="#e51813")+
  labs(x=NULL,y=NULL)
dev.off()

####最后用MEME的ASE做


###利用homer findMotifs.pl 查找启动子区域的motif
conda install -c bioconda homer
cd /public/home/chaohe/sorghum/chip/align/diff
cut -f 1 list.txt | while read i;
do
findMotifs.pl $i.bed.fa fasta "$i" -fasta /public/home/chaohe/db/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa  
done

###提取前十的motif，只用know_motif
cd /public/home/chaohe/sorghum/chip/align/diff
#cut -f 1 list.txt | while read i;
#do
#grep ">" "$i"/homerMotifs.all.motifs  >"$i".motif.txt
#done
#perl -p -i -e 's/,/\t/g' *.motif.txt
cut -f 1 list.txt | while read i;
do
cp "$i"/knownResults.txt  "$i".knownmotif.txt
done

###绘制p-value前10的motif
setwd("D:/高粱/ChIP/motif")
a = list.files(".")
mo <- data.frame()
for (i in 1:14) {
  fi <- read.csv(a[[i]],header=T)
  fi$motif <- paste(fi$Name,fi$Consensus,sep="_")
  fi$type <- rep(a[[i]],nrow(fi))
  fi$type <- gsub("R1_","",fi$type)
  fi$type <- gsub(".knownmotif.csv","",fi$type)
  fi <- fi[which(fi$P.value < 0.05),]
  mo <- rbind(mo,fi)
}
mo$P.value <- -log10(mo$P.value)
head(mo)
top<- mo %>% group_by(type) %>% top_n(n = 10, wt = P.value)
head(top)
write.csv(top,"top10_motif_up_down.csv")
top <- read.csv("top10_motif_up_down.csv")
#top$type <- factor(top$type,levels=c("CLKH1AZvsCRKH1AZ_up","CLK27acvsCRK27ac_up","CLK9acvsCRK9ac_up","CLK4me3vsCRK4me3_up",
#                                     "CLK27me3vsCRK27me3_up","CLK36me3vsCRK36me3_up","CLK4me2vsCRK4me2_up",
#                                     "CLKH1AZvsCRKH1AZ_down","CLK27acvsCRK27ac_down","CLK9acvsCRK9ac_down","CLK4me3vsCRK4me3_down",
#                                     "CLK27me3vsCRK27me3_down","CLK36me3vsCRK36me3_down","CLK4me2vsCRK4me2_down"),ordered=TRUE)
top$motif <- factor(top$motif,levels=unique(top$motif),ordered=TRUE)
top <- top[order(top$type),] 
pdf("top10_crcl_motif.pdf")
ggplot(top,aes(x=type,y=motif))+
  geom_point(aes(color=`P.value`),size=3)+
  theme_bw()+
  #guides(color=guide_legend(title = "-log10Pvalue")) + 
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient(low="#2193B0",high="#ff0080") +
  labs(x=NULL,y=NULL)
dev.off()

#############motif富集分析，MEME上的AME，数据库选择motif_db/CIS-BP_2.00/Sorghum_bicolor.meme
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/差异修饰的基因的顺式调控元件motif富集分析")
name=c("CLKH1AZvsCRKH1AZ_up","CLK27acvsCRK27ac_up","CLK9acvsCRK9ac_up","CLK4me3vsCRK4me3_up",
   "CLK27me3vsCRK27me3_up","CLK36me3vsCRK36me3_up","CLK4me2vsCRK4me2_up",
   "CLKH1AZvsCRKH1AZ_down","CLK27acvsCRK27ac_down","CLK9acvsCRK9ac_down","CLK4me3vsCRK4me3_down",
   "CLK27me3vsCRK27me3_down","CLK36me3vsCRK36me3_down","CLK4me2vsCRK4me2_down")
b <- data.frame()
for (i in name) {
  a <- read.table(paste0("R1_",i,".tsv"),header=T)
  a$type <- rep(i,nrow(a))
  b <- rbind(b,a)
}
head(b)
write.csv(b,"total_motif.csv",row.names=F)
####绘图
ge <- read.csv("total_motif_M.csv",header=T) 
head(ge)
topp<- ge %>% group_by(type) %>% top_n(n = 5, wt = X.TP)
top<- topp %>% group_by(type) %>% top_n(n = 5, wt = -log10(E.value))
top <- top[order(top$type,-log10(top$E.value)),]
top$symbol <- paste(top$Family_Name,top$consensus,top$DBID,sep="_")
top$symbol <- factor(top$symbol, levels=unique(top$symbol),ordered=TRUE)
head(top)
top[top == 0] <- 1e-320
top$E.value <- -log(top$E.value,10)
top <- top[which(top$E.value > 0),]
head(top)
coll=c("#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
       "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
       "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
       "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

pdf("top5_motif_barplot.pdf",width=14)
ggplot(top,mapping = aes(x=E.value,y=symbol,fill=consensus))+
  #geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  geom_bar(stat='identity',position = position_dodge(),show.legend = TRUE) +
  #labs(y = 'Number of peak region') +
  #scale_fill_manual(values=brewer.pal(10,"Paired")) +
  scale_fill_d3("category20c")+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),legend.position = "none") + 
  guides(fill = guide_legend(title = 'Type')) +
  xlab("-log10(E_value)") +
  geom_text(aes(label = X.TP), size = 3, hjust = 0.5, vjust = 1) +
  ylab(NULL) + 
  #coord_cartesian(xlim = c(0,800)) +
  facet_grid(~type)
dev.off()

pdf("top5_motif_气泡图.pdf")
ggplot(top,aes(x=type,y=symbol)) +
  geom_point(aes(size=`X.TP`,
                 color=`E.value`)) +
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  #scale_color_gradient(low="#c2e59c",high="#64b3f4") +
  #scale_color_distiller(palette = "Paired") +
  scale_color_gradient(low="#A8C0FF",high="#3F2B96")+
  labs(x=NULL,y=NULL)
dev.off()

ggplot(top,aes(x=type,y=symbol))+
  geom_point(aes(size=`X.TP`,
                 color=`E.value`))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#ffede3",high="#e51813")+
  labs(x=NULL,y=NULL)


ggplot(top,aes(x=type,y=symbol))+
  geom_point(aes(color=`E.value`),size=3)+
  theme_bw()+
  #guides(color=guide_legend(title = "-log10Pvalue")) + 
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  #scale_color_gradient(low="#c2e59c",high="#64b3f4") +
  #scale_color_distiller(palette = "Paired") +
  scale_color_stepsn(breaks = c(0,1e-316,1e-55),
                     colors = brewer.pal(12,"Paired")) +
labs(x=NULL,y=NULL)

#40E0D0
#40e0d0
#FF8C00
#ff8c00
#FF0080
#ff0080

#6A82FB
#6a82fb


#2193b0

###H3K27me3 peak定量,还是要合并1kb以内有重叠的peak，最终用合并相距1kb的peak
#setwd("D:/高粱/H3K27me3差异peak分析/只合并peak")
setwd("D:/高粱/H3K27me3差异peak分析/合并1kb以内的")
#构建tpm矩阵
a <- list.files(path="D:/高粱/H3K27me3差异peak分析/只合并peak", pattern="_counts_subread.txt", all.files=F, full.names=F)
dir = paste("./",a,sep="")     
n = length(dir)  
merge.data = read.table(file = dir[1],header=TRUE,dec = ".")
merge.data$count <- merge.data[,7]
merge.data <- merge.data[,-7]
merge.data$type <- rep(dir[1],nrow(merge.data))
for (i in 2:n){
  new.data = read.table(file = dir[i], header=TRUE, dec = ".",skip=1)
  new.data$count <- new.data[,7]
  new.data <- new.data[,-7]
  new.data$type <- rep(dir[i],nrow(new.data))
  merge.data = rbind(merge.data,new.data)
}
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
#计算TPM值
total_tpm <- data.frame()
for (i in 1:n)  {
  pea<-merge.data[which(merge.data$type == dir[i]),]
  head(pea)
  count <- as.data.frame(pea$count)
  Length <- as.data.frame(pea$Length)
  tpmm <- apply(count, 2, function(x) tpm(x, Length))
  peak_tpms <-cbind(pea,tpmm)
  peak_tpms$tpm <- peak_tpms[,9]
  peak_tpms <- peak_tpms[,-9]
  total_tpm <- rbind(total_tpm,peak_tpms)
}
#长变宽
total_tpmR1 <- total_tpm[,-7]
toal <- tidyr::spread(data=total_tpmR1,key=type,value=tpm)
head(toal)
write.csv(toal,"k273_tpm_matrix.csv")
###构建count矩阵
head(merge.data)
toal_count <- tidyr::spread(data=merge.data,key=type,value=count)
head(toal_count)
write.csv(toal_count,"k273_count_matrix.csv")
###差异peak分析，用Deseq2的方法做
library(DESeq2)
library(dplyr)
## 导入TPM数据矩阵
#countdata <- read.csv("H3K4me3_count.csv", row.names = 1)
countdata <- read.csv("k273_count_matrix.csv", row.names = 1)
row.names(countdata) <- countdata[,1]
countdata <- countdata[,-1:-6]
head(countdata)
tpmdata <- read.csv("k273_tpm_matrix.csv", row.names = 1)
row.names(tpmdata) <- tpmdata[,1]
tpmdata <- tpmdata[,-1:-6]
## 过滤在所有重复样本中小于1的基因
countdata = countdata[rowMax(as.matrix(countdata)) > 1,]
#导入样本注释信息
coldata  <- read.csv("../col_data.csv",row.names = 1)
coldata$sample <- coldata$condition
head(coldata)
#检查数据Counts文件与coldata数据是否匹配
all(rownames(coldata) %in% colnames(countdata))  
all(rownames(coldata) == colnames(countdata))
#PCA分析
##差异分析
# 制作差异矩阵
##CRK27me3vsCLK27me3
a1 <- countdata[,c(3:4,1:2)]
c1 <- tpmdata[,c(3:4,1:2)]
head(a1)
b1 <- coldata[c(3:4,1:2),]
head(b1)
dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
# 过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  
# 差异比较
dep <- DESeq(dds)
res <- results(dep,independentFiltering=TRUE) #避免padj为NA
diff = res
#diff <- na.omit(diff)  ## 去除NA
#加上tpm值
diff <- as.data.frame(diff)
diff_tpm <- merge(diff,c1,by="row.names",all.x=F)
#diff_tpm <- diff_tpm[which(diff_tpm$padj < 0.05),]
head(diff_tpm)
####调取peak的具体位置,excel中操作
#peak <- read.table("CRCL_H3K27me3_merged.peaks.bed")
#colnames(peak) <- c("chr","start","end","Row.names")
#head(peak)
#diff_tpm_peak <- merge(diff_tpm,peak,by="Row.names",all=T)
#dim(diff_tpm_peak)
write.csv(diff_tpm,"CRK27me3vsCLK27me3_diff_all.csv")  # 导出所有的差异文件
#加上peak的位置,有些位置名字太长，只能再excel中做合并
peak <- read.table("CRCL_H3K27me3_merged.peaks.bed",header=F)
peak$Row.names <- peak$V4
head(peak)
head(diff_tpm)
diff_tpm_peak <- merge(diff_tpm,peak,by="Row.names",all=F)
head(diff_tpm_peak)
##绘制差异peak长度的分布
#diff_tpm_peak <- read.csv("CRK27me3vsCLK27me3_diff_allR2.csv",header=T,row.names = 1)
diff_tpm_peak <- read.csv("CRK27me3vsCLK27me3_diff_allR1.csv",header=T,row.names = 1)
head(diff_tpm_peak)
diff_tpm_peak$length <- diff_tpm_peak$end - diff_tpm_peak$start
diff_tpm_peak$loglength <- log10(diff_tpm_peak$end - diff_tpm_peak$start)
ggplot(diff_tpm_peak, mapping = aes(x=loglength)) +
  geom_density(size=1.5,key_glyph = draw_key_path) +
  #scale_colour_manual(values=brewer.pal(9,"Paired")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(size=1,colour = "black"),
        axis.text.x = element_text(color="black", 
                                   size=12),
        axis.text.y = element_text(color="black", 
                                   size=12)) +
  annotation_logticks()
l00 <- diff_tpm_peak[which(diff_tpm_peak$length < 1000),]
l0 <- diff_tpm_peak[which(diff_tpm_peak$length >= 1000 & diff_tpm_peak$length < 5000),]
l1 <- diff_tpm_peak[which(diff_tpm_peak$length >= 5000 & diff_tpm_peak$length < 10000),]
l2 <- diff_tpm_peak[which(diff_tpm_peak$length >= 10000 & diff_tpm_peak$length < 20000),]
l3 <- diff_tpm_peak[which(diff_tpm_peak$length >= 20000 & diff_tpm_peak$length < 40000),]
l4 <- diff_tpm_peak[which(diff_tpm_peak$length >= 40000 & diff_tpm_peak$length < 60000),]
l5 <- diff_tpm_peak[which(diff_tpm_peak$length >= 60000 & diff_tpm_peak$length < 80000),]
l6 <- diff_tpm_peak[which(diff_tpm_peak$length >= 80000 & diff_tpm_peak$length < 100000),]
l7 <- diff_tpm_peak[which(diff_tpm_peak$length >= 100000),]
l00$group <- rep("<1",nrow(l00))
l0$group <- rep("1-5",nrow(l0))
l1$group <- rep("5-10",nrow(l1))
l2$group <- rep("10-20",nrow(l2))
l3$group <- rep("20-40",nrow(l3))
l4$group <- rep("40-60",nrow(l4))
l5$group <- rep("60-80",nrow(l5))
l6$group <- rep("80-100",nrow(l6))
l7$group <- rep(">100kb",nrow(l7))
length <- rbind(l00,l0,l1,l2,l3,l4,l5,l6,l7)
length <- length[which(length$padj < 0.05 & (length$log2FoldChange > 0.75 | length$log2FoldChange < -0.75 )),]
head(length)
total_number <- aggregate(length$chr,list(length$group),length)
total_number$Group.1 <- factor(total_number$Group.1,levels=c("<1","1-5",
                                                             "5-10","10-20","20-40",
                                                             "40-60","60-80","80-100",">100"),ordered=TRUE)
head(total_number)
pdf("H3K27me3-DEP.pdf",width=4.5,height=4,family="ArialMT")
ggplot(total_number,mapping = aes(x=Group.1,y=x,fill=Group.1))+
  geom_bar(stat='identity',position='dodge',show.legend = TRUE) +
  labs(x = 'Length (kb) of H3K27me3-DEPs', y = 'Number of regions') +
  #scale_fill_manual(values=brewer.pal(8,"Dark2")) +
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) + 
  ggtitle("H3K27me3-DEPs") +
  guides(fill = guide_legend(title = 'Type')) +
  geom_text(aes(label = x), size = 3, vjust =-0.2, position = position_dodge(0.9))
dev.off()
####将不差异变化的peak加上画堆积柱形图
write.csv(length,"DEP_H3K27me3.csv")
length <- length[,c(1:3,16)]
length$type <- rep("DMP",nrow(length))
head(length)
#在excel中将不变化的peak挑出来
non <- read.csv("notDEP_H3K27me3.csv",header=T,row.names=1)
non <- non[,-4]
colnames(non) <- c("chr","start","end")
non$length <- non$end - non$start
n00 <- non[which(non$length < 1000),]
n0 <- non[which(non$length >= 1000 & non$length < 5000),]
n1 <- non[which(non$length >= 5000 & non$length < 10000),]
n2 <- non[which(non$length >= 10000 & non$length < 20000),]
n3 <- non[which(non$length >= 20000 & non$length < 40000),]
n4 <- non[which(non$length >= 40000 & non$length < 60000),]
n5 <- non[which(non$length >= 60000 & non$length < 80000),]
n6 <- non[which(non$length >= 80000 & non$length < 100000),]
n7 <- non[which(non$length >= 100000),]
n00$group <- rep("<1",nrow(n00))
n0$group <- rep("1-5",nrow(n0))
n1$group <- rep("5-10",nrow(n1))
n2$group <- rep("10-20",nrow(n2))
n3$group <- rep("20-40",nrow(n3))
n4$group <- rep("40-60",nrow(n4))
n5$group <- rep("60-80",nrow(n5))
n6$group <- rep("80-100",nrow(n6))
n7$group <- rep(">100kb",nrow(n7))
non <- rbind(n00,n0,n1,n2,n3,n4,n5,n6,n7)
non$type <- rep("None",nrow(non))
non <- non[,-4]
colnames(non) <- colnames(length)
head(non)
head(length)
to <- rbind(length,non)
head(to)
total_number <- aggregate(to$chr,list(to$group,to$type),length)
total_number$Group.1 <- factor(total_number$Group.1,levels=c("<1","1-5",
                                                             "5-10","10-20","20-40",
                                                             "40-60","60-80","80-100",">100"),ordered=TRUE)
head(total_number)
pdf("H3K27me3-DEPR1.pdf",width=4.5,height=4,family="ArialMT")
ggplot(total_number,mapping = aes(x=Group.1,y=x,fill=Group.2))+
  geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  labs(x = 'Length (kb) of H3K27me3-DEPs', y = 'Number of regions') +
  #scale_fill_manual(values=brewer.pal(8,"Paired")[c(6,2)]) +
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) + 
  ggtitle("H3K27me3-DEPs") +
  guides(fill = guide_legend(title = 'Type')) +
  geom_text(aes(label = x), size = 2, vjust =-0.2, position = position_stack(vjust = 0.5))
dev.off()

#####绘制长peak的circos图
#定义长于10kb的peak为长DEP，提取出来
diff_tpm_peak <- read.csv("CRK27me3vsCLK27me3_diff_allR1.csv",header=T)
diff_tpm_peak <- diff_tpm_peak[which(diff_tpm_peak$padj < 0.05 & (diff_tpm_peak$log2FoldChange > 0.75 | diff_tpm_peak$log2FoldChange < -0.75 )),]
head(diff_tpm_peak)
diff_tpm_peak$length <- diff_tpm_peak$end - diff_tpm_peak$start
#ld <- diff_tpm_peak[which(diff_tpm_peak$length >= 10000),c(2:4,6)]
ld <- diff_tpm_peak[,c(2:4,6)]
head(ld)
write.table(ld,"h3k27me3_DEPs.bed", row.names = F)
#导出显著差异的长peak
ld <- diff_tpm_peak[which(diff_tpm_peak$length >= 10000),c(2:4,6)]
head(ld)
write.table(ld,"h3k27me3_DEPR2.bed", row.names = F,sep="\t")
###提取落在差异peak的区域的值
cd /public/home/chaohe/sorghum/chip/align/RPKM
perl -p -i -e 's/^M//g' h3k27me3_DEPs.bed
sort -k1,1 -k2n,2 h3k27me3_DEPs.bed | awk '{print $2"\t"$3"\t"$4}' >h3k27me3_DEP.bed
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEP.bed -b CLK27me3_rep1.bedGraph | awk '{if ($8 == 0) print $4"\t"$5"\t"$6"\t"$7}' >CLK27me3_rep1R1.bedGraph
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEP.bed -b CLK27me3_rep2.bedGraph | awk '{if ($8 == 0) print $4"\t"$5"\t"$6"\t"$7}' >CLK27me3_rep2R1.bedGraph
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEP.bed -b CRK27me3_rep1.bedGraph | awk '{if ($8 == 0) print $4"\t"$5"\t"$6"\t"$7}' >CRK27me3_rep1R1.bedGraph
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEP.bed -b CRK27me3_rep2.bedGraph | awk '{if ($8 == 0) print $4"\t"$5"\t"$6"\t"$7}' >CRK27me3_rep2R1.bedGraph
#qval < 0.05 log2FC大于0.75的长差异DEPs,length > 10kb
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEPR2.bed -b CLK27me3_rep1.bedGraph | awk '{if ($9 == 0) print $5"\t"$6"\t"$7"\t"$8}' >CLK27me3_rep1R2.bedGraph
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEPR2.bed -b CLK27me3_rep2.bedGraph | awk '{if ($9 == 0) print $5"\t"$6"\t"$7"\t"$8}' >CLK27me3_rep2R2.bedGraph
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEPR2.bed -b CRK27me3_rep1.bedGraph | awk '{if ($9 == 0) print $5"\t"$6"\t"$7"\t"$8}' >CRK27me3_rep1R2.bedGraph
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEPR2.bed -b CRK27me3_rep2.bedGraph | awk '{if ($9 == 0) print $5"\t"$6"\t"$7"\t"$8}' >CRK27me3_rep2R2.bedGraph
#调取基因的表达量
bedtools closest -D ref -t all -mdb all -a h3k27me3_DEPR2.bed -b /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed | awk '{if ($11 == 0) print $5"\t"$6"\t"$7"\t"$8}' > DEP_gene.bed
gee <- read.csv("DEP_gene.bed")
#circos.clear()
#画图
library(circlize)
library(stringr)
library(grid)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
#setwd("/public/home/chaohe/sorghum/chip/align/RPKM")
circos.clear()
setwd("D:/高粱/H3K27me3差异peak分析/合并1kb以内的")
seq_stat<-read.csv("D:/高粱/ChIP/CIRCOS/3chromosome_length.csv",header=T,row.names=1)
seq_stat$seq_start <- as.numeric(seq_stat$seq_start)
seq_stat$seq_end <- as.numeric(seq_stat$seq_end)
circle_size = unit(1, "snpc")
circos.par(gap.degree = 2)
circos.genomicInitialize(seq_stat, plotType = 'axis')
circos.track(
  ylim = c(0, 0.5), track.height = 0.08, bg.border = NA, bg.col =rep(c("#CAB2D6", "#B2DF8A"), 7),
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    seq_ID = CELL_META$sector.index
    circos.text(mean(xlim), mean(ylim), seq_ID, cex = 0.7, col = 'black', facing = 'inside', niceFacing = FALSE)
  } )
#基因密度
density<-read.table("D:/高粱/ChIP/CIRCOS/count_gene_100kb.txt",header=T)
color_assign <- colorRamp2(breaks = c(1, 10, 21), 
                           col = c('#FEE8C8', '#FC8D59', '#D7301F'))
circos.genomicTrackPlotRegion(
  density, track.height = 0.12, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )
#long DEPs
ld <- read.table("h3k27me3_DEPR2.bed",header=T)
head(ld)
colnames(ld) <- c("seq_chr","seq_start","seq_end","value")
circos.genomicTrack(ld, track.height = 0.08, bg.col = NA, bg.border = NA,
                    ylim=c(-6,6),
                    panel.fun = function(region,value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#377EB8',
                                         border = ifelse(value[[1]] > 0, "#ef5d77", "#377EB8"))
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    } )
#gene_expression
gee <- read.table("DEP_gene.bed")
gee <- gee[,c(1,2,3,5)]
head(gee)
circos.genomicTrack(
  gee,  
  #stack=TRUE,
  track.height = 0.08, bg.col = NA, bg.border = NA,
  ylim=c(-14,14),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#377EB8',
                       border = ifelse(value[[1]] > 0, "#ef3b36", "#4568dc"))
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
  } )
#CRK27me3_rep1
cr1 <- read.table("CRK27me3_rep1R2.bedGraph",sep="\t")
head(cr1)
colnames(cr1) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cr1,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#D95F02',border = '#D95F02')
  } )
cr2 <- read.table("CRK27me3_rep2R2.bedGraph")
colnames(cr2) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cr2,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#7570B3',border = '#7570B3')
  } )
cl1 <- read.table("CLK27me3_rep1R2.bedGraph")
colnames(cl1) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cl1,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#66A61E',border = '#E78AC3')
  } )
cl2 <- read.table("CLK27me3_rep2R2.bedGraph")
head(cl2)
colnames(cl2) <- c("seq_ID","seq_start","seq_end","value1")
circos.genomicTrack(
  cl2,  
  #stack=TRUE,
  track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  #ylim=c(2,90),
  panel.fun = function(region,value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = -1, lwd = 0.02, col ='#E7298A',border = '#FFD92F')
  } )
####选特殊的区域画图
conda activate py36
module load pyGenomeTracks/3.5
cd /public/home/chaohe/sorghum/chip/align/RPKM
pyGenomeTracks --tracks DEP.ini --region  8:53870358-54008184  -o  H3K27me3-DEP.pdf
pyGenomeTracks --tracks DEP.ini --region  5:70870000-71000000  -o  5H3K27me3-DEP.pdf

#####绘制chromMM生成的数据，热图，需要bar，用重跑后的结果来做，就是坐标倒过来，1-8变成8-1
setwd("D:/高粱/chromMM/newresult/")
state <- read.table("emissions_8.txt",header=T,sep="\t")
head(state)
row.names(state) <- state[,1]
state <- state[,-1]
head(state)
#绘制热图
library(circlize )
library(pheatmap)
library(ComplexHeatmap)
ms <- as.data.frame(colnames(state))
#ac = data.frame(group=str_split(ms,'', simplify = T)[,1])
rownames(ms) = colnames(state)
pdf("chromMMR1.pdf",family="ArialMT")
bk <- c(seq(0,0.19,by=0.01),seq(0.2,1,by=0.01))
p1 <- pheatmap::pheatmap(state, 
                   cluster_rows = F,
                   cluster_cols =  F,
                   #annotation_col = ms,
                   display_numbers=T,
                   color = c(colorRampPalette(colors = c("white","white"))(20), colorRampPalette(colors = c("#E68DB7","#922376"))(101)),
                   legend_breaks=seq(0,1,2),
                   breaks=bk)
dev.off()
#####统计每个state的基因的数量,重叠率大于40%
cd /public/home/chaohe/ChromHMM/OUTPUTSAMPLE
#sort -k1,1 -k2n,2 stage1_8_segments.bed | bedtools intersect -wa -wb -a /public/home/chaohe/sorghum/RNA-seq/db/geneR1.bed -b - -f 0.5 >stage1.txt
#sort -k1,1 -k2n,2 stage2_8_segments.bed | bedtools intersect -wa -wb -a /public/home/chaohe/sorghum/RNA-seq/db/geneR1.bed -b - -f 0.5 >stage2.txt
###新版
cd /public/home/chaohe/ChromHMM/OUTPUTSAMPLE
sort -k1,1 -k2n,2 stage1_8_segments.bed | bedtools closest -D ref -t all -mdb all -a - -b /public/home/chaohe/sorghum/RNA-seq/db/geneR1.bed | awk '{if($9 ==0) print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > stage1.txt
sort -k1,1 -k2n,2 stage2_8_segments.bed | bedtools closest -D ref -t all -mdb all -a - -b /public/home/chaohe/sorghum/RNA-seq/db/geneR1.bed | awk '{if($9 ==0) print $1"\t"$2"\t"$3"\t"$4"\t"$8}'  > stage2.txt

s1 <- read.table("stage1.txt")
s2 <- read.table("stage2.txt")
#s1_number <- aggregate(s1$V4,list(s1$V8),length)
s1_number <- aggregate(s1$V5,list(s1$V4),length)
s2_number <- aggregate(s2$V5,list(s2$V4),length)
head(s1_number)
head(s2_number)
#添加表达量
tpm <- read.csv("D:/高粱/RNA-seq/final_RNAseq/all_tpm_expression.csv")
tpm$leaf <- rowMeans(as.matrix(tpm[,2:4]))
tpm$root <- rowMeans(as.matrix(tpm[,5:7]))
tpm <- tpm[,c(1,14,15)]
tpm$V5 <- tpm[,1]
head(tpm)
s1_tpm <- merge(s1,tpm,by="V5",all.x=F)
s2_tpm <- merge(s2,tpm,by="V5",all.x=F)
###保存每个cluster的基因集及其表达量
write.csv(s1_tpm,"s1_tpm.csv")
write.csv(s2_tpm,"s2_tpm.csv")
#准备作图数据
s1_tpm <- s1_tpm[,c(5,7)]
s1_tpm$leaf <- log2(s1_tpm$leaf+1)
s2_tpm <- s2_tpm[,c(5,8)]
s2_tpm$root <- log2(s2_tpm$root+1)
s1_tpm$V4 <- gsub("E","",s1_tpm$V4)
s2_tpm$V4 <- gsub("E","",s2_tpm$V4)
s1_tpm$V4 <- factor(s1_tpm$V4 , levels=c("8","7","6","5",
                                         "4","3","2","1"), ordered=TRUE)
s2_tpm$V4 <- factor(s2_tpm$V4 , levels=c("8","7","6","5",
                                         "4","3","2","1"), ordered=TRUE)
#绘制柱形图可视化
p1 <- ggboxplot(s1_tpm, x = "V4", y ="leaf", fill = "V4",
                #add = "mean_sd", error.plot = "crossbar") +
                #add = "boxplot") +
                outlier.shape = NA,
                bxp.errorbar=T, 
                add = "none") +
  scale_fill_manual(values=c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a","#d62728","#ff9896","#6A3D9A","#CAB2D6")) +
  ggtitle("Leaf") + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) + ylab("log2 (TPM+1)")

p2 <- ggboxplot(s2_tpm, x = "V4", y ="root", fill = "V4",
                #add = "mean_sd", error.plot = "crossbar") +
                #add = "boxplot") +
                outlier.shape = NA,
                bxp.errorbar=T, 
                add = "none") +
  scale_fill_manual(values=c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a","#d62728","#ff9896","#6A3D9A","#CAB2D6")) +
  ggtitle("Root") + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) + ylab("log2 (TPM+1)")
grid.arrange(p1,p2,
             nrow=1,ncol=2)     %>%  ggsave("histone_cluterR1.pdf",.,width=100,height=100, units="mm")
####绘制overlap enrichmet热图
o1 <- read.table("stage1_8_overlap_new.txt",header=T,row.names=1,sep="\t")
o1 <- o1[-9,c(1,2,3)]
colnames(o1) <- c("genome","exon","gene")
head(o1)
library(ComplexHeatmap)
ms <- as.data.frame(colnames(o1))
#ac = data.frame(group=str_split(ms,'', simplify = T)[,1])
rownames(ms) = colnames(o1)
pdf("s1_overlap-1.pdf")
bk <- c(seq(0,0.19,by=0.01),seq(0.2,1,by=0.01))
q1 <- pheatmap::pheatmap(o1, 
                         cluster_rows = F,
                         cluster_cols =  F,
                         #annotation_col = ms,
                         display_numbers=T,
                         scale = "column",
                         color = colorRampPalette(c("#fffcdc", "#D9A7C7"))(100),
                         legend_breaks=seq(0,1,2),
                         breaks=bk)

dev.off()
pdf("s1_overlap-2.pdf")
q2 <- pheatmap::pheatmap(o1, 
                         cluster_rows = F,
                         cluster_cols =  F,
                         #scale = "column",
                         #annotation_col = ms,
                         #colorRampPalette(colors = c("blue","white","red"))(100)
                         color = colorRampPalette(c("#fffcdc", "#D9A7C7"))(100),
                         display_numbers=T,
                         breaks=bk)
dev.off()
o2 <- read.table("stage2_8_overlap_new.txt",header=T,row.names=1,sep="\t")
o2 <- o2[-9,c(1,2,3)]
colnames(o2) <- c("genome","exon","gene")
head(o2)
library(ComplexHeatmap)
ms <- as.data.frame(colnames(o2))
#ac = data.frame(group=str_split(ms,'', simplify = T)[,1])
rownames(ms) = colnames(o2)
pdf("s2_overlap-1.pdf")
q3 <- pheatmap::pheatmap(o2, 
                         cluster_rows = F,
                         cluster_cols =  F,
                         scale = "column",
                         #annotation_col = ms,
                         #colorRampPalette(colors = c("blue","white","red"))(100)
                         color = colorRampPalette(c("#fffcdc", "#D9A7C7"))(100),
                         legend_breaks=seq(0,1,2),
                         breaks=bk)
dev.off()
pdf("s2_overlap-2.pdf")
q4 <- pheatmap::pheatmap(o2, 
                         cluster_rows = F,
                         cluster_cols =  F,
                         #scale = "column",
                         #annotation_col = ms,
                         #colorRampPalette(colors = c("blue","white","red"))(100)
                         color = colorRampPalette(c("#fffcdc", "#D9A7C7"))(100),
                         display_numbers=T,
                         breaks=bk)
dev.off()


##绘制PCA聚类图
setwd("D:/高粱/ChIP/PCA")
matrix <- read.table("matrix.txt",header=T)
matrix <- matrix[,c(-1,-30)]
head(matrix)
colnames(matrix) <- c('CLK27ac_rep1','CLK27ac_rep2', 'CLK27me3_rep1', 'CLK27me3_rep2','CLK36me3_rep1', 'CLK36me3_rep2', 'CLK43_rep1',
                      'CLK43_rep2', 'CLK4me2_rep1', 'CLK4me2_rep2','CLK9ac_rep1', 'CLK9ac_rep2', 'CLKH2AZ_rep1','CLKH2AZ_rep2', 'CRK27ac_rep1', 'CRK27ac_rep2',
                      'CRK27me3_rep1', 'CRK27me3_rep2','CRK36me3_rep1', 'CRK36me3_rep2', 'CRK43_rep1','CRK43_rep2', 'CRK4me2_rep1', 'CRK4me2_rep2',
                      'CRK9ac_rep1', 'CRK9ac_rep2', 'CRKH2AZ_rep1','CRKH2AZ_rep2', 'PLK27ac_rep1', 'PLK27ac_rep2', 'PLK27me3_rep1',
                      'PLK27me3_rep2', 'PLK36me3_rep1','PLK36me3_rep2', 'PLK43_rep1', 'PLK43_rep2','PLK4me2_rep1', 'PLK4me2_rep2', 'PLK9ac_rep1',
                      'PLK9ac_rep2', 'PLKH2AZ_rep1', 'PLKH2AZ_rep2','PRK27ac_rep1', 'PRK27ac_rep2', 'PRK27me3_rep1','PRK27me3_rep2', 'PRK36me3_rep1',
                      'PRK36me3_rep2', 'PRK43_rep1', 'PRK43_rep2','PRK4me2_rep1', 'PRK4me2_rep2', 'PRK9ac_rep1','PRK9ac_rep2', 'PRKH2AZ_rep1', 'PRKH2AZ_rep2')
#绘制热图
M=cor(matrix,method = "pearson")
pdf("ChIP_correlation.pdf")
pheatmap::pheatmap(M,display_numbers=F)
dev.off()


##全基因组层面展示所有数据，用信号无法解决分辨率的问题，改用peak文件，代码在后面
cd /public/home/chaohe/sorghum/chip/align/circos
module load Perl/5.26.1
module load Circos/0.69-8
##准备数据
#基因密度文件
bedops --chop 10000 --stagger 10000 -x  /public/home/chaohe/sorghum/RNA-seq/db/Sorghum_bicolor.genome_table.txt > gene_10kb.bed
perl -p -i -e 's/chr//g' gene_10kb.bed
bedmap --echo --count --delim '\t'  gene_10kb.bed  /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed > gene_10kb_density.txt
#perl -p -i -e 's/chr//g' gene_10kb_density.txt
#准备基因表达量文件
setwd("D:/高粱/RNA-seq/final_RNAseq")
tpm <- read.csv("D:/高粱/RNA-seq/final_RNAseq/all_tpm_expression.csv")
head(tpm)
row.names(tpm) <- tpm[,1]
CL <- tpm[,2:4]
CR <- tpm[,5:7]
CL$average <- rowMeans(CL)
CR$average <- rowMeans(CR)
head(CL)
head(CR)
###CR CL添加基因的位置信息
pos <- read.table("geneR1.bed",sep="\t",row.names=5)
head(pos)
cl_pos <- merge(CL,pos,by="row.names",all=F)
cr_pos <- merge(CR,pos,by="row.names",all=F)
head(cl_pos)
head(cr_pos)
cl_pos <- cl_pos[,c(6,7,8,5)]
cr_pos <- cr_pos[,c(6,7,8,5)]
cl_pos$type <- rep("id=CL",nrow(cl_pos))
cr_pos$typpe <- rep("id=CR",nrow(cr_pos))
colnames(cl_pos) <- colnames(cr_pos)
rna <- rbind(cl_pos,cr_pos)
rna <- rna[order(rna$V1,rna$V2),]
head(rna)
write.table(rna,"RNA.txt",sep="\t",row.names=F)
perl -p -i -e 's/\"//g' RNA.txt
perl -p -i -e 's/^M//g' RNA.txt
###生成组蛋白修饰信号文件
conda install -c bioconda ucsc-bigwigtobedgraph
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"K4me2_rpkm.bw "$i"K4me2_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"K27me3_rpkm.bw "$i"K27me3_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"K27ac_rpkm.bw "$i"K27ac_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"K43_rpkm.bw "$i"K43_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"K36me3_rpkmR1.bw "$i"K36me3_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"KH2AZ_rpkm.bw "$i"KH2AZ_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "for i in CL CR; do bigWigToBedGraph ../"$i"K9ac_rpkm.bw "$i"K9ac_rpkm.bedGraph;done"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "bigWigToBedGraph /public/home/chaohe/sorghum/mnase/nucleosome/CRMN.like.bw CRMN.bedGraph"
bsub -J k42 -n 10 -o %J.out -e %J.err -q normal "bigWigToBedGraph /public/home/chaohe/sorghum/mnase/nucleosome/CLMN.like.bw CLMN.bedGraph"
#将基因组划分为10kb的区域文件
bedops --chop 10000 --stagger 10000 -x  /public/home/chaohe/sorghum/RNA-seq/db/Sorghum_bicolor.genome_table.txt > gene_10kb.bed
perl -p -i -e 's/chr//g' gene_10kb.bed
#计算10kb内的组蛋白修饰信号
awk '{if($1 != 10) print $0}' gene_10kb.bed > gene_10kb1.bed
awk '{if($1 == 10) print $0}' gene_10kb.bed > gene_10kb2.bed
awk '{print $1"\t"$2"\t"$3"\t.\t"$4}'  CRMN.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb.bed -  > CRMN_10k.bedgraph
awk '{print $1"\t"$2"\t"$3"\t.\t"$4}'  CLMN.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb.bed -  > CLMN_10k.bedgraph
i="CR"
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K4me2_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K4me2_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K27me3_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K27ac_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K43_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K43_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K36me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K36me3_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"KH2AZ_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"KH2AZ_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K9ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K9ac_10k.bedgraph1
i="CL"
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K4me2_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K4me2_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K27me3_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K27ac_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K43_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K43_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K36me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K36me3_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"KH2AZ_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"KH2AZ_10k.bedgraph1
awk '{if($1 != 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K9ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb1.bed -  > "$i"K9ac_10k.bedgraph1
i="CR"
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K4me2_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K4me2_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K27me3_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K27ac_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K43_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K43_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K36me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K36me3_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"KH2AZ_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"KH2AZ_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K9ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K9ac_10k.bedgraph2
i="CL"
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K4me2_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K4me2_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K27me3_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K27ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K27ac_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K43_rpkm.bedGraph | bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K43_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K36me3_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K36me3_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"KH2AZ_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"KH2AZ_10k.bedgraph2
awk '{if($1 == 10)print $1"\t"$2"\t"$3"\t.\t"$4}'  "$i"K9ac_rpkm.bedGraph |  bedmap --echo --sum --delim '\t'  gene_10kb2.bed -  > "$i"K9ac_10k.bedgraph2
for i in K4me2 K27me3 K27ac K43 K36me3 KH2AZ K9ac;
do
cat CR"$i"_10k.bedgraph1 CR"$i"_10k.bedgraph2 > CR"$i"_10k.bedgraph
cat CL"$i"_10k.bedgraph1 CL"$i"_10k.bedgraph2 > CL"$i"_10k.bedgraph
done
#rm "$i"K4me2_10k.bedgraph "$i"K27me3_10k.bedgraph "$i"K27ac_10k.bedgraph "$i"K43_10k.bedgraph "$i"K36me3_10k.bedgraph "$i"KH2AZ_10k.bedgraph "$i"K9ac_10k.bedgraph
####可视化
circos -conf circosR2.conf
####高保真压缩svg
module load nodejs/8.11.3
npm install -g svgo 
svgo circosR1.svg
#svg转pdf
module load Python/3.8.6
pip3 install cairosvg
python
import cairosvg
cairosvg.svg2pdf(url='circos.svg', write_to='circos_heatmap.pdf')
######################################################################无法解决分辨率的问题，分开做###################
####可视化
circos -conf circos-1.conf


###查看拟南芥的基因密度
cd /public/home/chaohe/sorghum/chip/align/circos
cp /public/home/chaohe/sorghum/other_species/ara_genome_table.txt ara_genome_table.txt
sed 's/^/chr&/g' -i ara_genome_table.txt
awk '{print $1"\t0\t"$2}' ara_genome_table.txt | bedops --chop 1000000 --stagger 1000000 -x  - > aragene_1M.bed
perl -p -i -e 's/chr//g' aragene_1M.bed
bedmap --echo --count --delim '\t'  aragene_1M.bed  /public/home/zhruan/ATAC-aradopisis/final/gene.bed > aragene_1M_density.txt
#高粱
bedops --chop 1000000 --stagger 1000000 -x  /public/home/chaohe/sorghum/RNA-seq/db/Sorghum_bicolor.genome_table.txt > gene_1M.bed
perl -p -i -e 's/chr//g' gene_1M.bed
bedmap --echo --count --delim '\t'  gene_1M.bed  /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed > gene_1M_density.txt
#水稻
cp /public/home/chaohe/sorghum/other_species/rice_genome_table.txt rice_genome_table.txt
sed 's/^/chr&/g' -i rice_genome_table.txt
awk '{print $1"\t0\t"$2}' rice_genome_table.txt | bedops --chop 1000000 --stagger 1000000 -x  - > ricegene_1M.bed
perl -p -i -e 's/chr//g' ricegene_1M.bed
bedmap --echo --count --delim '\t'  ricegene_1M.bed  /public/home/zhruan/ATAC-rice/final/gene.bed > ricegene_1M_density.txt






#####单倍型分析
#计算p-value
#提前准备好基因型文件，gene.txt 表型文件  NL18r2.txt
#准备NL18r2.txt
cd /public/home/tllu/enhancer_snp
awk '{print "SL\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >SL/NL18r0.txt
awk '{print "KW\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >KW/NL18r1.txt
awk '{print "KN\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >KN/NL18r2.txt
awk '{print "FSP_C_2019\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >FSP_C_2019/NL18r6.txt
awk '{print "FSPNP_C_2019\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >FSPNP_C_2019/NL18r5.txt
awk '{print "FSP_C_2020\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >FSP_C_2020/NL18r4.txt
awk '{print "FSPNP_C_2020\t.\t.\t.\t.\t.\t"$9"\t."}' gene.txt >FSPNP_C_2020/NL18r3.txt
#运行程序
cd /public/home/tllu/enhancer_snp
module load Python/3.8.6
module load plink/1.9
python NL18.py
bsub  -J sl -n 20 -o sl.out -e sl.err -q smp -R "rusage[mem=300GB]" "python NL18.py"
cd /public/home/tllu/enhancer_snp/KW
bsub  -J sl -n 20 -o sl.out -e sl.err -q smp -R "rusage[mem=300GB]" "python NL18r1.py"
cd /public/home/tllu/enhancer_snp/KN
bsub  -J KN -n 20 -o KN.out -e KN.err -q smp -R "rusage[mem=300GB]" "python NL18r2.py"
cd /public/home/tllu/enhancer_snp/FSPNP_C_2020
bsub  -J FSPNP_C_2020 -n 20 -o FSPNP_C_2020.out -e FSPNP_C_2020.err -q smp -R "rusage[mem=300GB]" "python NL18r3.py"
cd /public/home/tllu/enhancer_snp/FSP_C_2020
bsub  -J FSP_C_2020 -n 20 -o FSP_C_2020.out -e FSP_C_2020.err -q smp -R "rusage[mem=300GB]" "python NL18r4.py"
cd /public/home/tllu/enhancer_snp/FSPNP_C_2019
bsub  -J FSPNP_C_2019 -n 20 -o FSPNP_C_2019.out -e FSPNP_C_2019.err -q smp -R "rusage[mem=300GB]" "python NL18r5.py"
cd /public/home/tllu/enhancer_snp/FSP_C_2019
bsub  -J FSP_C_2019 -n 20 -o FSP_C_2019.out -e FSP_C_2019.err -q smp -R "rusage[mem=300GB]" "python NL18r6.py"


####看FigS5 A和C 两类基因在根和叶中其他修饰变化的数量（做个柱状图）
setwd("D:/高粱/H3K27me3Metaplot")
p456 <- read.table("P456_gene.bed",sep="\t",row.names=4)
p241 <- read.table("P241_gene.bed",sep="\t",row.names=4)
head(p241)
k9ac <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLK9acvsCRK9ac_DEG_DMG_clusterd.csv",header=T,row.names=1)
k27ac <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLK27acvsCRK27ac_DEG_DMG_clusterd.csv",header=T,row.names=1)
k4me2 <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLK4me2vsCRK4me2_DEG_DMG_clusterd.csv",header=T,row.names=1)
k43 <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLK4me3vsCRK4me3_DEG_DMG_clusterd.csv",header=T,row.names=1)
k273 <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLK27me3vsCRK27me3_DEG_DMG_clusterd.csv",header=T,row.names=1)
k36 <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLK36me3vsCRK36me3_DEG_DMG_clusterd.csv",header=T,row.names=1)
h1az <- read.csv("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CLKH1AZvsCRKH1AZ_DEG_DMG_clusterd.csv",header=T,row.names=1)
###P456
name <- c("k9ac","k27ac","K4me2","K4me3","K36me3","KH1AZ")
a <- data.frame()
for (i in name){
  b <- read.csv(paste0("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CL",i,"vsCR",i,"_DEG_DMG_clusterd.csv"),header=T,row.names=1)
  a1 <- merge(p456,b,by="row.names",all=F)
  a1$group <- rep(i,nrow(a1))
  a <- rbind(a,a1)
}
head(a)
#统计数量
a_number <- aggregate(a$Row.names,list(a$group,a$type),length)
a_number <- a_number[order(a_number$Group.1,-a_number$x),]
p1 <- ggplot(a_number,mapping = aes(x=Group.1,y=x,fill=Group.2))+
  geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  #geom_bar(stat='identity',position = position_dodge(),show.legend = TRUE) +
  #labs(y = 'Number of peak region') +
  #scale_fill_manual(values=brewer.pal(12,"Paired")[c(5,8)]) +
  scale_fill_manual(values=c("#c2e59c","#be93c5"))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),legend.position = "top") + 
  guides(fill = guide_legend(title = 'Type')) +
  ylab("Number") +
  ggtitle("P456") +
  geom_text(aes(label = x), size = 3, hjust = "inward", vjust = "inward", position = "stack") +
  xlab(NULL)
###P456
name <- c("k9ac","k27ac","K4me2","K4me3","K36me3","KH1AZ")
a <- data.frame()
for (i in name){
  b <- read.csv(paste0("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2/CL",i,"vsCR",i,"_DEG_DMG_clusterd.csv"),header=T,row.names=1)
  a1 <- merge(p241,b,by="row.names",all=F)
  a1$group <- rep(i,nrow(a1))
  a <- rbind(a,a1)
}
head(a)
#统计数量
a_number <- aggregate(a$Row.names,list(a$group,a$type),length)
a_number <- a_number[order(a_number$Group.1,-a_number$x),]
p2 <- ggplot(a_number,mapping = aes(x=Group.1,y=x,fill=Group.2))+
  geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  #geom_bar(stat='identity',position = position_dodge(),show.legend = TRUE) +
  #labs(y = 'Number of peak region') +
  #scale_fill_manual(values=brewer.pal(12,"Paired")[c(5,8)]) +
  scale_fill_manual(values=c("#c2e59c","#be93c5"))+
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),legend.position = "top") + 
  guides(fill = guide_legend(title = 'Type')) +
  ylab("Number") +
  ggtitle("P241") +
  geom_text(aes(label = x), size = 3, hjust = "inward", vjust = "inward", position = "stack") +
  xlab(NULL)
####组图
grid.arrange(p1,p2,
             nrow=1,ncol=2)     %>%  ggsave("456_241_histone.pdf",.,width=210,height=100, units="mm")

#####在crispr-ceral上可视化高粱基因组的可视化数据
###在mysql数据库中导入sorghum基因组
#基因组数据
##在mysql中创造数据集，去掉染色体名称中的多于的字符
cd /disk2/users/che/database/sorghum
python rename_fasta.py Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa > Sorghum_bicolor.fa
/usr/local/bin/bp_seqfeature_load.pl -u che -p che_123 -c -f -d sorghum_che Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gff3 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa



######将bw格式转为wiggle格式
cd /public/home/chaohe/sorghum/chip/align/RPKM
#conda install -c bioconda ucsc-bigwigtowig
cut -f 1 /public/home/chaohe/sorghum/chip/meta/frip0.list | uniq | while read i;
do
bigWigToWig "$i"_rep0.rpkm.bw "$i".wig
done
cd /public/home/chaohe/sorghum/mnase/align
for i in CL CR PL PR
do 
bigWigToWig "$i"MN_rep0_rpkm.bw "$i"MN.wig; done
cd /public/home/chaohe/sorghum/RNA-seq/fq
for i in CL CR PL PR ; do
bigWigToWig "$i"_RNA.bw "$i"_RNA.wig; done

####转化格式
perl -p -i -e 's/chr//g' *wig
PWD="/disk2/users/che/database/sorghum"
mkdir gbrowse_track_of_ChIP
wiggle2gff3.pl --source=CLK27ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A1 CLK27ac.wig >CLK27ac.gff3
wiggle2gff3.pl --source=CLK27me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A2 CLK27me3.wig >CLK27me3.gff3
wiggle2gff3.pl --source=CLK36me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A3 CLK36me3.wig >CLK36me3.gff3
wiggle2gff3.pl --source=CLK43 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A4 CLK43.wig >CLK43.gff3
wiggle2gff3.pl --source=CLK4me2 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A5 CLK4me2.wig >CLK4me2.gff3
wiggle2gff3.pl --source=CLK9ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A6 CLK9ac.wig >CLK9ac.gff3
wiggle2gff3.pl --source=CLKH2AZ --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A7 CLKH2AZ.wig >CLKH2AZ.gff3
wiggle2gff3.pl --source=CL_RNA --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A8 CL_RNA.wig >CL_RNA.gff3
wiggle2gff3.pl --source=CRK27ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A9 CRK27ac.wig >CRK27ac.gff3
wiggle2gff3.pl --source=CRK27me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A10 CRK27me3.wig >CRK27me3.gff3
wiggle2gff3.pl --source=CRK36me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A11 CRK36me3.wig >CRK36me3.gff3
wiggle2gff3.pl --source=CRK43 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A12 CRK43.wig >CRK43.gff3
wiggle2gff3.pl --source=CRK4me2 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A13 CRK4me2.wig >CRK4me2.gff3
wiggle2gff3.pl --source=CRK9ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A14 CRK9ac.wig >CRK9ac.gff3
wiggle2gff3.pl --source=CRKH2AZ --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A15 CRKH2AZ.wig >CRKH2AZ.gff3
wiggle2gff3.pl --source=CR_RNA --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A16 CR_RNA.wig >CR_RNA.gff3
wiggle2gff3.pl --source=PLK27ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A17 PLK27ac.wig >PLK27ac.gff3
wiggle2gff3.pl --source=PLK27me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A18 PLK27me3.wig >PLK27me3.gff3
wiggle2gff3.pl --source=PLK36me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A19 PLK36me3.wig >PLK36me3.gff3
wiggle2gff3.pl --source=PLK43 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A20 PLK43.wig >PLK43.gff3
wiggle2gff3.pl --source=PLK4me2 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A21 PLK4me2.wig >PLK4me2.gff3
wiggle2gff3.pl --source=PLK9ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A22 PLK9ac.wig >PLK9ac.gff3
wiggle2gff3.pl --source=PLKH2AZ --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A23 PLKH2AZ.wig >PLKH2AZ.gff3
wiggle2gff3.pl --source=PL_RNA --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A24 PL_RNA.wig >PL_RNA.gff3
wiggle2gff3.pl --source=PRK27ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A25 PRK27ac.wig >PRK27ac.gff3
wiggle2gff3.pl --source=PRK27me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A26 PRK27me3.wig >PRK27me3.gff3
wiggle2gff3.pl --source=PRK36me3 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A27 PRK36me3.wig >PRK36me3.gff3
wiggle2gff3.pl --source=PRK43 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A28 PRK43.wig >PRK43.gff3
wiggle2gff3.pl --source=PRK4me2 --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A29 PRK4me2.wig >PRK4me2.gff3
wiggle2gff3.pl --source=PRK9ac --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A30 PRK9ac.wig >PRK9ac.gff3
wiggle2gff3.pl --source=PRKH2AZ --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A31 PRKH2AZ.wig >PRKH2AZ.gff3
wiggle2gff3.pl --source=PR_RNA --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A32 PR_RNA.wig >PR_RNA.gff3
wiggle2gff3.pl --source=PRMN --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A33 PRMN.wig >PRMN.gff3
wiggle2gff3.pl --source=CRMN --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A34 CRMN.wig >CRMN.gff3
wiggle2gff3.pl --source=PLMN --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A35 PLMN.wig >PLMN.gff3
wiggle2gff3.pl --source=CLMN --method=ChIP --path=$PWD/gbrowse_track_of_ChIP --trackname=track_A36 CLMN.wig >CLMN.gff3
##导入数据
/usr/local/bin/bp_seqfeature_load.pl -u che -p che_123 -c -f -d test Sorghum.fa Sorghum_bicolor_NCBIv3.52.gff3
/usr/local/bin/bp_seqfeature_load.pl -u che -p che_123 -f -d test C*.gff3 P*.gff3 


SORBI_3002G146501
2:29886079-29888405


-rw-rw-r-- 1 che che 584257335 Apr  9 12:21 CLK27ac.wig
-rw-rw-r-- 1 che che 383094128 Apr  9 12:21 CLK27me3.wig
-rw-rw-r-- 1 che che 279708670 Apr  9 12:21 CLK36me3.wig
-rw-rw-r-- 1 che che 264896666 Apr  9 12:22 CLK43.wig
-rw-rw-r-- 1 che che 398418009 Apr  9 12:22 CLK4me2.wig
-rw-rw-r-- 1 che che 322475081 Apr  9 12:22 CLK9ac.wig
-rw-rw-r-- 1 che che 449225809 Apr  9 12:22 CLKH2AZ.wig
-rw-rw-r-- 1 che che  34670635 Apr  9 12:22 CL_RNA.wig
-rw-rw-r-- 1 che che 612357743 Apr  9 12:22 CRK27ac.wig
-rw-rw-r-- 1 che che 541538170 Apr  9 12:22 CRK27me3.wig
-rw-rw-r-- 1 che che 523397704 Apr  9 12:22 CRK36me3.wig
-rw-rw-r-- 1 che che 337923469 Apr  9 12:22 CRK43.wig
-rw-rw-r-- 1 che che 658028957 Apr  9 12:22 CRK4me2.wig
-rw-rw-r-- 1 che che 467919671 Apr  9 12:22 CRK9ac.wig
-rw-rw-r-- 1 che che 734886416 Apr  9 12:22 CRKH2AZ.wig
-rw-rw-r-- 1 che che  36850291 Apr  9 12:22 CR_RNA.wig
-rw-rw-r-- 1 che che 412689609 Apr  9 12:22 PLK27ac.wig
-rw-rw-r-- 1 che che 380883082 Apr  9 12:22 PLK27me3.wig
-rw-rw-r-- 1 che che 271381418 Apr  9 12:22 PLK36me3.wig
-rw-rw-r-- 1 che che 325651781 Apr  9 12:22 PLK43.wig
-rw-rw-r-- 1 che che 296541235 Apr  9 12:22 PLK4me2.wig
-rw-rw-r-- 1 che che 262896606 Apr  9 12:23 PLK9ac.wig
-rw-rw-r-- 1 che che 487558922 Apr  9 12:23 PLKH2AZ.wig
-rw-rw-r-- 1 che che  35019000 Apr  9 12:23 PL_RNA.wig
-rw-rw-r-- 1 che che 700912748 Apr  9 12:23 PRK27ac.wig
-rw-rw-r-- 1 che che 395841000 Apr  9 12:23 PRK27me3.wig
-rw-rw-r-- 1 che che 573272154 Apr  9 12:23 PRK36me3.wig
-rw-rw-r-- 1 che che 359786840 Apr  9 12:23 PRK43.wig
-rw-rw-r-- 1 che che 648500365 Apr  9 12:23 PRK4me2.wig
-rw-rw-r-- 1 che che 659582254 Apr  9 12:23 PRK9ac.wig
-rw-rw-r-- 1 che che 729294792 Apr  9 12:23 PRKH2AZ.wig
-rw-rw-r-- 1 che che  37595946 Apr  9 12:23 PR_RNA.wig


#####热图展示KEGG C4基因term内基因的组蛋白修饰变化情况
library(openxlsx)
library(pheatmap)
library(viridis)
library(RColorBrewer)
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("C4基因.xlsx")
head(c4)
cg4 <- c4[,-11]
row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
cg4 <- cg4[,c(4:10)]
head(cg4)
annotation_row = data.frame(
  CellType = rep(c("C4", "other"), c(9,18)))
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
ann_colors = list(CellType = c("C4" = brewer.pal(12,'Paired')[[7]], "other" = brewer.pal(12,'Set3')[[9]]),
  Modification = c(H3K9ac=brewer.pal(7,'Paired')[[5]],
                   H3K27ac=brewer.pal(7,'Accent')[[2]],
                   H3K4me3=brewer.pal(7,'Paired')[[1]],
                   H3K36me3=brewer.pal(8,'Paired')[[3]],
                   H2A.Z=brewer.pal(8,'Accent')[[8]],
                   H3K4me2=brewer.pal(7,'Set2')[[1]],
                   H3K27me3=brewer.pal(7,'Accent')[[7]]))
bk <- c(seq(-5,-0.1,by=0.01),seq(0.1,2,by=0.01))
p1 <- pheatmap(cg4, 
         cluster_cols = F, 
         cluster_rows = F, 
         color = colorRampPalette(c("#1F77B4", "#5F9DC9", "#7FB1D4", "#9FC4DE", "#BED8E9", "#DFEBF4",
                                    'white', "#F985A0", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
         #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
         #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         #color = viridis(8),
         #color = viridis(8, option = "G"),
         border_color = "black",
         annotation_col = annotation_col, 
         annotation_row = annotation_row, 
         #legend_breaks = c(-5, 0, 2),
         legend_labels = c("Down","","Up"), 
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         display_numbers =T)
rg4 <- c4[,c(1,2,3,11)]
rg4$nam <- paste(rg4$X3,rg4$X1,sep="_")
name <- rg4[,5]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "mRNA"
head(rg4)
annotation_col = data.frame(colnames(rg4))
rownames(annotation_col) = colnames(rg4)
colnames(annotation_col) = "Modification"
ann_colors = list(CellType = c("C4" = brewer.pal(12,'Paired')[[7]], "other" = brewer.pal(12,'Set3')[[9]]),
                  Modification = c(mRNA=brewer.pal(7,'Paired')[[5]]))
p2 <- pheatmap(rg4, 
         cluster_cols = F, 
         cluster_rows = F, 
         color = colorRampPalette(c("#1F77B4", "#5F9DC9", "#7FB1D4", "#9FC4DE", "#BED8E9", "#DFEBF4",
                                    'white', "#F985A0", "#F64971"))(80), #热图色块颜色是从蓝到红分为100个等级
         #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
         #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         #color = viridis(4),
         #color = viridis(8, option = "G"),
         #legend_breaks = c(-10, 0, 10),
         border_color = "black",
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         show_rownames = F,
         display_numbers =T)
library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("C4_heatmap.pdf",.,width=500,height=210, units="mm")

#####热图展示KEGG lignin biosynthesis gene root-leaf内基因的组蛋白修饰变化情况
library(openxlsx)
library(pheatmap)
library(viridis)
library(RColorBrewer)
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("lignin biosynthesis gene root-leaf.xlsx")
head(c4)
cg4 <- c4[,-10]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- cg4$X2
cg4 <- cg4[,c(3:9)]
head(cg4)
annotation_row = data.frame(c4$X1)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
library("scales")
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Summer")
max_n <- attr(pal, "max_n")
ann_colors = list(Modification = c(K9ac=brewer.pal(7,'Paired')[[5]],
                   K27ac=brewer.pal(7,'Accent')[[2]],
                   K4me3=brewer.pal(7,'Paired')[[1]],
                   K36me3=brewer.pal(8,'Paired')[[3]],
                   H2AZ=brewer.pal(8,'Accent')[[8]],
                   K4me2=brewer.pal(7,'Set2')[[1]],
                   K27me3=brewer.pal(7,'Accent')[[7]]))
#bk <- c(seq(-5,-0.1,by=0.01),seq(0.1,2,by=0.01))
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#3786BC", "#5095C4", "#69A4CD", "#82B3D5","#E6EFF6",'white',"#FEEAEF", "#FDD6DF", "#FCC2CF", "#FBAEBF", "#FA99B0", "#F985A0", "#F87190", "#F75D80", "#F64971"))(30), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,2,10)]
rg4$name <- rg4$X2
name <- rg4[,4]
rg4 <- data.frame(rg4[,c(-1:-2,-4)])
row.names(rg4) <- name
colnames(rg4) <- "mRNA"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#478FC1", "#5C9CC8","#3786BC", "#5095C4", "#69A4CD", "#82B3D5", "#9BC2DD", "#B4D1E6", "#CDE0EE", "#E6EFF6",
                                          "#FFFFFF", "#FEF5F7", "#FEEBF0", "#FDE2E8", "#FCC5D2", "#FBBBCA", "#FBB2C3",
                                          "#FAA8BB", "#FA9FB4", "#F995AC", "#F98CA5", "#F8829D",
                                          "#F87896", "#F76F8E", "#F76587", "#F65C7F", "#F65278", "#F64971"))(14), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               border_color = "black",
               show_rownames = F,
               display_numbers =T)

library(dplyr)
library(ggplot2)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("lignin biosynthesis gene root-leaf_heatmapr.pdf",.,width=500,height=210, units="mm")


#####热图展示KEGG lignin biosynthesis gene stress内基因的组蛋白修饰变化情况
library(openxlsx)
library(pheatmap)
library(viridis)
library(RColorBrewer)
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("lignin biosynthesis gene stress.xlsx")
head(c4)
cg4 <- c4[,-8]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- cg4$symbol
cg4 <- cg4[,c(3:7)]
head(cg4)
annotation_row = data.frame(c4$type)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
library("scales")
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Classic 20")
max_n <- attr(pal, "max_n")
ann_colors = list(
  c4.type = c("PAL" = pal(8)[1], "4CL" = pal(8)[2],
              "HCT" = pal(8)[3], "C3H" = pal(8)[4],
              "C4H" = pal(8)[5], "COMT" = pal(8)[6],
              "CCOAOMT" = pal(8)[7],"F5H" = pal(8)[8],
              "CCR" = pal(10)[9],"CAD" = pal(10)[10]),
  Modification = c(K9ac=brewer.pal(7,'Paired')[[5]],
                   K27ac=brewer.pal(7,'Accent')[[2]],
                   K4me3=brewer.pal(7,'Paired')[[1]],
                   H2A.Z=brewer.pal(8,'Accent')[[8]],
                   K4me2=brewer.pal(7,'Set2')[[1]]))
#bk <- c(seq(-5,-0.1,by=0.01),seq(0.1,2,by=0.01))
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#5095C4",  "#82B3D5","#E6EFF6",'white',"#FEEAEF", "#FDD6DF", "#FCC2CF", "#FBAEBF", "#FA99B0", "#F985A0", "#F87190", "#F75D80", "#F64971"))(30), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,2,8)]
rg4$name <- rg4$symbol
name <- rg4[,4]
rg4 <- data.frame(rg4[,c(-1:-2,-4)])
row.names(rg4) <- name
colnames(rg4) <- "mRNA"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4","#3383BA", "#69A4CD", "#82B3D5",  "#B4D1E6", "#CDE0EE", 
                                          "#FFFFFF", "#FEF5F7", "#FEEBF0", "#FDE2E8","#FDD8E1","#FCCFD9", "#FCC5D2", "#FBBBCA", "#FBB2C3",
                                          "#FAA8BB", "#FA9FB4", "#F995AC", "#F98CA5", "#F8829D",
                                          "#F87896", "#F76F8E", "#F76587", "#F65C7F", "#F65278", "#F64971"))(28), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               border_color = "black",
               show_rownames = F,
               display_numbers =T)

library(dplyr)
library(ggplot2)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("lignin_biosynthesis_gene_stress_heatmapr.pdf",.,width=500,height=210, units="mm")

####光合
setwd("D:/高粱/GO和KEGG/筛选后")
c4 <- read.csv("photosynthesis2.csv")
head(c4)
cg4 <- c4[,-11]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- paste(cg4$symbol,cg4$locus,sep="_")
cg4 <- cg4[,c(4:10)]
head(cg4)
annotation_row = data.frame(c4$type)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
unique(c4$type)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Summer")
max_n <- attr(pal, "max_n")
ann_colors = list(c4.type = c("Photosystem_I" = pal(8)[1], "Cytochrome_b6/f_complex" = pal(8)[2],
                              "Photosystem_II" = pal(8)[3], "ATP_synthase" = pal(8)[4],
                              "plastocyanin" = pal(8)[5], "cytc6" = pal(8)[6],
                              "chlorophyll_biosynthesis" = pal(8)[7]),
                  Modification = c(K9ac=brewer.pal(7,'Paired')[[5]],
                                   K27ac=brewer.pal(7,'Accent')[[2]],
                                   K4me3=brewer.pal(7,'Paired')[[1]],
                                   K36me3=brewer.pal(8,'Paired')[[3]],
                                   H2AZ=brewer.pal(8,'Accent')[[8]],
                                   K4me2=brewer.pal(7,'Set2')[[1]],
                                   K27me3=brewer.pal(7,'Accent')[[7]]))
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#2A7EB7", "#3685BB", "#428CBF", "#4E93C3", "#599AC7", "#65A1CB", "#71A9CF", 
                                          "#7DB0D3", "#89B7D7", "#94BEDB", "#A0C5DF", "#ACCCE3","#B8D4E7", "#C4DBEB", "#CFE2EF", 
                                          "#DBE9F3", "#E7F0F7", "#F3F7FB",
                                          "#FFFFFF", "#FEF5F7",  "#FDD8E1",
                                          "#FBBBCA", "#FBB2C3", "#FAA8BB", "#FA9FB4", "#F65C7F", "#F64971"))(20), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,2,3,11)]
rg4$name <- rg4$locus
name <- rg4[,5]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "mRNA"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#2A7EB7", "#3685BB", "#428CBF", "#4E93C3", "#599AC7", "#65A1CB", "#71A9CF", 
                                          "#7DB0D3", "#89B7D7", "#94BEDB", "#A0C5DF", "#ACCCE3","#B8D4E7", "#C4DBEB", "#CFE2EF", 
                                          "#DBE9F3", "#E7F0F7", "#F3F7FB",
                                          "#FFFFFF"))(50), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               border_color = "black",
               show_rownames = F,
               display_numbers =T)

library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("photosynthesis2_heatmapr.pdf",.,width=500,height=210, units="mm")

####POD_STRESS
setwd("D:/高粱/GO和KEGG/筛选后")
c4 <- read.xlsx("POD stress.xlsx")
head(c4)
cg4 <- c4[,-8]
row.names(cg4) <- cg4$lous
head(cg4)
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
#row.names(cg4) <- paste(cg4$symbol,cg4$locus,sep="_")
cg4 <- cg4[,c(3:7)]
head(cg4)
annotation_row = data.frame(c4$type)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
unique(c4$type)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Hue Circle")
max_n <- attr(pal, "max_n")
ann_colors = list(c4.type = c("POD" = pal(8)[1]),
                  Modification = c(K9ac=brewer.pal(7,'Paired')[[5]],
                                   K27ac=brewer.pal(7,'Accent')[[2]],
                                   K4me3=brewer.pal(7,'Paired')[[1]],
                                   H2A.Z=brewer.pal(8,'Accent')[[8]],
                                   K4me2=brewer.pal(7,'Set2')[[1]]))
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#2A7EB7", "#3685BB", "#428CBF", "#4E93C3", "#599AC7", "#65A1CB", "#71A9CF", 
                                          "#7DB0D3", "#89B7D7", "#94BEDB", "#A0C5DF", "#ACCCE3","#B8D4E7", "#C4DBEB", "#CFE2EF", 
                                          "#DBE9F3", "#E7F0F7", "#F3F7FB",
                                          "#FFFFFF", "#FDDDE5", "#FCCDD8", "#FBBCCB", "#FAACBE", "#FA9BB1", "#F98BA4", "#F87A97", "#F76A8A", "#F6597D", "#F64971"))(20), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               #annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,2,8)]
rg4$name <- rg4$lous
name <- rg4[,2]
rg4 <- data.frame(rg4[,c(-1:-2,-4)])
row.names(rg4) <- name
colnames(rg4) <- "mRNA"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#3786BC", "#448DC0", "#5095C4", "#5D9CC8", "#69A4CD", "#76ABD1", "#82B3D5",
                                          "#8FBBD9", "#9BC2DD", "#A7CAE1", "#B4D1E6", "#C0D9EA",
                                          "#CDE0EE", "#D9E8F2", "#E6EFF6", "#F2F7FA", "#FFFFFF",
                                          "#FEEEF2", "#FDDDE5", "#FCCDD8", "#FBBCCB", "#FAACBE", "#FA9BB1", 
                                          "#F98BA4", "#F87A97", "#F76A8A", "#F6597D", "#F64971"))(60), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               border_color = "black",
               show_rownames = F,
               display_numbers =T)

library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("POD_stress_heatmapr.pdf",.,width=500,height=210, units="mm")

#POD leaf root.xlsx
setwd("D:/高粱/GO和KEGG/筛选后")
c4 <- read.xlsx("POD leaf root.xlsx")
head(c4)
cg4 <- c4[,-9]
row.names(cg4) <- cg4$X1
head(cg4)
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
#row.names(cg4) <- paste(cg4$symbol,cg4$locus,sep="_")
cg4 <- cg4[,c(2:8)]
head(cg4)
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#3685BB", "#4E93C3", "#65A1CB","#7DB0D3", "#94BEDB", "#ACCCE3","#B8D4E7", "#C4DBEB", "#CFE2EF", 
                                          "#DBE9F3", "#E7F0F7", "#F3F7FB",
                                          "#FFFFFF", "#FEF3F6", "#FDE8ED", "#FDDCE4", "#FCD1DB", "#FCC6D2", "#FBBAC9", "#FBAFC0", "#FAA3B8", "#F998AF", "#F98DA6", "#F8819D", "#F87694",
                                          "#F76B8B", "#F75F82", "#F65479", "#F64971"))(16), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               #annotation_col = annotation_col, 
               #annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,9)]
rg4$name <- rg4$X1
name <- rg4[,1]
rg4 <- data.frame(rg4[,c(-1,-3)])
row.names(rg4) <- name
colnames(rg4) <- "mRNA"
head(rg4)
p2 <- pheatmap(rg4, 
         cluster_cols = F, 
         cluster_rows = F, 
         color = colorRampPalette(c("#1F77B4", "#448DC0", "#69A4CD", "#8FBBD9", "#B4D1E6", "#D9E8F2", "#FFFFFF",
                                    "#FEF3F6", "#FDE8ED", "#FDDCE4", "#FCD1DB", "#FCC6D2", "#FBBAC9", "#FBAFC0",
                                    "#FAA3B8", "#F998AF", "#F98DA6", "#F8819D", "#F87694",
                                    "#F76B8B", "#F75F82", "#F65479", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
         #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
         #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         #color = viridis(4),
         #color = viridis(8, option = "G"),
         #legend_breaks = c(-10, 0, 10),
         border_color = "black",
         show_rownames = F,
         display_numbers =T)

library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("POD_leaf_root_heatmapr.pdf",.,width=500,height=210, units="mm")


#photosynthesis stress
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("photosynthesis stress.xlsx")
head(c4)
cg4 <- c4[,-9]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- cg4$symbol
cg4 <- cg4[,c(4:8)]
head(cg4)
annotation_row = data.frame(c4$type)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
unique(c4$type)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Summer")
max_n <- attr(pal, "max_n")
ann_colors = list(c4.type = c("Photosystem_I" = pal(8)[1], "Cytochrome_b6/f_complex" = pal(8)[2],
                              "Photosystem_II" = pal(8)[3], "ATP_synthase" = pal(8)[4],
                              "plastocyanin" = pal(8)[5], "cytc6" = pal(8)[6],
                              "chlorophyll_biosynthesis" = pal(8)[7]),
                  Modification = c(K9ac=brewer.pal(7,'Paired')[[5]],
                                   K27ac=brewer.pal(7,'Accent')[[2]],
                                   K4me3=brewer.pal(7,'Paired')[[1]],
                                   K36me3=brewer.pal(8,'Paired')[[3]],
                                   H2AZ=brewer.pal(8,'Accent')[[8]],
                                   K4me2=brewer.pal(7,'Set2')[[1]],
                                   K27me3=brewer.pal(7,'Accent')[[7]]))
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4",   "#4E93C3",  "#7DB0D3",
                                          "#A0C5DF", "#B8D4E7", "#CFE2EF", 
                                           "#E7F0F7",
                                          "#FFFFFF", "#FEF4F6", "#FDE9EE", "#FDDEE5", "#FCD4DD", "#FCC9D5",
                                          "#FBBECC", "#FBB4C4", "#FAA9BC", "#FA9EB3", "#F993AB", "#F989A3",
                                          "#F87E9A", "#F87392", "#F7698A", "#F75E81", "#F65379","#F65378", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               #annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,2,3,9)]
rg4$name <- rg4$symbol
name <- rg4[,5]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "Expression"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#2A7EB7", "#3685BB", "#428CBF", "#4E93C3", "#599AC7", "#65A1CB", "#71A9CF", 
                                          "#7DB0D3", "#89B7D7", "#94BEDB", "#A0C5DF", "#ACCCE3","#B8D4E7", "#C4DBEB", "#CFE2EF", 
                                          "#DBE9F3", "#E7F0F7", "#F3F7FB",
                                          "#FFFFFF","#FDE3E9","#FCD4DE","#FCC7D3","#FBB9C8","#FAABBD","#FA9CB2","#F98FA7","#F8819C","#F87291","#F76486",
                                          "#F6567B", "#F64971"))(60), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               border_color = "black",
               show_rownames = F,
               display_numbers =T)

library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("photosynthesis_stress_heatmapr.pdf",.,width=500,height=210, units="mm")

#photosynthesis stress
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("C4_stress.xlsx")
head(c4)
cg4 <- c4[,-9]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- cg4$symbol
cg4 <- cg4[,c(4:8)]
head(cg4)
annotation_row = data.frame(c4$type)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
unique(c4$type)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Summer")
max_n <- attr(pal, "max_n")
ann_colors = list(CellType = c("C4" = brewer.pal(12,'Paired')[[7]], "other" = brewer.pal(12,'Set3')[[9]]),
                  Modification = c(H3K9ac=brewer.pal(7,'Paired')[[5]],
                                   H3K27ac=brewer.pal(7,'Accent')[[2]],
                                   H3K4me3=brewer.pal(7,'Paired')[[1]],
                                   H3K36me3=brewer.pal(8,'Paired')[[3]],
                                   H2A.Z=brewer.pal(8,'Accent')[[8]],
                                   H3K4me2=brewer.pal(7,'Set2')[[1]],
                                   H3K27me3=brewer.pal(7,'Accent')[[7]]))

p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4",   "#4E93C3",  "#7DB0D3",
                                          "#A0C5DF", "#B8D4E7", "#CFE2EF", 
                                          "#E7F0F7",
                                          "#FFFFFF", "#FEF4F6", "#FDE9EE", "#FDDEE5", "#FCD4DD", "#FCC9D5",
                                          "#FBBECC", "#FBB4C4", "#FAA9BC", "#FA9EB3", "#F993AB", "#F989A3",
                                          "#F87E9A", "#F87392", "#F7698A", "#F75E81", "#F65379","#F65378", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(8),
               #color = viridis(8, option = "G"),
               border_color = "black",
               #annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #legend_breaks = c(-5, 0, 2),
               legend_labels = c("Down","","Up"), 
               annotation_colors = ann_colors,
               show_rownames = TRUE,
               display_numbers =T)
rg4 <- c4[,c(1,2,3,9)]
rg4$name <- rg4$symbol
name <- rg4[,5]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "Expression"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#69A4CD", "#B4D1E6", "#FFFFFF","#FEF7F9", "#FEEFF3", "#FDE8ED", "#FDE0E7", "#FDD9E1",
                                          "#FCD1DB", "#FCC9D5", "#FCC2CF", "#FBBAC9", "#FBB3C3", "#FAABBD", "#FAA3B8", "#FA9CB2",
                                          "#F994AC", "#F98DA6", "#F985A0", "#F87E9A", "#F87694", "#F76E8E", "#F76788", "#F75F82",
                                          "#F6587C", "#F65076", "#F64971"))(60), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               #annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               border_color = "black",
               show_rownames = F,
               display_numbers =T)

library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("C4_stress_heatmapr.pdf",.,width=500,height=210, units="mm")


library(gridExtra)
library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ggthemes)
#ABA stress
#leaf
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("胁迫基因热图.xlsx",sheet=1)
head(c4)
cg4 <- c4[,-1:-2]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- cg4$Locus
cg4 <- cg4[,c(3:7)]
head(cg4)
annotation_row = data.frame(c4$symbol)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
unique(c4$type)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Summer")
max_n <- attr(pal, "max_n")
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#3786BC", "#5095C4", "#69A4CD", "#82B3D5", "#9BC2DD",
                                          "#B4D1E6", "#CDE0EE", "#E6EFF6", "#EAF2F8", "#FFFFFF", "#FEF4F6", 
                                          "#FDE9EE", "#FDDEE5", "#FCD4DD", "#FCC9D5",
                                          "#FBBECC", "#FBB4C4", "#FAA9BC", "#FA9EB3", "#F993AB", "#F989A3",
                                          "#F87E9A", "#F87392", "#F7698A", "#F75E81", "#F65379","#F65378", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
               border_color = "black",
               display_numbers = T,
               annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               show_rownames = TRUE)
rg4 <- c4[,c(1,2,3,4)]
rg4$name <- rg4$symbol
name <- rg4[,3]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "Expression"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#CDE0EE",  "#FFFFFF",
                                          "#FEF7F9", "#FEEFF3", "#FDE8ED", "#FDE0E7", "#FDD9E1",
                                          "#FCD1DB", "#FCC9D5", "#FCC2CF", "#FBBAC9", "#FBB3C3", "#FAABBD", "#FAA3B8", "#FA9CB2",
                                          "#F994AC", "#F98DA6", "#F985A0", "#F87E9A", "#F87694", "#F76E8E", "#F76788", "#F75F82",
                                          "#F6587C", "#F65076", "#F64971"))(60), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               #annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               border_color = "black",
               show_rownames = F,
               display_numbers =T)
library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("leaf_ABA_stress_related_heatmapr.pdf",.,width=500,height=210, units="mm")


#root
setwd("D:/高粱/GO和KEGG筛选后")
c4 <- read.xlsx("胁迫基因热图.xlsx",sheet=2)
head(c4)
cg4 <- c4[,-1:-2]
#row.names(cg4) <- paste(cg4$X3,cg4$X1,sep="_")
row.names(cg4) <- cg4$Locus
cg4 <- cg4[,c(3:7)]
head(cg4)
annotation_row = data.frame(c4$symbol)
rownames(annotation_row) = rownames(cg4)
annotation_col = data.frame(colnames(cg4))
rownames(annotation_col) = colnames(cg4)
colnames(annotation_col) = "Modification"
unique(c4$type)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal <- tableau_color_pal("Summer")
max_n <- attr(pal, "max_n")
p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#4B92C3", "#629FCA", "#78ADD2", "#8FBBD9",
                                          "#A5C8E1", "#BBD6E8", "#D2E3F0", "#E8F1F7", "#FFFFFF",
                                          "#FEF4F6", "#FDE2E8", "#FCD1DB", "#FBC0CE", "#FAAFC0",
                                          "#FA9EB3", "#F98DA6", "#F87C98", "#F76B8B", "#F65A7E", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
               border_color = "black",
               display_numbers = T,
               annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               show_rownames = TRUE)
rg4 <- c4[,c(1,2,3,4)]
rg4$name <- rg4$symbol
name <- rg4[,3]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "Expression"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#3786BC", "#5095C4", "#69A4CD", "#82B3D5", 
                                          "#9BC2DD", "#B4D1E6", "#CDE0EE", "#E6EFF6", "#FFFFFF",
                                          "#FEF7F9", "#FEEFF3", "#FDE8ED", "#FDE0E7", "#FDD9E1",
                                          "#FCD1DB", "#FCC9D5", "#FCC2CF", "#FBBAC9", "#FBB3C3", "#FAABBD", "#FAA3B8", "#FA9CB2",
                                          "#F994AC", "#F98DA6", "#F985A0", "#F87E9A", "#F87694", "#F76E8E", "#F76788", "#F75F82",
                                          "#F6587C", "#F65076", "#F64971"))(30), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               #annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               border_color = "black",
               show_rownames = F,
               display_numbers =T)
library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("root_ABA_stress_related_heatmapr.pdf",.,width=500,height=210, units="mm")


colorRampPalette(c("#1F77B4","#FFFFFF"))(10)





######差异修饰之间两两overlap
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2")
###up
#H3K9ac
k9ac <- read.csv("CLK9acvsCRK9ac_DEG_DMG_clusterd.csv",row.names=1)
head(k9ac)
k9ac <- k9ac[which(k9ac$CL != "Unexpress" | k9ac$CR != "Unexpress"),]
k9ac_up <- k9ac[which(k9ac$histone_log2FC >= 0.75 & k9ac$histone_padj < 0.05),]
#H3K4me3
K4me3 <- read.csv("CLK4me3vsCRK4me3_DEG_DMG_clusterd.csv",row.names=1)
head(K4me3)
K4me3 <- K4me3[which(K4me3$CL != "Unexpress" | K4me3$CR != "Unexpress"),]
K4me3_up <- K4me3[which(K4me3$histone_log2FC >= 0.75 & K4me3$histone_padj < 0.05),]
#K4me2
K4me2 <- read.csv("CLK4me2vsCRK4me2_DEG_DMG_clusterd.csv",row.names=1)
head(K4me2)
K4me2 <- K4me2[which(K4me2$CL != "Unexpress" | K4me2$CR != "Unexpress"),]
K4me2_up <- K4me2[which(K4me2$histone_log2FC >= 0.75 & K4me2$histone_padj < 0.05),]
#K27ac
K27ac <- read.csv("CLK27acvsCRK27ac_DEG_DMG_clusterd.csv",row.names=1)
head(K27ac)
K27ac <- K27ac[which(K27ac$CL != "Unexpress" | K27ac$CR != "Unexpress"),]
K27ac_up <- K27ac[which(K27ac$histone_log2FC >= 0.75 & K27ac$histone_padj < 0.05),]
#K27me3
K27me3 <- read.csv("CLK27me3vsCRK27me3_DEG_DMG_clusterd.csv",row.names=1)
head(K27me3)
K27me3 <- K27me3[which(K27me3$CL != "Unexpress" | K27me3$CR != "Unexpress"),]
K27me3_up <- K27me3[which(K27me3$histone_log2FC >= 0.75 & K27me3$histone_padj < 0.05),]
#K36me3
K36me3 <- read.csv("CLK36me3vsCRK36me3_DEG_DMG_clusterd.csv",row.names=1)
head(K36me3)
K36me3 <- K36me3[which(K36me3$CL != "Unexpress" | K36me3$CR != "Unexpress"),]
K36me3_up <- K36me3[which(K36me3$histone_log2FC >= 0.75 & K36me3$histone_padj < 0.05),]
#KH1AZ
KH1AZ <- read.csv("CLKH1AZvsCRKH1AZ_DEG_DMG_clusterd.csv",row.names=1)
head(KH1AZ)
KH1AZ <- KH1AZ[which(KH1AZ$CL != "Unexpress" | KH1AZ$CR != "Unexpress"),]
KH1AZ_up <- KH1AZ[which(KH1AZ$histone_log2FC >= 0.75 & KH1AZ$histone_padj < 0.05),]
###down
#H3K9ac
k9ac <- read.csv("CLK9acvsCRK9ac_DEG_DMG_clusterd.csv",row.names=1)
head(k9ac)
k9ac <- k9ac[which(k9ac$CL != "Unexpress" | k9ac$CR != "Unexpress"),]
k9ac_down <- k9ac[which(k9ac$histone_log2FC <= -0.75 & k9ac$histone_padj < 0.05),]
#H3K4me3
K4me3 <- read.csv("CLK4me3vsCRK4me3_DEG_DMG_clusterd.csv",row.names=1)
head(K4me3)
K4me3 <- K4me3[which(K4me3$CL != "Unexpress" | K4me3$CR != "Unexpress"),]
K4me3_down <- K4me3[which(K4me3$histone_log2FC <= -0.75 & K4me3$histone_padj < 0.05),]
#K4me2
K4me2 <- read.csv("CLK4me2vsCRK4me2_DEG_DMG_clusterd.csv",row.names=1)
head(K4me2)
K4me2 <- K4me2[which(K4me2$CL != "Unexpress" | K4me2$CR != "Unexpress"),]
K4me2_down <- K4me2[which(K4me2$histone_log2FC <= -0.75 & K4me2$histone_padj < 0.05),]
#K27ac
K27ac <- read.csv("CLK27acvsCRK27ac_DEG_DMG_clusterd.csv",row.names=1)
head(K27ac)
K27ac <- K27ac[which(K27ac$CL != "Unexpress" | K27ac$CR != "Unexpress"),]
K27ac_down <- K27ac[which(K27ac$histone_log2FC <= -0.75 & K27ac$histone_padj < 0.05),]
#K27me3
K27me3 <- read.csv("CLK27me3vsCRK27me3_DEG_DMG_clusterd.csv",row.names=1)
head(K27me3)
K27me3 <- K27me3[which(K27me3$CL != "Unexpress" | K27me3$CR != "Unexpress"),]
K27me3_down <- K27me3[which(K27me3$histone_log2FC <= -0.75 & K27me3$histone_padj < 0.05),]
#K36me3
K36me3 <- read.csv("CLK36me3vsCRK36me3_DEG_DMG_clusterd.csv",row.names=1)
head(K36me3)
K36me3 <- K36me3[which(K36me3$CL != "Unexpress" | K36me3$CR != "Unexpress"),]
K36me3_down <- K36me3[which(K36me3$histone_log2FC <= -0.75 & K36me3$histone_padj < 0.05),]
#KH1AZ
KH1AZ <- read.csv("CLKH1AZvsCRKH1AZ_DEG_DMG_clusterd.csv",row.names=1)
head(KH1AZ)
KH1AZ <- KH1AZ[which(KH1AZ$CL != "Unexpress" | KH1AZ$CR != "Unexpress"),]
KH1AZ_down <- KH1AZ[which(KH1AZ$histone_log2FC <= -0.75 & KH1AZ$histone_padj < 0.05),]

#####以k27ac_up为参考与其它修饰进行overlap
same <- length(intersect(x=k9ac_up$Row.names, y = K4me3_up$Row.names))
former <- length(setdiff(k9ac_up$Row.names, K4me3_up$Row.names))
latter <- length(setdiff(K4me3_up$Row.names, k9ac_up$Row.names))
k9ac_up_vs_K4me3_up <- data.frame(same,former,latter)
k9ac_up_vs_K4me3_up$type <- "k9ac_up_vs_K4me3_up"
same <- length(intersect(x=k9ac_up$Row.names, y = K4me2_down$Row.names))
former <- length(setdiff(k9ac_up$Row.names, K4me2_down$Row.names))
latter <- length(setdiff(K4me2_down$Row.names, k9ac_up$Row.names))
k9ac_up_vs_K4me2_down <- data.frame(same,former,latter)
k9ac_up_vs_K4me2_down$type <- "k9ac_up_vs_K4me2_down"
same <- length(intersect(x=k9ac_up$Row.names, y = K27ac_up$Row.names))
former <- length(setdiff(k9ac_up$Row.names, K27ac_up$Row.names))
latter <- length(setdiff(K27ac_up$Row.names, k9ac_up$Row.names))
k9ac_up_vs_K27ac_up <- data.frame(same,former,latter)
k9ac_up_vs_K27ac_up$type <- "k9ac_up_vs_K27ac_up"
same <- length(intersect(x=k9ac_up$Row.names, y = K27me3_down$Row.names))
former <- length(setdiff(k9ac_up$Row.names, K27me3_down$Row.names))
latter <- length(setdiff(K27me3_down$Row.names, k9ac_up$Row.names))
k9ac_up_vs_K27me3_down <- data.frame(same,former,latter)
k9ac_up_vs_K27me3_down$type <- "k9ac_up_vs_K27me3_down"
same <- length(intersect(x=k9ac_up$Row.names, y = K27me3_down$Row.names))
former <- length(setdiff(k9ac_up$Row.names, K27me3_down$Row.names))
latter <- length(setdiff(K27me3_down$Row.names, k9ac_up$Row.names))
k9ac_up_vs_K27me3_down <- data.frame(same,former,latter)
k9ac_up_vs_K27me3_down$type <- "k9ac_up_vs_K27me3_down"
same <- length(intersect(x=k9ac_up$Row.names, y = K36me3_up$Row.names))
former <- length(setdiff(k9ac_up$Row.names, K36me3_up$Row.names))
latter <- length(setdiff(K36me3_up$Row.names, k9ac_up$Row.names))
k9ac_up_vs_K36me3_up <- data.frame(same,former,latter)
k9ac_up_vs_K36me3_up$type <- "k9ac_up_vs_K36me3_up"
same <- length(intersect(x=k9ac_up$Row.names, y = KH1AZ_down$Row.names))
former <- length(setdiff(k9ac_up$Row.names, KH1AZ_down$Row.names))
latter <- length(setdiff(KH1AZ_down$Row.names, k9ac_up$Row.names))
k9ac_up_vs_KH1AZ_down <- data.frame(same,former,latter)
k9ac_up_vs_KH1AZ_down$type <- "k9ac_up_vs_KH1AZ_down"
#####以K4me3_up为参考与其它修饰进行overlap
same <- length(intersect(x=K4me3_up$Row.names, y = K4me2_down$Row.names))
former <- length(setdiff(K4me3_up$Row.names, K4me2_down$Row.names))
latter <- length(setdiff(K4me2_down$Row.names, K4me3_up$Row.names))
K4me3_up_vs_K4me2_down <- data.frame(same,former,latter)
K4me3_up_vs_K4me2_down$type <- "K4me3_up_vs_K4me2_down"
same <- length(intersect(x=K4me3_up$Row.names, y = K27ac_up$Row.names))
former <- length(setdiff(K4me3_up$Row.names, K27ac_up$Row.names))
latter <- length(setdiff(K27ac_up$Row.names, K4me3_up$Row.names))
K4me3_up_vs_K27ac_up <- data.frame(same,former,latter)
K4me3_up_vs_K27ac_up$type <- "K4me3_up_vs_K27ac_up"
same <- length(intersect(x=K4me3_up$Row.names, y = K27me3_down$Row.names))
former <- length(setdiff(K4me3_up$Row.names, K27me3_down$Row.names))
latter <- length(setdiff(K27me3_down$Row.names, K4me3_up$Row.names))
K4me3_up_vs_K27me3_down <- data.frame(same,former,latter)
K4me3_up_vs_K27me3_down$type <- "K4me3_up_vs_K27me3_down"
same <- length(intersect(x=K4me3_up$Row.names, y = K36me3_up$Row.names))
former <- length(setdiff(K4me3_up$Row.names, K36me3_up$Row.names))
latter <- length(setdiff(K36me3_up$Row.names, K4me3_up$Row.names))
K4me3_up_vs_K36me3_up <- data.frame(same,former,latter)
K4me3_up_vs_K36me3_up$type <- "K4me3_up_vs_K36me3_up"
same <- length(intersect(x=K4me3_up$Row.names, y = KH1AZ_down$Row.names))
former <- length(setdiff(K4me3_up$Row.names, KH1AZ_down$Row.names))
latter <- length(setdiff(KH1AZ_down$Row.names, K4me3_up$Row.names))
K4me3_up_vs_KH1AZ_down <- data.frame(same,former,latter)
K4me3_up_vs_KH1AZ_down$type <- "K4me3_up_vs_KH1AZ_down"
#####以K4me2_down为参考与其它修饰进行overlap
same <- length(intersect(x=K4me2_down$Row.names, y = K27ac_up$Row.names))
former <- length(setdiff(K4me2_down$Row.names, K27ac_up$Row.names))
latter <- length(setdiff(K27ac_up$Row.names, K4me2_down$Row.names))
K4me2_down_vs_K27ac_up <- data.frame(same,former,latter)
K4me2_down_vs_K27ac_up$type <- "K4me2_down_vs_K27ac_up"
same <- length(intersect(x=K4me2_down$Row.names, y = K27me3_down$Row.names))
former <- length(setdiff(K4me2_down$Row.names, K27me3_down$Row.names))
latter <- length(setdiff(K27me3_down$Row.names, K4me2_down$Row.names))
K4me2_down_vs_K27me3_down <- data.frame(same,former,latter)
K4me2_down_vs_K27me3_down$type <- "K4me2_down_vs_K27me3_down"
same <- length(intersect(x=K4me2_down$Row.names, y = K36me3_up$Row.names))
former <- length(setdiff(K4me2_down$Row.names, K36me3_up$Row.names))
latter <- length(setdiff(K36me3_up$Row.names, K4me2_down$Row.names))
K4me2_down_vs_K36me3_up <- data.frame(same,former,latter)
K4me2_down_vs_K36me3_up$type <- "K4me2_down_vs_K36me3_up"
same <- length(intersect(x=K4me2_down$Row.names, y = KH1AZ_down$Row.names))
former <- length(setdiff(K4me2_down$Row.names, KH1AZ_down$Row.names))
latter <- length(setdiff(KH1AZ_down$Row.names, K4me2_down$Row.names))
K4me2_down_vs_KH1AZ_down <- data.frame(same,former,latter)
K4me2_down_vs_KH1AZ_down$type <- "K4me2_down_vs_KH1AZ_down"
#####以K27ac_up为参考与其它修饰进行overlap
same <- length(intersect(x=K27ac_up$Row.names, y = K27me3_down$Row.names))
former <- length(setdiff(K27ac_up$Row.names, K27me3_down$Row.names))
latter <- length(setdiff(K27me3_down$Row.names, K27ac_up$Row.names))
K27ac_up_vs_K27me3_down <- data.frame(same,former,latter)
K27ac_up_vs_K27me3_down$type <- "K27ac_up_vs_K27me3_down"
same <- length(intersect(x=K27ac_up$Row.names, y = K36me3_up$Row.names))
former <- length(setdiff(K27ac_up$Row.names, K36me3_up$Row.names))
latter <- length(setdiff(K36me3_up$Row.names, K27ac_up$Row.names))
K27ac_up_vs_K36me3_up <- data.frame(same,former,latter)
K27ac_up_vs_K36me3_up$type <- "K27ac_up_vs_K36me3_up"
same <- length(intersect(x=K27ac_up$Row.names, y = KH1AZ_down$Row.names))
former <- length(setdiff(K27ac_up$Row.names, KH1AZ_down$Row.names))
latter <- length(setdiff(KH1AZ_down$Row.names, K27ac_up$Row.names))
K27ac_up_vs_KH1AZ_down <- data.frame(same,former,latter)
K27ac_up_vs_KH1AZ_down$type <- "K27ac_up_vs_KH1AZ_down"
#####以K27me3_down为参考与其它修饰进行overlap
same <- length(intersect(x=K27me3_down$Row.names, y = K36me3_up$Row.names))
former <- length(setdiff(K27me3_down$Row.names, K36me3_up$Row.names))
latter <- length(setdiff(K36me3_up$Row.names, K27me3_down$Row.names))
K27me3_down_vs_K36me3_up <- data.frame(same,former,latter)
K27me3_down_vs_K36me3_up$type <- "K27me3_down_vs_K36me3_up"
same <- length(intersect(x=K27me3_down$Row.names, y = KH1AZ_down$Row.names))
former <- length(setdiff(K27me3_down$Row.names, KH1AZ_down$Row.names))
latter <- length(setdiff(KH1AZ_down$Row.names, K27me3_down$Row.names))
K27me3_down_vs_KH1AZ_down <- data.frame(same,former,latter)
K27me3_down_vs_KH1AZ_down$type <- "K27me3_down_vs_KH1AZ_down"
#####以K36me3_up为参考与其它修饰进行overlap
same <- length(intersect(x=K36me3_up$Row.names, y = KH1AZ_down$Row.names))
former <- length(setdiff(K36me3_up$Row.names, KH1AZ_down$Row.names))
latter <- length(setdiff(KH1AZ_down$Row.names, K36me3_up$Row.names))
K36me3_up_vs_KH1AZ_down <- data.frame(same,former,latter)
K36me3_up_vs_KH1AZ_down$type <- "K36me3_up_vs_KH1AZ_down"

####绘制堆积柱形图
up <- rbind(k9ac_up_vs_K4me3_up, k9ac_up_vs_K4me2_down, k9ac_up_vs_K27ac_up, k9ac_up_vs_K27me3_down, k9ac_up_vs_K36me3_up, k9ac_up_vs_KH1AZ_down,
            K4me3_up_vs_K4me2_down, K4me3_up_vs_K27ac_up, K4me3_up_vs_K27me3_down, K4me3_up_vs_K36me3_up, K4me3_up_vs_KH1AZ_down,
            K4me2_down_vs_K27ac_up, 
            K4me2_down_vs_K27me3_down, K4me2_down_vs_K36me3_up, K4me2_down_vs_KH1AZ_down,
            K27ac_up_vs_K27me3_down, K27ac_up_vs_K36me3_up, K27ac_up_vs_KH1AZ_down,  K27me3_down_vs_K36me3_up, 
            K27me3_down_vs_KH1AZ_down, K36me3_up_vs_KH1AZ_down)
#宽变长
up<- up %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer( cols =  c("same":"latter"),
                names_to = 'stage',
                values_to = 'expr')
head(up)
up <- as.data.frame(up)
tg <- up %>%
  group_by(type) %>%
  summarise(N=expr) %>%
  mutate(freq = (N/sum(N)*100), Cumsum = cumsum(freq))
tg$fz <- paste(tg$type,tg$N,sep="_")
up$fz <-  paste(up$type,up$expr,sep="_")
mt <- merge(up,tg,"fz",all=TRUE)
mt <- mt[,c(-1,-2,-6,-7,-9)]
mt$pertage <- paste0(round(mt$freq,2),"%")
colnames(mt) <- c("type","stage","Number","freq","pertage")
mt$stage <- factor(mt$stage,levels=c("former","same","latter"),order=TRUE)
head(mt)
xx <- mt[which(mt$stage == "same"),]
xx <- xx[order(xx$stage,xx$freq),]
mt$type <- factor(mt$type, levels=unique(xx$type), ordered=TRUE)
pdf("激活类up抑制类down两两之间的交集.pdf",width=12)
ggplot(mt,aes(type,Number,fill=as.factor(stage))) +
  #geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  geom_bar(stat='identity',position='fill',show.legend = TRUE) +
  labs(y = 'Peak number') +
  scale_fill_manual(values=brewer.pal(8,"Set2")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1.2),
        axis.ticks.x=element_blank(),
        legend.position="top")
  #geom_text(aes(label = pertage), size = 3, hjust = 0.5, vjust = 1, position = "fill")
dev.off()

######保存数据
write.csv(mt,"统计数据_激活类up抑制类down两两之间的交集.csv")
xx <- mt[which(mt$stage == "same"),]
xx <- xx[order(xx$freq),]
mt$type <- factor(mt$type, levels=unique(xx$type), ordered=TRUE)
pdf("激活类up抑制类down两两之间的交集.pdf",width=12)
ggplot(mt,aes(type,Number,fill=as.factor(stage))) +
  #geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  geom_bar(stat='identity',position='fill',show.legend = TRUE) +
  labs(y = 'Peak number') +
  scale_fill_manual(values=brewer.pal(8,"Set2")) +
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()


#交集数据
#####以k27ac_up为参考与其它修饰进行overlap
same <- data.frame(intersect(x=k9ac_up$Row.names, y = K4me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, K4me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K4me3_up$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_K4me3_up <- rbind(same,former,latter)
k9ac_up_vs_K4me3_up$type <- "k9ac_up_vs_K4me3_up"
same <- data.frame(intersect(x=k9ac_up$Row.names, y = K4me2_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, K4me2_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K4me2_down$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_K4me2_down <- rbind(same,former,latter)
k9ac_up_vs_K4me2_down$type <- "k9ac_up_vs_K4me2_down"
same <- data.frame(intersect(x=k9ac_up$Row.names, y = K27ac_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, K27ac_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27ac_up$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_K27ac_up <- rbind(same,former,latter)
k9ac_up_vs_K27ac_up$type <- "k9ac_up_vs_K27ac_up"
same <- data.frame(intersect(x=k9ac_up$Row.names, y = K27me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, K27me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_down$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_K27me3_down <- rbind(same,former,latter)
k9ac_up_vs_K27me3_down$type <- "k9ac_up_vs_K27me3_down"
same <- data.frame(intersect(x=k9ac_up$Row.names, y = K27me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, K27me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_down$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_K27me3_down <- rbind(same,former,latter)
k9ac_up_vs_K27me3_down$type <- "k9ac_up_vs_K27me3_down"
same <- data.frame(intersect(x=k9ac_up$Row.names, y = K36me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, K36me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_up$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_K36me3_up <- rbind(same,former,latter)
k9ac_up_vs_K36me3_up$type <- "k9ac_up_vs_K36me3_up"
same <- data.frame(intersect(x=k9ac_up$Row.names, y = KH1AZ_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_up$Row.names, KH1AZ_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_down$Row.names, k9ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_up_vs_KH1AZ_down <- rbind(same,former,latter)
k9ac_up_vs_KH1AZ_down$type <- "k9ac_up_vs_KH1AZ_down"
#####以K4me3_up为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K4me3_up$Row.names, y = K4me2_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_up$Row.names, K4me2_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K4me2_down$Row.names, K4me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_up_vs_K4me2_down <- rbind(same,former,latter)
K4me3_up_vs_K4me2_down$type <- "K4me3_up_vs_K4me2_down"
same <- data.frame(intersect(x=K4me3_up$Row.names, y = K27ac_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_up$Row.names, K27ac_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27ac_up$Row.names, K4me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_up_vs_K27ac_up <- rbind(same,former,latter)
K4me3_up_vs_K27ac_up$type <- "K4me3_up_vs_K27ac_up"
same <- data.frame(intersect(x=K4me3_up$Row.names, y = K27me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_up$Row.names, K27me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_down$Row.names, K4me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_up_vs_K27me3_down <- rbind(same,former,latter)
K4me3_up_vs_K27me3_down$type <- "K4me3_up_vs_K27me3_down"
same <- data.frame(intersect(x=K4me3_up$Row.names, y = K36me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_up$Row.names, K36me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_up$Row.names, K4me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_up_vs_K36me3_up <- rbind(same,former,latter)
K4me3_up_vs_K36me3_up$type <- "K4me3_up_vs_K36me3_up"
same <- data.frame(intersect(x=K4me3_up$Row.names, y = KH1AZ_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_up$Row.names, KH1AZ_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_down$Row.names, K4me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_up_vs_KH1AZ_down <- rbind(same,former,latter)
K4me3_up_vs_KH1AZ_down$type <- "K4me3_up_vs_KH1AZ_down"
#####以K4me2_up为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K4me2_down$Row.names, y = K27ac_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_down$Row.names, K27ac_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
former$group <- "former_unique"
colnames(same) <- c("gene","group")
same$group <- "same"
latter <- data.frame(setdiff(K27ac_up$Row.names, K4me2_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
same$group <- "same"
colnames(same) <- c("gene","group")
K4me2_down_vs_K27ac_up <- rbind(same,former,latter)
K4me2_down_vs_K27ac_up$type <- "K4me2_down_vs_K27ac_up"
same <- data.frame(intersect(x=K4me2_down$Row.names, y = K27me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_down$Row.names, K27me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_down$Row.names, K4me2_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me2_down_vs_K27me3_down <- rbind(same,former,latter)
K4me2_down_vs_K27me3_down$type <- "K4me2_down_vs_K27me3_down"
same <- data.frame(intersect(x=K4me2_down$Row.names, y = K36me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_down$Row.names, K36me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_up$Row.names, K4me2_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me2_down_vs_K36me3_up <- rbind(same,former,latter)
K4me2_down_vs_K36me3_up$type <- "K4me2_down_vs_K36me3_up"
same <- data.frame(intersect(x=K4me2_down$Row.names, y = KH1AZ_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_down$Row.names, KH1AZ_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_down$Row.names, K4me2_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me2_down_vs_KH1AZ_down <- rbind(same,former,latter)
K4me2_down_vs_KH1AZ_down$type <- "K4me2_down_vs_KH1AZ_down"
#####以K27ac_up为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K27ac_up$Row.names, y = K27me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27ac_up$Row.names, K27me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_down$Row.names, K27ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27ac_up_vs_K27me3_down <- rbind(same,former,latter)
K27ac_up_vs_K27me3_down$type <- "K27ac_up_vs_K27me3_down"
same <- data.frame(intersect(x=K27ac_up$Row.names, y = K36me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27ac_up$Row.names, K36me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_up$Row.names, K27ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27ac_up_vs_K36me3_up <- rbind(same,former,latter)
K27ac_up_vs_K36me3_up$type <- "K27ac_up_vs_K36me3_up"
same <- data.frame(intersect(x=K27ac_up$Row.names, y = KH1AZ_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27ac_up$Row.names, KH1AZ_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_down$Row.names, K27ac_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27ac_up_vs_KH1AZ_down <- rbind(same,former,latter)
K27ac_up_vs_KH1AZ_down$type <- "K27ac_up_vs_KH1AZ_down"
#####以K27me3_down为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K27me3_down$Row.names, y = K36me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27me3_down$Row.names, K36me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_up$Row.names, K27me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27me3_down_vs_K36me3_up <- rbind(same,former,latter)
K27me3_down_vs_K36me3_up$type <- "K27me3_down_vs_K36me3_up"
same <- data.frame(intersect(x=K27me3_down$Row.names, y = KH1AZ_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27me3_down$Row.names, KH1AZ_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_down$Row.names, K27me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27me3_down_vs_KH1AZ_down <- rbind(same,former,latter)
K27me3_down_vs_KH1AZ_down$type <- "K27me3_down_vs_KH1AZ_down"
#####以K36me3_up为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K36me3_up$Row.names, y = KH1AZ_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K36me3_up$Row.names, KH1AZ_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_down$Row.names, K36me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K36me3_up_vs_KH1AZ_down <- rbind(same,former,latter)
K36me3_up_vs_KH1AZ_down$type <- "K36me3_up_vs_KH1AZ_down"
up <- rbind(k9ac_up_vs_K4me3_up, k9ac_up_vs_K4me2_down, k9ac_up_vs_K27ac_up, k9ac_up_vs_K27me3_down, k9ac_up_vs_K36me3_up, k9ac_up_vs_KH1AZ_down,
            K4me3_up_vs_K4me2_down, K4me3_up_vs_K27ac_up, K4me3_up_vs_K27me3_down, K4me3_up_vs_K36me3_up, K4me3_up_vs_KH1AZ_down,
            K4me2_down_vs_K27ac_up, 
            K4me2_down_vs_K27me3_down, K4me2_down_vs_K36me3_up, K4me2_down_vs_KH1AZ_down,
            K27ac_up_vs_K27me3_down, K27ac_up_vs_K36me3_up, K27ac_up_vs_KH1AZ_down,  K27me3_down_vs_K36me3_up, 
            K27me3_down_vs_KH1AZ_down, K36me3_up_vs_KH1AZ_down)
write.csv(up,"激活类up抑制类down两两之间的交集.csv")

#######激活类down抑制类up两两之间交集
#####以k27ac_down为参考与其它修饰进行overlap
same <- length(intersect(x=k9ac_down$Row.names, y = K4me3_down$Row.names))
former <- length(setdiff(k9ac_down$Row.names, K4me3_down$Row.names))
latter <- length(setdiff(K4me3_down$Row.names, k9ac_down$Row.names))
k9ac_down_vs_K4me3_down <- data.frame(same,former,latter)
k9ac_down_vs_K4me3_down$type <- "k9ac_down_vs_K4me3_down"
same <- length(intersect(x=k9ac_down$Row.names, y = K4me2_up$Row.names))
former <- length(setdiff(k9ac_down$Row.names, K4me2_up$Row.names))
latter <- length(setdiff(K4me2_up$Row.names, k9ac_down$Row.names))
k9ac_down_vs_K4me2_up <- data.frame(same,former,latter)
k9ac_down_vs_K4me2_up$type <- "k9ac_down_vs_K4me2_up"
same <- length(intersect(x=k9ac_down$Row.names, y = K27ac_down$Row.names))
former <- length(setdiff(k9ac_down$Row.names, K27ac_down$Row.names))
latter <- length(setdiff(K27ac_down$Row.names, k9ac_down$Row.names))
k9ac_down_vs_K27ac_down <- data.frame(same,former,latter)
k9ac_down_vs_K27ac_down$type <- "k9ac_down_vs_K27ac_down"
same <- length(intersect(x=k9ac_down$Row.names, y = K27me3_up$Row.names))
former <- length(setdiff(k9ac_down$Row.names, K27me3_up$Row.names))
latter <- length(setdiff(K27me3_up$Row.names, k9ac_down$Row.names))
k9ac_down_vs_K27me3_up <- data.frame(same,former,latter)
k9ac_down_vs_K27me3_up$type <- "k9ac_down_vs_K27me3_up"
same <- length(intersect(x=k9ac_down$Row.names, y = K27me3_up$Row.names))
former <- length(setdiff(k9ac_down$Row.names, K27me3_up$Row.names))
latter <- length(setdiff(K27me3_up$Row.names, k9ac_down$Row.names))
k9ac_down_vs_K27me3_up <- data.frame(same,former,latter)
k9ac_down_vs_K27me3_up$type <- "k9ac_down_vs_K27me3_up"
same <- length(intersect(x=k9ac_down$Row.names, y = K36me3_down$Row.names))
former <- length(setdiff(k9ac_down$Row.names, K36me3_down$Row.names))
latter <- length(setdiff(K36me3_down$Row.names, k9ac_down$Row.names))
k9ac_down_vs_K36me3_down <- data.frame(same,former,latter)
k9ac_down_vs_K36me3_down$type <- "k9ac_down_vs_K36me3_down"
same <- length(intersect(x=k9ac_down$Row.names, y = KH1AZ_up$Row.names))
former <- length(setdiff(k9ac_down$Row.names, KH1AZ_up$Row.names))
latter <- length(setdiff(KH1AZ_up$Row.names, k9ac_down$Row.names))
k9ac_down_vs_KH1AZ_up <- data.frame(same,former,latter)
k9ac_down_vs_KH1AZ_up$type <- "k9ac_down_vs_KH1AZ_up"
#####以K4me3_down为参考与其它修饰进行overlap
same <- length(intersect(x=K4me3_down$Row.names, y = K4me2_up$Row.names))
former <- length(setdiff(K4me3_down$Row.names, K4me2_up$Row.names))
latter <- length(setdiff(K4me2_up$Row.names, K4me3_down$Row.names))
K4me3_down_vs_K4me2_up <- data.frame(same,former,latter)
K4me3_down_vs_K4me2_up$type <- "K4me3_down_vs_K4me2_up"
same <- length(intersect(x=K4me3_down$Row.names, y = K27ac_down$Row.names))
former <- length(setdiff(K4me3_down$Row.names, K27ac_down$Row.names))
latter <- length(setdiff(K27ac_down$Row.names, K4me3_down$Row.names))
K4me3_down_vs_K27ac_down <- data.frame(same,former,latter)
K4me3_down_vs_K27ac_down$type <- "K4me3_down_vs_K27ac_down"
same <- length(intersect(x=K4me3_down$Row.names, y = K27me3_up$Row.names))
former <- length(setdiff(K4me3_down$Row.names, K27me3_up$Row.names))
latter <- length(setdiff(K27me3_up$Row.names, K4me3_down$Row.names))
K4me3_down_vs_K27me3_up <- data.frame(same,former,latter)
K4me3_down_vs_K27me3_up$type <- "K4me3_down_vs_K27me3_up"
same <- length(intersect(x=K4me3_down$Row.names, y = K36me3_down$Row.names))
former <- length(setdiff(K4me3_down$Row.names, K36me3_down$Row.names))
latter <- length(setdiff(K36me3_down$Row.names, K4me3_down$Row.names))
K4me3_down_vs_K36me3_down <- data.frame(same,former,latter)
K4me3_down_vs_K36me3_down$type <- "K4me3_down_vs_K36me3_down"
same <- length(intersect(x=K4me3_down$Row.names, y = KH1AZ_up$Row.names))
former <- length(setdiff(K4me3_down$Row.names, KH1AZ_up$Row.names))
latter <- length(setdiff(KH1AZ_up$Row.names, K4me3_down$Row.names))
K4me3_down_vs_KH1AZ_up <- data.frame(same,former,latter)
K4me3_down_vs_KH1AZ_up$type <- "K4me3_down_vs_KH1AZ_up"
#####以K4me2_up为参考与其它修饰进行overlap
same <- length(intersect(x=K4me2_up$Row.names, y = K27ac_down$Row.names))
former <- length(setdiff(K4me2_up$Row.names, K27ac_down$Row.names))
latter <- length(setdiff(K27ac_down$Row.names, K4me2_up$Row.names))
K4me2_up_vs_K27ac_down <- data.frame(same,former,latter)
K4me2_up_vs_K27ac_down$type <- "K4me2_up_vs_K27ac_down"
same <- length(intersect(x=K4me2_up$Row.names, y = K27me3_up$Row.names))
former <- length(setdiff(K4me2_up$Row.names, K27me3_up$Row.names))
latter <- length(setdiff(K27me3_up$Row.names, K4me2_up$Row.names))
K4me2_up_vs_K27me3_up <- data.frame(same,former,latter)
K4me2_up_vs_K27me3_up$type <- "K4me2_up_vs_K27me3_up"
same <- length(intersect(x=K4me2_up$Row.names, y = K36me3_down$Row.names))
former <- length(setdiff(K4me2_up$Row.names, K36me3_down$Row.names))
latter <- length(setdiff(K36me3_down$Row.names, K4me2_up$Row.names))
K4me2_up_vs_K36me3_down <- data.frame(same,former,latter)
K4me2_up_vs_K36me3_down$type <- "K4me2_up_vs_K36me3_down"
same <- length(intersect(x=K4me2_up$Row.names, y = KH1AZ_up$Row.names))
former <- length(setdiff(K4me2_up$Row.names, KH1AZ_up$Row.names))
latter <- length(setdiff(KH1AZ_up$Row.names, K4me2_up$Row.names))
K4me2_up_vs_KH1AZ_up <- data.frame(same,former,latter)
K4me2_up_vs_KH1AZ_up$type <- "K4me2_up_vs_KH1AZ_up"
#####以K27ac_down为参考与其它修饰进行overlap
same <- length(intersect(x=K27ac_down$Row.names, y = K27me3_up$Row.names))
former <- length(setdiff(K27ac_down$Row.names, K27me3_up$Row.names))
latter <- length(setdiff(K27me3_up$Row.names, K27ac_down$Row.names))
K27ac_down_vs_K27me3_up <- data.frame(same,former,latter)
K27ac_down_vs_K27me3_up$type <- "K27ac_down_vs_K27me3_up"
same <- length(intersect(x=K27ac_down$Row.names, y = K36me3_down$Row.names))
former <- length(setdiff(K27ac_down$Row.names, K36me3_down$Row.names))
latter <- length(setdiff(K36me3_down$Row.names, K27ac_down$Row.names))
K27ac_down_vs_K36me3_down <- data.frame(same,former,latter)
K27ac_down_vs_K36me3_down$type <- "K27ac_down_vs_K36me3_down"
same <- length(intersect(x=K27ac_down$Row.names, y = KH1AZ_up$Row.names))
former <- length(setdiff(K27ac_down$Row.names, KH1AZ_up$Row.names))
latter <- length(setdiff(KH1AZ_up$Row.names, K27ac_down$Row.names))
K27ac_down_vs_KH1AZ_up <- data.frame(same,former,latter)
K27ac_down_vs_KH1AZ_up$type <- "K27ac_down_vs_KH1AZ_up"
#####以K27me3_up为参考与其它修饰进行overlap
same <- length(intersect(x=K27me3_up$Row.names, y = K36me3_down$Row.names))
former <- length(setdiff(K27me3_up$Row.names, K36me3_down$Row.names))
latter <- length(setdiff(K36me3_down$Row.names, K27me3_up$Row.names))
K27me3_up_vs_K36me3_down <- data.frame(same,former,latter)
K27me3_up_vs_K36me3_down$type <- "K27me3_up_vs_K36me3_down"
same <- length(intersect(x=K27me3_up$Row.names, y = KH1AZ_up$Row.names))
former <- length(setdiff(K27me3_up$Row.names, KH1AZ_up$Row.names))
latter <- length(setdiff(KH1AZ_up$Row.names, K27me3_up$Row.names))
K27me3_up_vs_KH1AZ_up <- data.frame(same,former,latter)
K27me3_up_vs_KH1AZ_up$type <- "K27me3_up_vs_KH1AZ_up"
#####以K36me3_down为参考与其它修饰进行overlap
same <- length(intersect(x=K36me3_down$Row.names, y = KH1AZ_up$Row.names))
former <- length(setdiff(K36me3_down$Row.names, KH1AZ_up$Row.names))
latter <- length(setdiff(KH1AZ_up$Row.names, K36me3_down$Row.names))
K36me3_down_vs_KH1AZ_up <- data.frame(same,former,latter)
K36me3_down_vs_KH1AZ_up$type <- "K36me3_down_vs_KH1AZ_up"
####绘制堆积柱形图
down <- rbind(k9ac_down_vs_K4me3_down, k9ac_down_vs_K4me2_up, k9ac_down_vs_K27ac_down, k9ac_down_vs_K27me3_up, k9ac_down_vs_K36me3_down, k9ac_down_vs_KH1AZ_up,
              K4me3_down_vs_K4me2_up, K4me3_down_vs_K27ac_down, K4me3_down_vs_K27me3_up, K4me3_down_vs_K36me3_down, K4me3_down_vs_KH1AZ_up,
              K4me2_up_vs_K27ac_down, 
              K4me2_up_vs_K27me3_up, K4me2_up_vs_K36me3_down, K4me2_up_vs_KH1AZ_up,
              K27ac_down_vs_K27me3_up, K27ac_down_vs_K36me3_down, K27ac_down_vs_KH1AZ_up,  K27me3_up_vs_K36me3_down, 
              K27me3_up_vs_KH1AZ_up, K36me3_down_vs_KH1AZ_up)
#宽变长
down<- down %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer( cols =  c("same":"latter"),
                names_to = 'stage',
                values_to = 'expr')
head(down)
down <- as.data.frame(down)
tg <- down %>%
  group_by(type) %>%
  summarise(N=expr) %>%
  mutate(freq = (N/sum(N)*100), Cumsum = cumsum(freq))
tg$fz <- paste(tg$type,tg$N,sep="_")
down$fz <-  paste(down$type,down$expr,sep="_")
mt <- merge(down,tg,"fz",all=TRUE)
mt <- mt[,c(-1,-2,-6,-7,-9)]
mt$pertage <- paste0(round(mt$freq,2),"%")
colnames(mt) <- c("type","stage","Number","freq","pertage")
xx <- mt[which(mt$stage == "same"),]
xx <- xx[order(xx$stage,xx$freq),]
mt$type <- factor(mt$type, levels=unique(xx$type), ordered=TRUE)
mt$stage <- factor(mt$stage,levels=c("former","same","latter"),order=TRUE)
pdf("激活类down抑制类up两两之间的交集.pdf",width=12)
ggplot(mt,aes(type,Number,fill=as.factor(stage))) +
  #geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  geom_bar(stat='identity',position='fill',show.legend = TRUE) +
  labs(y = 'Peak number') +
  scale_fill_manual(values=brewer.pal(8,"Set2")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1.2),
        axis.ticks.x=element_blank(),
        legend.position="top")
#geom_text(aes(label = pertage), size = 3, hjust = 0.5, vjust = 1, position = "fill")
dev.off()
######保存数据
write.csv(mt,"统计数据_激活类down抑制类up两两之间的交集.csv")
#交集数据
#####以k27ac_down为参考与其它修饰进行overlap
same <- data.frame(intersect(x=k9ac_down$Row.names, y = K4me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, K4me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K4me3_down$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_K4me3_down <- rbind(same,former,latter)
k9ac_down_vs_K4me3_down$type <- "k9ac_down_vs_K4me3_down"
same <- data.frame(intersect(x=k9ac_down$Row.names, y = K4me2_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, K4me2_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K4me2_up$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_K4me2_up <- rbind(same,former,latter)
k9ac_down_vs_K4me2_up$type <- "k9ac_down_vs_K4me2_up"
same <- data.frame(intersect(x=k9ac_down$Row.names, y = K27ac_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, K27ac_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27ac_down$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_K27ac_down <- rbind(same,former,latter)
k9ac_down_vs_K27ac_down$type <- "k9ac_down_vs_K27ac_down"
same <- data.frame(intersect(x=k9ac_down$Row.names, y = K27me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, K27me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_up$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_K27me3_up <- rbind(same,former,latter)
k9ac_down_vs_K27me3_up$type <- "k9ac_down_vs_K27me3_up"
same <- data.frame(intersect(x=k9ac_down$Row.names, y = K27me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, K27me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_up$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_K27me3_up <- rbind(same,former,latter)
k9ac_down_vs_K27me3_up$type <- "k9ac_down_vs_K27me3_up"
same <- data.frame(intersect(x=k9ac_down$Row.names, y = K36me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, K36me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_down$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_K36me3_down <- rbind(same,former,latter)
k9ac_down_vs_K36me3_down$type <- "k9ac_down_vs_K36me3_down"
same <- data.frame(intersect(x=k9ac_down$Row.names, y = KH1AZ_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(k9ac_down$Row.names, KH1AZ_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_up$Row.names, k9ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
k9ac_down_vs_KH1AZ_up <- rbind(same,former,latter)
k9ac_down_vs_KH1AZ_up$type <- "k9ac_down_vs_KH1AZ_up"
#####以K4me3_down为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K4me3_down$Row.names, y = K4me2_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_down$Row.names, K4me2_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K4me2_up$Row.names, K4me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_down_vs_K4me2_up <- rbind(same,former,latter)
K4me3_down_vs_K4me2_up$type <- "K4me3_down_vs_K4me2_up"
same <- data.frame(intersect(x=K4me3_down$Row.names, y = K27ac_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_down$Row.names, K27ac_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27ac_down$Row.names, K4me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_down_vs_K27ac_down <- rbind(same,former,latter)
K4me3_down_vs_K27ac_down$type <- "K4me3_down_vs_K27ac_down"
same <- data.frame(intersect(x=K4me3_down$Row.names, y = K27me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_down$Row.names, K27me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_up$Row.names, K4me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_down_vs_K27me3_up <- rbind(same,former,latter)
K4me3_down_vs_K27me3_up$type <- "K4me3_down_vs_K27me3_up"
same <- data.frame(intersect(x=K4me3_down$Row.names, y = K36me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_down$Row.names, K36me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_down$Row.names, K4me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_down_vs_K36me3_down <- rbind(same,former,latter)
K4me3_down_vs_K36me3_down$type <- "K4me3_down_vs_K36me3_down"
same <- data.frame(intersect(x=K4me3_down$Row.names, y = KH1AZ_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me3_down$Row.names, KH1AZ_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_up$Row.names, K4me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me3_down_vs_KH1AZ_up <- rbind(same,former,latter)
K4me3_down_vs_KH1AZ_up$type <- "K4me3_down_vs_KH1AZ_up"
#####以K4me2_down为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K4me2_up$Row.names, y = K27ac_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_up$Row.names, K27ac_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
former$group <- "former_unique"
colnames(same) <- c("gene","group")
same$group <- "same"
latter <- data.frame(setdiff(K27ac_down$Row.names, K4me2_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
same$group <- "same"
colnames(same) <- c("gene","group")
K4me2_up_vs_K27ac_down <- rbind(same,former,latter)
K4me2_up_vs_K27ac_down$type <- "K4me2_up_vs_K27ac_down"
same <- data.frame(intersect(x=K4me2_up$Row.names, y = K27me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_up$Row.names, K27me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_up$Row.names, K4me2_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me2_up_vs_K27me3_up <- rbind(same,former,latter)
K4me2_up_vs_K27me3_up$type <- "K4me2_up_vs_K27me3_up"
same <- data.frame(intersect(x=K4me2_up$Row.names, y = K36me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_up$Row.names, K36me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_down$Row.names, K4me2_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me2_up_vs_K36me3_down <- rbind(same,former,latter)
K4me2_up_vs_K36me3_down$type <- "K4me2_up_vs_K36me3_down"
same <- data.frame(intersect(x=K4me2_up$Row.names, y = KH1AZ_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K4me2_up$Row.names, KH1AZ_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_up$Row.names, K4me2_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K4me2_up_vs_KH1AZ_up <- rbind(same,former,latter)
K4me2_up_vs_KH1AZ_up$type <- "K4me2_up_vs_KH1AZ_up"
#####以K27ac_down为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K27ac_down$Row.names, y = K27me3_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27ac_down$Row.names, K27me3_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K27me3_up$Row.names, K27ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27ac_down_vs_K27me3_up <- rbind(same,former,latter)
K27ac_down_vs_K27me3_up$type <- "K27ac_down_vs_K27me3_up"
same <- data.frame(intersect(x=K27ac_down$Row.names, y = K36me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27ac_down$Row.names, K36me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_down$Row.names, K27ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27ac_down_vs_K36me3_down <- rbind(same,former,latter)
K27ac_down_vs_K36me3_down$type <- "K27ac_down_vs_K36me3_down"
same <- data.frame(intersect(x=K27ac_down$Row.names, y = KH1AZ_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27ac_down$Row.names, KH1AZ_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_up$Row.names, K27ac_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27ac_down_vs_KH1AZ_up <- rbind(same,former,latter)
K27ac_down_vs_KH1AZ_up$type <- "K27ac_down_vs_KH1AZ_up"
#####以K27me3_up为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K27me3_up$Row.names, y = K36me3_down$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27me3_up$Row.names, K36me3_down$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(K36me3_down$Row.names, K27me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27me3_up_vs_K36me3_down <- rbind(same,former,latter)
K27me3_up_vs_K36me3_down$type <- "K27me3_up_vs_K36me3_down"
same <- data.frame(intersect(x=K27me3_up$Row.names, y = KH1AZ_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K27me3_up$Row.names, KH1AZ_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_up$Row.names, K27me3_up$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K27me3_up_vs_KH1AZ_up <- rbind(same,former,latter)
K27me3_up_vs_KH1AZ_up$type <- "K27me3_up_vs_KH1AZ_up"
#####以K36me3_down为参考与其它修饰进行overlap
same <- data.frame(intersect(x=K36me3_down$Row.names, y = KH1AZ_up$Row.names))
same$group <- "same"
colnames(same) <- c("gene","group")
former <- data.frame(setdiff(K36me3_down$Row.names, KH1AZ_up$Row.names))
former$group <- "former_unique"
colnames(former) <- c("gene","group")
latter <- data.frame(setdiff(KH1AZ_up$Row.names, K36me3_down$Row.names))
latter$group <- "latter_unique"
colnames(latter) <- c("gene","group")
K36me3_down_vs_KH1AZ_up <- rbind(same,former,latter)
K36me3_down_vs_KH1AZ_up$type <- "K36me3_down_vs_KH1AZ_up"
down <- rbind(k9ac_down_vs_K4me3_down, k9ac_down_vs_K4me2_up, k9ac_down_vs_K27ac_down, k9ac_down_vs_K27me3_up, k9ac_down_vs_K36me3_down, k9ac_down_vs_KH1AZ_up,
              K4me3_down_vs_K4me2_up, K4me3_down_vs_K27ac_down, K4me3_down_vs_K27me3_up, K4me3_down_vs_K36me3_down, K4me3_down_vs_KH1AZ_up,
              K4me2_up_vs_K27ac_down, 
              K4me2_up_vs_K27me3_up, K4me2_up_vs_K36me3_down, K4me2_up_vs_KH1AZ_up,
              K27ac_down_vs_K27me3_up, K27ac_down_vs_K36me3_down, K27ac_down_vs_KH1AZ_up,  K27me3_up_vs_K36me3_down, 
              K27me3_up_vs_KH1AZ_up, K36me3_down_vs_KH1AZ_up)
write.csv(down,"激活类down抑制类up两两之间的交集.csv")


######以基因为目标，构建不同变化类型的矩阵
install.packages("UpSetR")
#载入包
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
#构建数据
setwd("D:/高粱/ChIP/定量/CLvsCR/DEPR2/DEG_DMG_correlation/DEG_DMG_correlationR2")
#绘制环状柱形图


#####绘制胁迫相关基因圆形热图
library(circlize )
setwd("D:/高粱/GO和KEGG筛选后")
mat1 <- read.xlsx("D:/高粱/文章撰写/胁迫基因热图R1.xlsx",sheet=1)
mat2 <- read.xlsx("D:/高粱/文章撰写/胁迫基因热图R1.xlsx",sheet=2)
head(mat1)
row.names(mat1) <- paste(mat1[,2],mat1[,3],sep="-")
row.names(mat2) <- paste(mat2[,2],mat1[,3],sep="-")
at1 <- mat1[,c(-1:-4)]
at2 <- mat2[,c(-1:-4)]
head(at1)
head(at2)
st1 <- as.factor(mat1[,2])
st2 <- as.factor(mat2[,2])
head(st1)
head(st2)


circos.clear()
#circos.par(gap.after = c(10))##为添加列名留出空间
circos.heatmap.initialize(at1, split = st1)
#col_fun <- colorRamp2(c(-2, 0, 14), c("#af8dc3", "#f7f7f7", "#7fbf7b"))
col_fun <- colorRamp2(c(-1.3, 0, 14.08), c("#1F77B4", "#FFFFFF", "#F64971"))
circos.heatmap(mat1[,4], col = col_fun,dend.track.height = 0.2,rownames.side="inside")
col_fun2 <- colorRamp2(c(-2, 0, 4), c("#c6ffdd", "#fbd786", "#f7797d"))
circos.heatmap(at1, col = col_fun2,rownames.side="inside")



circos.track(ylim=c(0,1), track.index = 2, ##将列名添加在第二个轨道（就是热图所在的环形轨道）
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 cn = "mRNA" ##取得列名
                 n = length(mat2)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), ##x轴坐标
                             1:n - convert_y(0.5, "mm"), ##y轴坐标
                             cn, ##输入要展示的列名
                             cex = 0.25, ##列名的大小
                             adj = c(0, 0.5),
                             facing = "inside")
               }
             }, bg.border = NA)



p1 <- pheatmap(cg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#4B92C3", "#629FCA", "#78ADD2", "#8FBBD9",
                                          "#A5C8E1", "#BBD6E8", "#D2E3F0", "#E8F1F7", "#FFFFFF",
                                          "#FEF4F6", "#FDE2E8", "#FCD1DB", "#FBC0CE", "#FAAFC0",
                                          "#FA9EB3", "#F98DA6", "#F87C98", "#F76B8B", "#F65A7E", "#F64971"))(50), #热图色块颜色是从蓝到红分为100个等级
               border_color = "black",
               display_numbers = T,
               annotation_col = annotation_col, 
               annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               show_rownames = TRUE)
rg4 <- c4[,c(1,2,3,4)]
rg4$name <- rg4$symbol
name <- rg4[,3]
rg4 <- data.frame(rg4[,c(-1:-3,-5)])
row.names(rg4) <- name
colnames(rg4) <- "Expression"
head(rg4)
p2 <- pheatmap(rg4, 
               cluster_cols = F, 
               cluster_rows = F, 
               color = colorRampPalette(c("#1F77B4", "#3786BC", "#5095C4", "#69A4CD", "#82B3D5", 
                                          "#9BC2DD", "#B4D1E6", "#CDE0EE", "#E6EFF6", "#FFFFFF",
                                          "#FEF7F9", "#FEEFF3", "#FDE8ED", "#FDE0E7", "#FDD9E1",
                                          "#FCD1DB", "#FCC9D5", "#FCC2CF", "#FBBAC9", "#FBB3C3", "#FAABBD", "#FAA3B8", "#FA9CB2",
                                          "#F994AC", "#F98DA6", "#F985A0", "#F87E9A", "#F87694", "#F76E8E", "#F76788", "#F75F82",
                                          "#F6587C", "#F65076", "#F64971"))(30), #热图色块颜色是从蓝到红分为100个等级
               #color = colorRampPalette(brewer.pal(11, "PiYG"))(50),
               #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
               #color = viridis(4),
               #color = viridis(8, option = "G"),
               #legend_breaks = c(-10, 0, 10),
               #annotation_row = annotation_row, 
               #annotation_colors = ann_colors,
               border_color = "black",
               show_rownames = F,
               display_numbers =T)
library(dplyr)
cowplot::plot_grid(p1$gtable, p2$gtable, labels=c("modification","expression")) %>%  ggsave("root_ABA_stress_related_heatmapr.pdf",.,width=500,height=210, units="mm")





circos.clear()
circos.heatmap.initialize(at2, split = st2)
col_fun <- colorRamp2(c(-4, 0, 9), c("#c6ffdd", "#fbd786", "#f7797d"))
circos.heatmap(mat2[,4], col = col_fun)
col_fun2 <- colorRamp2(c(-3, 0, 3), c("#c6ffdd", "#fbd786", "#f7797d"))
circos.heatmap(at2, col = col_fun2,rownames.side="inside")


CELL_META$row_order(simply CELL_META$order)

circos.track(
  track.index = get.current.track.index(), 
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$cell.ylim[2] + convert_y(2, "mm"),
      paste0("this is group ", CELL_META$sector.index),
      facing = "bending.inside", cex = 0.8, 
      adj = c(0.5, 0), niceFacing = TRUE
    )
  }, bg.border = NA)




col_fun = colorRamp2(c(-2, 0, 14), c("#26B9CB", "#FFFFFF", "#B72865"))##设置热图颜色

circos.par(gap.after = c(10))##为添加列名留出空间

circos.heatmap(mat1[, column_od], ##将列聚类后重新排序的矩阵
               
               col = col_fun1, ##设置颜色
               
               dend.side = "inside",##树状图在圈内
               
               rownames.side = "outside",##行名在圈外
               
               dend.track.height = 0.2,
               
               dend.callback = function(dend, m, si) {
                 
                 # when k = 1, it renders one same color for the whole dendrogram
                 
                 color_branches(dend, k = 4, col = 2:5)##对树状图进行着色
                 
               }
               
)

circos.track(track.index = 2, ##将列名添加在第二个轨道（就是热图所在的环形轨道）
             
             panel.fun = function(x, y) {
               
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 
                 cn = colnames(mat1[, column_od])##取得列名
                 
                 n = length(cn)
                 
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), ##x轴坐标
                             
                             1:n - convert_y(0.5, "mm"), ##y轴坐标
                             
                             cn, ##输入要展示的列名
                             
                             cex = 0.25, ##列名的大小
                             
                             adj = c(0, 0.5),
                             
                             facing = "inside")
                 
               }
               
             }, bg.border = NA)



circos.clear()





col_fun <- colorRamp2(c(-2, 0, 2), c("#fc8d59", "#ffffbf", "#91bfdb"))
circos.heatmap(
  mat, split = split, col = col_fun, 
  rownames.side = "inside"
)
circos.track(
  track.index = get.current.track.index(), 
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$cell.ylim[2] + convert_y(2, "mm"),
      paste0("this is group ", CELL_META$sector.index),
      facing = "bending.inside", cex = 0.8, 
      adj = c(0.5, 0), niceFacing = TRUE
    )
  }, bg.border = NA)
circos.clear()


####提取表型数据
setwd("D:/高粱/单倍型分析")
phe <- read.xlsx("breed_information.xlsx")
head(phe)
t <- data.frame()
for (i in 1:289) {
  a <- read.xlsx("phenotype.xlsx",sheet=i)
  b <- a[-1:-20,2:5]
  b$type <- rep(a[11,1],nrow(b))
  colnames(b) <- b[1,]
  b <- b[-1,]
  colnames(b) <- c("category","descriptor","value","environment","type")
  t <- rbind(t,b)
}
head(t)
write.csv(t,"total_phenotype.csv")
##分别导出每个性状的表型情况
#excel中删除了一些非数字表型
t <- read.csv("total_phenotype.csv")
name <- unique(t$descriptor)
for (i in name) {
  x <- t[which(t$descriptor == i),]
  x <- x[,c(6,6,7)]
  x$type <- gsub("Phenotype\\/","",x$type)
  x$type.1 <- gsub("Phenotype\\/","",x$type.1)
  write.table(x,paste0(i,".txt"),sep="\t",row.names=F,quote = FALSE)
}
head(x)
###选择short_day_anthesis vigor infloresence_exsertion panicle_length rust stalk basil_tiller nodal_tiller
#plant_color sprouting_tendency yield_potential anthracnose height_uniformity mid-rid_color stalk_juiciness panicle_erectness
#plant_height seed_weight
#kernel_color flowering_rating 待定
setwd("D:/高粱/单倍型分析/单个表型文件/")
name <- list.files(path = ".") 
###获得所有材料的名称
name <- list.files(path = ".") 
x <- read.table(name[[1]],header=T)
head(x)
for(i in 2:20) {
  y <- read.table(name[[i]],header=T)
  x <- merge(x,y,by="type",all=T)
  i=i+1
}
head(x)
#x <- unique(x)
#x <- x[!duplicated(x$type),]
mm <- na.omit(x)
mm <- mm[,1:2]
mm$V2 <- mm$type
yz <- read.table("../output.snp.filtered-temporary.fam")
yz <- merge(mm,yz,by="V2",all=F)
yz <- yz[,1:2]
head(yz)
name <- list.files(path = ".") 
for (i in name) {
  sam <- read.table(i,header=T)
  xyz <- merge(yz,sam,by="type",all=F)
  xyz <- xyz[,c(1,2,4)]
  write.table(xyz,paste0(i,"R1.txt"),row.names=FALSE,na="",col.names=FALSE,sep="\t",quote=FALSE)
}

#####单倍型分析结果的可视化
setwd("D:/高粱/单倍型分析/单倍型分析结果")
library(openxlsx)
name <- c("C4单倍型分析情况","PP2C单倍型分析结果","POD单倍型分析结果","LEA单倍型分析结果")
hat <- data.frame()
for (i in name) {
  ha <- read.xlsx("单倍型分析结果总表.xlsx",sheet=i,colNames = TRUE)
  head(ha)
  ha <-ha[order(ha[,1],ha[,2],ha[,4]),]
  #ha$fz <- paste(ha$trait,ha$gene,sep="-")
  ha$type <- rep(i,nrow(ha))
  #ha1 <-ha[!duplicated(ha$fz),]
  head(ha)
  hat <- rbind(hat,ha)
}
head(hat)
write.csv(hat,"select_type_gene_hypotype.csv")
list <- unique(hat$SNP_id)
write.table(list,"list.txt",row.names=F,quote = FALSE,col.names=FALSE)
###集群上的R
cd /public/home/tllu/sorghum
R
bim <- read.table("output.snp.filtered.bim")
list <- read.table("list.txt")
list$V2 <- list[,1]
snp <- merge(bim,list,by="V2",all=F)
snp$start <- snp$V4-1
snp$end <- snp$V4+1
snp <- snp[,c(2,8,9,1,5,6)]
write.table(snp,"selected_snp.bed",row.names=F,quote = FALSE,col.names=FALSE)
###对SNP进行注释
perl -p -i -e 's/ /\t/g' selected_snp.bed
sort -k1,1 -k2n,2 selected_snp.bed | bedtools closest -D ref -t all -mdb all -a - -b /public/home/chaohe/sorghum/chip/align/Sorghum_geneR1.bed >selected_snp.txt
awk '{if($12=="+" && $13 >= 0) print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$12"\t"$13"\t"$3-$8; else if($12 =="+" && $13 < 0) print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$12"\t"$13"\t"$2-$8}' selected_snp.txt | uniq  >minus.txt
awk '{if($12=="-" && $13 >= 0) print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$12"\t"$13"\t"$9-$3; else if($12 =="-" && $13 < 0) print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$12"\t"$13"\t"$9-$2}' selected_snp.txt | uniq  >plus.txt
cat minus.txt plus.txt >selected_snp_gene.peaks.txt
#注释
setwd("D:/高粱/单倍型分析/单倍型分析结果")
narrow <- read.table("selected_snp_gene.peaks.txt")
head(narrow)
narrow$fz<-abs(narrow[,8])
narrow<-narrow[order(narrow[,1],narrow[,4],narrow[,9]),]
#narrow<-tapply(narrow[,4],INDEnarrow=abs(narrow[,6]),FUN=min)   #求相同基因的平均peak value
narrow<-narrow[!duplicated(narrow[,4]),]
narrow<-narrow[,-9]
narrow<-narrow[order(narrow[,1],narrow[,2],narrow[,8]),]
genebody<-unique(narrow[which((narrow[,7]==0) & (narrow[,8]>= 0)),])
gene_downstream<-unique(narrow[which((narrow[,7] != 0) & (narrow[,8]>0)),])
promoter <- unique(narrow[which(((narrow[,7]!=0) & (narrow[,8]>= -2000)  & (narrow[,8] < 0))),])
intergeric <- unique(narrow[which((narrow[,7]!=0) & (narrow[,8]< -2000)),])
genebody$type<-c(rep("genebody",nrow(genebody)))
gene_downstream$type<-c(rep("gene_downstream",nrow(gene_downstream)))
promoter$type<-c(rep("promoter",nrow(promoter)))  
intergeric$type<-c(rep("intergeric",nrow(intergeric)))  
T<-rbind(genebody, gene_downstream, promoter, intergeric)
head(T)
colnames(T) <- c("chr","start","end","snp_id","gene","strand","dis","dis_TSS","element")
write.csv(T,"select_SNP_annotation.csv",row.names = F)

#####调取相关基因的注释信息
head(hat)
ha <- read.xlsx("单倍型分析结果总表.xlsx",sheet="2kb启动子区域的SNP",colNames = TRUE)
head(ha)
hata <- merge(hat,ha,by="SNP_id",all=F)
head(hata)
hata <- hata[order(hata[,2],hata[,9],hata[,4]),]
hata$fz <- paste(hata$gene.x,hata$trait,sep="-")
ha1 <-hata[!duplicated(hata$fz),]
head(ha1)
write.csv(ha1,"单倍型绘图SNP.csv",row.names=F)

#####18个SNP绘图数据准备
#order1.txt
chr_15749743	panicle_length	SORBI_3003G149200	0.005950476	LEA单倍型分析结果
chr_9912327	panicle_length	SORBI_3004G105100	8.13E-06	POD单倍型分析结果
chr_58971762	panicle_length	SORBI_3009G255900	0.004327919	PP2C单倍型分析结果
#order2.txt
chr_66991701	plant_color	SORBI_3004G338000	0.004951091	C4单倍型分析情况
#order3.txt
chr_66991701	plant_height	SORBI_3004G338000	3.91E-06	C4单倍型分析情况
chr_15748670	plant_height	SORBI_3003G149200	2.99E-06	LEA单倍型分析结果
chr_9913319	plant_height	SORBI_3004G105100	0.000103366	POD单倍型分析结果
chr_58971098	plant_height	SORBI_3009G255900	1.25E-12	PP2C单倍型分析结果
#order4.txt
chr_66991701	short_day_anthesis	SORBI_3004G338000	0.000671105	C4单倍型分析情况
chr_15748683	short_day_anthesis	SORBI_3003G149200	0.004782215	LEA单倍型分析结果
chr_9913689	short_day_anthesis	SORBI_3004G105100	0.003260309	POD单倍型分析结果
#order5.txt
chr_66990742	stalk_juiciness	SORBI_3004G338000	1.41E-05	C4单倍型分析情况
chr_15749834	stalk_juiciness	SORBI_3003G149200	8.17E-11	LEA单倍型分析结果
chr_9912327	stalk_juiciness	SORBI_3004G105100	8.82E-13	POD单倍型分析结果
chr_58971972	stalk_juiciness	SORBI_3009G255900	2.68E-07	PP2C单倍型分析结果
#order6.txt
chr_15750069	stalk_waxiness	SORBI_3003G149200	0.004440549	LEA单倍型分析结果
#order7.txt
chr_66990742	yield_potential	SORBI_3004G338000	2.09E-08	C4单倍型分析情况
chr_9912327	yield_potential	SORBI_3004G105100	8.14E-17	POD单倍型分析结果
chr_58971762	yield_potential	SORBI_3009G255900	1.07E-09	PP2C单倍型分析结果
####调取这些SNP的基因型和表型
cd /public/home/tllu/sorghum/panicle_length
plink --bfile target_snp --extract ../order1.txt  --recode --out SNP_C4_panicle_length
paste -d'\t' SNP_C4_panicle_length.ped panicle_length.txt >SNP-C4_panicle_length.txt
cd /public/home/tllu/sorghum/plant_color
plink --bfile target_snp --extract ../order2.txt  --recode --out SNP_C4_plant_color
paste -d'\t' SNP_C4_plant_color.ped plant_color.txt >SNP-C4_plant_color.txt
cd /public/home/tllu/sorghum/plant_height
plink --bfile target_snp --extract ../order3.txt  --recode --out SNP_C4_plant_height
paste -d'\t' SNP_C4_plant_height.ped plant_height.txt >SNP-C4_plant_height.txt
cd /public/home/tllu/sorghum/short_day_anthesis
plink --bfile target_snp --extract ../order4.txt  --recode --out SNP_C4_short_day_anthesis
paste -d'\t' SNP_C4_short_day_anthesis.ped short_day_anthesis.txt >SNP-C4_short_day_anthesis.txt
cd /public/home/tllu/sorghum/stalk_juiciness
plink --bfile target_snp --extract ../order5.txt  --recode --out SNP_C4_stalk_juiciness
paste -d'\t' SNP_C4_stalk_juiciness.ped stalk_juiciness.txt >SNP-C4_stalk_juiciness.txt
cd /public/home/tllu/sorghum/stalk_waxiness
plink --bfile target_snp --extract ../order6.txt  --recode --out SNP_C4_stalk_waxiness
paste -d'\t' SNP_C4_stalk_waxiness.ped stalk_waxiness.txt >SNP-C4_stalk_waxiness.txt
cd /public/home/tllu/sorghum/yield_potential
plink --bfile target_snp --extract ../order7.txt  --recode --out SNP_C4_yield_potential
paste -d'\t' SNP_C4_yield_potential.ped yield_potential.txt >SNP-C4_yield_potential.txt

#绘图
name <- c("panicle_length","plant_color","plant_height","short_day_anthesis","stalk_juiciness","stalk_waxiness","yield_potential")
i=name[[1]]
a1 <- read.table(paste0("SNP-C4_",i,".txt"))
b1 <- a1[,c(1,7,15)]
b1$type <- rep("chr_15749743",nrow(b1))
b1$cluster <- "LEA"
b2 <- a1[,c(1,9,15)]
b2$type <- rep("chr_9912327",nrow(b2))
b2$cluster <- "POD"
b3 <- a1[,c(1,11,15)]
b3$type <- rep("chr_9912327",nrow(b3))
b3$cluster <- "PP2C"
colnames(b1) <- colnames(b2) <- colnames(b3) <- c("sample","Haplotype","value","type","cluster")
a1 <- rbind(b1,b2,b3)
a1$pattern <- "Paniche_Length"

i=name[[2]]
a2 <- read.table(paste0("SNP-C4_",i,".txt"))
head(a2)
a2 <- a2[,c(1,7,11)]
a2$type <- rep("chr_66991701",nrow(a2))
a2$cluster <- "C4"
colnames(a2) <- c("sample","Haplotype","value","type","cluster")
a2$pattern <- "Plant_Color"

i=name[[3]]
a3 <- read.table(paste0("SNP-C4_",i,".txt"))
head(a3)
b1 <- a3[,c(1,7,17)]
b1$type <- rep("chr_66991701",nrow(b1))
b1$cluster <- "C4"
b2 <- a3[,c(1,9,17)]
b2$type <- rep("chr_15748670",nrow(b2))
b2$cluster <- "LEA"
b3 <- a3[,c(1,11,17)]
b3$type <- rep("chr_9913319",nrow(b3))
b3$cluster <- "POD"
b4 <- a3[,c(1,13,17)]
b4$type <- rep("chr_58971098",nrow(b4))
b4$cluster <- "PP2C"
colnames(b1) <- colnames(b2) <- colnames(b3) <- colnames(b4) <- c("sample","Haplotype","value","type","cluster")
a3 <- rbind(b1,b2,b3,b4)
a3$pattern <- "Plant_Height"

i=name[[4]]
a4 <- read.table(paste0("SNP-C4_",i,".txt"))
b1 <- a4[,c(1,7,15)]
b1$type <- rep("chr_66991701",nrow(b1))
b1$cluster <- "C4"
b2 <- a4[,c(1,9,15)]
b2$type <- rep("chr_15748683",nrow(b2))
b2$cluster <- "LEA"
b3 <- a4[,c(1,11,15)]
b3$type <- rep("chr_9913689",nrow(b3))
b3$cluster <- "POD"
colnames(b1) <- colnames(b2) <- colnames(b3) <- c("sample","Haplotype","value","type","cluster")
a4 <- rbind(b1,b2,b3)
a4$pattern <- "Short_Day_Anthesis"

i=name[[5]]
a5 <- read.table(paste0("SNP-C4_",i,".txt"))
head(a5)
b1 <- a5[,c(1,7,17)]
b1$type <- rep("chr_66991701",nrow(b1))
b1$cluster <- "C4"
b2 <- a5[,c(1,9,17)]
b2$type <- rep("chr_15748670",nrow(b2))
b2$cluster <- "LEA"
b3 <- a5[,c(1,11,17)]
b3$type <- rep("chr_9913319",nrow(b3))
b3$cluster <- "POD"
b4 <- a5[,c(1,13,17)]
b4$type <- rep("chr_58971098",nrow(b4))
b4$cluster <- "PP2C"
colnames(b1) <- colnames(b2) <- colnames(b3) <- colnames(b4) <- c("sample","Haplotype","value","type","cluster")
a5 <- rbind(b1,b2,b3,b4)
a5$pattern <- "Stalk_Juiciness"

i=name[[6]]
a6 <- read.table(paste0("SNP-C4_",i,".txt"))
head(a6)
a6 <- a6[,c(1,7,11)]
a6$type <- rep("chr_15750069",nrow(a6))
a6$cluster <- "LEA"
colnames(a6) <- c("sample","Haplotype","value","type","cluster")
a6$pattern <- "Stalk_Waxiness"

i=name[[7]]
a7 <- read.table(paste0("SNP-C4_",i,".txt"))
b1 <- a7[,c(1,7,15)]
b1$type <- rep("chr_66990742",nrow(b1))
b1$cluster <- "C4"
b2 <- a7[,c(1,9,15)]
b2$type <- rep("chr_9912327",nrow(b2))
b2$cluster <- "POD"
b3 <- a7[,c(1,11,15)]
b3$type <- rep("chr_58971762",nrow(b3))
b3$cluster <- "PP2C"
colnames(b1) <- colnames(b2) <- colnames(b3) <- c("sample","Haplotype","value","type","cluster")
a7 <- rbind(b1,b2,b3)
a7$pattern <- "Yield_Potential"
#合并
at <- rbind(a1,a2,a3,a4,a5,a6,a7)
###箱式图
head(at)
at <- at[which(at$Haplotype != "0"),]
at <- at[which(at$Haplotype != "*"),]
#a1$Haplotype <- gsub("A","H1",a1$Haplotype)
#a1$Haplotype <- gsub("G","H2",a1$Haplotype)
plot_list=list()
att <- at[which(at$pattern != "Stalk_Waxiness"),]
name <- unique(att$type)
for (i in name) {
  x <- att[which(att$type == i),]
  if (length(unique(x$Haplotype)) == 2) {
    my_comparisons <- list(c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[2]]))
  } else if (length(unique(x$Haplotype)) == 3) {
    my_comparisons <- list(c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[2]]),c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[3]]),
                           c(unique(x$Haplotype)[[2]],unique(x$Haplotype)[[3]]))
  } else if (length(unique(x$Haplotype)) == 4) {
    my_comparisons <- list(c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[2]]),c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[3]]),
                           c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[4]]),c(unique(x$Haplotype)[[2]],unique(x$Haplotype)[[3]]),
                           c(unique(x$Haplotype)[[2]],unique(x$Haplotype)[[4]]),c(unique(x$Haplotype)[[3]],unique(x$Haplotype)[[4]]))
  }
  p1 <- ggboxplot(x, x = "Haplotype", y = "value",
                  color = "Haplotype", palette = "jco",
                  add = c("jitter"),
                  bxp.errorbar=TRUE,width = 0.5,
                  bxp.errorbar.width = 0.1,
                  add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
    #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
    #scale_fill_manual(values=brewer.pal(12,"Paired")) +
    scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
    ggtitle(x$cluster) + 
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    xlab(x$type) + ylab(x$pattern) +
    stat_compare_means(comparisons = my_comparisons,
                       #aes(group = Haplotype),
                       #label = "p.signif",
                       method = "t.test")
  plot_list[[i]] <- p1
}
###组图
library(gridExtra)
grid.arrange(plot_list[[1]],plot_list[[4]], 
             plot_list[[7]],plot_list[[3]], 
             plot_list[[9]],plot_list[[2]],
             plot_list[[5]],plot_list[[8]],
             plot_list[[6]], plot_list[[10]],
             nrow=2,ncol=5)     %>%  ggsave("单倍型分析结果.pdf",.,width=210,height=130, units="mm")


#####可视化单倍型分析结果 seed weight 和 产量潜能 相关的POD C4 LEA PP2C各一个基因做箱式图
#order1.txt
cd /public/home/tllu/sorghum/seed_weight
chr_3316945	seed_weight.txt	SORBI_3003G036200	4.03E-09	C4单倍型分析情况
chr_9912327	seed_weight.txt	SORBI_3004G105100	9.68E-20	POD单倍型分析结果
chr_39654419	seed_weight.txt	SORBI_3007G109800	1.05E-07	LEA单倍型分析结果
chr_57655596	seed_weight.txt	SORBI_3009G238600	9.20E-05	PP2C单倍型分析结果
####调取这些SNP的基因型和表型
cd /public/home/tllu/sorghum/seed_weight
plink --bfile target_snp --extract order1.txt  --recode --out order1_seed_weight
paste -d'\t' order1_seed_weight.ped seed_weight.txt >order1_seed_weight.txt

#order2.txt
cd /public/home/tllu/sorghum/yield_potential
chr_3315210	yield_potential.txt	SORBI_3003G036200	5.53E-10	C4单倍型分析情况
chr_9912327	yield_potential.txt	SORBI_3004G105100	8.14E-17	POD单倍型分析结果
chr_39654419	yield_potential.txt	SORBI_3007G109800	1.13E-08	LEA单倍型分析结果
chr_57657057	yield_potential.txt	SORBI_3009G238600	0.006326004	PP2C单倍型分析结果
####调取这些SNP的基因型和表型
cd /public/home/tllu/sorghum/yield_potential
plink --bfile target_snp --extract order2.txt  --recode --out order2_yield_potential
paste -d'\t' order2_yield_potential.ped yield_potential.txt >order2_yield_potential.txt

#绘图
setwd("D:/高粱/单倍型分析/单倍型分析结果")
name <- c("seed_weight","yield_potential")
i=name[[1]]
a1 <- read.table(paste0("order1_",i,".txt"))
b1 <- a1[,c(1,7,17)]
b1$type <- rep("chr_3316945",nrow(b1))
b1$cluster <- "C4"
b1$gene <- "SORBI_3003G036200"
b2 <- a1[,c(1,9,17)]
b2$type <- rep("chr_9912327",nrow(b2))
b2$cluster <- "LEA"
b2$gene <- "SORBI_3004G105100"
b3 <- a1[,c(1,11,17)]
b3$type <- rep("chr_39654419",nrow(b3))
b3$cluster <- "POD"
b3$gene <- "SORBI_3007G109800"
b4 <- a1[,c(1,13,17)]
b4$type <- rep("chr_57655596",nrow(b4))
b4$cluster <- "PP2C"
b4$gene <- "SORBI_3009G238600"
colnames(b1) <- colnames(b2) <- colnames(b3) <- colnames(b4) <- c("sample","Haplotype","value","type","cluster","gene")
a1 <- rbind(b1,b2,b3,b4)
a1$pattern <- i

i=name[[2]]
a2 <- read.table(paste0("order2_",i,".txt"))
b1 <- a2[,c(1,7,17)]
b1$type <- rep("chr_3315210",nrow(b1))
b1$cluster <- "C4"
b1$gene <- "SORBI_3003G036200"
b2 <- a2[,c(1,9,17)]
b2$type <- rep("chr_9912327",nrow(b2))
b2$cluster <- "LEA"
b2$gene <- "SORBI_3004G105100"
b3 <- a2[,c(1,11,17)]
b3$type <- rep("chr_39654419",nrow(b3))
b3$cluster <- "POD"
b3$gene <- "SORBI_3007G109800"
b4 <- a2[,c(1,13,17)]
b4$type <- rep("chr_57657057",nrow(b4))
b4$cluster <- "PP2C"
b4$gene <- "SORBI_3009G238600"
colnames(b1) <- colnames(b2) <- colnames(b3) <- colnames(b4) <- c("sample","Haplotype","value","type","cluster","gene")
a2 <- rbind(b1,b2,b3,b4)
a2$pattern <- i
#合并
at <- rbind(a1,a2)
###箱式图
head(at)
at <- at[which(at$Haplotype != "0"),]
at <- at[which(at$Haplotype != "*"),]
#a1$Haplotype <- gsub("A","H1",a1$Haplotype)
#a1$Haplotype <- gsub("G","H2",a1$Haplotype)
plot_list=list()
m <- 1
name <- c("seed_weight","yield_potential")
for (i in name) {
  y <- at[which(at$pattern == i),]
  list <- unique(y$type)
  for (j in list) {
    x <- y[which(y$type == j),]
    if (length(unique(x$Haplotype)) == 2) {
      my_comparisons <- list(c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[2]]))
    } else if (length(unique(x$Haplotype)) == 3) {
      my_comparisons <- list(c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[2]]),c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[3]]),
                             c(unique(x$Haplotype)[[2]],unique(x$Haplotype)[[3]]))
    } else if (length(unique(x$Haplotype)) == 4) {
      my_comparisons <- list(c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[2]]),c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[3]]),
                             c(unique(x$Haplotype)[[1]],unique(x$Haplotype)[[4]]),c(unique(x$Haplotype)[[2]],unique(x$Haplotype)[[3]]),
                             c(unique(x$Haplotype)[[2]],unique(x$Haplotype)[[4]]),c(unique(x$Haplotype)[[3]],unique(x$Haplotype)[[4]]))
    }
    p1 <- ggboxplot(x, x = "Haplotype", y = "value",
                    color = "Haplotype", palette = "jco",
                    add = c("jitter"),
                    bxp.errorbar=TRUE,width = 0.5,
                    bxp.errorbar.width = 0.1,
                    add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
      #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
      #scale_fill_manual(values=brewer.pal(12,"Paired")) +
      scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
      ggtitle(paste(x$gene,x$cluster,sep="-")) + 
      theme(legend.position = "none") +
      theme(plot.title = element_text(hjust = 0.5,size=18),
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
      xlab(x$type) + ylab(x$pattern) +
      stat_compare_means(comparisons = my_comparisons,
                         #aes(group = Haplotype),
                         #label = "p.signif",
                         method = "t.test")
    plot_list[[m]] <- p1
    m <- m+1
  }
}
###组图
library(gridExtra)
grid.arrange(plot_list[[1]],plot_list[[2]], 
             plot_list[[3]],plot_list[[4]], 
             plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],
             nrow=2,ncol=4)     %>%  ggsave("单倍型分析结果-seedweight_yieldPotential.pdf",.,width=210,height=130, units="mm")

#####调取SNP上下游100bp序列做motif discovery分析
chr3_3316945
chr4_9912327
chr7_39654419
chr9_57655596
chr3_3315210
chr4_9912327
chr7_39654419
chr9_57657057
3	3316895	3316995
3	3315160	3315260
4	9912277	9912377
4	9912277	9912377
7	39654369	39654469
7	39654369	39654469
9	57655546	57655646
9	57657007	57657107
cd /public/home/tllu/sorghum
bedtools getfasta -fi /public/home/chaohe/db/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -bed sor.bed -fo sor.fa
#####结果 CDF3
motif_ID	motif_ALT_ID	seq_ID	site_Start	site_End	site_Strand	site_Score	site_Sequence
#1-GCACTT	STREME-1	4	9912321	9912326	+	9.78	GCACTT
1-GCACTT	STREME-1	3	3315207	3315212	+	5.25	GCACAT
#1-GCACTT	STREME-1	3	3316939	3316944	-	0.18	GGACTG
1-GCACTT	STREME-1	7	39654415	39654420	+	0.02	GAATTT
2-CAAGGT	STREME-2	3	3316945	3316950	+	11.36	CAAGGT
3-GCCACCGGCGAT	STREME-3	9	57655590	57655601	-	18.89	GCCACCGGCGAT
4-CATTAATTTTTA	STREME-4	9	57657051	57657062	+	16.83	CATTAATTTTTA
5-GATGTA	STREME-5	4	9912327	9912332	+	11.32	GATGTA







