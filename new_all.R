library( "DESeq2" )
library(ggplot2)
library(readr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(sjmisc)
library(purrr)
library(carData)
library(car)
library(fdrtool)
# 需要计算的数据
# countData <- read.csv("/data1/DYY/data/288sample_ensg_int_tpm_by_tissue/all_gene.csv",header = T, sep = '\t')
# countData <- read.csv("/data1/DYY/bambu/counts_gene.csv",header = T, sep = '\t',row.names = 1)
# countData <- read.csv("/data1/DYY/bambu/bambu_data/gene_counts.csv",header = T, sep = '\t',row.names = 1)
# countData <- read.csv("/data1/DYY/bambu/288count_gene.csv",header = T, sep = '\t',row.names = 1)
countData <- read.csv("/data1/DYY/bambu/288count_transcript.csv",header = T, sep = '\t',row.names = 1)
head(data, 4)
data[1:5,1:5]
nrow(data)
# rownames(countData) <- countData[,1]
data <- as.matrix(countData)
ncol(data)
# data <- data[-((nrow(data)-1):nrow(data)),c(1:ncol(data))]
data <- apply(data, c(1, 2), as.integer)

# 
# countdata <- countData[,c(1:ncol(countData))]
# countData <- countData[-c(-2,-1),]
# tail(countData)
# nrow(countData)
# 数据的分组列表
# countdata <- countData[3:nrow(countData),]
# countdata[2:ncol(countdata)] <- as.data.frame(lapply(countdata[,2:ncol(countdata)],as.numeric))
# metadata <- t(countData[0:2,])
# colnames(metadata) <- metadata[1,]
# metaData <- read.csv("/data1/DYY/data/288sample_ensg_int_tpm_by_tissue/group.csv",header = T, sep = '\t',row.names = 1)
# metaData <- read.csv("/data1/DYY/bambu/350_samples_group.csv",header = T, sep = '\t')
metaData <- read.csv("/data1/DYY/bambu/288_group.csv",header = T, sep = '\t')
head(metaData)
# # metaData <- read.csv("/data1/DYY/bambu/288_group.csv",header = T, sep = '\t')
# # group <- metaData[,c('ID','tissue','sample')]
# rownames(metaData) <- metaData[,1]
# head(countData)
# metaData <- metaData[colnames(data), ]
# # rownames(metaData) <- 1:nrow(metaData)
# metaData <- metaData[,-1]
metaData <- data.frame(metaData)

# metadata <- read.csv("/data1/DYY/data/ensg_int_tpm_by_tissue/group.csv",header = T, sep = '\t',row.names = 1)
# rownames(metadata) <- metadata$index
# metadata <- metadata[-1,]
# colnames(metaData) <- c("ID","Type")
# 
# countData <- countData[,-1]
# rownames(metaData) <- c(1:nrow(metaData))

# countData <- countData[3:nrow(countData),]
# metaData <- t(metaData)
# metaData <- metaData[-1,]
# countData<- as.data.frame(lapply(countData,as.integer))
# row.names(countData) = countData$ensg
# countData = countData[,2:ncol(countData)]
# countData <- t(countData)
###############
####DESeq2#####
###############
# countData <- subset(countData,select=-c(group))
# countData[(1:20),(1:20)]
path <- '/data1/DYY/bambu/288_count_matrix_deseq2/'
metaData$sample <- substr(metaData$Sample.ID,1,3)
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=metaData, 
                              design=~ sample + Type) #dds建库
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table
summary(res,alpha=0.05)
r <- results(dds, tidy=TRUE)
filepath <- paste(path,"all_trans_result.csv",sep = '') #输出文件路径
write.table(r,file = filepath,row.names=FALSE,col.names=TRUE,sep="\t") #输出
dds <- estimateSizeFactors(dds)
filepath <- paste(path,"all_trans_norm.counts.fn",sep="")
write.table(counts(dds,normalized=TRUE),file = filepath,sep = '\t',col.names = T)

##################
###正态分布检验###
##################
# filelist <- paste(path,filename,sep = '/')
# for(i in 1:length(filelist)){
#   csvname <- filelist[i] #取出数据文件
#   countdata <- read.csv(csvname,header=T,sep = '\t',row.names = 1) #读文件时要把第一列读为索引
#   for (c in seq(1,nrow(countdata),5000)) {
#     a <- as.vector(t(countdata[c,]))
#     qqnorm(a,main="Normal Q-Q Plot",col = "blue")
#     qqline(a, col = "red")
#   }
# }
################
###leveneTest###
################
P_value <- c()
file_path <- paste(path,'all_trans_norm.counts.fn',sep = '')
countdata <- read.csv(file_path,header=T,sep = '\t',row.names = 1)
# metadata <- read.csv(datname,header=T,sep = '\t',row.names = 1)
group <- as.factor(as.vector(c(metaData[2]))$Type) #取出分组信息为factor
# group <- as.factor(as.vector(c(metaData[2]))$tissue)
std <- apply(countdata,1,sd) #求标准差
countdata <- cbind(countdata,std) #添加一列标准差
# countdata <- filter(countdata,std>1) #筛选标准差>1的列
countdata <- subset(countdata,select=-c(std)) #去掉std列
for (c in 1:nrow(countdata)) {
  data <- as.vector(t(countdata[c,])) #将data转为向量
  Pr <- format(leveneTest(data, group)[[3]][1],scientific=T) #方差齐性检验
  P_value <- append(P_value,Pr) #将结果加入新的vector
  #     # p_value <- matrix(P_value) #list2matrix
}
countdata <- cbind(countdata,P_value) #给data添加新的一列
filepath <- "/data1/DYY/bambu/288_LEVENE/all_trans.csv"
print(filepath)
write.table(countdata,file = filepath,sep = '\t',col.names = T)

##################
###BH校正并筛选###
##################
# i = 1
data <- read.csv(filepath,sep = '\t',header = T)
dat <- data[order(data$P_value),] #按p值排序
p_adjust <- p.adjust(
  dat[,ncol(dat)],
  method = "BH" #使用BH校正
)
dat$p_adjust <- p_adjust
filter_data <- filter(.data=dat,p_adjust<0.05)
files <- "/data1/DYY/bambu/288_FDR/all_trans.csv"
write.table(filter_data,file=files,sep = '\t',col.names = T)

