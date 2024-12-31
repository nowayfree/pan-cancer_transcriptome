library( "DESeq2" )
library(ggplot2)
library(readr)
library(tidyverse)
# BiocManager::install("tidyverse")
library(tidyr)
library(dplyr)
library(ggplot2)
library(sjmisc)
library(purrr)
library(carData)
library(car)
library(fdrtool)
###############
####DESeq2#####
###############
path <- "/data1/DYY/bambu/" #文件路径
fileName <- dir(paste(path,"288_count_gene_by_tissue/",sep = "")) #取出全部csv文件名
filename = unique(substr(fileName,1,3)) #取出全部癌症
# files <- paste(path,"count_gene_by_tissue/",filename,'.csv',sep = '') #数据路径
files <- paste(path,"288_count_gene_by_tissue/",filename,'.csv',sep='')
files <- as.matrix(files)
datas <- paste(path,"group_by_tissue/",filename,'.csv',sep = '') #分组信息路径
datas <- as.matrix(datas)
# a <- array(1,dim = c(1,1))
# a[1,1]<-b
# b<-array(c(1,2,3,4), dim = c(2,2))
# b
# i = 1
a <- NULL
for(i in 1:length(filename)){
  # i=1
  csvname <- files[i] #取出数据文件
  datname <- datas[i] #取出分组文件
  countdata <- read.csv(csvname,header=T,sep = '\t',stringsAsFactor = FALSE)
  countdata <- as.data.frame(countdata)
  metadata <- read.csv(datname,header=T,sep='\t')
  metadata <- as.data.frame(metadata)
  # metadata <- t(countdata[(1:2),])
  # colnames(metadata) = metadata[1,]
  # metadata <- as.data.frame(metadata[-1,])
  # countdata <- countdata[(3:nrow(countdata)),]
  # rownames(countdata) <- 1:nrow(countdata)
  # countdata <- data[3:nrow(data),]
  # countdata[2:ncol(countdata)] <- as.data.frame(lapply(countdata[,2:ncol(countdata)],as.numeric))
  # metadata <- t(data[0:2,])
  # colnames(metadata) <- metadata[1,]
  # metadata <- metadata[-1,]
  # colnames(metadata) <- c("Type","ID",'tissue')
  metadata$Type <- factor(metadata$Type)
  # metadata$tissue <- factor(metadata$tissue)
  # countdata <- countdata[, rownames(metadata)]
  dds <- DESeqDataSetFromMatrix(countData=countdata, 
                                colData=metadata, 
                                design=~Type, tidy = TRUE) #dds建库
  dds <- DESeq(dds)
  res <- results(dds)
  # head(results(dds, tidy=TRUE)) #let's look at the results table
  summary(res,alpha=0.05)
  r <- results(dds, tidy=TRUE)
  filepath <- paste(path,"288_count_matrix_deseq2/",filename[i],"_result",".csv",sep = '') #输出文件路径
  write.table(r,file = filepath,row.names=FALSE,col.names=TRUE,sep="\t") #输出
  dds <- estimateSizeFactors(dds)
  filepath <- paste(path,"288_count_matrix_deseq2/",filename[i],"_norm.counts.fn",sep="")
  write.table(counts(dds,normalized=TRUE),file = filepath,sep = '\t',col.names = T)
}
################
###leveneTest###
################
file_des <- paste(path,"288_count_matrix_deseq2/",filename,"_norm.counts.fn",sep = '') #数据路径
file_des <- as.matrix(file_des)
for(i in 1:length(filename)){
  # i = 1
  P_value <- c()
  descsvname <- file_des[i] #取出数据文件
  datname <- datas[i] #取出分组文件
  # countdata <- head(read.csv(csvname,header=T,sep = '\t'))
  desdata <- read.csv(descsvname,header=T,sep = '\t',row.names = 1)
  # countdata <- countdata[1:10,1:10]
  # csvname <- files[i] #取出数据文件
  # datname <- datas[i] #取出分组文件
  # data <- read.csv(csvname,header=T,sep = '\t',stringsAsFactor = FALSE)
  # metadata <- t(countdata[(1:2),])
  # colnames(metadata) = metadata[1,]
  # metadata <- as.data.frame(metadata[-1,])
  # countdata <- countdata[(3:nrow(countdata)),]
  # rownames(countdata) <- 1:nrow(countdata)
  # countdata <- data[3:nrow(data),]
  # countdata[2:ncol(countdata)] <- as.data.frame(lapply(countdata[,2:ncol(countdata)],as.numeric))
  # metadata <- t(data[0:2,])
  # colnames(metadata) <- metadata[1,]
  # metadata <- metadata[-1,]
  # metadata <- read.csv('',header=T,sep = '\t')
  metadata <- read.csv(datname,header=T,sep='\t')
  metadata <- as.data.frame(metadata)
  group <- as.factor(as.vector(c(metadata[,2]))) #取出分组信息为factor
  std <- apply(desdata[,1:ncol(desdata)],1,sd) #求标准差
  desdata <- cbind(desdata,std) #添加一列标准差
  desdata <- filter(.data=desdata,std>1) #筛选标准差>1的列
  desdata <- subset(desdata,select=-c(std)) #去掉std列
  # metadata[metadata$DATA_ID%in%colnames(countdata)]
  for (c in 1:nrow(desdata)) {
    data <- as.vector(t(desdata[c,])) #将data转为向量
    #   std <- sd(data) #标准差
    #   if (is.na(std) | std < 1) {
    #     countdata <- countdata[-c(c),]
    #   }else {
    Pr <- format(leveneTest(data, group)[[3]][1],scientific=T) #方差齐性检验
    P_value <- append(P_value,Pr) #将结果加入新的vector
    #     # p_value <- matrix(P_value) #list2matrix
  }
  desdata <- cbind(desdata,P_value) #给data添加新的一列
  filepath <- paste(path,"288_LEVENE","/",filename[i],".csv",sep = "")
  # print(filepath)
  write.table(desdata,file = filepath,sep = '\t')
}

##################
###BH校正并筛选###
##################
for(i in 1:length(filename)){
  # i=1
  csvname <- paste(path,"288_LEVENE/",filename[i],".csv",sep = "")
  data <- read.csv(csvname,sep = '\t',header = T)
  dat <- data[order(data$P_value),] #按p值排序
  p_adjust <- p.adjust(
    dat[,ncol(dat)],
    method = "BH" #使用BH检验
  )
  dat$p_adjust <- p_adjust
  filter_data <- filter(.data=dat,p_adjust<0.05)
  files <- paste(path,"288_FDR/",filename[i],".csv",sep="")
  write.table(filter_data,file=files,sep = '\t',col.names = T) 
}

