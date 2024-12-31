library('clusterProfiler') 
library('org.Hs.eg.db')
library(stringr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library("ggsci")
library(svglite)
data <- read.table('/home/bwu4/go_kegg/Supplementary_Table_2.csv',sep="\t",header = TRUE)
head(data)

df <- data[,c("tissue",'Description','down_gene_ratio','down_p.adjust')]

#df$Description <- factor(df$Description,levels = rev(df$Description))
mytheme <- theme(axis.title=element_text(size=14,colour = 'black'), #坐标轴标题
                 axis.text=element_text(size=14,colour = 'black'), #坐标轴标签
                 axis.line = element_line(linewidth = 0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank() #关闭图例边框
)
df <- df %>%
  mutate(down_gene_ratio = sapply(down_gene_ratio, function(x) {
    parts <- unlist(strsplit(x, "/"))
    numerator <- as.numeric(parts[1])
    denominator <- as.numeric(parts[2])
    return(numerator / denominator)
  }))
plot_data <- df[order(df$down_p.adjust,decreasing = TRUE),]


df <- df[df$Description %in% c("anion transmembrane transport", "glucose metabolic process", "sodium ion homeostasis"), ]

p <- ggplot(head(plot_data,60),aes(x=tissue,y=Description,colour=-1*log10(down_p.adjust),size=-log10(down_gene_ratio)))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO Pathway Terms")+
  xlab("")+
  labs(color=expression(-log[10](p.adjust)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
p
ggsave('/data1/DYY/plot/supplementary_figure_3.pdf',width = 10, height = 10,plot = p)
