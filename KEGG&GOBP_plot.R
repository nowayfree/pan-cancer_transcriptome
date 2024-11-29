library('clusterProfiler') 
library('org.Hs.eg.db')
library(stringr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library("ggsci")
library(svglite)
library("enrichplot")
tis = 'tumor_spec_gene'
tis = 'ir_sig_up'
tis = 'tumor_>_normal'
tis = 'tumor_<_normal'
tis = 'liv_kid_share'
tis = 'down_reg_shared'
tis = 'magoh'

gene_name_list <- read.table('/data1/DYY/data/ir_sig_deg_barplot/up_gene.txt')
gene_name_list <- read.table('/data1/DYY/data/all_tumor_spec_gene.txt')
gene_name_list <- read.table('/data1/DYY/bambu/ir_sig_up_gene.csv')
gene_name_list <- read.table('/data1/DYY/bambu/tumor_specific_gene_in_2_tissues.csv',header = F)
gene_name_list <- read.table('/data1/DYY/bambu/all_ts_gene.csv')
head(gene_name_list)
ENTREZID<- bitr(gene_name_list$V1, fromType = "ENSEMBL", toType=c("SYMBOL","ENTREZID"),OrgDb = org.Hs.eg.db)
ego_BP <- enrichGO(gene = ENTREZID$ENTREZID, #universe = names(geneList),
                   OrgDb = org.Hs.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.5,qvalueCutoff = 0.5,readable = TRUE)
write.table(ego_BP,paste('~/bambu_data_kegg_go/',tis,'_go_bp.csv',sep = ''), sep = '\t', quote = FALSE, row.names = FALSE)
# test <- read.table(paste('~/tmp_data_kegg_go/',tis,'_go_bp.csv',sep = ''), sep = '\t')
plotpath <- '/data1/DYY/bambu/plot/enrichment_plot/'
kegg <- enrichKEGG(
  gene = ENTREZID$ENTREZID,  #基因列表文件中的基因名称
  keyType = 'kegg',  #KEGG 富集
  organism = 'human',  
  pAdjustMethod = 'fdr',  #指定 p 值校正方法
  pvalueCutoff = 0.5,  #指定 p 值阈值（可指定 1 以输出全部）
  qvalueCutoff = 0.5)  #指定 q 值阈值（可指定 1 以输出全部）
if (length(kegg)!=0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

  write.table(kegg,paste('~/bambu_data_kegg_go/',tis,'_kegg.csv',sep = ''), sep = '\t', quote = FALSE, row.names = FALSE)
  KEGG_dataset <- read.table(file = paste('~/bambu_data_kegg_go/',tis,'_kegg.csv',sep = ''),header = TRUE, sep = "\t")
  ##对kegg数据进行按照p值排序
  KEGG_dataset <- arrange(KEGG_dataset,KEGG_dataset[,5])
  KEGG_dataset$enrichment_fold <- apply(KEGG_dataset,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    enrichment_fold
  })
  #Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
  #将Pathway列转化为因子型
  KEGG_dataset$Description <- factor(KEGG_dataset$Description,levels = rev(KEGG_dataset$Description))
  mytheme <- theme(axis.title=element_text(size=14,colour = 'black'), #坐标轴标题
                   axis.text=element_text(size=14,colour = 'black'), #坐标轴标签
                   axis.line = element_line(linewidth = 0.5, colour = 'black'), #轴线
                   panel.background = element_rect(color='black'), #绘图区边框
                   legend.key = element_blank() #关闭图例边框
  )
  p <- ggplot(head(KEGG_dataset[order(KEGG_dataset$enrichment_fold, decreasing = TRUE), ],10),aes(x=enrichment_fold,y=Description,colour=-1*log10(p.adjust),size=Count))+
    geom_point()+
    scale_size(range=c(2, 10))+
    # scale_colour_gradient(low = "#6d8fbf", high = "#e4e4ea")+
    # scale_colour_gradient(low = "#843F0B", high = "#EE822F")+
    scale_fill_gradient(low = "#4399ab", high = "#e5eca2")+
    theme_bw()+
    ylab("KEGG Pathway Terms")+
    xlab("Enrichment fold")+
    labs(color=expression(-log[10](PValue)))+
    theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
    theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
    theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
  kegg_plot <- p+mytheme
  kegg_plot
  kegg_plottitle <- paste(tis,"_kegg_dotq_0.1",sep = '_')
  ggsave(paste(plotpath,kegg_plottitle,'.pdf',sep = ''),width = 10, height = 10,plot = go_plot)
}

###############
#####go_bp#####
###############
go_bp = read.table(paste('~/bambu_data_kegg_go/',tis,'_go_bp.csv',sep = ''), sep = '\t',header = T)
go_bp$enrichment_fold <- apply(go_bp,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"]))
  enrichment_fold=round(GeneRatio/BgRatio,2)
  enrichment_fold
})
##对go数据进行按照p值排序
go_bp <- arrange(go_bp,go_bp[,10])
#Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
#将Pathway列转化为因子型
go_bp$Description <- factor(go_bp$Description,levels = rev(go_bp$Description))
p <- ggplot(go_bp[1:10,],aes(y=as.numeric(enrichment_fold),x=Description,fill=-1*log10(p.adjust))) + 
  geom_bar(stat="identity",position = "dodge") +
  # facet_grid(Category~.,scales = "free",space = "free") + 
  coord_flip() + 
  theme_bw() + 
  scale_fill_gradient(low = "#4399ab", high = "#e5eca2")+
  # scale_fill_gradient(low = "#843F0B", high = "#EE822F")+
  # scale_color_gsea()+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size = 14),
        legend.position="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        # axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        # axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
go_plot <- p  #+ scale_x_discrete(labels=function(x) str_wrap(x, width=10))
go_plot
#对y轴标签分行
# p + scale_x_discrete(labels=function(x) str_wrap(x, width=10))
# unlink(file.path('~/tmp_data_kegg_go/', '*'), recursive = TRUE)
# go_plottitle <- paste(tis,"_go_bp_plot",sep = '_')
# mf_title <- paste(tis,"EnrichmentGO_MF_dot",sep = '_')
# cc_title <- paste(tis,"EnrichmentGO_CC_dot",sep = '_')
bp_title <- paste(tis,"EnrichmentGO_BPplot",sep = '_')
# mf <- dotplot(ego_MF,title=mf_title)
# cc <- dotplot(ego_CC,title=cc_title)
# bp <- dotplot(ego_BP,title=bp_title)
# ggsave(paste(plotpath,mf_title,'.pdf',sep = ''), plot = mf)
# ggsave(paste(plotpath,cc_title,'.pdf',sep = ''), plot = cc)
ggsave(paste(plotpath,bp_title,'.pdf',sep = ''),width = 10, height = 10, plot = go_plot)

