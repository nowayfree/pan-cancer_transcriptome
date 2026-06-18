library(transPlotR)
library("rtracklayer")


args <- commandArgs(trailingOnly = TRUE)

gtf_data = import(args[1])
gtf_data = as.data.frame(gtf_data)
gtf_data$gene_name <- args[2]
gtf_data$transcript_name <- 'test'
trancriptVis(gtfFile = gtf_data,
             gene = args[2],
             arrowType = 'open',myTranscript=c(args[3:length(args)]))
ggsave(paste0('/data1/DYY/bambu/plot/gene_structure/',args[2],'.pdf',sep=''),height = 8)
