library(ggplot2)
library(reshape2)

dtt<-read.table("/data1/DYY/bambu/3f_data_generatio.csv",sep ="\t",header= T)
dtt<-melt(dtt)

p_value <- read.table("/data1/DYY/bambu/3f_data_p.csv",sep ="\t",header= T)
p_value <- melt(p_value)
p_value$logp <- -log10(p_value$value)
p_value
dtt$'-logp' <- p_value$logp
# dtn$'-logp' <- p_value$logp
p1 <- ggplot() +
  # 绘制小于 1.3 的点
  geom_point(
    data = subset(dtt, `-logp` < 1.3),
    aes(x = variable, y = Description, size = value),
    color = "lightgray"  # 直接固定灰色
  ) +
  # 绘制大于等于 1.3 的点
  geom_point(
    data = subset(dtt, `-logp` >= 1.3),
    aes(x = variable, y = Description, size = value, color = `-logp`)
  ) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('darkgreen', 'orange'),
    limits = c(1.3, max(dtt$`-logp`))
  ) +
  theme_bw() +
  scale_size_continuous(range = c(1, 10)) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
  ) +
  xlab(NULL) +
  ylab(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  geom_vline(xintercept = c(20.5, 39.5, 58.5, 78.5, 95.5, 114.5, 134.5), size = 0.5) +
  ggtitle('normal')

p1
ggsave("/data1/DYY/bambu/plot/figure3/3f.pdf", p1, device = "pdf", width = 10, height = 5, dpi = 300)
