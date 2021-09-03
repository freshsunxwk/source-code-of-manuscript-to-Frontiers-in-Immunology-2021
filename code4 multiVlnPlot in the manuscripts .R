##############xwk###############
#multiVlnPlot in the manuscripts 
library(ggplot2)
xxxx<- function(obj,feature, pt.size = 0.1,  plot.margin = unit(c(-0.75, 0.5, -0.75, 0.5), "cm"), ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +xlab("") + ylab(feature) + ggtitle("      ") +theme(legend.position = "none",plot.margin = plot.margin )
  return(p)}
multiVlnPlot<- function(obj, features,pt.size = 0.1, plot.margin = unit(c(-0.75, 0.5, -0.75, 0.5), "cm"),...) {
  plot_list<- purrr::map(features, function(x) xxxx(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)
  return(p)}
multiVlnPlot(immune.combined, markers.to.plot, pt.size=0.1, cols=col)
#################volcano with  label
library(TCGAbiolinks)
b.interferon.response$Gene=rownames(b.interferon.response)
TCGAVisualize_volcano(b.interferon.response$avg_log2FC, b.interferon.response$p_val_adj,title = "",
                      filename = volFig, 
                      xlab = "logFC",
                      names = b.interferon.response$Gene, show.names = "highlighted",
                      x.cut = 0.25, y.cut = 0.25, 
                      highlight = b.interferon.response$Gene[which(abs(b.interferon.response$avg_log2FC) >=0.25)],
                      highlight.color = "orange") 
theme(legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.direction = "horizontal",
      legend.position = c(0.2,0.93),
      panel.background = element_rect(fill = "transparent",colour = "black"),
      plot.background = element_rect(fill = "transparent",colour = "black"))
dev.off()