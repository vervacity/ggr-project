



# script to plot datasets


library(gplots)
library(RColorBrewer)

library(ggplot2)
library(reshape)


my_palette <- colorRampPalette(brewer.pal(9, 'Blues'))(49)


args <- commandArgs(trailingOnly=TRUE)
in_table_file <- args[1]
out_pdf <- args[2]

# load table
data <- read.table(in_table_file, sep='\t', header=TRUE, row.names=1)
data[is.na(data)] <- 0

# plot table
pdf(out_pdf)
heatmap.2(as.matrix(data),
          Rowv=FALSE, Colv=FALSE, dendrogram='none',
          trace='none',
          sepwidth=c(0.02, 0.02),
          col=my_palette)
dev.off()



# ggplot style
data$assay <- factor(rownames(data), levels=rownames(data)[rev(order(rownames(data)))])
#data$assay <- with(data, reorder(assay, assay))
#data_rev <- data[rev(1:nrow(data)),]
data.m <- melt(data)
data.m$timepoint <- sub('day.', '', data.m$variable)

p <- ggplot(data.m, aes(timepoint, assay)) + geom_tile(aes(fill=value), colour="black") + scale_fill_gradient(low="white", high="steelblue")
base_size <- 9
p + theme_bw() + theme_grey(base_size=base_size) +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))  + theme(legend.position="none") +
        labs(x="", y="") + theme(axis.ticks=element_blank())
        
#+
#    labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) +
#        scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none", axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size *0.8, angle = 330, hjust = 0, colour = "grey50"))

ggsave("testing.pdf", width=6, height=1.5)

