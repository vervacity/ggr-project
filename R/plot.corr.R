#!/usr/bin/env Rscript

# description: plot simple correlation matrix
library(gplots)
library(RColorBrewer)

# args
args <- commandArgs(trailingOnly=TRUE)
mat_file <- args[1]
plot_file <- args[2]

# load data
data <- read.table(gzfile(mat_file), sep="\t", header=TRUE, row.names=1)
print(head(data))


# grid
mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
mylwid = c(0.05,4,0.05)
mylhei = c(0.05,4,0.5)

my_palette <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(49))
    
# plot
pdf(plot_file, height=6, width=5)
heatmap.2(
    as.matrix(data),
    Rowv=FALSE,
    Colv=FALSE,
    dendrogram="none",
    trace="none",
    density.info="none",
    #labCol=c(
    #    rep("", floor(ncol(plot_data)/2)),
    #    label[i],
    #    rep("", ceiling(ncol(plot_data)/2)-1)),
    #labRow=rep("", nrow(plot_data)),
    keysize=0.1,
    key.title=NA,
    key.xlab=NA,
    key.par=list(pin=c(4,0.1),
        mar=c(2.1,0,2.1,0),
        mgp=c(3,1,0),
        cex.axis=1.0,
        font.axis=2),
    srtCol=45,
    cexCol=1.25,
    margins=c(1,0),
    lmat=mylmat,
    lwid=mylwid,
    lhei=mylhei,
    col=my_palette,
    #breaks=my_breaks,
    #rowsep=rowsep,
    sepcolor="black")
dev.off()


quit()

# plot
ggplot(pca_data, aes(x=x, y=y, colour=group)) +
    geom_point(size=3) + xlim(-250,250) + ylim(-250,250) + 
    theme_bw()
ggsave(plot_file)
