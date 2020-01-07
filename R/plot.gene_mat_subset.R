#!/usr/bin/env Rscript

library(gplots)
library(RColorBrewer)

library(ggplot2)
library(reshape2)
library(scales)

# args
args <- commandArgs(trailingOnly=TRUE)
mat_file <- args[1]
plot_file <- args[2]
#plot_file <- paste(
#    unlist(strsplit(mat_file, ".mat.txt", fixed=TRUE))[1],
#    ".pdf", sep="")

# read data, normalize to d0
data <- read.table(mat_file, sep="\t", header=TRUE)
data[,3:ncol(data)] <- data[,3:ncol(data)] - data[,3]
max_val <- max(abs(data[,3:ncol(data)]))
max_val <- 6

# adjust
data$index <- NULL
rownames(data) <- data$hgnc_symbol
data$hgnc_symbol <- NULL
data[data > max_val] <- max_val
data[data < -max_val] <- -max_val


# heatmap2 grid
mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
#mylwid = c(0.02,0.5,0.02)
mylwid = c(0.02,0.5,0.5)
mylhei = c(0.1,2,0.27)

my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

pdf(plot_file, height=3, width=1.25, family="ArialMT")
heatmap.2(
    as.matrix(data),
    Rowv=TRUE,
    Colv=FALSE,
    dendrogram="none",
    trace='none',
    density.info="none",
    
    srtCol=60,
    cexCol=0.5,
    offsetCol=-0.5,
    #labRow="",
    rowsep=0:(nrow(data)+1),
    colsep=c(0, ncol(data)+1),
    sepcolor="black",
    sepwidth=c(0.0001, 0.0001),
    
    keysize=0.1,
    key.title=NA,
    key.xlab=NA,
    key.par=list(
        mar=c(0.9,1,1,1),
        mgp=c(0,-0.1,0),
        tcl=-0,1,
        lend=2,
        cex.axis=0.6,
        bty="n"),
    key.xtickfun=function() {
        breaks <- pretty(parent.frame()$breaks)
        breaks <- breaks[c(1,length(breaks))]
        list(at = parent.frame()$scale01(breaks),
             labels = breaks)},

    margins=c(0,0),
    lmat=mylmat,
    lwid=mylwid,
    lhei=mylhei,
    
    col=my_palette)
dev.off()



quit()




print(data)
data_melt <- melt(data, id.vars=c("index", "hgnc_symbol"))
print(head(data_melt))

my.scale_fill_distiller <- function(
    ...,
    type = "seq",
    palette = 1,
    direction = -1,
    values = NULL,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") {
    # warn about using a qualitative brewer palette to generate the gradient
    type <- match.arg(type, c("seq", "div", "qual"))
    if (type == "qual") {
        warning("Using a discrete colour palette in a continuous scale.\n  Consider using type = \"seq\" or type = \"div\" instead", call. = FALSE)
    }
    continuous_scale(
        aesthetics, "distiller",
        gradient_n_pal(brewer_pal(type, palette, direction)(7), values, space), na.value = na.value, guide = guide, ...)
    # NB: 6 colours per palette gives nice gradients; more results in more saturated colours which do not look as good
}

# plot
ggplot(data_melt, aes(x=variable, y=hgnc_symbol, fill=value)) +
    geom_tile() +
    theme_bw() +
    #scale_fill_distiller(palette="RdBu", limits=c(-max_val,max_val))
    my.scale_fill_distiller(palette="RdBu", limits=c(-max_val,max_val))
ggsave(plot_file)
