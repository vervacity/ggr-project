#!/usr/bin/env Rscript

library(gplots)
library(RColorBrewer)
library(reshape2)

set.seed(42)

# functions
smoothen <- function(data) {
    # given data, smoothen the data
    data_tmp <- data[2:(length(data)-1)]
    data_tmp <- data_tmp + data[1:(length(data)-2)] + data[3:length(data)]
    data_tmp <- data_tmp / 3.
    data_new <- c(data[1], data_tmp, data[length(data)])

    return(data_new)
}


# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# load data
data <- read.table(
    gzfile(data_file), sep="\t", header=TRUE, row.names=1)
#data[,159:162] <- 0
data[,161] <- NA
data <- data[,130:190]
#data <- data[,110:210]
print(dim(data))

# clean up names
pwm_names <- gsub("HCLUST-\\d+_", "", rownames(data))
pwm_names <- gsub(".UNK.0.A", "", pwm_names)
rownames(data) <- pwm_names

# smoothen
#data <- apply(data, 2, smoothen)

# clean up data
data_melt <- melt(data)
min_val <- quantile(data_melt$value, 0.10, na.rm=TRUE)
data[data < min_val] <- min_val

# clip the top?
max_val <- quantile(data_melt$value, 0.99, na.rm=TRUE) # 0.993
data[data > max_val] <- max_val

# ordering
# order on a smoothed profile?
data_smooth <- apply(data, 2, smoothen)

#data_smooth <- data
#mid_idx <- as.integer(dim(data_smooth)[2] / 2.)
#data_smooth <- data_smooth[,mid_idx:ncol(data_smooth)]

hc <- hclust(dist(data_smooth), method="ward.D2")
ordering <- rev(hc$order)
data <- data[ordering,]


plot_heatmap <- function(data, plot_title) {

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.2,0.7,0.07)
    mylhei = c(0.09,2,0.20)

    # color
    my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(49)
    my_palette <- colorRampPalette(brewer.pal(9, "BuGn"))(49)
    
    #my_breaks <- c(
    #    seq(0.005, 0.035, length=30),
    #    seq(0.036, 0.040, length=20))
    
    # labCol
    labCol <- rep("", ncol(data))
    for (i in 1:ncol(data)) {
        if (i%%10-1 == 0) {
            labCol[i] <- i - 1 - as.integer(ncol(data) / 2)
        }
    }
    
    # heatmap
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        colsep=c(0,(ncol(data)+1)),
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001,0.001),

        cexRow=0.5,
        offsetRow=-9.5,

        labCol=labCol,
        cexCol=0.5,
        offsetCol=-0.25,
        srtCol=0,
        adjCol=c(1,0),
        
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            #mar=c(2,6,1.4,6),
            mar=c(1,1,0.5,1),
            #mar=c(0.5,1,0.5,1),
            mgp=c(0,-0.1,0),
            tcl=-0.1,
            lend=2,
            cex.axis=0.6,
            bty="n"),
        key.xtickfun=function() {
            breaks <- pretty(parent.frame()$breaks)
            breaks <- breaks[c(1,length(breaks))]
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)},

        margins=c(1,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        na.color="grey",
        col=my_palette)
        #breaks=my_breaks)

    title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
}


pdf(plot_file, height=3.5, width=1.3, useDingbats=FALSE, family="ArialMT")
plot_heatmap(data, "Summary")
dev.off()




