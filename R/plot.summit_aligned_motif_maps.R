#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(gplots)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

library(grid)
library(gridGraphics)
library(gridExtra)

load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

color_granularity <- 50
my_palette <- get_ggr_assay_palette("Blues", color_granularity)


# heatmap fn
plot_profile_heatmap <- function(plot_data, i, plot_key, my_breaks) {

    #label <- c("d0", "d3", "d6")

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.02,0.5,0.02)
    mylhei = c(0.1,2,0.27)

    # heatmap
    heatmap.2(
        as.matrix(plot_data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        srtCol=60,
        cexCol=0.5,
        offsetCol=-0.5,
        labCol=c(
            rep("", floor(ncol(plot_data)/2)),
            #label[i],
            rep("", ceiling(ncol(plot_data)/2)-1)),
        labRow=rep("", nrow(plot_data)),

        colsep=c(0, ncol(plot_data)),
        #rowsep=rowsep,
        sepcolor="black",
        sepwidth=c(0.0001, 0.0001),

        key=plot_key,
        keysize=0.1,
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(0.9,1,1,1) / (par("cex") * 0.66),
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

        margins=c(0,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        
        col=my_palette,
        breaks=my_breaks,
        #add.expr=add_borders(i),
        useRaster=FALSE, # CHANGE BACK IF NEEDED
        par(xpd=FALSE)) # xpd - where plotting shows up

    if (plot_key) {
        title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
    }
    
}

set.seed(42)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]
title <- args[3]

# read in data
data <- read.table(data_file, header=TRUE)
data$region <- NULL
data$start_summit <- NULL
data$extra <- NULL

#data[data > 0] <- 1
data[data == -1] <- NA

data <- data[rowSums(data, na.rm=TRUE) != 0,]

#data_binary <- data
data_present <- data
data_present[!is.na(data_present)] <- 1

#data_binary[data_binary != 0] = 1

#data_sum <- colSums(data, na.rm=TRUE) / as.numeric(nrow(data))
data_sum <- colSums(data, na.rm=TRUE) / colSums(data_present, na.rm=TRUE)
data_sum[is.na(data_sum)] <- 0
data_sum <- data.frame(pos=(1:272 - 136), vals=data_sum)

test_plot_file <- paste(plot_file, "sum.pdf", sep="")
ggplot(data_sum, aes(x=pos, y=vals)) +
    geom_line(size=0.230) +
    labs(x="Motif position relative to ATAC summit",
         y="Motif affinity score",
         title=title) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_text(vjust=-1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0))
    #scale_y_continuous(limits=c(0,0.002), expand=c(0,0))
        
ggsave(test_plot_file, height=1, width=2)

q()

# subsample
#data <- data[sample(1:dim(data)[1], 1000, replace=FALSE),]
data <- data[1:1000,86:186]


data_melted <- melt(data)
my_breaks <- seq(
    quantile(data_melted$value, 0),
    quantile(data_melted$value, 1),
    length.out=color_granularity)
if (sum(my_breaks) == 0) {
    my_breaks <- seq(
        quantile(data_melted$value, 0),
        quantile(data_melted$value, 1),
        length.out=color_granularity)
}



# plot
pdf(plot_file, height=2, width=1.4)
plot_profile_heatmap(data, "", FALSE, my_breaks)
dev.off()


