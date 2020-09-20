#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
#library(ggsci)

plot_fn <- function(data, plot_file, xlab) {

    # get R vals (spearman and pearson)
    r_s <- cor(data$NN, data$MPRA, method="spearman")
    r_s <- format(round(r_s, 3), nsmall=3)
    r_s <- paste("R[S]==", r_s, sep="")
    
    r_p <- cor(data$NN, data$MPRA, method="pearson")
    r_p <- format(round(r_p, 3), nsmall=3)
    r_p <- paste("R[P]==", r_p, sep="")

    
    p <- ggplot(data, aes(x=NN, y=MPRA)) +
        geom_hex(bins=30, show.legend=FALSE) +
        #geom_point(
        #    shape=20, size=0.5, alpha=0.5, show.legend=FALSE) +
        geom_smooth(method = "lm", se=TRUE, size=0.5, colour="black") +
        annotate(
            "text", x=1.1, y=7.3, size=1.5, vjust="inward", hjust="inward",
            label=r_s, parse=TRUE) +
        annotate(
            "text", x=1.1, y=6.8, size=1.5, vjust="inward", hjust="inward",
            label=r_p, parse=TRUE) +
        labs(
            x=xlab, y="MPRA",
            title=paste("Correlating MPRA expression \nto ", xlab, sep="")) + 
           #coord_fixed() +
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
            axis.title.y=element_text(vjust=1),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.text=element_text(size=6),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.key.size=unit(0.1, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.spacing.x=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=5),
            legend.position="bottom") +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0), limits=c(1,4))
    
    ggsave(plot_file, height=1.75, width=1.75, useDingbats=FALSE)

}

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# load data
data <- read.table(data_file, sep="\t", header=TRUE)

# corr val
r <- cor.test(data$NN, data$MPRA, method="spearman")
print(r)

xlab <- "NN"

# plot
plot_fn(data, plot_file, xlab)
