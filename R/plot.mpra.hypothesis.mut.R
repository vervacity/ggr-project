#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(ggsci)

plot_fn <- function(data, plot_file) {

    p <- ggplot(data, aes(x=pwms, y=value, colour=pwms)) +
        geom_hline(yintercept=0, size=0.115, linetype="dashed") +
        #geom_line(alpha=0.7, size=0.115, colour="gray", aes(group=example_id)) +
        #geom_violin(
        #    position=position_dodge(), size=0.230, trim=FALSE, scale="width", show.legend=FALSE) +
        geom_boxplot(
            position=position_dodge(), size=0.230, width=0.3, outlier.shape=NA, show.legend=FALSE) +
        geom_point(
            position=position_jitterdodge(jitter.width=0.01), size=0.5, alpha=0.5, show.legend=FALSE) +
        ylab("log2(FC)") +
        xlab("Combinations") +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.margin=margin(5,10,6,5),
            plot.title=element_text(size=8, margin=margin(b=1)),
            
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            
            axis.title=element_text(size=6, margin=margin(0,0,0,0)),
            axis.title.x=element_text(vjust=1.5),
            axis.title.y=element_text(vjust=1),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            #axis.text.x=element_text(size=6, angle=30, vjust=0, hjust=1),
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
        scale_colour_brewer(palette="Dark2") +
        coord_flip()
            
    ggsave(plot_file, height=1.5, width=2, useDingbats=FALSE)

}


# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
day <- args[2]
plot_prefix <- args[3]

# load data
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)

# set levels to order results correctly
combo_order <- c("1,0", "0,1", "0,0")
pwm_order <- c()
for (i in 1:length(combo_order)) {
    combo_val <- combo_order[i]
    pwm <- unique(as.character(data[data$combo == combo_val, "pwms"]))[1]
    pwm_order <- c(pwm_order, pwm)
}

# setup
keep_headers <- c("example_id", "pwms", "combos", day)
data <- data[keep_headers]

# take out 1,1 and -1
plot_file <- paste(plot_prefix, ".normalized.pdf", sep="")
data_filt <- data[(data$combos != "1,1") & (data$combos != "-1"),]
data_filt$pwms <- factor(as.character(data_filt$pwms), levels=pwm_order) # order
colnames(data_filt) <- c("example_id", "pwms", "combos", "value")
plot_fn(data_filt, plot_file)



