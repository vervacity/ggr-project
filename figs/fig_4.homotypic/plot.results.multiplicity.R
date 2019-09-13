#!/usr/bin/env Rscript

library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data, adjust
data <- read.table(gzfile(data_file), header=TRUE)
data$sample_id <- factor(data$sample_id)

# get pwm name
pwm_name <- unlist(strsplit(basename(data_file), ".", fixed=TRUE))
pwm_name <- pwm_name[1]

# plot
ggplot(data, aes(x=pwm_count, y=prediction, group=sample_id)) +
    geom_line(alpha=0.3, size=0.230, colour="black") +
    stat_summary(
        data=data,
        aes(x=pwm_count, y=prediction, group=1), fun.y=median, geom='line',
        size=0.460, colour="red") +
    labs(title=pwm_name, x="Motif count", y="Predicted log2(FC)") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        #panel.grid.major.y=element_blank(),
        #panel.grid.major.x=element_line(size=0.115),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in")) +
    scale_x_continuous(limits=c(1,6), expand=c(0,0))

ggsave(plot_file, height=1.5, width=1, useDingbats=FALSE)






