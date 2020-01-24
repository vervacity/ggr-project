#!/usr/bin/env Rscript

library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(gzfile(data_file), header=TRUE)

# params
dist_limits <- c(7, 25)

# make factor
data$sample_id <- factor(data$sample_id)

# get pwm name
pwm_name <- unlist(strsplit(basename(data_file), ".", fixed=TRUE))
pwm_name <- pwm_name[1]

# determine desired span
unique_distances <- unique(data$pwm_dist)
span <- 3. / length(unique_distances)

# plot
ggplot(data, aes(x=pwm_dist, y=prediction, group=sample_id)) +
    #geom_line(alpha=0.5, size=0.230, colour="black") +
    stat_smooth(
        geom="line", se=FALSE,
        span=span,
        alpha=0.5, size=0.230, colour="black") +
    stat_summary(
        data=data,
        aes(x=pwm_dist, y=prediction, group=1), fun.y=median, geom='line',
        size=0.460, colour="red") +
    labs(title=pwm_name, x="Distance between motifs (bp)", y="Predicted log2(FC)") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(size=0.115),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.background=element_blank(),
        legend.box.background=element_rect(colour="black", fill="white", size=0.115),
        legend.margin=margin(1,1,1,1),
        legend.key=element_blank(),
        legend.key.size=unit(0.1, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5),
        legend.position=c(0.9, 0.9),
        strip.background=element_blank(),
        strip.text=element_blank()) +
    scale_x_continuous(limits=dist_limits, expand=c(0,0))
ggsave(plot_file, height=1, width=2, useDingbats=FALSE)






