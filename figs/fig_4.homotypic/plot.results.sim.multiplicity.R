#!/usr/bin/env Rscript

library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(gzfile(data_file), header=TRUE)
#print(head(data))

# make factor
data$sample_id <- factor(data$sample_id)

# get pwm name
pwm_name <- unlist(strsplit(basename(data_file), ".", fixed=TRUE))
pwm_name <- pwm_name[1]

# plot
ggplot(data, aes(x=pwm_count, y=prediction, group=sample_id)) +
    geom_line(alpha=0.5, size=0.230, colour="black") +
    stat_summary(
        data=data,
        aes(x=pwm_count, y=prediction, group=1), fun.y=median, geom='line',
        size=0.460, colour="red") +
    labs(title=pwm_name, x="Motif count", y="log2(Fold change)") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
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
    scale_x_continuous(limits=c(1,6), expand=c(0,0))
ggsave(plot_file, height=1.5, width=1, useDingbats=FALSE)






