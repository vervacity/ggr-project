#!/usr/bin/env Rscript

library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
logits_file <- args[1]
plot_file <- args[2]

# read data
data <- read.table(gzfile(logits_file), sep="\t", header=TRUE)
data$mut_type <- factor(data$mut_type, levels=c("endog", "motif,scr", "scr,motif", "scr,scr"))
background_val <- data$logit[data$mut_type == "scr,scr"]

# plot
ggplot(data, aes(x=mut_type, y=logits, fill=mut_type, colour=mut_type)) +
    geom_hline(yintercept=background_val, size=0.115, linetype="dashed") +
    geom_col(alpha=0.5, size=0.115, width=0.6, show.legend=FALSE) +
    labs(title="Predicted\naccessibility", x="", y="Log2(FC)") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(hjust=0.5, size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, vjust=1, hjust=1, angle=60),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in")) +
    scale_fill_brewer(palette="Dark2") +
    scale_colour_brewer(palette="Dark2") +
    scale_y_continuous(expand=c(0,0))

ggsave(plot_file, height=1.2, width=0.75, units="in", useDingbats=FALSE)

