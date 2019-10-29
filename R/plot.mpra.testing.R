#!/usr/bin/env Rscript


library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read data
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)

# set up
data$labels <- paste(data$pwm1_clean, data$pwm2_clean, sep=",")
data <- data[data$interaction != "FAILED.NEGATIVE_ENDOG",]
data <- data[data$interaction != "FAILED.TRAJ",]

# plot
ggplot(data, aes(x=expected, y=actual, colour=interaction)) +
    geom_point() +
    geom_abline(intercept=0, slope=1) +
    geom_text(aes(label=labels, vjust=1, hjust=1), size=1) +
    coord_fixed() +
    theme_bw()
ggsave(plot_file)


