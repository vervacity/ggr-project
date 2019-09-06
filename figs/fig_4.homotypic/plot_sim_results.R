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

# plot
ggplot(data, aes(x=pwm_count, y=prediction, color=sample_id)) +
    geom_line()
ggsave(plot_file)






