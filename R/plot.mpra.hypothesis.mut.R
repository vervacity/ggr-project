#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]

# params
days <- c("d0", "d3", "d6")

# load data
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)
data_melt <- melt(data, id.vars=c("example_id", "combos"), na.rm=TRUE)
data_melt$day <- gsub("_.+", "", data_melt$variable)
data_melt$rep <- gsub(".+_", "", data_melt$variable)
#print(head(data_melt))

# aggregate reps
if (FALSE) {
    data_melt <- aggregate(
        data_melt$value,
        by=list(data_melt$example_id, data_melt$combos, data_melt$day),
        mean)
    colnames(data_melt) <- c("example_id", "combos", "day", "value")
}

# plot by day
for (day_idx in 1:length(days)) {
    day <- days[day_idx]
    day_prefix <- paste(plot_prefix, ".", day, sep="")
    day_data <- data_melt[data_melt$day == day,]

    if (median(day_data$value[day_data$combos=="0,0"]) < 0) {
        next
    }

    # figure out median of background seq (1,1) and normalize to that
    correction_val <- median(day_data$value[day_data$combos=="1,1"])
    day_data$value <- day_data$value - correction_val
    
    plot_file <- paste(day_prefix, ".pdf", sep="")
    ggplot(day_data, aes(x=combos, y=value, colour=rep)) +
        #geom_line(alpha=0.7, size=0.115, colour="gray", aes(group=example_id)) +
        geom_violin(position=position_dodge(), trim=FALSE, scale="width") +
        geom_boxplot(position=position_dodge(), alpha=0.5, width=0.3, outlier.shape=NA) +
        geom_point(position=position_jitterdodge(jitter.width=0), alpha=0.5) +
        theme_bw()        
    ggsave(plot_file)
    

}






