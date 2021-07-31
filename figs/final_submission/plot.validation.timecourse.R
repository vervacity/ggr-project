#!/usr/bin/env Rscript

library(ggplot2)
library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]

# read in data
data <- read.table(data_file, sep="\t", header=TRUE)
plots <- unique(data$plot)

# plot color
palette <- brewer.pal(9, "Paired")
my_color <- palette[2]
my_color <- "black"

# make each plot
for (plot_i in 1:length(plots)) {
    plot <- plots[plot_i]
    print(plot)

    # set up title
    ids <- strsplit(plot, ".", fixed=TRUE)[[1]]
    title <- paste(ids[2], "/", ids[3], " motifs (", ids[4], ")", sep="")

    # get subset and set y_lim
    data_subset <- data[data$plot == plot,]
    data_subset <- data_subset[!is.na(data_subset$value),]
    y_max <- ceiling(max(data_subset$value) / 10000) * 10000
    y_max <- max(c(y_max, 25000))

    # plot
    ggplot(data_subset, aes(x=day, y=value)) +
        geom_bar(
            colour=my_color, size=0.230, width=0.5, fill=NA,
            stat="summary", fun="mean",
            show.legend=FALSE) +
        stat_summary(
            geom="errorbar", fun.data=mean_se,
            colour=my_color, size=0.230, width= 0.40) +
        geom_jitter(
            shape=21, size=0.5, colour=my_color,
            width=0.15, height=0,
            show.legend=FALSE) +
        labs(x="", y="RLU/ug lysate/lentiviral copy", title=title) +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.margin=margin(5,10,1,5),
            plot.title=element_text(size=8, hjust=0.5, margin=margin(b=5)), #1
            
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            
            axis.title=element_text(size=6, margin=margin(0,0,0,0)),
            axis.title.x=element_text(vjust=1.5),
            axis.title.y=element_text(vjust=1),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.text.x=element_text(size=6, angle=45, hjust=1, vjust=1),
            axis.text.y=element_text(size=6),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in")) +
    scale_y_continuous(expand=c(0,0), limits=c(0, y_max))

    # save
    plot_file <- paste(plot_prefix, ".", plot, ".pdf", sep="")
    ggsave(plot_file, height=1.5, width=1.5, useDingbats=FALSE)

}

