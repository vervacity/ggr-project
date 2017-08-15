#!/usr/bin/env Rscript

# Description: plot nice trajectories

library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
trajectory_files <- args[2:length(args)]


for (i in 1:length(trajectory_files)) {

    # read in file
    data_z <- read.table(gzfile(trajectory_files[i]), header=FALSE, row.names=1, stringsAsFactors=FALSE)
    colnames(data_z) <- c("d0.0", "d0.5", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0")
    data_z$region <- rownames(data_z)
    prefix <- unlist(strsplit(basename(trajectory_files[i]), ".txt", fixed=TRUE))[1]
    
    # plot with ggplot2
    data_melted <- melt(data_z, id="region")
    data_melted$variable <- sub("d", "", data_melted$variable)
    data_melted$variable <- as.numeric(as.character(data_melted$variable))

    plot_file <- paste(out_dir, "/", prefix, ".plot.png", sep="")
    if (!file.exists(plot_file)) {

        ggplot(data_melted, aes(x=variable, y=value, group=region)) + geom_line(colour="gray44") +
            theme_bw() +
            theme(
                axis.line=element_line(colour="black"),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank()) +
            ggtitle(paste("Soft trajectory", prefix)) + 
            theme(plot.title=element_text(hjust=0.5)) + 
            xlab("Time (days)") + 
            ylab("Openness (z-score)") + 
            scale_x_continuous(limits=c(0,6), expand=c(0,0)) + 
            scale_y_continuous(expand=c(0,0)) + 
            stat_summary(aes(y=value, group=1), fun.y=mean, colour="black", geom="line", group=1)
                
        ggsave(plot_file)

    }

}






