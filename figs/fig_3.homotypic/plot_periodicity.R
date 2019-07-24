#!/usr/bin/env Rscript

library(rhdf5)
library(reshape2)
library(ggplot2)
library(viridis)

data_file <- "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.background/ggr.scanmotifs.h5"
data_file <- "/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-17.footprinting/motifs.input_x_grad.lite/ggr.scanmotifs.h5"
#data_file <- "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.late/ggr.scanmotifs.h5"

# densities
#densities <- h5read(data_file, "sequence.active.pwm-hits.densities.max")[,1,]
activity <- h5read(data_file, "ATAC_SIGNALS.NORM")
#activity <- h5read(data_file, "H3K27ac_SIGNALS.NORM")
#pool_window <- 20

# pwm names
features <- h5read(data_file, "sequence-weighted.active.pwm-scores.thresh.max.idx", read.attributes=TRUE)
#features <- h5read(data_file, "sequence.active.pwm-scores.thresh.max.idx", read.attributes=TRUE)
pwm_names <- attr(features, "pwm_names")
#print(pwm_names)

# keep desired pwms
# test KLF, ZNF
pwm_indices <- c(169, 185) # ZNF, KLF GOOD
#pwm_indices <- c(171, 107) # TP53, ATF4
#pwm_indices <- c(169, 134) # ZNF750, CEBPD
#pwm_indices <- c(175, 1) # CEBPA, GATA
#pwm_indices <- c(146, 168) # REL, NFKB
#pwm_indices <- c(177, 146) # RUNX, REL
#pwm_indices <- c(177, 149) # RUNX, SMAD3 GOOD
#pwm_indices <- c(168, 149) # NFKB, SMAD3 GOOD
features <- features[1,pwm_indices,]

# make a plot for each time point
for (task_idx in 1:nrow(activity)) {
    print(task_idx)
    task_activity <- activity[task_idx,]

    # make a new dataframe with pwm and activity
    task_data <- data.frame(t(features))

    # attach activity
    task_data$activity <- task_activity
    
    # clean up
    task_data <- task_data[task_data$X1 != 432,]
    task_data <- task_data[task_data$X2 != 432,]
    task_data <- task_data[abs(task_data$X1 - task_data$X2) > 12,]
    task_data <- task_data[task_data$activity > 1,]
    task_data$present <- 1
    
    # center on one of the pwms
    task_data$dist <- task_data$X1 - task_data$X2
    print(dim(task_data))
    #print(head(task_data))
    agg_data <- aggregate(
        task_data$present,
        by=list(task_data$dist),
        sum)
    
    # and plot
    plot_file <- paste("task-", task_idx, ".periodicity.pdf", sep="")

    #ggplot(task_data, aes(x=dist, stat(count))) + # color=x, y=activity
    #ggplot(task_data, aes(x=Group.1, y=x)) + # color=x, y=activity
    ggplot(task_data, aes(x=dist, y=activity)) +
        #geom_density(adjust=1/10) +
        geom_point(stroke=0, shape=20, alpha=0.05, size=5) +
        #geom_line(data=agg_data, aes(x=Group.1, y=x)) +
        #stat_density_2d(aes(fill=..density..), geom="raster", contour=FALSE) +
        #stat_binhex(bins=50) +
        #scale_fill_distiller(palette= "Spectral", direction=1) +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            #panel.grid=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.text.y=element_text(size=5),
            axis.text.x=element_text(size=5),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.key.size=unit(0.1, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=6),
            legend.position=c(1.1, 0.9),
            strip.background=element_blank(),
            strip.text=element_blank()) #+
        #scale_color_viridis()
        #scale_x_continuous(limits=c(0,8), expand=c(0,0))#
        #scale_size_continuous(limits=limits) +
        #scale_colour_distiller(palette="Blues", direction=1, na.value="white") +
        #scale_fill_distiller(limits=limits, palette=palette, direction=1, na.value="white")
    ggsave(plot_file, height=1, width=10)

}



