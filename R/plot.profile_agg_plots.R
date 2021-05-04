#!/usr/bin/env Rscript

# description: take in agg values and plot
library(ggplot2)
library(RColorBrewer)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]

# read in data
data <- read.table(data_file, sep="\t", header=TRUE)
data$assay_day <- paste(data$assay, data$day, sep="_")

# plot all together
# TODO adjust colors
# TODO add title
atac_palette <- colorRampPalette(brewer.pal(9, "Blues"))(9)[5:7]
H3K27ac_palette <- colorRampPalette(brewer.pal(9, "Reds"))(9)[5:7]
H3K4me1_palette <- colorRampPalette(brewer.pal(9, "Oranges"))(9)[3:5]
H3K27me3_palette <- colorRampPalette(brewer.pal(9, "Greens"))(9)[5:7]
my_colors <- c(atac_palette, H3K27ac_palette, H3K27me3_palette, H3K4me1_palette)

plot_file <- paste(plot_prefix, ".ALL.pdf", sep="")
ggplot(data, aes(x=position, y=value, color=assay_day)) +
    geom_line(size=0.230) +
    labs(
        x="Position relative to TSS",
        y="log10 Fold Change",
        title="ALL") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        #axis.title=element_blank(),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_text(vjust=-1),
        #axis.line=element_blank(),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        #legend.position="bottom",
        #legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0)) +
    scale_colour_manual(values=my_colors)
                                        #+
    #scale_x_continuous(limits=c(-2500, 2500))
        
ggsave(plot_file, height=1, width=2.5)

# and also do this per day?
days <- c("d0", "d3", "d6")
day_colors <- c(atac_palette[3], H3K27ac_palette[3], H3K27me3_palette[3], H3K4me1_palette[3])
max_val <- max(data$value)
for (day_idx in 1:length(days)) {
    day <- days[day_idx]
    day_data <- data[data$day == day,]

    plot_file <- paste(plot_prefix, ".", day, ".pdf", sep="")
    ggplot(day_data, aes(x=position, y=value, color=assay)) +
        geom_line(size=0.230) +
        labs(
            x="Position relative to TSS",
            y="log10 Fold Change",
            title=day) +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.margin=margin(5,10,1,5),
            plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title=element_text(size=6, margin=margin(0,0,0,0)),
            axis.title.x=element_text(vjust=1.5),
            axis.title.y=element_text(vjust=-1),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.text=element_text(size=6),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            legend.position="bottom",
            legend.text=element_text(size=5),
            legend.title=element_blank(),
            legend.key.size=unit(0.05, "in"),
            legend.spacing.x=unit(0.05, "in"),
            legend.margin=margin(0,0,0,0)) +
        scale_colour_manual(values=day_colors) +
        scale_y_continuous(limit=c(0, max_val))
            
    ggsave(plot_file, height=1.5, width=1.5)
}








