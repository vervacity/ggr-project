

# description - plot out the datasets

library(ggplot2)
library(extrafont)

# args
args <- commandArgs(trailingOnly=TRUE)
dataset_file <- args[1]

# read in file
data <- read.table(dataset_file, sep="\t", header=TRUE)
data$timepoint_val <- as.numeric(sub("d", "", data$timepoint))

# test colors
library(RColorBrewer)
total_timepoints <- length(unique(data$timepoint_val))
my_palette <- colorRampPalette(brewer.pal(9, "Reds"))(total_timepoints+4)[5:(total_timepoints+4)]

# order in the order seen in the matrix
data$assay <- factor(data$assay, levels=rev(unique(data$assay)), ordered=TRUE)

# adjust timepoints
data$timepoint <- sub("d", "", data$timepoint)

# TODO consider plotting lines with specific hex colors

# plot out as sets
p <- ggplot() + geom_line(
    data=data,
    aes_string(
        group="assay",
        x="timepoint",
        y="assay"),
    size=1)
p <- p + geom_point(
    data=data,
    aes_string(
        group="timepoint",
        x="timepoint",
        y="assay",
        colour="timepoint"),
    shape=21,
    fill="white",
    #colour="black",
    size=3,
    stroke=1) +
    scale_colour_manual(values=my_palette)
p <- p + #scale_x_discrete(breaks=0:max(data$timepoint_val), labels=0:max(data$timepoint_val)) +
    ylab(NULL) + 
    xlab("Timepoint (days)") + 
    theme_bw() + 
    theme(
        text=element_text(family="ArialMT", size=10),
        axis.text.x=element_text(angle=30, hjust=1, colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.title.x=element_text(size=8),
        panel.grid.minor=element_blank(),
        legend.position="none")

#scale_x_discrete(breaks=0:max(data$timepoint_val), labels=0:max(data$timepoint_val)) +


ggsave("testing.pdf", height=1.75, width=4)
embed_fonts("testing.pdf")
