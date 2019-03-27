#!/usr/bin/env/ Rscript

# description: maintain consistent styling in GGR figures
library(ggplot2)
library(RColorBrewer)
library(ggsci) # useful palettes: d3, igv, simpsons
library(viridis)

get_ggr_timepoint_colors <- function() {
    # set up color range
    num_timepoints <- 13
    #timepoint_palette <- colorRampPalette(brewer.pal(9, "YlOrBr"))(num_timepoints+4)[5:(num_timepoints+1)]
    timepoint_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(num_timepoints+4)[5:(num_timepoints+1)]
    #timepoint_palette <- colorRampPalette(brewer.pal(9, "RdPu"))(num_timepoints+4)[5:(num_timepoints+1)]
    #timepoint_palette <- colorRampPalette(brewer.pal(9, "PuBuGn"))(num_timepoints+4)[5:(num_timepoints+1)]
    #timepoint_palette <- colorRampPalette(brewer.pal(9, "PuBu"))(num_timepoints+4)[5:(num_timepoints+1)]
    #timepoint_palette <- colorRampPalette(brewer.pal(9, "BuPu"))(num_timepoints+4)[5:(num_timepoints+1)]
    return(timepoint_palette)
}


get_ggr_assay_palette <- function(color_string, granularity) {
    if (granularity != 50) {
        stop("GGR has fixed color granularity for this, do not adjust!")
    }
    if (color_string == "Blues") {
        assay_palette <- colorRampPalette(
            c("white", brewer.pal(9, color_string)))(granularity+3)[1:granularity-1]
    } else if (color_string == "Oranges") {
        assay_palette <- colorRampPalette(
            c("white", brewer.pal(9, color_string)))(granularity+15)[1:granularity-1]
    } else if (color_string == "Reds") {
        assay_palette <- colorRampPalette(
            c("white", brewer.pal(9, color_string)))(granularity+10)[1:granularity-1]
    } else if (color_string == "Greens") {
        assay_palette <- colorRampPalette(
            c("white", brewer.pal(9, color_string)))(granularity+10)[1:granularity-1]
    } else if (color_string == "Purples") {
        assay_palette <- colorRampPalette(
            c("white", brewer.pal(9, color_string)))(granularity-1)
    } else {
        assay_palette <- colorRampPalette(
            c("white", brewer.pal(9, color_string)))(granularity-1)
    }
    return(assay_palette)
}


ggplot_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    palette <- hcl(h=hues, l=65, c=100)[1:n]
    return(palette)
}


get_trajectory_palette <- function(num_trajectories) {
    if (num_trajectories == 15) {
        palette <- rev(viridis(num_trajectories))
    } else {
        palette <- rev(plasma(num_trajectories))
    }
    return(palette)
}


