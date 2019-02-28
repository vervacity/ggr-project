#!/usr/bin/env/ Rscript

# description: maintain consistent styling in GGR figures
library(ggplot2)
library(RColorBrewer)


get_ggr_timepoint_colors <- function() {
    # set up color range
    num_timepoints <- 13
    timepoint_palette <- colorRampPalette(brewer.pal(9, "Reds"))(num_timepoints+4)[5:(num_timepoints+1)]
    return(timepoint_palette)
}







