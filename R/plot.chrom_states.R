#!/usr/bin/env Rscript

# description: plot chrom states
library(gplots)
library(RColorBrewer)

library(grid)
library(gridGraphics)
library(gridExtra)

# load style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
chrom_states_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(chrom_states_file, sep="\t", header=TRUE)
print(head(data))

# heatmap fn
plot_profile_heatmap <- function(plot_data, color_set) {

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.05,1,0.05)
    mylhei = c(0.25,4,0.5)

    # get GGR colors
    my_palette <- get_ggr_assay_palette(color_set, 50)

    # adjust color range
    if (color_set != "Purples") {
        my_breaks <- seq(6, max(plot_data), length.out=50)
        print(my_breaks)
    } else {
        my_breaks=NULL
    }
    
    # heatmap
    heatmap.2(
        as.matrix(plot_data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        #labCol=c(
        #    rep("", floor(ncol(plot_data)/2)),
        #    label[i],
        #    rep("", ceiling(ncol(plot_data)/2)-1)),
        labRow=rep("", nrow(plot_data)),
        keysize=0.1,
        key.title=NA,
        key.xlab=NA,
        key.par=list(pin=c(4,0.1),
            mar=c(2.1,0,2.1,0),
            mgp=c(3,1,0),
            cex.axis=1.0,
            font.axis=2),
        srtCol=45,
        cexCol=1.25,
        margins=c(1,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        col=my_palette,
        breaks=my_breaks,
        #rowsep=rowsep,
        sepcolor="black")
        #add.expr=add_borders(i),
        #useRaster=TRUE,
       # par(xpd=FALSE)) # xpd - where plotting shows up

}


# grab grob fn
grab_grob <- function(fn) {
    grid.grabExpr(grid.echo(fn))
}

# sets
start_stop_sets <- c("2,2", "3,17", "18,27", "28,30", "32,34", "36,38")
colors <- c("Blues", "Purples", "Blues", "Reds", "Oranges", "Greens")

# get grobs
grob_list <- lapply(
    1:length(start_stop_sets),
    function(i) {
        print(i)
        
        # get coords
        coords <- as.numeric(unlist(strsplit(start_stop_sets[i], ",")))
        start <- coords[1]
        end <- coords[2]

        # get data subset
        data_subset <- data[,start:end]
        print(head(data_subset))

        if (is.null(dim(data_subset))) {
            data_subset <- data.frame(x=data_subset, y=data_subset)
        }
        
        # plot heatmap and save out
        fn <- function() plot_profile_heatmap(data_subset, colors[i])

        # grab grob
        return(grab_grob(fn))
    })

# plot joint heatmap
pdf(file=plot_file, height=7, width=3, onefile=FALSE, family="ArialMT")
grid.newpage()
grid.arrange(grobs=grob_list, nrow=1, ncol=length(start_stop_sets), clip=TRUE)
dev.off()
