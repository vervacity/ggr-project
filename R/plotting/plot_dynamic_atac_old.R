#!/usr/bin/env Rscript


# Description: given soft ATAC clusters and hard histone clusters
# sort by those clusterings and output the ordered list and
# heatmap


library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(fastcluster)

args <- commandArgs(trailingOnly=TRUE)
atac_timepoints_file <- args[1]
deeptools_ordered_bed <- args[2]
out_plot_file <- args[3]
out_subsample_file <- args[4]

atac_timepoints <- read.table(gzfile(atac_timepoints_file), header=TRUE, row.names=1, stringsAsFactors=FALSE)
deeptools_ordering <- read.table(gzfile(deeptools_ordered_bed), header=FALSE, skip=1)


# adjust files
atac_timepoints <- data.frame(t(scale(t(atac_timepoints))))
atac_timepoints$region <- rownames(atac_timepoints)
deeptools_ordering <- data.frame(region=deeptools_ordering$V4, order=rownames(deeptools_ordering))

# merge and order by deeptools ordering
atac_timepoints_ordered <- merge(deeptools_ordering, atac_timepoints, by="region", all.x=TRUE, all.y=FALSE, sort=FALSE)

# and now plot out
rownames(atac_timepoints_ordered) <- atac_timepoints_ordered$region
atac_timepoints_ordered$region <- NULL
atac_timepoints_ordered$order <- NULL

evenly_spaced_subsample <- atac_timepoints_ordered[seq(1, nrow(atac_timepoints_ordered), 20), ]

# TODO save out subsample for histone plotting
write.table(rownames(evenly_spaced_subsample), file=gzfile(out_subsample_file), quote=FALSE, row.names=FALSE, col.names=FALSE)

# plotting
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
mylwid = c(2,6,2)
mylhei = c(0.5,12,1.5)

png(out_plot_file, height=18, width=6, units="in", res=200)
heatmap.2(
    as.matrix(evenly_spaced_subsample),
    Rowv=FALSE,
    Colv=FALSE,
    dendrogram="none",
    trace='none',
    density.info="none",
    keysize=0.1,
    key.title=NA,
    key.xlab=NA,
    key.par=list(pin=c(4,0.1),
        mar=c(6.1,0,5.1,0),
        mgp=c(3,2,0),
        cex.axis=2.0,
        font.axis=2),
    srtCol=45,
    cexCol=3.0,
    labRow="",
    margins=c(3,0),
    col=my_palette,
    lmat=mylmat,
    lwid=mylwid,
    lhei=mylhei)
dev.off()
