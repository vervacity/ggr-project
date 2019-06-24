#!/usr/bin/env Rscript

# description: code to plot out PWM scores and RNA scores
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim
library(gplots)
library(RColorBrewer)
library(reshape2)

library(grid)
library(gridGraphics)
library(gridExtra)

library(rhdf5)

# args
args <- commandArgs(trailingOnly=TRUE)
h5_file <- args[1]
main_group <- args[2]
col_labels <- c("d0.0", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0")
out_dir <- "."
debug <- FALSE

make_heatmap <- function(data, data_max, use_colsep, use_labrow, rowlab_left) {

    # grid
    mylmat = rbind(c(0,0,3,0),c(0,1,2,0),c(0,4,5,0))
    mylwid = c(0,0.07,0.7,0.50) # 1.5
    mylhei = c(0.10,2,0.25) # 3
    
    if (!is.null(data_max)) {
        # for the rowside colors
        # make sure the values lie between 0-1
        pal <- colorRamp(brewer.pal(9, "RdBu"))
        data_max <- 1 - data_max
        rowsidecolors <- rgb(pal(data_max), maxColorVal=256)
        key <- TRUE
    } else {
        rowsidecolors <- NULL
        key <- FALSE
    }

    # adjust this for multi-column or single column
    if (use_colsep) {
        colsep <- 0:(ncol(data)+1)
        colsepwidth <- 0.001
        labCol <- colnames(data)
    } else {
        colsep <- NULL
        colsepwidth <- 0.0
        labCol <- c(colnames(data)[1], "")
    }

    # adjust for row labels
    if (use_labrow) {
        labRow <- rownames(data)
    } else {
        labRow <- rep("", nrow(data))
    }
    
    # adjust for rowlab_left
    if (rowlab_left) {
        mylwid = c(0.5,0.07,0.7,0.05)
        offsetRow <- -10
        my_palette <- colorRampPalette(brewer.pal(9, "Oranges"))(49)
    } else {
        offsetRow <- -0.2
        my_palette <- colorRampPalette(brewer.pal(9, "Purples"))(49)
    }
    offsetCol <- 0
    
    # adjust color again
    if (!key) {
        my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(49)
        if (colsepwidth > 0) {
            mylhei = c(0.05,1,0.9)
            offsetRow <- -0.2
            offsetCol <- 10
        }
    }
    
    # plot
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        colsep=colsep,
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001,colsepwidth),

        labRow=labRow,
        cexRow=0.5,
        offsetRow=offsetRow,
        
        labCol=labCol,
        cexCol=0.5,
        srtCol=45,
        offsetCol=offsetCol,

        key=key,
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(2,6,1.7,6),
            mgp=c(0,-0.1,0),
            tcl=-0.1,
            lend=2,
            cex.axis=0.6,
            bty="n"),
        key.xtickfun=function() {
            breaks <- pretty(parent.frame()$breaks)
            breaks <- breaks[c(1,length(breaks))]
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)},

        margins=c(1,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        col=my_palette,
        RowSideColors=rowsidecolors)

}

# grab grob fn
grab_grob <- function(fn) {
    #grid.echo(fn)
    grid.grabExpr(grid.echo(fn))
    #grid.grab()
}


# find groups
metadata <- h5ls(h5_file)
groups <- metadata$name[metadata$group == paste("/", main_group, sep="")]

# run for each group (label set)
for (group_i in 1:length(groups)) {
    group <- paste("/", main_group, "/", groups[group_i], sep="")
    print(group)

        # extract relevant data
        group_data <- h5read(h5_file, group, read.attributes=TRUE)
        
        # ===================
        # 1) correlation plot
        # ===================
        correlation_scores <- aperm(group_data$correlations)
        rownames(correlation_scores) <- attr(group_data, "hgnc_ids")

        # filter function for good correlation as needed
        good_cor <- correlation_scores > -2
        if (!any(good_cor)) {
            good_cor <- correlation_scores > -2
        }
        
        if (sum(good_cor, na.rm = TRUE) == 1) {
            good_cor <- correlation_scores > -2
        }

        # adjust correlation score from [-1,1] to [0,1]
        correlation_scores <- (correlation_scores + 1) / 2

        # set up placeholder
        placeholder <- data.frame(correlation=correlation_scores, cor.dup=correlation_scores)

        # filter
        correlation_scores <- correlation_scores[good_cor]
        placeholder <- placeholder[good_cor,]

       # ===================
        # 2) pwm scores (w max vector)
        # ===================
        pwm_scores <- aperm(group_data$pwm_patterns)

        # adjust pwm names and put as rownames
        pwm_names <- attr(group_data, "pwm_names")
        pwm_names <- gsub("HCLUST-\\d+_", "", pwm_names)
        pwm_names <- gsub(".UNK.*", "", pwm_names)
        pwm_names <- paste(pwm_names, "MOTIF", sep=" ")
        rownames(pwm_scores) <- pwm_names
        colnames(pwm_scores) <- col_labels

        # min/max norm
        pwm_scores_max <- apply(pwm_scores, 1, max)
        pwm_scores_min <- apply(pwm_scores, 1, min)
        pwm_scores <- (pwm_scores - pwm_scores_min) / (pwm_scores_max - pwm_scores_min)

        # adjust max vector
        pwm_unit_max <- pwm_scores_max / max(pwm_scores_max)
        pwm_unit_max <- pwm_unit_max / 2 + 0.5

        # filter
        pwm_scores <- pwm_scores[good_cor,]
        pwm_unit_max <- pwm_unit_max[good_cor]

        # reorder rows by max point across pattern
        pwm_scores_peaks <- apply(pwm_scores, 1, which.max)
        reorder_indices <- order(pwm_scores_peaks)
    
        if (TRUE) {
            correlation_scores <- correlation_scores[reorder_indices]
            placeholder <- placeholder[reorder_indices,]
            pwm_unit_max <- pwm_unit_max[reorder_indices]
            pwm_scores <- pwm_scores[reorder_indices,]
        }

        # ===================
        # 3) RNA scores (w max vector)
        # ===================
        rna_scores <- aperm(group_data$rna_patterns)
        rownames(rna_scores) <- attr(group_data, "hgnc_ids")
        colnames(rna_scores) <- col_labels
        
        # min/max norm
        rna_scores_max <- apply(rna_scores, 1, max)
        rna_scores_min <- apply(rna_scores, 1, min)
        rna_scores <- (rna_scores - rna_scores_min) / (rna_scores_max - rna_scores_min)

        # adjust max vector
        rna_unit_max <- rna_scores_max / max(rna_scores_max)
        rna_unit_max <- rna_unit_max / 2 + 0.5

        # filter
        rna_scores <- rna_scores[good_cor,]
        rna_unit_max <- rna_unit_max[good_cor]

        # reorder rows
        if (TRUE) {
            rna_unit_max <- rna_unit_max[reorder_indices]
            rna_scores <- rna_scores[reorder_indices,]
        }

        # ===================
        # 4) other outputs
        # ===================
        keys <- metadata$name[metadata$group == paste(group, "/other", sep="")]
        for (key_i in 1:length(keys)) {
            key_data <- unlist(group_data$other[keys[key_i]])
            if (key_i == 1) {
                outputs <- key_data
                output_names <- keys[key_i]
            } else {
                outputs <- rbind(outputs, key_data)
                output_names <- c(output_names, keys[key_i])
            }
        }
        rownames(outputs) <- output_names
        colnames(outputs) <- col_labels
        
        fn4 <- function() make_heatmap(outputs, NULL, TRUE, TRUE, TRUE)

        # ===================
        # 5) plots
        # ===================

        plot_prefix <- paste(out_dir, "/fig_3-c.", groups[group_i], sep="")

    # plot pwm scores
    fn1 <- function() make_heatmap(pwm_scores, pwm_unit_max, TRUE, TRUE, TRUE)
    if (debug) {
        pdf(paste(plot_prefix, ".pwm_scores.pdf", sep=""))
        fn1()
        dev.off()
    }
    
    # plot correlation
    fn3 <- function() make_heatmap(placeholder, NULL, FALSE, FALSE, FALSE)
    if (debug) {
        pdf(paste(plot_prefix, ".pwm_rna_cor.pdf", sep=""))
        fn3()
        dev.off()
    }
    
    # plot rna
    fn2 <- function() make_heatmap(rna_scores, rna_unit_max, TRUE, TRUE, FALSE)
    if (debug) {
        pdf(paste(plot_prefix, ".rna_scores.pdf", sep=""))
        fn2()
        dev.off()
    }

    # now try plot all together
    grob_list <- list(
        "1"=grab_grob(fn1),
        "2"=grab_grob(fn3),
        "3"=grab_grob(fn2),
        "4"=grab_grob(fn4))
    
    pdf(
        file=paste(plot_prefix, ".all.pdf", sep=""),
        height=3.5, width=2.75, onefile=FALSE, useDingbats=FALSE, family="ArialMT")
    plot.new()
    grid.newpage()
    grid.arrange(grobs=grob_list, nrow=2, ncol=3, heights=c(3, 0.25), widths=c(1, 0.08, 1), clip=FALSE)
    title(groups[group_i], adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
    dev.off()

}











