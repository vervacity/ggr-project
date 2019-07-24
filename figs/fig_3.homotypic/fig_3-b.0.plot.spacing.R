#!/usr/bin/env Rscript

library(rhdf5)
library(reshape2)
library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
sig_pwms_file <- args[1]
data_file <- args[2]
max_val_key <- "sequence-weighted.active.pwm-scores.thresh.max.val"
max_idx_key <- "sequence-weighted.active.pwm-scores.thresh.max.idx"
pwm_scores_key <- "sequence-weighted.active.pwm-scores.thresh"
raw_pwm_scores_key <- "sequence.active.pwm-scores.thresh"
signal_key <- "ATAC_SIGNALS.NORM"
left_clip <- 420 + 12
num_positions <- 200
final_extend_len <- 160
out_prefix <- "fig_5-b"

# read in sig pwms
sig_pwms <- read.table(sig_pwms_file, sep="\t", header=TRUE, row.names=1)
sig_pwms <- rownames(sig_pwms)
print(sig_pwms)

# get pwm names
max_vals <- h5read(data_file, max_val_key, read.attributes=TRUE)
max_indices <- h5read(data_file, max_idx_key, read.attributes=TRUE)
pwm_names <- attr(max_vals, "pwm_names")

# analyze each pwm
for (pwm_idx in 1:length(sig_pwms)) {
    
    # name and global index
    pwm_name <- sig_pwms[pwm_idx]
    global_pwm_idx <- which(pwm_name == pwm_names)[[1]]
    global_pwm_idx <- global_pwm_idx[1]
    print(pwm_name)
    
    # check to see which have this pwm and also adjust indices
    pwm_max_vals <- aperm(max_vals)[,global_pwm_idx,1]
    example_indices <- which(pwm_max_vals != 0)
    pwm_max_indices <- aperm(max_indices)[example_indices,global_pwm_idx,1]
    
    # select out this slice
    pwm_scores <- h5read(data_file, pwm_scores_key, index=list(global_pwm_idx, NULL, NULL, example_indices))
    raw_pwm_scores <- h5read(data_file, raw_pwm_scores_key, index=list(global_pwm_idx, NULL, NULL, example_indices))
    score_length <- dim(pwm_scores)[2]
    
    # make zeros arrays
    aligned_array <- matrix(0, nrow=length(example_indices), ncol=2*(num_positions))
    aligned_raw_array <- matrix(0, nrow=length(example_indices), ncol=2*(num_positions))

    # insert each row with offset
    for (example_idx in 1:nrow(aligned_array)) {

        # get offset, and scores
        offset <- num_positions - (pwm_max_indices[example_idx] - left_clip)
        scores <- pwm_scores[1,,,example_idx]
        scores <- apply(scores, 1, max)
        aligned_array[example_idx,offset:(offset+score_length-1)] <- scores

        # raw
        scores <- raw_pwm_scores[1,,1,example_idx]
        aligned_raw_array[example_idx,offset:(offset+score_length-1)] <- scores
    }

    # center on max peak
    means <- apply(aligned_array, 2, mean)
    max_idx <- which(means == max(means))[1]
    clip_start <- max_idx - final_extend_len
    clip_stop <- max_idx + final_extend_len
    aligned_array <- aligned_array[,clip_start:clip_stop]
    aligned_raw_array <- aligned_raw_array[,clip_start:clip_stop]
    
    # and now build mean array
    max_pos <- (dim(aligned_array)[2] - 1) / 2
    
    means <- apply(aligned_array, 2, mean)
    means <- data.frame(mean=means, pos=-max_pos:max_pos)
    means$scores <- "weighted"
    means$mean[155:165] <- 0
    
    raw_means <- apply(aligned_raw_array, 2, mean)
    raw_means <- data.frame(mean=raw_means, pos=-max_pos:max_pos)
    raw_means$scores <- "original"
    raw_means$mean[155:165] <- 0
    
    # pull together
    means$mean <- means$mean / max(means$mean)
    raw_means$mean <- raw_means$mean / max(raw_means$mean)
    means <- rbind(means, raw_means)
    
    # plot
    plot_file <- paste(out_prefix, ".", pwm_name, ".spacings.pdf", sep="")
    ggplot(means, aes(x=pos, y=mean, colour=scores)) +
        geom_line(size=0.230) + 
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.minor=element_blank(),
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
            legend.position=c(0.9, 0.9),
            strip.background=element_blank(),
            strip.text=element_blank()) +
        scale_x_continuous(limits=c(-100,100), breaks=seq(-100, 100, 10), expand=c(0,0)) +
        #scale_size_continuous(limits=limits) +
        #scale_colour_distiller(limits=limits, palette=palette, direction=1, na.value="white") +
        #scale_fill_distiller(limits=limits, palette=palette, direction=1, na.value="white")
    ggsave(plot_file, height=1.5, width=5.5)


    # TODO collect results and plot heatmap
    
}
