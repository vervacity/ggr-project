#!/usr/bin/env Rscript

# fig 3b
# description: plot importance scores in SNPs and delta output
# data: /srv/scratch/dskim89/ggr/ggr.tronn.2019-06-11.ase

library(ggplot2)
library(reshape2)
library(ggsci)

# args
args <- commandArgs(trailingOnly=TRUE)
imptscore_file <- args[1]
sig_file <- args[2]

# plot the impt score results
data <- read.table(imptscore_file, header=TRUE, stringsAsFactors=FALSE)

# extra filter?
#data <- data[!(data$sig == 1 & data$variants.sig == 0),]
data <- data[data$variant.impt_score.present > 0.10,]
print(dim(data))

# totals
num_sig <- sum(data$sig == 1)
sig_name <- paste("asATAC SNP\n(n=", num_sig, ")", sep="")

num_non_sig <- sum(data$sig == 0)
non_sig_name <- paste("non-asATAC SNP\n(n=", num_non_sig, ")", sep="")

# clean up
data$delta_logit <- abs(data$delta_logit)
data$sig[data$sig == 1] <- sig_name
data$sig[data$sig == 0] <- non_sig_name
data$sig <- factor(data$sig, levels=c(non_sig_name, sig_name))

plot_file <- "fig_3.asATAC.imptscore_present.pdf"
dodge <- position_dodge(width=0.9)
ggplot(data, aes(x=sig, y=delta_logit, colour=sig)) +
    geom_violin(size=0.115, width=0.8, position=dodge, scale="width") +
    geom_boxplot(size=0.115, width=0.2, position=dodge, outlier.shape=NA, show.legend=FALSE) +
    #geom_point(alpha=0.5, shape=16, stroke=0, size=0.5, position=position_jitterdodge(jitter.width=0.15, dodge.width=0.9)) +
    labs(title="SNPs with impt scores", y="Predicted difference |ref - alt|") +        
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=3)),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, angle=45, hjust=1, vjust=1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="none",
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=4),
        legend.text=element_text(size=4)) +
    scale_y_continuous(limits=c(-0.05, 1), expand=c(0,0)) + 
    scale_color_jama()
ggsave(plot_file, width=1.25, height=1.75)


# plot the beta max (sig only) 
data <- read.table(sig_file, header=TRUE, stringsAsFactors=FALSE)

# extra filters?
data <- data[abs(data$variant.impt_score.present) > 0.10 | (data$variant.impt_score.present == 0),]
print(dim(data))

# totals
num_impt <- sum(data$variant.impt_score.present != 0)
impt_name <- paste("+Impt score\n(n=", num_impt, ")", sep="")

num_non_impt <- sum(data$variant.impt_score.present == 0)
non_impt_name <- paste("-Impt score\n(n=", num_non_impt, ")", sep="")

# clean up
data$beta_max <- abs(data$beta_max)
data$variant.impt_score.present[data$variant.impt_score.present != 0 ] <- impt_name
data$variant.impt_score.present[data$variant.impt_score.present == 0 ] <- non_impt_name
data$variant.impt_score.present <- factor(
    data$variant.impt_score.present,
    levels=c(non_impt_name, impt_name))

plot_file <- "fig_3.asATAC.sig_only.pdf"
dodge <- position_dodge(width=0.9)
ggplot(data, aes(x=variant.impt_score.present, y=beta_max, colour=variant.impt_score.present)) +
    geom_violin(size=0.115, width=0.8, position=dodge, scale="width") +
        geom_boxplot(size=0.115, width=0.2, position=dodge, outlier.shape=NA, show.legend=FALSE) +
    geom_point(alpha=0.5, shape=16, stroke=0, size=0.5, position=position_jitterdodge(jitter.width=0.15, dodge.width=0.9)) +
    labs(title="SNPs with asATAC", y="Effect size (Beta)") +        
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=3)),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title.y=element_text(size=6),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, angle=45, hjust=1, vjust=1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="none",
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=4),
        legend.text=element_text(size=4)) +
    scale_y_continuous(limits=c(1.35, 3), expand=c(0,0)) +
    scale_color_jama()
ggsave(plot_file, width=1.25, height=1.75)
