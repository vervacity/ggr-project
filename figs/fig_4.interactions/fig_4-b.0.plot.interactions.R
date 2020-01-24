#!/usr/bin/env Rscript


library(ggplot2)
library(reshape2)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

ggr_colors <- get_ggr_timepoint_colors()
print(ggr_colors)
unique_indices <- c(0,1,2,3,4,5,6,9,10,12)
category_levels <- c("buffer", "synergy", "additive")

# args
args <- commandArgs(trailingOnly=TRUE)
endo_file <- args[1]
sims_file <- args[2]


# plot endogenous summary
plot_file <- "fig_4-b.0.endogenous_seq.pdf"
data_endo <- read.table(endo_file, header=TRUE)
data_endo$best_task_index <- as.factor(data_endo$best_task_index)
data_endo$pval <- -log10(data_endo$pval)
endo_task_indices <- as.numeric(as.character(data_endo$best_task_index))
data_endo$best_task_index <- factor(
    endo_task_indices,
    levels=unique_indices)
endo_colors <- ggr_colors[unique_indices %in% unique(endo_task_indices)]
data_endo$category <- factor(data_endo$category, levels=category_levels)
ggplot(data_endo, aes(x=expected, y=actual)) +
    geom_point(shape=21, aes(size=pval, colour=category, fill=best_task_index)) +
    geom_abline() +
    xlim(0,1.5) +
    ylim(0,1.5) +
    theme_bw() +
    theme(
        aspect.ratio=1) +
    scale_fill_manual(values=endo_colors)
ggsave(plot_file)

# different plot for endogenous - rank order?
data_endo <- data_endo[order(data_endo$diff),]
data_endo$diff_order <- nrow(data_endo):1
plot_file <- "fig_4.b.0.endogenous_seq.TEST.pdf"
#ggplot(data_endo, aes(x=diff_order, y=diff)) +
ggplot(data_endo, aes(x=expected, y=diff)) +
    geom_hline(size=0.115, yintercept=0, linetype="dashed") +
    #geom_point(shape=21, aes(size=pval, fill=best_task_index)) +
    geom_point(shape=21, stroke=0.115, aes(size=pval, fill=best_task_index)) +
    #geom_bar(stat="identity", aes(fill=best_task_index)) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
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
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +
    scale_fill_manual(values=endo_colors) +
    scale_size_continuous(range=c(0,2)) +
    scale_x_continuous(limits=c(0, 1.5), expand=c(0,0)) +
    scale_y_continuous(limits=c(-0.2, 0.3), expand=c(0,0))
ggsave(plot_file, height=1, width=2)
data_endo$diff_order <- NULL



# plot sims summary
plot_file <- "fig_4-b.1.simulations_seq.pdf"
data_sims <- read.table(sims_file, header=TRUE)
data_sims$best_task_index <- as.factor(data_sims$best_task_index)
data_sims$pval <- -log10(data_sims$pval)
data_sims$best_task_index <- factor(
    as.numeric(as.character(data_sims$best_task_index)),
    levels=unique_indices)
data_sims$category <- factor(data_sims$category, levels=category_levels)
ggplot(data_sims, aes(x=expected, y=actual)) +
    geom_point(shape=21, aes(size=pval, colour=category, fill=best_task_index)) +
    geom_abline() +
    xlim(0,1.5) +
    ylim(0,1.5) +
    theme_bw() +
    theme(
        aspect.ratio=1) +
    scale_fill_manual(values=ggr_colors)
ggsave(plot_file)


data_sims <- data_sims[order(data_sims$diff),]
data_sims$diff_order <- nrow(data_sims):1
plot_file <- "fig_4.b.1.simulations_seq.TEST.pdf"
#ggplot(data_sims, aes(x=diff_order, y=diff)) +
ggplot(data_sims, aes(x=expected, y=diff)) +
    geom_hline(size=0.115, yintercept=0, linetype="dashed") +
    geom_point(shape=21, stroke=0.115, aes(size=pval, fill=best_task_index)) +
    #geom_point(shape=21, aes(fill=best_task_index)) +
    #geom_bar(stat="identity", aes(fill=best_task_index)) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
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
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +
    scale_fill_manual(values=ggr_colors) +
    scale_size_continuous(range=c(0,2)) +
    scale_x_continuous(limits=c(0, 1.5), expand=c(0,0)) +
    scale_y_continuous(limits=c(-0.2, 0.3), expand=c(0,0))
ggsave(plot_file, height=1, width=2)
data_sims$diff_order <- NULL


# plot them together 
pwm1_endo <- data_endo$pwm1
pwm2_endo <- data_endo$pwm2
data_endo$pwm2 <- pwm1_endo
data_endo$pwm1 <- pwm2_endo

data <- rbind(data_endo, data_sims)
print(dim(data))

# set up adjacency matrix
data_unmelted <- dcast(data, pwm1 ~ pwm2, value.var="actual")
data_unmelted[is.na(data_unmelted)] <- 0
rownames(data_unmelted) <- data_unmelted$pwm1
data_unmelted$pwm1 <- NULL
my_hc <- hclust(dist(data_unmelted), method="ward.D2")
pwm_order <- rownames(data_unmelted)[my_hc$order]
write.table(pwm_order, "ordering.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


# set as ordered factor
data$pwm1 <- factor(data$pwm1, levels=pwm_order)
data$pwm2 <- factor(data$pwm2, levels=rev(pwm_order))

plot_file <- "fig_4-c.adjacency.pdf"
#ggplot(data, aes(x=pwm1, y=pwm2, fill=best_task_index, colour=category)) +
ggplot(data, aes(x=pwm1, y=pwm2, fill=category, colour=best_task_index)) +
#ggplot(data, aes(x=pwm1, y=pwm2, fill=category)) +
    geom_point(shape=21, stroke=0.75, aes(size=pval)) +
    theme_bw() +
    theme(
        axis.text.x=element_text(angle=90),
        aspect.ratio=1) +
    scale_colour_manual(values=ggr_colors)
    #scale_fill_manual(values=ggr_colors)
ggsave(plot_file)
