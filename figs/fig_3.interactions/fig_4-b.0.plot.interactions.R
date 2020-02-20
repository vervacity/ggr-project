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
category_levels <- c("additive", "buffer", "synergy")

# args
args <- commandArgs(trailingOnly=TRUE)
endo_file <- args[1]
sims_file <- args[2]

# get palette
colors <- brewer.pal(8, "Set2")
#colors <- brewer.pal(8, "Dark2")
colors <- rev(colors[4:6])

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
    geom_abline(size=0.115, linetype="dashed") +
    geom_point(shape=21, stroke=0.115, aes(fill=category)) +
    labs(x="Linear expectation", y="Observed", title="Interaction scores") +
    theme_bw() +
    theme(
        aspect.ratio=1,
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +
    scale_fill_manual(values=colors, drop=FALSE) +
    scale_x_continuous(limits=c(0.5, 1.3), expand=c(0,0)) +
    scale_y_continuous(limits=c(0.5, 1.3), expand=c(0,0))

ggsave(plot_file, height=2, width=2, units="in", useDingbats=FALSE)


# different plot for endogenous - rank order?
data_endo <- data_endo[order(data_endo$diff),]
data_endo$diff_order <- nrow(data_endo):1
plot_file <- "fig_4.b.0.endogenous_seq.DIFF.pdf"
#ggplot(data_endo, aes(x=diff_order, y=diff)) +
ggplot(data_endo, aes(x=expected, y=diff)) +
    geom_hline(size=0.115, yintercept=0, linetype="dashed") +
    geom_point(shape=21, stroke=0.115, aes(fill=category)) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +
    scale_fill_manual(values=colors, drop=FALSE) +
    scale_size_continuous(range=c(0,2)) +
    scale_x_continuous(limits=c(0.5, 1.3), expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

ggsave(plot_file, height=1, width=2, units="in", useDingbats=FALSE)
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
    geom_abline(size=0.115, linetype="dashed") +
    geom_point(shape=21, stroke=0.115, aes(fill=category)) +
    labs(x="Linear expectation", y="Observed", title="Interaction scores") +
    theme_bw() +
    theme(
        aspect.ratio=1,
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +
    scale_fill_manual(values=colors, drop=FALSE) +
    scale_x_continuous(limits=c(0, 1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, 1), expand=c(0,0))

ggsave(plot_file, height=2, width=2, units="in", useDingbats=FALSE)

# plot differently
data_sims <- data_sims[order(data_sims$diff),]
data_sims$diff_order <- nrow(data_sims):1
plot_file <- "fig_4.b.1.simulations_seq.DIFF.pdf"
#ggplot(data_sims, aes(x=diff_order, y=diff)) +
ggplot(data_sims, aes(x=expected, y=diff)) +
    geom_hline(size=0.115, yintercept=0, linetype="dashed") +
    geom_point(shape=21, stroke=0.115, aes(fill=category)) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +
    scale_fill_manual(values=colors, drop=FALSE) +
    scale_size_continuous(range=c(0,2)) +
    scale_x_continuous(limits=c(0, 1), expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

ggsave(plot_file, height=1, width=2, units="in", useDingbats=FALSE)
data_sims$diff_order <- NULL

# --------------------
# adjacency matrix
# --------------------

# build a joint dataframe
pwm1_endo <- data_endo$pwm1
pwm2_endo <- data_endo$pwm2
data_endo$pwm2 <- pwm1_endo
data_endo$pwm1 <- pwm2_endo
data_endo$seq_type <- "endogenous"
data_sims$seq_type <- "simulation"
data <- rbind(data_endo, data_sims)
print(dim(data))

# get an ordering
data_unmelted <- dcast(data, pwm1 ~ pwm2, value.var="actual")
data_unmelted[is.na(data_unmelted)] <- 0
rownames(data_unmelted) <- data_unmelted$pwm1
data_unmelted$pwm1 <- NULL
my_hc <- hclust(dist(data_unmelted), method="ward.D2")
pwm_order <- rownames(data_unmelted)[my_hc$order]
#write.table(pwm_order, "ordering.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# then, for each row, check which comes first in ordering
for (row_i in 1:nrow(data)) {
    pwm1 <- data$pwm1[row_i]
    pwm1_pos <- match(pwm1, pwm_order)
    #print(pwm1)
    #print(pwm1_pos)
    
    pwm2 <- data$pwm2[row_i]
    pwm2_pos <- match(pwm2, pwm_order)
    #print(pwm2)
    #print(pwm2_pos)

    # check order
    if (data$seq_type[row_i] == "endogenous") {
        # make sure pwm1_pos > pwm2_pos
        if (pwm2_pos < pwm1_pos) {
            data$pwm1[row_i] <- pwm2
            data$pwm2[row_i] <- pwm1
        }
    } else {
        # make sure pwm1_pos < pwm2_pos
        if (pwm2_pos > pwm1_pos) {
            data$pwm1[row_i] <- pwm2
            data$pwm2[row_i] <- pwm1
        }
    }
}

# set as ordered factor
data$pwm1 <- factor(data$pwm1, levels=pwm_order)
data$pwm2 <- factor(data$pwm2, levels=rev(pwm_order))

# thresh pval
thresh <- quantile(data$pval, 0.80)
print(thresh)
data$pval[data$pval > thresh] <- thresh

# plot abs val of diff
data$diff_abs <- abs(data$diff)

plot_file <- "fig_4-c.adjacency.pdf"
#ggplot(data, aes(x=pwm1, y=pwm2, fill=best_task_index, colour=category)) +
#ggplot(data, aes(x=pwm1, y=pwm2, fill=category, colour=best_task_index)) +
ggplot(data, aes(x=pwm1, y=pwm2, fill=category)) +
    #geom_point(shape=21, stroke=0.115, aes(size=pval), show.legend=FALSE) +
    geom_point(shape=21, stroke=0.115, aes(size=diff_abs)) +
    labs(x="Motif", y="Motif", title="Significant pairwise interactions") +
    theme_bw() +
    theme(
        aspect.ratio=1,
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=0)),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.115),
        panel.grid=element_line(size=0.115),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5, angle=90, vjust=0.5),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5)) +
    scale_colour_manual(values=ggr_colors) +
    scale_size_continuous(range=c(0.3,3)) +
    scale_fill_manual(values=colors, drop=FALSE)

ggsave(plot_file, height=3.25, width=3.25, units="in", useDingbats=FALSE)