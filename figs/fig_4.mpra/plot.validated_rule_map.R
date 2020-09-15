#!/usr/bin/env Rscript


library(ggplot2)
library(RColorBrewer)
library(reshape2)

library(scales)

library(gridGraphics)
library(gridExtra)

# better adjust scale fill distiller
my.scale_fill_distiller <- function(
    ...,
    type = "seq",
    palette = 1,
    direction = -1,
    values = NULL,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") {
    # warn about using a qualitative brewer palette to generate the gradient
    type <- match.arg(type, c("seq", "div", "qual"))
    if (type == "qual") {
        warning("Using a discrete colour palette in a continuous scale.\n  Consider using type = \"seq\" or type = \"div\" instead", call. = FALSE)
    }
    continuous_scale(
        aesthetics,
        "distiller",
        gradient_n_pal(brewer_pal(type, palette, direction)(7), values, space),
        na.value = na.value, guide = guide, ...)
    # NB: 6 colours per palette gives nice gradients; more results in more saturated colours which do not look as good
}


# args
args <- commandArgs(trailingOnly=TRUE)
mpra_summary_file <- args[1]
predicted_summary_file <- args[2] # NOTE: requires that all other predicted files are ordered same as this one
atac_mat_file <- args[3]
rna_mat_file <- args[4]
go_terms_file <- args[5]
plot_file <- args[6]

# params
bottom_adjust <- 85
motif_bottom_adjust <- 70
interact_bottom_adjust <- 60
exp_v_obs_bottom_adjust <- 75
tile_height <- 1 # 0.75

# read in data
data <- read.table(gzfile(mpra_summary_file), sep="\t", header=TRUE)

# remove failed
data <- data[!grepl("FAILED", data$interaction),]

# adjust levels
interaction_levels <- rev(c("SYNERGY", "ADDITIVE", "BUFFER"))
data$interaction <- factor(as.character(data$interaction), levels=interaction_levels)
#data$day_empirical <- factor(as.character(data$day_empirical), levels=c("d6_AVG", "d3_AVG", "d0_AVG"))
data$day <- factor(as.character(data$day), levels=c("d6", "d3", "d0"))

# order by: day, interaction, pwm
#data <- data[with(data, order(day_empirical, interaction, pwm1)),]
data <- data[with(data, order(day, interaction, pwm1)),]

# set row order (make factor)
data$grammar <- factor(as.character(data$grammar), levels=data$grammar)

# interaction results bar
interaction_colors <- brewer.pal(8, "Dark2")
interaction_colors <- c(interaction_colors[5], interaction_colors[6], interaction_colors[4])

# plot interaction results
p_interaction_type <- ggplot(data, aes(x="interaction type", y=grammar, fill=interaction)) +
    geom_tile(colour="black", height=tile_height, show.legend=FALSE) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(t=4,r=1,b=interact_bottom_adjust,l=2),
        plot.title=element_text(size=8, margin=margin(b=1)),
        
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
            
        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=6, angle=60, hjust=1),
        axis.ticks=element_blank()) +
    scale_fill_manual(values=interaction_colors, drop=FALSE) 


# expected vs actual
keep_columns <- c("grammar", "expected", "actual")
expected_v_actual <- data[keep_columns]
expected_v_actual <- melt(expected_v_actual, id.vars="grammar")
limit <- max(abs(expected_v_actual$value))
p_expected_v_observed <- ggplot(expected_v_actual, aes(x=variable, y=grammar, fill=value)) +
    geom_tile(colour="black", height=tile_height, show.legend=FALSE) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(t=4,r=1,b=exp_v_obs_bottom_adjust,l=1),
        plot.title=element_text(size=8, margin=margin(b=1)),
        
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
            
        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=6, angle=60, hjust=1),
        axis.ticks=element_blank()) +
    my.scale_fill_distiller(palette="RdBu", direction=-1, limits=c(-limit, limit))


# ----------------
# pwm map
# ----------------

# then reduce
keep_columns <- c("grammar", "pwm1_clean", "X0.1", "day_empirical")
pwm1_data  <- data[keep_columns]
colnames(pwm1_data) <- c("grammar", "pwm", "log2FC", "day")

keep_columns <- c("grammar", "pwm2_clean", "X1.0", "day_empirical")
pwm2_data <- data[keep_columns]
colnames(pwm2_data) <- c("grammar", "pwm", "log2FC", "day")

pwm_data <- rbind(pwm1_data, pwm2_data)
limit <- max(abs(pwm_data$log2FC))

# get ordering by weighting by day
pwm_data$day <- gsub("d", "", pwm_data$day)
pwm_data$day <- as.numeric(gsub("_AVG", "", pwm_data$day))
pwm_ordering <- aggregate(pwm_data$day, by=list(pwm_data$pwm), mean)
pwm_ordering <- pwm_ordering[order(pwm_ordering$x),]
pwm_order <- pwm_ordering$Group.1

pwm_data$pwm <- factor(as.character(pwm_data$pwm), levels=pwm_order)

# plot
p_motifs <- ggplot(pwm_data, aes(x=pwm, y=grammar, fill=log2FC)) +
    geom_line(size=0.230, aes(group=grammar), colour="black") +
    #geom_point(shape=21, size=1, stroke=0.230, color="black", aes(fill=log2FC), show.legend=FALSE) +
    geom_point(shape=21, size=1.3, stroke=0.345, color="black", fill="white", show.legend=FALSE) +
    theme_bw() +    
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(t=4,r=1,b=motif_bottom_adjust,l=1),
        plot.title=element_text(size=8, margin=margin(b=1)),
            
        panel.background=element_blank(),
        panel.border=element_rect(size=0.115),
        panel.grid.major=element_line(size=0.115),
        panel.grid.minor=element_blank(),
            
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_blank(), #element_text(vjust=1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text.x=element_text(size=6, angle=60, hjust=1),
        axis.text.y=element_blank(),
        axis.ticks.x=element_line(size=0.115),
        axis.ticks.y=element_blank(),
        axis.ticks.length=unit(0.01, "in")) +
    my.scale_fill_distiller(palette="RdBu", direction=-1, limits=c(-limit, limit))

# trajectory pattern
keep_columns <- c("grammar", "traj.medians.d0", "traj.medians.d3", "traj.medians.d6")
data_traj <- data[keep_columns]
data_traj <- melt(data_traj, id.vars="grammar")
data_traj$variable <- gsub("traj.medians.", "", data_traj$variable, fixed=TRUE)

# thresh
thresh <- quantile(abs(data_traj$value), 0.90)
data_traj$value[data_traj$value > thresh] <- thresh
data_traj$value[data_traj$value < -thresh] <- -thresh

limit <- max(abs(data_traj$value))
p_mpra_traj <- ggplot(data_traj, aes(x=variable, y=grammar, fill=value)) +
    geom_tile(colour="black", height=tile_height) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=1)),
        plot.margin=margin(t=4,r=5,b=bottom_adjust,l=1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),

        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5)) +
    my.scale_fill_distiller(palette="RdBu", direction=-1, limits=c(-limit, limit)) +
    guides(fill=guide_colorbar(ticks=FALSE))

# ==========================
# MERGE WITH PREDICTED INFO
# ==========================

# pull in summary so know how to filter
pred_summary <- read.table(predicted_summary_file, sep="\t", header=TRUE, row.names=1)
pred_summary$grammar <- gsub("^.+/", "", pred_summary$filename)
pred_summary$grammar <- gsub(".gml", "", pred_summary$grammar)

# get ordering
ordering <- match(data$grammar, pred_summary$grammar)

# ATAC
atac_data <- read.table(atac_mat_file, header=TRUE, sep="\t", row.names=1)
atac_data <- atac_data[ordering,]
atac_data <- atac_data - atac_data[,1]
atac_data$ids <- rownames(atac_data)
atac_melt <- melt(atac_data, id.vars="ids")
atac_melt$ids <- factor(atac_melt$ids, levels=rownames(atac_data))
limit <- max(abs(atac_melt$value))
p_atac <- ggplot(atac_melt, aes(x=variable, y=ids, fill=value)) +
    geom_tile(colour="black", height=tile_height) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=1)),
        plot.margin=margin(t=4,r=1,bottom_adjust,l=1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),

        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5)) +
        #legend.spacing.x=unit(0.05, "in")) +
    my.scale_fill_distiller(palette="RdBu", direction=-1, limits=c(-limit, limit)) +
    guides(fill=guide_colorbar(ticks=FALSE))

# RNA
rna_data <- read.table(rna_mat_file, header=TRUE, sep="\t", row.names=1)
rna_data <- rna_data[ordering,]
rna_data <- rna_data - rna_data[,1]
rna_data$ids <- rownames(rna_data)
rna_melt <- melt(rna_data, id.vars="ids")
rna_melt$ids <- factor(rna_melt$ids, levels=rownames(rna_data))
limit <- max(abs(rna_melt$value))
p_rna <- ggplot(rna_melt, aes(x=variable, y=ids, fill=value)) +
    geom_tile(colour="black", height=tile_height) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=1)),
        plot.margin=margin(t=4,r=1,bottom_adjust,l=1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),

        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5)) +
    my.scale_fill_distiller(palette="RdBu", direction=-1, limits=c(-limit, limit)) +
    guides(fill=guide_colorbar(ticks=FALSE))


# GO term map, reorder, plot
go_data <- read.table(go_terms_file, header=TRUE, sep="\t", row.names=1)
go_data <- go_data[ordering,]
print(dim(go_data[apply(go_data, 1, sum) != 0,]))
print(dim(go_data))
go_data$ids <- rownames(go_data)
go_data_melted <- melt(go_data, id.vars="ids")
go_data_melted <- go_data_melted[go_data_melted$value != 0,]
go_data_melted$ids <- factor(go_data_melted$ids, levels=rownames(go_data))
go_term_levels <- gsub(".", " ", colnames(go_data), fixed=TRUE)
go_terms <- gsub(".", " ", go_data_melted$variable, fixed=TRUE)
go_data_melted$variable <- factor(go_terms, levels=go_term_levels)

p_go <- ggplot(go_data_melted, aes(x=variable, y=ids)) +
    #geom_point(shape=21, fill="white", aes(size=value)) +
    geom_point(shape=20, aes(size=value)) +
    labs(size="-log10(p)", y=NULL, x=NULL) +
    theme_bw() +
    theme(
        text=element_text("ArialMT"),
        plot.margin=margin(t=4, r=0, l=1),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.115),
        panel.grid.major=element_line(size=0.115),
        axis.text.x=element_text(size=6, angle=60, hjust=1),
        axis.text.y=element_blank(),
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
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        legend.spacing.x=unit(0.05, "in")) +
    scale_size_continuous(range=c(0.25,3)) +
    scale_y_discrete(drop=FALSE)

# plot all together
pdf(plot_file, height=4.25, width=7.25, family="ArialMT", useDingbats=FALSE)
grid.arrange(
    p_atac,
    p_motifs,
    p_rna,
    p_interaction_type,
    p_expected_v_observed,
    p_mpra_traj,
    p_go,
    ncol=7,
    widths=c(0.45, 3.1, 0.45, 0.15, 0.2, 0.3, 2.2))
dev.off()

