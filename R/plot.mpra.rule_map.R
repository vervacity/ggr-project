#!/usr/bin/env Rscript


library(ggplot2)
library(RColorBrewer)
library(reshape2)

library(scales)

library(gridGraphics)
library(gridExtra)

# betterr adjust scale fill distiller
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
summary_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(gzfile(summary_file), sep="\t", header=TRUE)

# remove failed
data <- data[!grepl("FAILED", data$interaction),]

# adjust levels
interaction_levels <- rev(c("SYNERGY", "ADDITIVE", "BUFFER"))
data$interaction <- factor(as.character(data$interaction), levels=interaction_levels)
data$day_empirical <- factor(as.character(data$day_empirical), levels=c("d6_AVG", "d3_AVG", "d0_AVG"))

# order by: day, interaction, pwm
data <- data[with(data, order(day_empirical, interaction, pwm1)),]

# set row order (make factor)
data$grammar <- factor(as.character(data$grammar), levels=data$grammar)

# interaction results bar
interaction_colors <- brewer.pal(8, "Dark2")
interaction_colors <- c(interaction_colors[5], interaction_colors[6], interaction_colors[4]) #rev(interaction_colors[4:6])

p1 <- ggplot(data, aes(x=0, y=grammar, fill=interaction)) +
    geom_tile(colour="black", show.legend=FALSE) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,40,5),
        plot.title=element_text(size=8, margin=margin(b=1)),
        
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
            
        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
    scale_fill_manual(values=interaction_colors, drop=FALSE) 


# expected vs actual
keep_columns <- c("grammar", "expected", "actual")
expected_v_actual <- data[keep_columns]
expected_v_actual <- melt(expected_v_actual, id.vars="grammar")
limit <- max(abs(expected_v_actual$value))
p2 <- ggplot(expected_v_actual, aes(x=variable, y=grammar, fill=value)) +
    geom_tile(colour="black", show.legend=FALSE) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,10,5),
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
p3 <- ggplot(pwm_data, aes(x=pwm, y=grammar, fill=log2FC)) +
    geom_line(aes(group=grammar), colour="black") +
    geom_point(shape=21, color="black", aes(fill=log2FC), show.legend=FALSE) +
    theme_bw() +    
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,10,5),
        plot.title=element_text(size=8, margin=margin(b=1)),
            
        panel.background=element_blank(),
        panel.border=element_rect(size=0.115),
        panel.grid.major.y=element_blank(),
        #panel.grid.major=element_blank(),
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
limit <- max(abs(data_traj$value))
p4 <- ggplot(data_traj, aes(x=variable, y=grammar, fill=value)) +
    geom_tile(colour="black", show.legend=FALSE) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,40,5),
        plot.title=element_text(size=8, margin=margin(b=1)),
        
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),

        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) + 
    my.scale_fill_distiller(palette="RdBu", direction=-1, limits=c(-limit, limit)) 
       


# plot all together
pdf(plot_file, height=4, width=7, family="ArialMT", useDingbats=FALSE)
grid.arrange(p1, p2, p3, p4, ncol=4, widths=c(0.4,0.4,5,0.75))
dev.off()

