#!/usr/bin/env Rscript


library(ggplot2)
library(RColorBrewer)
library(scales)


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
summary_file <- args[1]
prefix <- args[2]

# read in data
data <- read.table(gzfile(summary_file), sep="\t", header=TRUE)

# remove unused
data <- data[grepl("ggr", data$unique_id),]
#data <- data[grepl("combo-00", data$unique_id),]
print(head(data))

# average barcodes and THEN collapse replicates

# now collapse down on barcodes
keep_cols <- c(
    "d0_b1",
    "d0_b2",
    "d0_b3",
    "d0_b4",
    "d3_b1",
    "d3_b2",
    "d3_b3",
    "d3_b4",
    "d6_b1",
    "d6_b2",
    "d6_b3",
    "d6_b4")

d0_reps <- c(
    "d0_b1",
    "d0_b2",
    "d0_b3",
    "d0_b4")

d3_reps <- c(
    "d3_b1",
    "d3_b2",
    "d3_b3",
    "d3_b4")

d6_reps <- c(
    "d6_b1",
    "d6_b2",
    "d6_b3",
    "d6_b4")

# remove rows of zeros and set remaining zeros to NA
data_max <- apply(data[,keep_cols], 1, max)
data <- data[data_max > 0,]
data[data == 0] <- NA

# ===================================
# ATAC
# ===================================

data_atac <- strsplit(as.character(data$ATAC_SIGNALS.NORM), ",", fixed=TRUE)
data_atac <- do.call(rbind, data_atac)
data_atac <- apply(data_atac, 2, as.numeric)

# atac days
data$atac_d0 <- data_atac[,1]
data$atac_d3 <- data_atac[,7]
data$atac_d6 <- data_atac[,13]
atac_cols <- c(keep_cols, c("atac_d0", "atac_d3", "atac_d6"))

# average barcodes
data_agg <- aggregate(
    data[,atac_cols],
    by=list(id=data$example_combo_id),
    FUN=mean,
    na.rm=TRUE)


# average reps
data_agg$mpra_d0 <- apply(data_agg[, d0_reps], 1, mean, na.rm=TRUE)
data_agg$mpra_d3 <- apply(data_agg[, d3_reps], 1, mean, na.rm=TRUE)
data_agg$mpra_d6 <- apply(data_agg[, d6_reps], 1, mean, na.rm=TRUE)


cor_val <- cor(data_agg$atac_d0, data_agg$mpra_d0, use="complete.obs")
print(cor_val)

# plot
data_agg <- data_agg[data_agg$atac_d0 > 0,]
data_agg <- data_agg[data_agg$mpra_d0 > 0,]

ggplot(data_agg, aes(x=atac_d0, y=mpra_d0)) +
    geom_hex(bins=100) +
    geom_smooth(method="lm") +
    #geom_point(alpha=0.5) +
    theme_bw() +
    my.scale_fill_distiller(palette="Blues", direction=1)

ggsave("test_atac.pdf")


# ===================================
# H3K27ac
# ===================================

data_h3k27ac <- strsplit(as.character(data$H3K27ac_SIGNALS.NORM), ",", fixed=TRUE)
data_h3k27ac <- do.call(rbind, data_h3k27ac)
data_h3k27ac <- apply(data_h3k27ac, 2, as.numeric)

# atac days
data$h3k27ac_d0 <- data_h3k27ac[,1]
data$h3k27ac_d3 <- data_h3k27ac[,2]
data$h3k27ac_d6 <- data_h3k27ac[,3]
h3k27ac_cols <- c(keep_cols, c("h3k27ac_d0", "h3k27ac_d3", "h3k27ac_d6"))

# average barcodes
data_agg <- aggregate(
    data[,h3k27ac_cols],
    by=list(id=data$example_combo_id),
    FUN=mean,
    na.rm=TRUE)

# average reps
data_agg$mpra_d0 <- apply(data_agg[, d0_reps], 1, mean, na.rm=TRUE)
data_agg$mpra_d3 <- apply(data_agg[, d3_reps], 1, mean, na.rm=TRUE)
data_agg$mpra_d6 <- apply(data_agg[, d6_reps], 1, mean, na.rm=TRUE)

print(head(data_agg))

# plot
data_agg <- data_agg[data_agg$h3k27ac_d0 > 0,]
data_agg <- data_agg[data_agg$mpra_d0 > 0,]

ggplot(data_agg, aes(x=h3k27ac_d0, y=mpra_d0)) +
    geom_hex(bins=100) +
    geom_smooth(method="lm") +
    #geom_point(alpha=0.5) +
    theme_bw() +
    my.scale_fill_distiller(palette="Blues", direction=1)

ggsave("test_h3k27ac.pdf")

# logits
data_logits <- strsplit(as.character(data$logits.motif_mut), ",", fixed=TRUE)
data_logits <- do.call(rbind, data_logits)
data_logits <- apply(data_logits, 2, as.numeric)

print(head(data_logits))

# atac days
data$logits_d0 <- data_logits[,1]
data$logits_d3 <- data_logits[,7]
data$logits_d6 <- data_logits[,13]
logits_cols <- c(keep_cols, c("logits_d0", "logits_d3", "logits_d6"))

# average barcodes
data_agg <- aggregate(
    data[,logits_cols],
    by=list(id=data$example_combo_id),
    FUN=mean,
    na.rm=TRUE)

# average reps
data_agg$mpra_d0 <- apply(data_agg[, d0_reps], 1, mean, na.rm=TRUE)
data_agg$mpra_d3 <- apply(data_agg[, d3_reps], 1, mean, na.rm=TRUE)
data_agg$mpra_d6 <- apply(data_agg[, d6_reps], 1, mean, na.rm=TRUE)

print(head(data_agg))
print(dim(data_agg))

# plot
data_agg <- data_agg[data_agg$logits_d0 > 0,]
data_agg <- data_agg[data_agg$mpra_d0 > 0,]

ggplot(data_agg, aes(x=logits_d0, y=mpra_d0)) +
    geom_hex(bins=100) +
    geom_smooth(method="lm") +
    #geom_point(alpha=0.5) +
    theme_bw() +
    my.scale_fill_distiller(palette="Blues", direction=1)

ggsave("test_logits.pdf")


# plot mpra reps
ggplot(data_agg, aes(x=d0_b1, y=d0_b2)) + geom_point(alpha=0.5) + theme_bw()
ggsave("test_mpra_reps.pdf")


