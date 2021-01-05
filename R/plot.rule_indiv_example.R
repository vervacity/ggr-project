

library(ggplot2)
library(reshape2)
library(gridExtra)


args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]

# out file
prefix <- strsplit(data_file, ".txt")[[1]]
prefix <- prefix[1]
out_file <- paste(prefix, ".pdf", sep="")

# headers for mpra signal
days <- c("d0", "d3", "d6")
reps <- c("b1", "b2", "b3", "b4")
val_headers <- c()
for (day_i in 1:length(days)) {
    for (rep_i in 1:length(reps)) {
        val_headers <- c(val_headers, paste(days[day_i], reps[rep_i], sep="_"))
    }
}

val_headers <- c(val_headers, "unique_id", "example_fileidx")

# load data
data <- read.table(data_file, sep="\t", header=TRUE, row.names=1)
motifs <- as.character(unique(data$motifs))
data_vals <- data[,val_headers]
data_melt <- melt(data_vals, id.vars=c("unique_id", "example_fileidx"))
data_melt$day <- gsub("_.+$", "", data_melt$variable)
data_melt$rep <- gsub("^.+_", "", data_melt$variable)


# plot
p1 <- ggplot(data_melt, aes(x=day, y=value, group=unique_id)) +
    geom_line(aes(colour=unique_id), show.legend=FALSE) +
    geom_point(size=0.5, aes(colour=unique_id), show.legend=FALSE) +
    labs(title=motifs, x="Day", y="MPRA signal") +  
    theme_bw() +
    facet_grid(example_fileidx ~ rep, switch="y") +
    theme(
        strip.background=element_blank(),
        strip.placement="outside") +
    scale_x_discrete(expand=c(0.1,0.1)) #+
    #scale_y_continuous(position="right")
#ggsave(out_file, height=10, width=4)

# ATAC data
atac_headers <- c("ATAC_SIGNALS.NORM", "example_fileidx")
data_atac <- data[,atac_headers]
data_atac <- data_atac[!duplicated(data_atac),]
data_atac_vals <- do.call(
    rbind,
    strsplit(as.character(data_atac$ATAC_SIGNALS.NORM), ","))
class(data_atac_vals) <- "numeric"
keep_cols <- c(1, 2, 3, 4, 5, 6, 7, 10, 11, 13)
days <- c("d0.0", "d0.5", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0")

data_atac_vals <- data.frame(data_atac_vals[, keep_cols])
colnames(data_atac_vals) <- days
data_atac_vals$example_fileidx <- data_atac$example_fileidx
data_atac_melt <- melt(data_atac_vals, id.vars=c("example_fileidx"))
data_atac_melt$atac <- "ATAC"

p2 <- ggplot(data_atac_melt, aes(x=variable, y=value, group=example_fileidx)) +
    geom_line() +
    geom_point(size=0.5) +
    labs(title="", x="Day", y="ATAC signal") +  
    theme_bw() +
    facet_grid(example_fileidx ~ atac, switch="y") +
    theme(
        strip.background=element_blank(),
        strip.placement="outside") +
    scale_x_discrete(expand=c(0.1,0.1), labels=c("d0", "", "", "", "", "", "d3", "", "", "d6")) #+
    #scale_y_continuous(position="right")
#ggsave(out_file, height=10, width=4)

# H3K27ac data
h3k27ac_headers <- c("H3K27ac_SIGNALS.NORM", "example_fileidx")
data_h3k27ac <- data[,h3k27ac_headers]
data_h3k27ac <- data_h3k27ac[!duplicated(data_h3k27ac),]
data_h3k27ac_vals <- do.call(
    rbind,
    strsplit(as.character(data_h3k27ac$H3K27ac_SIGNALS.NORM), ","))
class(data_h3k27ac_vals) <- "numeric"
data_h3k27ac_vals <- data.frame(data_h3k27ac_vals)
colnames(data_h3k27ac_vals) <- c("d0", "d3", "d6")
data_h3k27ac_vals$example_fileidx <- data_h3k27ac$example_fileidx
data_h3k27ac_melt <- melt(data_h3k27ac_vals, id.vars=c("example_fileidx"))
data_h3k27ac_melt$h3k27ac <- "H3K27ac"

p3 <- ggplot(data_h3k27ac_melt, aes(x=variable, y=value, group=example_fileidx)) +
    geom_line() +
    geom_point(size=0.5) +
    labs(title="", x="Day", y="H3K27ac signal") +  
    theme_bw() +
    facet_grid(example_fileidx ~ h3k27ac, switch="y") +
    theme(
        strip.background=element_blank(),
        strip.placement="outside") +
    scale_x_discrete(expand=c(0.1,0.1)) #+
    #scale_y_continuous(position="right")

# NN predictions
nn_headers <- c("logits.norm", "example_fileidx")
data_nn <- data[,nn_headers]
data_nn <- data_nn[!duplicated(data_nn),]
data_nn_vals <- do.call(
    rbind,
    strsplit(as.character(data_nn$logits.norm), ","))
class(data_nn_vals) <- "numeric"
data_nn_vals <- data_nn_vals[,keep_cols]
data_nn_vals <- data.frame(data_nn_vals)
colnames(data_nn_vals) <- days
data_nn_vals$example_fileidx <- data_nn$example_fileidx

data_nn_melt <- melt(data_nn_vals, id.vars=c("example_fileidx"))
data_nn_melt$nn <- "NN predictions"

p4 <- ggplot(data_nn_melt, aes(x=variable, y=value, group=example_fileidx)) +
    geom_line() +
    geom_point(size=0.5) +
    labs(title="", x="Day", y="Logits") +  
    theme_bw() +
    facet_grid(example_fileidx ~ nn, switch="y") +
    theme(
        strip.background=element_blank(),
        strip.placement="outside") +
    scale_x_discrete(expand=c(0.1,0.1), labels=c("d0", "", "", "", "", "", "d3", "", "", "d6")) #+
    #scale_y_continuous(position="right")

# merge all
pdf(out_file, height=10, width=9)
grid.arrange(p1, p2, p3, p4, ncol=4, widths=c(3, 2, 1.5, 2))
dev.off()




