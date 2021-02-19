#!/usr/bin/env Rscript

library(VennDiagram)
library(grid)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(data_file, sep="\t", header=TRUE, row.names=1)
data <- data * (as.numeric(rownames(data))+1)
data <- as.list(data)

# remove nulls
for (l in names(data)) {
    elems <- data[[l]]
    elems <- elems[elems != 0]
    data[[l]] <- elems
}

# plot
temp <- venn.diagram(
    data,
    fill=c("brown1", "cornflowerblue"),
    main="Enriched motif identification:\nNeural Net Motifs vs HOMER+ATAC",
    main.fontfamily="ArialMT",
    main.cex=0.5,
    alpha=c(0.5, 0.5),
    cex=1,
    lty=1,
    lwd=1,
    fontfamily="ArialMT",
    cat.fontfamily="ArialMT",
    cat.cex=0.5,
    filename=NULL)

# save out
pdf(file=plot_file, height=2, width=2)
    grid.draw(temp)
dev.off()
