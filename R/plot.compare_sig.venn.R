#!/usr/bin/env Rscript

library(VennDiagram)
library(grid)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]


data <- read.table(data_file, sep="\t", header=TRUE, row.names=1)
data <- data * (as.numeric(rownames(data))+1)
data <- as.list(data)

for (l in names(data)) {
    elems <- data[[l]]
    elems <- elems[elems != 0]
    data[[l]] <- elems
}

# plot
temp <- venn.diagram(
    data,
    fill = c("red", "green"),
    alpha = c(0.5, 0.5),
    cex = 2,
    #cat.fontface = 4,
    lty =2,
    #fontfamily =3,
    filename = NULL)

# save out
pdf(file=plot_file)
    grid.draw(temp)
dev.off()
