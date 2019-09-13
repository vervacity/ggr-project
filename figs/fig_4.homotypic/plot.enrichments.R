#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

# args
args <- commandArgs(trailingOnly=TRUE)
multiplicity_file <- args[1]
spacing_file <- args[2]
plot_file <- args[3]

# read in (check if file exists)
if (file.exists(multiplicity_file)) {
    multiplicity <- read.table(gzfile(multiplicity_file), sep="\t", header=TRUE)
} else {
    multiplicity <- data.frame(
        term.id=character(),
        term.name=character(),
        count_threshold=double(),
        log10pval=double(),
        stringsAsFactors=FALSE)
}

multiplicity$count_threshold <- factor(multiplicity$count_threshold)

print(head(multiplicity))

if (file.exists(spacing_file)) {
    spacing <- read.table(gzfile(spacing_file), sep="\t", header=TRUE)
} else {
    spacing <- data.frame(
        term.id=character(),
        term.name=character(),
        range=double(),
        log10pval=double(),
        stringsAsFactors=FALSE)
}

# reconcile term names
term_name_union = union(multiplicity$term.name, spacing$term.name)
multiplicity$term.name <- factor(multiplicity$term.name, levels=term_name_union)
spacing$term.name <- factor(spacing$term.name, levels=term_name_union)

# and plot both in grid
p1 <- ggplot(multiplicity, aes(x=count_threshold, y=term.name)) +
    geom_point(aes(size=log10pval))

p2 <- ggplot(spacing, aes(x=range, y=term.name)) +
    geom_point(aes(size=log10pval))

ggsave("test.pdf", width=14, height=7, arrangeGrob(p1, p2, nrow=1, ncol=2))


# check terms







