#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

# args
args <- commandArgs(trailingOnly=TRUE)
multiplicity_file <- args[1]
spacing_file <- args[2]
plot_file <- args[3]

count_thresholds <- c(1, 2, 3, 4)
ranges <- c("5-15", "15-50", "50-101")


# read in multiplicity, clean up
if (file.exists(multiplicity_file)) {
    multiplicity <- read.table(gzfile(multiplicity_file), sep="\t", header=TRUE)
    multiplicity$range <- NA
    #multiplicity$count_threshold <- factor(multiplicity$count_threshold)
    multiplicity$category <- "multiplicity"
    multiplicity$show <- 1
    count_thresholds <- unique(multiplicity$count_threshold)
    multiplicity$count_threshold <- as.character(multiplicity$count_threshold)
} else {
    multiplicity <- data.frame(
        term.id=character(),
        term.name=character(),
        range=double(),
        count_threshold=character(),
        log10pval=double(),
        category=character(),
        stringsAsFactors=FALSE)
    #multiplicity$count_threshold <- factor(
    #    multiplicity$count_threshold, levels=c(1,2))
}
multiplicity$x <- multiplicity$count_threshold
multiplicity <- multiplicity[,order(colnames(multiplicity))]


# read in spacing, clean up
if (file.exists(spacing_file)) {
    # read in file, adjust ranges
    spacing <- read.table(gzfile(spacing_file), sep="\t", header=TRUE)
    spacing$count_threshold <- NA
    spacing$category <- "spacing"
    spacing$show <- 1
    ranges <- unique(spacing$range)
} else {
    # make empty frame
    spacing <- data.frame(
        term.id=character(),
        term.name=character(),
        range=double(),
        count_threshold=double(),
        log10pval=double(),
        category=character(),
        show=double(),
        stringsAsFactors=FALSE)
    # add in placeholder data (to hold open ggplot facet)

}


    

spacing$range <- factor(spacing$range, levels=ranges)
spacing$x <- spacing$range
spacing <- spacing[,order(colnames(spacing))]

# bind and factorize
data <- rbind(multiplicity, spacing)
#levels <- c(count_thresholds, ranges)
#data$x <- factor(data$x, levels=levels)
data$category <- factor(data$category, levels=c("multiplicity", "spacing"))

# placeholder data as needed for plotting
all_term_names <- unique(data$term.name)

# TODO same for multiplicity

if (nrow(data[data$category == "spacing",]) == 0) {
    
    for (range_idx in 1:length(ranges)) {
        for (term_idx in 1:length(all_term_names)) {
            spacing_row <- data.frame(
                category="spacing",
                count_threshold=NA,
                log10pval=NA,
                range=ranges[range_idx],
                term.id=NA,
                term.name=all_term_names[term_idx],
                x=ranges[range_idx],
                show=0,
                stringsAsFactors=FALSE)
            data <- rbind(data, spacing_row)
        }
    }
}


# plot
ggplot(data, aes(y=term.name)) +
    geom_point(
        shape=21,
        fill="white",
        aes(x=count_threshold, size=log10pval)) +
    geom_point(
        shape=21,
        fill="white",
        aes(x=range, size=log10pval)) +
    facet_grid(. ~ category, scales="free_x", space="free_x") +
    theme_bw()        
ggsave("test.pdf", width=6, height=6)
