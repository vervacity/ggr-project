#!/usr/bin/env Rscript

library(VennDiagram)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(grid)

library(DESeq2)


my_hclust <- function(data) {
    hc <- hclust(data, method="ward.D2")
    return(hc)
}

plot_heatmap <- function(data) {

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.5,1,0.5)
    mylhei = c(0.09,1,0.7)

    # color
    my_palette <- colorRampPalette(brewer.pal(7, "BuGn"))(49)
    #my_breaks <- seq(min(data), max(data), length=50)
    my_breaks <- seq(min(data), 0.9, length=50)
    
    labRow <- gsub("_b", ", rep ", rownames(data))
    labRow <- gsub("d", "day ", labRow)
    labCol <- gsub("_b", ", rep ", colnames(data))
    labCol <- gsub("d", "day ", labCol)

    
    # heatmap
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        hclustfun=my_hclust,
        dendrogram="none",
        trace="none",
        density.info="none",

        labRow=labRow,
        labCol=labCol,
        srtCol=30,
        
        colsep=0:(ncol(data)+1),
        rowsep=0:(nrow(data)+1),
        sepcolor="white",
        #sepwidth=c(0.001,0.001),
        sepwidth=c(0.001,0.5),

        cexRow=0.5,
        cexCol=0.5,
        
        offsetRow=0,
        offsetCol=0,
        
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            #mar=c(4,6,1.4,6),
            mar=c(1.5,1,3,1),
            mgp=c(0,-0.1,0),
            tcl=-0.1,
            lend=2,
            cex.axis=0.6,
            bty="n"),
        key.xtickfun=function() {
            breaks <- pretty(parent.frame()$breaks)
            breaks <- breaks[c(1,length(breaks))]
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)},

        margins=c(0,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        na.color="white",
        col=my_palette,
        breaks=my_breaks)

}


args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
out_dir <- args[2]

# read in data
data <- read.table(data_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
days <- c("d0", "d3", "d6")
sc_days <- c("sc.d0", "sc.d3", "sc.d6")
bulk_days <- c("bulk.d0", "bulk.d3", "bulk.d6")

# comparison plot of scATAC vs bulk ATAC across timepoints
data_cor <- data[!is.na(data$mapped),]
data_cor <- data_cor[data_cor$mapped == 1.0,]

# extra filtering
data_cor <- data_cor[c(sc_days, bulk_days)]
data_rlog <- data.frame(rlog(as.matrix(data_cor)))


data_rlog[sc_days][data_rlog[sc_days] < 5] <- 3
data_rlog <- data_rlog[apply(data_rlog[sc_days] > 5, 1, any),]

data_cor <- cor(data_rlog, method="pearson", use="complete.obs")
data_cor <- data_cor[sc_days, bulk_days]
print(data_cor)

plot_file <- paste(out_dir, "/sc_vs_bulk.corr.pdf", sep="")
pdf(plot_file, height=2, width=2, useDingbats=FALSE)
plot_heatmap(data_cor)
dev.off()

quit()

for (day_idx in 1:length(days)) {

    day <- days[day_idx]

    # pull out data for scatter plot
    data_scatter <- data[!is.na(data$mapped),]
    data_scatter <- data_scatter[data_scatter$mapped == 1.0,]

    # test
    data_scatter <- data_scatter[c(sc_days, bulk_days)]
    data_rlog <- data.frame(rlog(as.matrix(data_scatter)))
    
    sc_header <- paste("sc.", day, sep="")
    bulk_header <- paste("bulk.", day, sep="")
    headers <- c(sc_header, bulk_header)
    data_rlog <- data_rlog[headers]
    colnames(data_rlog) <- c("sc", "bulk")

    data_scatter <- data_scatter[headers]
    colnames(data_scatter) <- c("sc", "bulk")
    
    data_rlog <- data_rlog[data_scatter$sc > 0 & data_scatter$bulk > 0,]
    data_scatter <- data_rlog
    
    # log transform
    #data_scatter <- data.frame(rlog(as.matrix(data_scatter)))
    
    # plot scatter
    ggplot(data_scatter, aes(x=sc, y=bulk)) +
        #geom_point(alpha=0.1) +
        geom_hex() +
        theme_bw() +
        theme(
            aspect.ratio=1) +
        scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") +
        scale_x_continuous(limits=c(0,11), expand=c(0,0)) +
        scale_y_continuous(limits=c(0,11), expand=c(0,0))
    
    plot_file <- paste(out_dir, "/sc_vs_bulk.", day, ".scatter.pdf", sep="")
    ggsave(plot_file)


    # and also make venn diagrams
    headers <- c(sc_header, bulk_header, "sc_atac_ids")
    data_venn <- data[headers]
    colnames(data_venn) <- c("sc", "bulk", "ids")
    sc <- rownames(data_venn)[data_venn$sc > 0 & !is.na(data_venn$sc)]
    bulk <- rownames(data_venn)[data_venn$bulk > 0 & !is.na(data_venn$bulk)]
    results <- list("sc"=sc, "bulk"=bulk)

    temp <- venn.diagram(
        results,
        fill=c("brown1", "cornflowerblue"),
        main="scATAC-seq vs bulk ATAC-seq",
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
    plot_file <- paste(out_dir, "/sc_vs_bulk.", day, ".venn.pdf", sep="")
    pdf(file=plot_file, height=2, width=2)
    grid.draw(temp)
    dev.off()

    
}
