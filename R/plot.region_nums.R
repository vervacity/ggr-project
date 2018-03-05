#!/usr/bin/env Rscript

# description: plot out region nums
#library(extrafont)
library(ggplot2)
library(reshape2)

library(grid)
library(gridGraphics)

library(extrafont)

# read in args
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
out_plot <- args[2]

# read data, adjust as needed
data <- read.table(filename, sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
data$timepoint <- as.numeric(sub("d", "", data$timepoint))

# unmelt and remelt to make things match
#data <- dcast(data = data,formula = timepoint~type,fun.aggregate = sum,value.var = "count")
#data_melted <- melt(data, id.vars="timepoint")
data_melted <- data

data_melted$variable <- data$type
data_melted$value <- data$count
data_melted$timepoint <- factor(data_melted$timepoint)

#assays <- unique(data_melted$variable)


p <- ggplot(data_melted, aes(x=timepoint, y=value, fill=variable)) +
    #geom_line() + 
    geom_bar(
        stat="identity",
        width=0.2) +
    facet_grid(. ~ variable, scales="free", space="free") + 
    theme_bw() + 
    theme(
        text=element_text(family="ArialMT"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave(out_plot, height=3, width=10)
embed_fonts(out_plot)

quit()


















make_barplot <- function(data_melted) {
    p <- ggplot(data_melted, aes(x=timepoint, y=value, fill=variable)) +
        geom_bar(
            stat="identity",
            position=position_dodge()) +
        theme_bw()
    return(p)
}

# set up viewportsx
grab_grob <- function(fn) {
    grid.echo(fn)
    grid.grab()
}

gl <- lapply(1:length(assays), function(i) {
                 data_tmp <- data_melted[data_melted$variable == assays[i],]
                 ggplot(data_tmp, aes(x=timepoint, y=value, fill=variable)) +
                     geom_bar(
                         stat="identity",
                         position=position_dodge()) +
                         theme_bw()
             })

grid.newpage()
library(gridExtra)

pdf(out_plot)
grid.arrange(grobs=gl, nrow=length(assays), ncol=1)

dev.off()





quit()

grobs <- c()
for (i in 1:length(assays)) {
    
    data_tmp <- data_melted[data_melted$variable == assays[i],]
    
    fn <- function() make_barplot(data_tmp)
    grobs <- c(grobs, fn())
}




#print("here")
#print(grobs)

library(gridExtra)
print(class(grobs))

pdf(out_plot)
grid.arrange(grobs, nrow=length(assays), ncol=1)

dev.off()

quit()

grid.newpage()
lay <- grid.layout(nrow=length(assays), ncol=1)
pushViewport(viewport(layout=lay))

print(length(assays))

for (i in 1:length(grobs)){
    #print(i)
    #print(grobs[i],
    #      vp=viewport(
    #          layout.pos.row=i,
    #          layout.pos.col=1, clip=FALSE))
    grid.draw(editGrob(
        grobs[[i]],
        vp=viewport(
            layout.pos.row=i,
            layout.pos.col=1, clip=FALSE)))
    
}

upViewport()

dev.off()


quit()




# plot
p <- ggplot(data_melted, aes(x=timepoint, y=value, fill=variable)) +
    geom_bar(
        stat="identity",
        position=position_dodge()) +
    theme_bw()

ggsave(out_plot)





