#!/usr/bin/env Rscript

# run DAVID GO term analysis
# requires gene set and background gene set

library("RDAVIDWebService")

# folders files etc
args <- commandArgs(trailingOnly=TRUE)
gene_list_file <- args[1]
background_list_file <- args[2]
out_dir <- args[3]

# read in gene list and background gene list
gene_list <- read.table(gzfile(gene_list_file), stringsAsFactors=FALSE)$V1
background_list <- read.table(gzfile(background_list_file), stringsAsFactors=FALSE)$V1

prefix <- sub('\\.txt.gz$', "", basename(gene_list_file))
padj_cutoff <- 0.1

# open david
david<-DAVIDWebService(
    email="danielskim@stanford.edu",
    url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

# add lists
result <- addList(david, background_list,
                  idType="ENSEMBL_GENE_ID",
                  listName='background',
                  listType="Background")
result <- addList(david, gene_list,
                  idType="ENSEMBL_GENE_ID",
                  listName=prefix,
                  listType="Gene")

# TODO figure out what column is actually being thresholded
annotation_categories <- c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL")
for (i in 1:length(annotation_categories)) {

    # choose annotations and get clusters
    setAnnotationCategories(david, annotation_categories[i])
    termCluster<-getClusterReport(david, type="Term")
    
    # and save out
    out_file <- paste(out_dir, "/", prefix, ".go_rdavid.", annotation_categories[i],".txt", sep="")
    getFunctionalAnnotationChartFile(object=david, fileName=out_file, threshold=padj_cutoff)

    # TODO - make plots?
    
}
