#!/usr/bin/env Rscript

# run DAVID GO term analysis
# requires gene set and background gene set

library("RDAVIDWebService")

# folders files etc
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
gene_list_file <- args[2]
background_list_file <- args[3]

gene_list <- read.table(gene_list_file, stringsAsFactors=FALSE)$V1
background_list <- read.table(background_list_file, stringsAsFactors=FALSE)$V1

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

# choose annotations and get clusters
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
termCluster<-getClusterReport(david, type="Term")
out_file <- paste(prefix, '.david_go_enrich.txt', sep='')
getFunctionalAnnotationChartFile(object=david, fileName=out_file, threshold=padj_cutoff)
