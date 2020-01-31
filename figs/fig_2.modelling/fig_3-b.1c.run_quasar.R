#!/usr/bin/env Rscript

# description: script to run quasar analysis
library(QuASAR)
library(qvalue)

args <- commandArgs(trailingOnly=TRUE)
results_prefix <- args[1]
fileNames <- args[2:length(args)]

# params
min_genotype_coverage <- 50
min_infer_coverage <- 50
fdr <- 0.10

# set up sample data
ase.dat <- UnionExtractFields(fileNames, combine=TRUE)
ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=min_genotype_coverage)
sample.names <- colnames(ase.dat.gt$ref)

# genotype from reads
ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))

# save out genotypes
genotypes_file <- paste(results_prefix,  ".genotypes.txt", sep="")
print(genotypes_file)
out_dat <- data.frame(ase.dat.gt$annotations[, -5], map=ase.joint$gt)
write.table(out_dat, file=genotypes_file, row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")

# infer asATAC
ourInferenceData <- aseInference(
    gts=ase.joint$gt,
    eps.vect=ase.joint$eps,
    priors=ase.dat.gt$gmat,
    ref.mat=ase.dat.gt$ref,
    alt.mat=ase.dat.gt$alt,
    min.cov=min_infer_coverage,
    sample.names=sample.names,
    annos=ase.dat.gt$annotations)

# for each sample, save out with qval
start <- 0
for (sample_name in names(ourInferenceData)) {

    # get results and get qvals
    res <- ourInferenceData[[sample_name]]$dat
    qobj <- qvalue(res$pval2.het.ind.)
    res$qval <- qobj$qvalues

    # write all
    prefix <- strsplit(sample_name, ".in.gz")[[1]][1]
    out_file <- paste(prefix, "results_all.txt.gz", sep=".")
    write.table(res, gzfile(out_file), sep="\t", quote=FALSE, row.names=FALSE)

    # write sig (at FDR)
    res_sig <- res[res$qval < fdr,]
    out_file <- paste(prefix, "results_sig.txt.gz", sep=".")
    write.table(res_sig, gzfile(out_file), sep="\t", quote=FALSE, row.names=FALSE)

    # adjust prefix and colnames
    prefix <- basename(prefix)
    prefix <- strsplit(prefix, ".", fixed=TRUE)[[1]][1]
    prefix <- strsplit(prefix, "-", fixed=TRUE)[[1]][2]
    res$id <- paste(res$annotations.rsID, res$annotations.chr, res$annotations.pos0, sep="_")
    simple_results <- res[, c("id", "betas", "qval")]
    new_colnames <- colnames(simple_results)
    new_colnames[2:length(new_colnames)] <- paste(
        prefix,
        results_prefix,
        new_colnames[2:length(new_colnames)],
        sep=".")
    colnames(simple_results) <- new_colnames

    # merge into master file
    if (start == 0) {
        all_results <- simple_results
        start <- 1
    } else {
        all_results <- merge(
            all_results,
            simple_results,
            by="id",
            all=TRUE)
    }

    print(prefix)
    print(dim(all_results))
    
}

# clean up
all_results[is.na(all_results)] <- 0

# make an overlap file with all variants in set
all_results_file <- paste(results_prefix,  ".results_ALL.txt", sep="")
write.table(all_results, all_results_file, sep="\t", quote=FALSE, row.names=FALSE)
