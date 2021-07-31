#!/usr/bin/env Rscript

library(GEOquery)

# description: use GEOquery to get GSMs for each GSE

args <- commandArgs(trailingOnly=TRUE)
metadata_file <- args[1]
out_file <- args[2]

# read metadata
metadata <- read.table(metadata_file, sep="\t", skip=1, header=TRUE)
print(ncol(metadata))

# set up output dataframe
updated_metadata <- data.frame(matrix(ncol=ncol(metadata)+3, nrow=0))
new_cols <- c("GSM", "bio_rep", "tech_rep")
colnames(updated_metadata) <- c(colnames(metadata), new_cols)

# get GSMs for each GSE
for (row_idx in 1:nrow(metadata)) {
#for (row_idx in 1:2) {
    print(metadata$Description[row_idx])
    row_metadata <- metadata[row_idx,]

    # get GSE
    gse_id <- strsplit(metadata$Dbxrefs[row_idx], ":")[[1]]
    gse_id <- gse_id[2]
    gse <- getGEO(gse_id, GSEMatrix=FALSE)

    # extract GSMs
    gsm_ids <- names(GSMList(gse))
    gsms_str <- paste(gsm_ids, collapse=";")

    # for each GSM, figure out what rep?
    for (gsm_idx in 1:length(gsm_ids)) {
        gsm_id <- gsm_ids[gsm_idx]
        gsm <- getGEO(gsm_id)

        # get bio rep, tech rep
        bio_rep <- NULL
        tech_rep <- NULL
        description <- Meta(gsm)$description
        for (description_idx in 1:length(description)) {

            if (grepl("biological rep", description[description_idx])) {
                bio_rep <- description[description_idx]
                bio_rep <- strsplit(bio_rep, ": ")[[1]]
                bio_rep <- bio_rep[2]
            }
            if (grepl("technical rep", description[description_idx])) {
                tech_rep <- description[description_idx]
                tech_rep <- strsplit(tech_rep, ": ")[[1]]
                tech_rep <- tech_rep[2]
            }
        }

        # update row metadata
        row_metadata$GSM <- gsm_id
        row_metadata$bio_rep <- bio_rep
        row_metadata$tech_rep <- tech_rep
        
        # write out new line
        updated_metadata <- rbind(updated_metadata, row_metadata)
        
    }
    
}

# timepoint
updated_metadata$timepoint <- gsub("^.+in ", "", updated_metadata$Description)
updated_metadata$timepoint <- gsub(" of.+$", "", updated_metadata$timepoint)
updated_metadata$timepoint <- gsub(" ", "-", updated_metadata$timepoint)

# add a title column
updated_metadata$title <- paste(
    gsub(" .+$", "", updated_metadata$Description),
    gsub(".", "", updated_metadata$timepoint, fixed=TRUE),
    updated_metadata$bio_rep,
    updated_metadata$tech_rep, sep="_")

# sort
updated_metadata <- updated_metadata[with(updated_metadata, order(Description, bio_rep, tech_rep)),]

# and save out to updated file
write.table(updated_metadata, out_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

