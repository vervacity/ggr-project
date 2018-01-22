"""Description: histone analyses
"""



def merge_overlap_regions(
        histone_overlap_files,
        out_master_bed,
        # end here

        args, histone, method="atac_midpoint", filter_w_overlap=True):
    """modularize the master region setup for easy switching
    1) use naive overlap set
    2) go atac-centric with fixed extend length from atac midpoint
    """
    # set up naive overlap based master regions
    logging.info("HISTONE: {}: Generating master regions...".format(histone))
    args.chipseq["histones"][histone]["overlap_master_regions"] = "{0}/ggr.{1}.overlap.master.bed.gz".format(
        args.folders["data_dir"], histone)
    if not os.path.isfile(args.chipseq["histones"][histone]["overlap_master_regions"]):
        histone_overlap_files = sorted(
            glob.glob("{0}/{1}".format(
                args.chipseq["data_dir"], args.chipseq["histones"][histone]["overlap_glob"])))
        logging.info("Master regions using: {}".format(" ".join(histone_overlap_files)))
        merge_regions(histone_overlap_files, args.chipseq["histones"][histone]["overlap_master_regions"])
    
    if method == "atac_midpoint":
        # If centered on ATAC region midpoint, then extract the midpoint and then extend out with bedtools slop
        args.atac["master_slop_bed"] = "{}.slop_{}bp.bed.gz".format(
            args.atac["master_bed"].split(".bed")[0],
            args.params["histones"][histone]["overlap_extend_len"])
        # this bit actually belongs in ATAC? integrative?
        if not os.path.isfile(args.atac["master_slop_bed"]):
            slop_bed = (
                "zcat {0} | "
                "awk -F '\t' 'BEGIN{{OFS=\"\t\"}} "
                "{{ midpoint=$2+int(($3-$2)/2); "
                "$2=midpoint; $3=midpoint+1; print }}' | "
                "bedtools slop -i stdin -g {1} -b {2} | "
                "gzip -c > {3}").format(
                    args.atac["master_bed"],
                    args.annot["chromsizes"],
                    args.params["histones"][histone]["overlap_extend_len"],
                    args.atac["master_slop_bed"])
            print slop_bed
            run_shell_cmd(slop_bed)

        if filter_w_overlap:
            # now intersect ATAC with the naive overlap files and only keep region if has an overlap
            args.chipseq["histones"][histone]["master_slop_marked_bed"] = "{}.{}-marked.bed.gz".format(
                args.atac["master_slop_bed"].split(".bed")[0],
                histone)
            if not os.path.isfile(args.chipseq["histones"][histone]["master_slop_marked_bed"]):
                keep_marked = (
                    "bedtools intersect -u -a {0} -b {1} | "
                    "gzip -c > {2}").format(
                        args.atac["master_slop_bed"],
                        args.chipseq["histones"][histone]["overlap_master_regions"],
                        args.chipseq["histones"][histone]["master_slop_marked_bed"])
                print keep_marked
                run_shell_cmd(keep_marked)
            master_regions = args.chipseq["histones"][histone]["master_slop_marked_bed"]
            
        else:
            master_regions = args.atac["master_slop_bed"]

    elif method == "naive_overlap":
        # If naive overlap, don't do anything extra - already generated the master file
        master_regions = args.chipseq["histones"][histone]["overlap_master_regions"]

    else:
        raise Exception("non existent master regions method!")
            
    return master_regions
