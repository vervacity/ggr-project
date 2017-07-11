"""Contains functions to run ATAC-specific analyses
"""


def generate_master_regions(inputs_bed, out_bed):
    """Given a list of bed files, make into master file
    """
    merge_all = ("zcat {0} | "
                 "sort -k1,1 -k2,2n | "
                 "bedtools merge -i stdin | "
                 "gzip -c "
                 "> {1}").format(' '.join(inputs), out_bed)
    if not os.path.isfile(all_atac_bed):
        os.system(merge_all)

    return None





def run(args):
    """Run everything
    """
    # make folders
    DATA_DIR = args.folders['data_dir']
    ATAC_DIR = args.folders['atac_dir']
    os.system('mkdir -p {0} {1}'.format(DATA_DIR, ATAC_DIR))
    atac_prefix = 'ggr.atac'
    
    # generate a master regions file
    args.atac['master_bed'] = '{0}/{1}.idr.master.bed.gz'.format(DATA_DIR, atac_prefix)
    atac_peak_files = glob.glob('{0}/{1}'.format(args.atac['data_dir'], args.atac['idr_peak_glob']))
    generate_master_regions(atac_peak_files, args.atac['master_bed'])
        
    # get read counts in these regions
    atac_counts_mat_file = '{0}/{1}.cuts.counts.mat.txt.gz'.format(DATA_DIR, atac_prefix)
    if not os.path.isfile(atac_counts_mat_file):
        ggr_preprocess.make_count_matrix(all_atac_bed, atac_tagalign_files, atac_counts_mat_file, tmp_dir=RESULTS_ATAC_DIR)
        
        # remove media influenced timepoints first
        
    # use DESeq to run rlog normalization
    atac_rlog_norm_file = '{0}/{1}.cuts.rlog.mat.txt.gz'.format(DATA_DIR, atac_prefix)
    if not os.path.isfile(atac_rlog_norm_file):
        run_rlog = 'normalize_count_mat.R {0} {1}'.format(atac_counts_mat_file, atac_rlog_norm_file)
        os.system(run_rlog)
    

    
    
    
    

    


    return None
