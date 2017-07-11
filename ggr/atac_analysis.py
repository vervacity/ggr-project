"""Contains functions to run ATAC-specific analyses
"""

import os
import glob
import signal

import pandas as pd

from ggr.util.bed_utils import merge_regions
from ggr.util.parallelize import *
from ggr.util.filtering import remove_media_timepoints, filter_for_ids

def make_ends_file(tagAlign_file, out_file):
    '''
    Just get cut sites instead of the whole read. Use for ATAC-seq
    '''
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    
    # awk is significantly faster than opening in python
    make_ends = ("zcat {0} | "
                 #"head -n 10000 | "
                 "awk -F '\t' '{{ printf $1\"\t\"$2\"\t\"$2+1\"\tN\t1000\t\"$6\"\\n"
                 "\"$1\"\t\"$3-1\"\t\"$3\"\tN\t1000\t\"$6\"\\n\" }}' | "
                 "sort -k1,1 -k2,2n | "
                 "gzip -c "
                 "> {1}").format(tagAlign_file, out_file)
    print make_ends
    os.system(make_ends)

    return None


def make_atac_count_matrix(master_regions_file, tagalign_file_list, out_mat_file, tmp_dir='.'):
    '''
    For the regions in the master region file, gets cut site counts (5' and 3' sites of PE ATAC)
    and returns a matrix file with this information
    '''

    # Generate (sorted) cuts files from tagAlign files
    cuts_dir = '{}/cuts'.format(tmp_dir)
    os.system('mkdir -p {}'.format(cuts_dir))
    ends_queue = setup_multiprocessing_queue()
    for tagalign_file in tagalign_file_list:
        cuts_file = '{0}/{1}.cuts.tagAlign.gz'.format('{}'.format(cuts_dir),
                                                      os.path.basename(tagalign_file).split('.tagAlign')[0])
        if not os.path.isfile(cuts_file):
            print cuts_file
            ends_queue.put([make_ends_file, 
                [tagalign_file, cuts_file]])
    run_in_parallel(ends_queue)

    # Overlap with master regions and save into a master file
    # and then run bedtools coverage
    cuts_files = sorted(glob.glob('{}/cuts/*.tagAlign.gz'.format(tmp_dir)))
    count_cuts_queue = setup_multiprocessing_queue()
    counts_dir = '{}/counts'.format(tmp_dir)
    os.system('mkdir -p {}'.format(counts_dir))
    for cuts_file in cuts_files:
        cuts_in_regions_file = '{0}/counts/{1}.counts.txt.gz'.format(tmp_dir, os.path.basename(cuts_file).split('.tagAlign')[0])
        if not os.path.isfile(cuts_in_regions_file):
            get_coverage = ("bedtools coverage "
                            "-sorted "
                            "-counts "
                            "-a {0} -b {1} | "
                            "gzip -c > {2}").format(master_regions_file,
                                                    cuts_file,
                                                    cuts_in_regions_file)
            print get_coverage
            count_cuts_queue.put([os.system, [get_coverage]])
    run_in_parallel(count_cuts_queue)

    # Now combine into cuts for all timepoints file
    counts_files = sorted(glob.glob('{}/counts/*counts.txt.gz'.format(tmp_dir)))
    coverage_df = pd.DataFrame()
    if not os.path.isfile(out_mat_file):
        for counts_file in counts_files:
            filename_fields = os.path.basename(counts_file).split('.')
            header = '{0}_{1}'.format(filename_fields[0].split('-')[1], filename_fields[4])

            counts = pd.read_csv(counts_file, sep='\t', header=None)
            counts.columns = ['chrom', 'start', 'stop', 'counts']

            coverage_df[header] = counts['counts']
        coverage_df = coverage_df.set_index(counts['chrom'] + ':' + counts['start'].map(str) + '-' + counts['stop'].map(str))
        coverage_df.to_csv(out_mat_file, sep='\t', index_label=False, compression='gzip')

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
    atac_peak_files = glob.glob('{0}/{1}'.format(args.atac['data_dir'],
                                                 args.atac['idr_peak_glob']))
    merge_regions(atac_peak_files, args.atac['master_bed'])

    # TODO any plots here?
    # chromatin state annotations? enh/prom/etc?
    # ChromHMM enrichments
    
    
    # get read counts in these regions
    args.atac['counts_mat'] = '{0}/{1}.cuts.counts.mat.txt.gz'.format(DATA_DIR, atac_prefix)
    atac_tagalign_files = glob.glob('{0}/{1}'.format(args.atac['data_dir'],
                                                     args.atac['tagalign_glob']))
    if not os.path.isfile(args.atac['counts_mat']):
        make_atac_count_matrix(args.atac['master_bed'],
                               atac_tagalign_files,
                               args.atac['counts_mat'],
                               tmp_dir=ATAC_DIR)

    # use DESeq to run rlog normalization (for visualization and timeseries)
    args.atac['counts_rlog_mat'] = '{0}/{1}.cuts.rlog.mat.txt.gz'.format(DATA_DIR, atac_prefix)
    if not os.path.isfile(args.atac['counts_rlog_mat']):
        run_rlog = 'normalize_count_mat.R {0} {1}'.format(args.atac['counts_mat'], 
        	                                          args.atac['counts_rlog_mat'])
        print run_rlog
        os.system(run_rlog)
        
    # Remove media influenced timepoints first (both counts and rlog norm matrices)
    args.atac['counts_nomedia_mat'] = '{}.nomedia.mat.txt.gz'.format(args.atac['counts_mat'].split('.mat')[0])
    if not os.path.isfile(args.atac['counts_nomedia_mat']):
	remove_media_timepoints(args.atac['counts_mat'], args, args.atac['counts_nomedia_mat'])
    args.atac['counts_rlog_nomedia_mat'] = '{}.nomedia.mat.txt.gz'.format(args.atac['counts_rlog_mat'].split('.mat')[0])
    if not os.path.isfile(args.atac['counts_rlog_nomedia_mat']):
	remove_media_timepoints(args.atac['counts_rlog_mat'], args, args.atac['counts_rlog_nomedia_mat'])
    
    # with read counts run DESeq2 as timeseries mode
    args.atac['dynamic_region_ids'] = '{0}/{1}.dynamic.ids.txt.gz'.format(args.folders['atac_timeseries_dir'],
                                                                       atac_prefix)
    if not os.path.isfile(args.atac['dynamic_region_ids']):
	run_timeseries_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3}".format(args.atac['counts_nomedia_mat'],
                                                                               '{0}/deseq2/{1}'.format(args.folders['atac_timeseries_dir'], atac_prefix),
                                                                               args.params['deseq2_fdr'],
                                                                               args.atac['dynamic_region_ids'])
        print run_timeseries_deseq2
        os.system(run_timeseries_deseq2)
        
    # filter matrix for the dynamic and stable ones
    args.atac['dynamic_mat'] = '{0}/{1}.dynamic.mat.txt.gz'.format(args.folders['atac_dynamic_dir'], atac_prefix)
    os.system('mkdir -p {}'.format(args.folders['atac_dynamic_dir']))
    if not os.path.isfile(args.atac['dynamic_mat']):
        filter_for_ids(args.atac['counts_rlog_nomedia_mat'],
                       args.atac['dynamic_region_ids'],
                       args.atac['dynamic_mat'])

    args.atac['stable_mat'] = '{0}/{1}.stable.mat.txt.gz'.format(args.folders['atac_stable_dir'], atac_prefix)
    os.system('mkdir -p {}'.format(args.folders['atac_stable_dir']))
    if not os.path.isfile(args.atac['stable_mat']):
        filter_for_ids(args.atac['counts_rlog_nomedia_mat'],
                       args.atac['dynamic_region_ids'],
                       args.atac['stable_mat'],
                       opposite=True)

    # quick plot of numbers (dynamic vs stable)
    
    
    quit()
        
    # Run trajectories using DP_GP clustering
    if not os.path.isdir(args.folders['atac_dp-gp_dir']):
        os.system('mkdir -p {}'.format(args.folders['atac_dp-gp_dir']))
        cluster_dp_gp = ("DP_GP_cluster.py ",
                         "-i {0} ",
                         "-o {1}").format(args.atac['dynamic_mat'],
                                          '{0}/{1}'.format(args.atac['atac_dp-gp_dir'], atac_prefix))
        print cluster_dp_gp
        os.system(cluster_dp_gp)
        
        
    # Make plotting functions to plot trajectories nicely
    


    # Run homer on each group


    # as well as GREAT

    

    # overlap with histone information
    # TODO - actually use peak calls to determine overlaps and marks
    
    

    


    return None
