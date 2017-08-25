"""Contains wrappers to quickly work with bioinformatics tools
"""

import os

def run_homer(
        positives_bed,
        background_bed,
        out_dir,
        parallel=12):
    """Generic wrapper to run homer on a set of regions (vs background)
    """
    run_homer = ("findMotifsGenome.pl <(zcat {0}) hg19 {1} "
                 "-bg <(zcat {2}) "
                 "-p {3} "
                 "-nomotif").format(positives_bed,
                                    out_dir,
                                    background_bed,
                                    parallel)
    print run_homer
    os.system('GREPDB="{}"; /bin/bash -c "$GREPDB"'.format(run_homer))

    return None


def make_deeptools_heatmap(
        point_file,
        bigwig_files,
        prefix,
        sort=False,
        kval=4,
        referencepoint='TSS',
        extend_dist=2000, # 1000
        color="Blues"):
    '''
    Uses deeptools to make a profile heatmap
    '''

    # TODO(dk) try sorting using a BED with groups.
    
    # do the same with TSS
    point_matrix = '{}.point.mat.gz'.format(prefix)
    deeptools_compute_matrix = (
        "computeMatrix reference-point "
        "--referencePoint {0} "
        "-b {1} -a {1} "
        "-R {2} "
        "-S {3} "
        #"--skipZeros "
        "-o {4} ").format(
            referencepoint,
            extend_dist,
            point_file,
            ' '.join(bigwig_files),
            point_matrix)
    if not os.path.isfile(point_matrix):
        print deeptools_compute_matrix
        os.system(deeptools_compute_matrix)
        
    # set up sample labels as needed
    sample_labels = []
    for bigwig_file in bigwig_files:
        fields = os.path.basename(bigwig_file).split('.')
        sample_labels.append('{0}_{1}_{2}'.format(fields[3].split('-')[0],
                                                  fields[0].split('-')[1],
                                                  fields[4]))
        
    # make plot
    point_plot = '{}.heatmap.profile.png'.format(prefix)
    point_sorted_file = '{}.point.sorted.bed'.format(prefix)
    if sort == False:
        sorting = '--sortRegions=no'
    elif kval == 0:
        sorting = '--sortRegions=descend'
    else:
        sorting = '--kmeans {0} --regionsLabel {1}'.format(kval, ' '.join([str(i) for i in range(kval)]))
    deeptools_plot_heatmap = (
        "plotHeatmap -m {0} "
        "-out {1} "
        "--outFileSortedRegions {2} "
        "--colorMap {3} "
        "{4} "
        "--samplesLabel {5} "
        #"--zMax 450 " # TESTING
        "--xAxisLabel '' "
        "--refPointLabel Summit "
        "--legendLocation none "
        "--heatmapHeight 50 "
        "--boxAroundHeatmaps no "
        "--whatToShow 'heatmap and colorbar'").format(
            point_matrix,
            point_plot,
            point_sorted_file,
            color,
            sorting,
            ' '.join(sample_labels))
    if not os.path.isfile(point_plot):
        print deeptools_plot_heatmap
        os.system(deeptools_plot_heatmap)
    
        
    return None


