"""Contains wrappers to quickly work with bioinformatics tools
"""

import os

from ggr.util.utils import run_shell_cmd

def run_homer(
        bed_file,
        background_bed_file,
        out_dir,
        mknown=None,
        parallel=24):
    """Generic wrapper to run homer on a set of regions (vs background)
    """
    assert bed_file.endswith(".gz")
    assert background_bed_file.endswith(".gz")
    
    run_homer_cmd = (
        "findMotifsGenome.pl <(zcat {0}) hg19 {2} "
        "-bg <(zcat {1}) "
        "-p {3} "
        "-nomotif").format(
            bed_file,
            background_bed_file,
            out_dir,
            parallel)

    if mknown is not None:
        run_homer_cmd += " -mknown {}".format(mknown)
    
    print run_homer_cmd
    os.system('GREPDB="{}"; /bin/bash -c "$GREPDB"'.format(run_homer_cmd))

    return None


def run_great(
        bed_file,
        background_bed_file,
        out_dir):
    """Generic wrapper to run GREAT from R (rGREAT package)
    """
    assert bed_file.endswith(".gz")
    assert background_bed_file.endswith(".gz")
    
    prefix = os.path.basename(bed_file).split(".bed")[0].split(".narrowPeak")[0]
    
    run_rgreat = (
        "bioinformatics.go.rgreat.R {0} {1} {2}").format(
            bed_file,
            background_bed_file,
            "{}/{}".format(out_dir, prefix))
    run_shell_cmd(run_rgreat)
    
    return None

def run_gprofiler(
        gene_set_file,
        background_gene_set_file,
        out_dir,
        ordered=False,
        header=True):
    """
    """
    if ordered:
        ordered_val = 1
    else:
        ordered_val = 0
    
    if header:
        header_val = 1
    else:
        header_val = 0
        
    #gprofiler_cmd = "bioinformatics.go.gProfileR.R {} {} {} {}".format(
    gprofiler_cmd = "Rscript /users/dskim89/git/ggr-project/R/bioinformatics.go.gProfileR.R {} {} {} {} {}".format(
        gene_set_file,
        background_gene_set_file,
        out_dir,
        header_val,
        ordered_val)
    run_shell_cmd(gprofiler_cmd)

    return


def run_gsea(rank_file, gene_sets=[]):
    """take gene set and run gsea
    make sure the file has HGNC ids and rank values
    """
    
    
    return


def run_bioinformatics_on_bed(bed_file, background_bed_file, out_dir,
                              mknown=None, mknown_name="CUSTOM"):
    """Given a bed file and background bed file, 
    run HOMER/GREAT
    mknown - homer style motifs file to run with different motifs
    """
    assert bed_file.endswith(".gz")
    assert background_bed_file.endswith(".gz")
    
    prefix = os.path.basename(bed_file).split(".bed")[0].split(".narrowPeak")[0]
    
    # run homer
    homer_dir = "{}/homer/{}".format(out_dir, prefix)
    run_shell_cmd("mkdir -p {}".format(homer_dir))
    run_homer(
        bed_file,
        background_bed_file,
        homer_dir)

    # run homer with custom motifs if needed
    if mknown is not None:
        homer_dir = "{}/homer_{}/{}".format(out_dir, mknown_name, prefix)
        run_shell_cmd("mkdir -p {}".format(homer_dir))
        run_homer(
            bed_file,
            background_bed_file,
            homer_dir,
            mknown=mknown)
    
    # run GREAT
    great_dir = "{}/great/{}".format(out_dir, prefix)
    run_shell_cmd("mkdir -p {}".format(great_dir))
    run_great(
        bed_file,
        background_bed_file,
        great_dir)
        
    return homer_dir, great_dir


def make_deeptools_heatmap(
        point_file,
        bigwig_files,
        prefix,
        sort=False,
        kval=4,
        referencepoint='TSS',
        extend_dist=2000, # 1000
        bin_total=100,
        color="Blues"):
    """Uses deeptools to make a profile heatmap
    """
    # set up bin size
    bin_size = extend_dist * 2 / bin_total
    
    # compute matrix - extract from bigwigs
    point_matrix = '{}.point.mat.gz'.format(prefix)
    deeptools_compute_matrix = (
        "computeMatrix reference-point "
        "--referencePoint {0} "
        "-b {1} -a {1} -bs {2} "
        "-R {3} "
        "-S {4} "
        #"--skipZeros "
        "-o {5} ").format(
            referencepoint,
            extend_dist,
            bin_size,
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
    point_plot = '{}.heatmap.profile.pdf'.format(prefix)
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


