#!/bin/bash

# code
GIT_DIR=~/git/ggr-project
R_DIR=$GIT_DIR/R
SCRIPT_DIR=$GIT_DIR/figs/fig_1.datasets

# dirs
GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a

# --------------------
# MAIN
# --------------------

# dataset matrix
#$SCRIPT_DIR/plot.dataset_matrix.R $GIT_DIR/ggr/data/ggr_datasets.txt

# PCA plots
atac_rep1_mat_file=$GGR_DIR/results/atac/timeseries/plots/ggr.atac.ends.counts.rep1.rlog.dynamic.filt.mat.txt.gz
atac_rep2_mat_file=$GGR_DIR/results/atac/timeseries/plots/ggr.atac.ends.counts.rep2.rlog.dynamic.filt.mat.txt.gz
$R_DIR/plot.pca.R fig_1.pca.atac.pdf $atac_rep1_mat_file $atac_rep2_mat_file

rna_rep1_mat_file=$GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.rep1.rlog.dynamic.mat.txt.gz
rna_rep2_mat_file=$GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.rep2.rlog.dynamic.mat.txt.gz
$R_DIR/plot.pca.R fig_1.pca.rna.pdf $rna_rep1_mat_file $rna_rep2_mat_file

# gene set enrichment results
#$SCRIPT_DIR/plot.gsea.R $GGR_DIR/results/rna/timeseries/gsea/ggr.rna.gsea_results.txt.gz

# vignettes - pulled from vis files
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000170423 gene_vals.KRT78.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000245848 gene_vals.CEBPA.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000132470 gene_vals.ITGB4.pdf

# heatmaps
# NOTE: this includes stable
#$SCRIPT_DIR/plot_trajs.py $GGR_DIR

# linking plots
#$SCRIPT_DIR/plot_linking.py $GGR_DIR


# --------------------
# SUPPLEMENTS
# --------------------

# biomarker panel code
biomarker_mat_file=$GGR_DIR/results/rna/plots/biomarkers.mat.txt.gz
#$R_DIR/plot.gene_mat_subset.R $biomarker_mat_file fig_1.biomarker_panel.pdf

# global summary stats - ATAC
atac_static_counts_file=$GGR_DIR/results/atac/plots/ggr.atac.idr.timepoint_count_summary.static.txt
#$R_DIR/plot.counts.global.R $atac_static_counts_file fig_1.SUPPL.atac.timepoint_count_summary.static.pdf static
atac_dynamic_counts_file=$GGR_DIR/results/atac/timeseries/deseq2/differential_summary.d0_baseline_only.txt.gz
#$R_DIR/plot.counts.global.R $atac_dynamic_counts_file fig_1.SUPPL.atac.timepoint_count_summary.dynamic.pdf dynamic

# global summary stats - RNA
rna_static_counts_file=$GGR_DIR/results/rna/plots/ggr.rna.gene_count_summary.static.txt
#$R_DIR/plot.counts.global.R $rna_static_counts_file fig_1.SUPPL.rna.gene_count_summary.static.pdf static
rna_dynamic_counts_file=$GGR_DIR/results/rna/timeseries/deseq2/differential_summary.d0_baseline_only.txt.gz
#$R_DIR/plot.counts.global.R $rna_dynamic_counts_file fig_1.SUPPL.rna.gene_count_summary.dynamic.pdf dynamic

# PCA for other datasets
h3k27ac_rep1_mat_file=$GGR_DIR/results/histones/H3K27ac/plots/ggr.histone.H3K27ac.midpoints.counts.rep1.rlog.filt.mat.txt.gz
h3k27ac_rep2_mat_file=$GGR_DIR/results/histones/H3K27ac/plots/ggr.histone.H3K27ac.midpoints.counts.rep2.rlog.filt.mat.txt.gz
$R_DIR/plot.pca.R fig_1.pca.histone.H3K27ac.pdf $h3k27ac_rep1_mat_file $h3k27ac_rep2_mat_file

h3k4me1_rep1_mat_file=$GGR_DIR/results/histones/H3K4me1/plots/ggr.histone.H3K4me1.midpoints.counts.rep1.rlog.filt.mat.txt.gz
h3k4me1_rep2_mat_file=$GGR_DIR/results/histones/H3K4me1/plots/ggr.histone.H3K4me1.midpoints.counts.rep2.rlog.filt.mat.txt.gz
$R_DIR/plot.pca.R fig_1.pca.histone.H3K4me1.pdf $h3k4me1_rep1_mat_file $h3k4me1_rep2_mat_file

h3k27me3_rep1_mat_file=$GGR_DIR/results/histones/H3K27me3/plots/ggr.histone.H3K27me3.midpoints.counts.rep1.rlog.filt.mat.txt.gz
h3k27me3_rep2_mat_file=$GGR_DIR/results/histones/H3K27me3/plots/ggr.histone.H3K27me3.midpoints.counts.rep2.rlog.filt.mat.txt.gz
$R_DIR/plot.pca.R fig_1.pca.histone.H3K27me3.pdf $h3k27me3_rep1_mat_file $h3k27me3_rep2_mat_file

hichip_mat_file=$GGR_DIR/results/linking/hichip/ggr.linking.ALL.reps.mat.txt.gz
$R_DIR/plot.pca.R fig_1.pca.hichip.pdf $hichip_mat_file


# --------------------
# MISC
# --------------------

# other vignettes?
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000091409 gene_vals.ITGA6.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000065618 gene_vals.COL17A1.pdf
