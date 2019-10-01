#!/bin/bash

# code
GIT_DIR=~/git/ggr-project
SCRIPT_DIR=$GIT_DIR/figs/fig_1.datasets

# dirs
GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a

# dataset matrix
#$SCRIPT_DIR/plot.dataset_matrix.R $GIT_DIR/ggr/data/ggr_datasets.txt

# PCA plots
#$SCRIPT_DIR/plot.pca.R fig_1.pca.atac.pdf $GGR_DIR/results/atac/timeseries/plots/ggr.atac.ends.counts.rep1.rlog.dynamic.filt.mat.txt.gz $GGR_DIR/results/atac/timeseries/plots/ggr.atac.ends.counts.rep2.rlog.dynamic.filt.mat.txt.gz
$SCRIPT_DIR/plot.pca.R fig_1.pca.rna.pdf $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.rep1.rlog.dynamic.mat.txt.gz $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.rep2.rlog.dynamic.mat.txt.gz

# gene set enrichment results
$SCRIPT_DIR/plot.gsea.R $GGR_DIR/results/rna/timeseries/gsea/ggr.rna.gsea_results.txt.gz

# vignettes - pulled from vis files


# heatmaps
#$SCRIPT_DIR/plot_trajs.py $GGR_DIR


# linking plots
