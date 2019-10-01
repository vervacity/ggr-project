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
#$SCRIPT_DIR/plot.pca.R fig_1.pca.rna.pdf $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.rep1.rlog.dynamic.mat.txt.gz $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.rep2.rlog.dynamic.mat.txt.gz

# gene set enrichment results
#$SCRIPT_DIR/plot.gsea.R $GGR_DIR/results/rna/timeseries/gsea/ggr.rna.gsea_results.txt.gz

# vignettes - pulled from vis files
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000091409 gene_vals.ITGA6.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000170423 gene_vals.KRT78.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000245848 gene_vals.CEBPA.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000065618 gene_vals.COL17A1.pdf
#$SCRIPT_DIR/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000132470 gene_vals.ITGB4.pdf

# heatmaps
#$SCRIPT_DIR/plot_trajs.py $GGR_DIR

# linking plots
