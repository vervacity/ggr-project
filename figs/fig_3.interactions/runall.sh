#!/bin/bash

SCRIPT_DIR=/users/dskim89/git/ggr-project/figs/fig_3.interactions
TRONN_R=/users/dskim89/git/tronn/R

GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a
ANNOT_DIR=/mnt/lab_data3/dskim89/ggr/annotations

# -------------
# A
# -------------

# schematic

# -------------
# B
# -------------

# MAIN: pairwise interactions
WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle.complete
endogenous_file=$WORK_DIR/synergy_calcs.endog/interactions.summary.txt
simulations_file=$WORK_DIR/synergy_calcs.sims/interactions.summary.txt
#$SCRIPT_DIR/fig_4-b.0.plot.interactions.R $endogenous_file $simulations_file

# SUPPL: pairwise combinations, using various signal sources
ordering_file=ordering.txt # from prev script run, apply ONLY to mut matrices
summary_files=/mnt/lab_data3/dskim89/ggr/nn/2020-01-13/grammars*/grammars.annotated/grammar_summary.txt
#$TRONN_R/plot.combinations.adjacencies.R fig_S3.unfiltered FALSE $summary_files
#$TRONN_R/plot.combinations.adjacencies.R fig_S3.go_filt TRUE $summary_files

# -------------
# C
# -------------

# MAIN: some nice vignettes <- actually, only pull from those that succeeded?
synergy_dir=$WORK_DIR/synergy
summary_file=$WORK_DIR/grammars.annotated.manual_filt.merged.final/grammars_summary.txt
summary_file=/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results/results/combinatorial_rules/summary.txt.gz
#links_file=$GGR_DIR/results/linking/hichip/ggr.linking.ALL.overlap.interactions.txt.gz
links_file=$GGR_DIR/results/linking/proximity/ggr.linking.ALL.overlap.interactions.txt.gz
tss_file=$GGR_DIR/annotations/ggr.rna.tss.expressed.bed.gz
filter_genes_file=$GGR_DIR/annotations/hg19.genodermatoses.ensembl_ids.txt.gz
#filter_genes_fil=$GGR_DIR/data/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat.txt.gz
chromsizes=$ANNOT_DIR/hg19.chrom.sizes

# genodermatoses filt
rm -r vignettes.genodermatoses # REMOVE THIS LATER
$SCRIPT_DIR/select_vignettes.py fig_3.vignettes $summary_file $synergy_dir $links_file $tss_file $filter_genes_file $chromsizes vignettes.genodermatoses


# SUPPL: predicted rule map
WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle/grammars.annotated.manual_filt.merged.final

# SUPPL: vignette of region with GWAS hit?


# just copy over
