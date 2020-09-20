#!/bin/bash

# setup
GIT_DIR=/users/dskim89/git/ggr-project
MPRA_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results
SCRIPT_DIR=$GIT_DIR/figs/fig_4.mpra

# analysis happens through script
#$GIT_DIR/scripts/analyze.mpra.py


# ----------------------------------
# VIGNETTES
# ----------------------------------

# Main

# HOXA1,GLI3 (early)
#rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/ggr.TRAJ_LABELS-0.grammar-22.annot-15*pdf ./

# GRHL, ATF (mid)
#rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/ggr.TRAJ_LABELS-1.grammar-99.annot-167*pdf ./

# FOXA, CEBP (late)
#rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/ggr.TRAJ_LABELS-1.grammar-71.annot-137*pdf ./

# others of interest
#ggr.TRAJ_LABELS-2.grammar-74
#ggr.TRAJ_LABELS-2.grammar-55
#TRAJ_LABELS-1.grammar-99

# ----------------------------------
# SUMMARY
# ----------------------------------

# fig 5c - summary
#rsync -avz --progress $MPRA_DIR/summary.expected_v_actual.pdf ./

# ----------------------------------
# FULL SUMMARY
# ----------------------------------

# fig 5d - summary of rules
#rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/summary.rule_map.pdf ./

# build new summary
PREDICTED_RULE_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle.complete/grammars.annotated.manual_filt.merged.final
#$SCRIPT_DIR/plot.validated_rule_map.R $MPRA_DIR/results/combinatorial_rules/summary.txt.gz $PREDICTED_RULE_DIR/grammars_summary.txt $PREDICTED_RULE_DIR/grammars.filt.atac.mat.txt $PREDICTED_RULE_DIR/grammars.filt.rna.mat.txt $PREDICTED_RULE_DIR/grammars.filt.go_terms.mat.txt fig_4.validated_rule_map.pdf

# ----------------------------------
# SUPPL
# ----------------------------------

# QC
# plasmid lib - distribution of representation across whole lib
PLASMID_FILE=$MPRA_DIR/results/plasmid.dedup.counts.txt.gz
COUNTS_FILE=$MPRA_DIR/results/ggr.counts.raw.mat.txt.gz
#SIGNAL_FILE=$MPRA_DIR/results/ggr.signal.w_metadata.mat.txt.gz
SIGNAL_FILE=$MPRA_DIR/results/ggr.signal.mat.txt.gz
#$GIT_DIR/R/qc.mpra.R $PLASMID_FILE $COUNTS_FILE $SIGNAL_FILE

# mpra vs epigenome signals and fitting
# NOTE that epigenome signals are first plotted using analyze mpra script analyze.mpra.py
# and then the nn ones are in the other analysis folder
fit_dir=/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2020-06-29.fit_nn
