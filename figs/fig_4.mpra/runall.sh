#!/bin/bash

# setup
GIT_DIR=/users/dskim89/git/ggr-project
MPRA_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results/results/combinatorial_rules
SCRIPT_DIR=$GIT_DIR/figs/fig_4.mpra

# analysis happens through script
#$GIT_DIR/scripts/analyze.mpra.py

# fig 5b - vignettes

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

# fig 5c - summary
#rsync -avz --progress $MPRA_DIR/summary.expected_v_actual.pdf ./

# fig 5d - summary of rules
#rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/summary.rule_map.pdf ./
PREDICTED_RULE_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle.complete/grammars.annotated.manual_filt.merged.final

# build new summary
$SCRIPT_DIR/plot.validated_rule_map.R $MPRA_DIR/summary.txt.gz $PREDICTED_RULE_DIR/grammars_summary.txt $PREDICTED_RULE_DIR/grammars.filt.atac.mat.txt $PREDICTED_RULE_DIR/grammars.filt.rna.mat.txt $PREDICTED_RULE_DIR/grammars.filt.go_terms.mat.txt fig_4.validated_rule_map.pdf

