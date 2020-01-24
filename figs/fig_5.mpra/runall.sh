#!/bin/bash

# setup
GIT_DIR=/users/dskim89/git/ggr-project
MPRA_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results

# analysis happens through script
#$GIT_DIR/scripts/analyze.mpra.py

# fig 5b - vignettes

# HOXA1,GLI3 (early)
rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/ggr.TRAJ_LABELS-0.grammar-22.annot-15*pdf ./

# GRHL, ATF (mid)
rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/ggr.TRAJ_LABELS-1.grammar-99.annot-167*pdf ./

# FOXA, CEBP (late)
rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/ggr.TRAJ_LABELS-1.grammar-71.annot-137*pdf ./

# others of interest
#ggr.TRAJ_LABELS-2.grammar-74
#ggr.TRAJ_LABELS-2.grammar-55
#TRAJ_LABELS-1.grammar-99

# fig 5c - summary
rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/summary.expected_v_actual.pdf ./

# fig 5d - summary of rules
rsync -avz --progress $MPRA_DIR/results/combinatorial_rules/summary.rule_map.pdf ./


