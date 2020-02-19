#!/bin/bash

SCRIPT_DIR=/users/dskim89/git/ggr-project/figs/fig_3.interactions
TRONN_R=/users/dskim89/git/tronn/R

# -------------
# A
# -------------


# -------------
# B
# -------------

# MAIN: pairwise interactions
WORK_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/ggr.tronn.2019-06-19.interactions
endogenous_file=$WORK_DIR/interactions.summary.endogenous.2.txt
simulations_file=$WORK_DIR/interactions.summary.sims.txt

WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle.complete
endogenous_file=$WORK_DIR/synergy_calcs.endog/interactions.summary.txt
simulations_file=$WORK_DIR/synergy_calcs.sims/interactions.summary.txt

$SCRIPT_DIR/fig_4-b.0.plot.interactions.R $endogenous_file $simulations_file

# SUPPL: pairwise combinations, using various signal sources
summary_files=/mnt/lab_data3/dskim89/ggr/nn/2020-01-13/grammars*/grammars.annotated/grammar_summary.txt
#$TRONN_R/plot.combinations.adjacencies.R fig_S3.unfiltered FALSE $summary_files
#$TRONN_R/plot.combinations.adjacencies.R fig_S3.go_filt TRUE $summary_files

# -------------
# C
# -------------

# MAIN: some nice vignettes



# SUPPL: predicted rule map
WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle/grammars.annotated.manual_filt.merged.final
# just copy over
