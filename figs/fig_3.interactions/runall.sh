#!/bin/bash

SCRIPT_DIR=/users/dskim89/git/ggr-project/figs/fig_4.interactions

# fig 4b
WORK_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/ggr.tronn.2019-06-19.interactions
endogenous_file=$WORK_DIR/interactions.summary.endogenous.2.txt
simulations_file=$WORK_DIR/interactions.summary.sims.txt
#$SCRIPT_DIR/fig_4-b.0.plot.interactions.R $endogenous_file $simulations_file

# 4c - predicted rule map
WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle/grammars.annotated.manual_filt.merged.final
# just copy over from relevant dir

#$SCRIPT_DIR/plot.wrapper.rule_map.py $WORK_DIR/grammars_summary.txt .
