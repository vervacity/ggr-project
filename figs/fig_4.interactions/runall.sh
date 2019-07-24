WORK_DIR=/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-19.interactions
endogenous_file=$WORK_DIR/interactions.summary.endogenous.2.txt
simulations_file=$WORK_DIR/interactions.summary.sims.txt

SCRIPT_DIR=/users/dskim89/git/ggr-project/figs/fig_4.interactions

# fig 4b
$SCRIPT_DIR/fig_4-b.0.plot.interactions.R $endogenous_file $simulations_file
