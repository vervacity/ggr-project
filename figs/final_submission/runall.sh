# final adjustments
WORK_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/final_submission/2021-07-23.fig_adjustments
CODE_DIR=/users/dskim89/git/ggr-project/figs/final_submission

# adjust timecourse plots
data_file=$WORK_DIR/validation.timecourse.melt.tsv
$CODE_DIR/plot.validation.timecourse.R $data_file $WORK_DIR/fig_5d

# adjust KO plots
data_file=$WORK_DIR/validation.KO.melt.tsv
$CODE_DIR/plot.validation.KO.R $data_file $WORK_DIR/fig_10c

# adjust ChIP-reChIP plots
data_file=$WORK_DIR/validation.chip_rechip.melt.tsv
$CODE_DIR/plot.validation.chip_rechip.R $data_file $WORK_DIR/fig_5e



