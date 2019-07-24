

SCRIPT_DIR=~/git/ggr-project/figs/fig_3.homotypic
MOTIFS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.background/ggr.scanmotifs.h5
MOTIFS_EXTENDED_FILE=/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-17.footprinting/motifs.input_x_grad.lite/ggr.scanmotifs.h5
SIG_PWMS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_patterns_summary.txt

# fig 3 - homotypic spacing
$SCRIPT_DIR/fig_3-b.0.plot.spacing.R $SIG_PWMS_FILE $MOTIFS_FILE
