

SCRIPT_DIR=~/git/ggr-project/figs/fig_3.motifs_and_tfs
SIG_MOTIFS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim

# fig 3 b


# fig 3 c
PVALS_FILE=$SIG_MOTIFS_DIR/pvals.rna_filt.corr_filt.h5
#$SCRIPT_DIR/fig_3-c.0.plot.pwm_x_rna.R $PVALS_FILE pvals

# fig 3d, 3f
SUMMARY_DIR=$SIG_MOTIFS_DIR/summary
$SCRIPT_DIR/fig_3-d.0.plot.motifs_and_tfs.R $SUMMARY_DIR/ggr.pwms_present_summary.txt $SUMMARY_DIR/ggr.pwms_patterns_summary.txt $SUMMARY_DIR/ggr.tfs_corr_summary.txt $SUMMARY_DIR/ggr.tfs_patterns_summary.txt
