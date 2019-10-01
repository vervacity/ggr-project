

SCRIPT_DIR=~/git/ggr-project/figs/fig_3.motifs_and_tfs
SIG_MOTIFS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim

GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a
PWM_FILE=$GGR_DIR/annotations/HOCOMOCOv11_core_pwms_HUMAN_mono.renamed.nonredundant.txt

# vignette - visualize region, and plot matching pwms
# help choosing
#$SCRIPT_DIR/choose_vignette.py
#~/git/ggr-project/figs/fig_1.datasets/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000189182 gene_vals.KRT77.pdf
#$SCRIPT_DIR/plot_pwm.py $PWM_FILE TP53
#$SCRIPT_DIR/plot_pwm.py $PWM_FILE CEBP

# allelic imbalance ATAC to validate importance scores
#$SCRIPT_DIR/fig_3-b.3.prep.py
#$SCRIPT_DIR/fig_3-b.4.plot.R results.w_imptscore.txt results.sig_only.txt

# pwm x rna (supplement)
PVALS_FILE=$SIG_MOTIFS_DIR/pvals.rna_filt.corr_filt.h5
#$SCRIPT_DIR/fig_3-c.0.plot.pwm_x_rna.R $PVALS_FILE pvals

# motifs and tfs
SUMMARY_DIR=$SIG_MOTIFS_DIR/summary
#$SCRIPT_DIR/fig_3-d.0.plot.motifs_and_tfs.R $SUMMARY_DIR/ggr.pwms_present_summary.txt $SUMMARY_DIR/ggr.pwms_patterns_summary.txt $SUMMARY_DIR/ggr.tfs_corr_summary.txt $SUMMARY_DIR/ggr.tfs_patterns_summary.txt
