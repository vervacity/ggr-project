

SCRIPT_DIR=~/git/ggr-project/figs/fig_3.motifs_and_tfs
SIG_MOTIFS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim

GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a
PWM_FILE=$GGR_DIR/annotations/HOCOMOCOv11_core_pwms_HUMAN_mono.renamed.nonredundant.txt

# vignette - visualize region, and plot matching pwms
#$SCRIPT_DIR/choose_vignette.py
#~/git/ggr-project/figs/fig_1.datasets/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000189182 gene_vals.KRT77.pdf
#$SCRIPT_DIR/plot_pwm.py $PWM_FILE TP53
#$SCRIPT_DIR/plot_pwm.py $PWM_FILE CEBP
#~/git/ggr-project/figs/fig_1.datasets/plot.rna_bar.R $GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz ENSG00000073282 gene_vals.TP63.pdf

# allelic imbalance ATAC to validate importance scores
#$SCRIPT_DIR/fig_3-b.3.prep.py
#$SCRIPT_DIR/fig_3-b.4.plot.R results.w_imptscore.txt results.sig_only.txt

# chipseq plots: TP63, KLF4, ZNF750
# DON'T FORGET UCSC TOOLS
MOTIFS="/srv/scratch/dskim89/ggr/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_present_summary.txt"
SCANMOTIFS="/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-17.footprinting/motifs.input_x_grad.lite/ggr.scanmotifs.h5"
MAPCHAIN_FILE="/srv/scratch/dskim89/ggr/ggr.tronn.2019-07-15.chipseq/hg19ToHg38.over.chain.gz"

BIGWIG_FILE="/srv/scratch/dskim89/ggr/ggr.tronn.2019-04-30.cistrome/cistrome.44759.TP63.hg38.bw"
python /users/dskim89/git/tronn/scripts/ggr/ggr_agg_chipseq_signal.py \
       --motif_list $MOTIFS \
       --data_file $SCANMOTIFS \
       --motif_string TP53 \
       --bigwig_files $BIGWIG_FILE \
       --mapchain $MAPCHAIN_FILE \
       -o chipseq \
       --prefix chipseq

BIGWIG_FILE="/srv/scratch/dskim89/ggr/ggr.tronn.2019-04-30.cistrome/cistrome.48800.KLF4.hg38.bw"
python /users/dskim89/git/tronn/scripts/ggr/ggr_agg_chipseq_signal.py \
       --motif_list $MOTIFS \
       --data_file $SCANMOTIFS \
       --motif_string KLF12 \
       --bigwig_files $BIGWIG_FILE \
       --mapchain $MAPCHAIN_FILE \
       -o chipseq \
       --prefix chipseq


BIGWIG_FILE="/srv/scratch/dskim89/ggr/ggr.tronn.2019-04-30.cistrome/cistrome.48799.ZNF750.hg38.bw"
python /users/dskim89/git/tronn/scripts/ggr/ggr_agg_chipseq_signal.py \
       --motif_list $MOTIFS \
       --data_file $SCANMOTIFS \
       --motif_string TFAP2A \
       --bigwig_files $BIGWIG_FILE \
       --mapchain $MAPCHAIN_FILE \
       -o chipseq \
       --prefix chipseq

# footprinting plots


# pwm x rna (supplement)
PVALS_FILE=$SIG_MOTIFS_DIR/pvals.rna_filt.corr_filt.h5
#$SCRIPT_DIR/fig_3-c.0.plot.pwm_x_rna.R $PVALS_FILE pvals

# motifs and tfs
SUMMARY_DIR=$SIG_MOTIFS_DIR/summary
#$SCRIPT_DIR/fig_3-d.0.plot.motifs_and_tfs.R $SUMMARY_DIR/ggr.pwms_present_summary.txt $SUMMARY_DIR/ggr.pwms_patterns_summary.txt $SUMMARY_DIR/ggr.tfs_corr_summary.txt $SUMMARY_DIR/ggr.tfs_patterns_summary.txt
