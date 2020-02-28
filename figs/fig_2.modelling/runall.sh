#!/bin/bash

# code
GIT_DIR=~/git/ggr-project
R_DIR=$GIT_DIR/R
SCRIPT_DIR=$GIT_DIR/figs/fig_2.modelling

# dirs
GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a
EVAL_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03
#PREDICT_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05
PREDICT_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-10-08

SIG_MOTIFS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim
PWM_FILE=$GGR_DIR/annotations/HOCOMOCOv11_core_pwms_HUMAN_mono.renamed.nonredundant.txt
ANNOT_DIR=/mnt/lab_data3/dskim89/ggr/annotations

# -----------
# A
# -----------

# MAIN: evaluation metrics for GGR regression
#$SCRIPT_DIR/fig_2-b.2.plot.ggr_regr.R fig_2-b.1 $EVAL_DIR/ggr.basset.regr.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.regr.pretrained.folds.testfold-*/by_task/*txt

# SUPPL: evaluation for pretrain -> classification
#$SCRIPT_DIR/fig_2-b.0.plot.encode_clf.R fig_2-b.0 $EVAL_DIR/encode-roadmap.basset.clf.testfold-*/by_task/*txt
#$SCRIPT_DIR/fig_2-b.1.plot.ggr_clf.R fig_2-b.1 $EVAL_DIR/ggr.basset.clf.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.clf.pretrained.folds.testfold-*/by_task/*txt

# -----------
# B
# -----------

# MAIN: scatter plots for predicted/actual d0,3,6
#$SCRIPT_DIR/fig_2-c.0.plot.scatter.R $EVAL_DIR/ggr.basset.regr.pretrained.folds.testfold-4/by_task/correlation

# SUPPL: evaluation for prediction - heatmaps
#$SCRIPT_DIR/fig_2-d.0.plot.heatmaps.R $PREDICT_DIR/motifs.input_x_grad.early/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.mid/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.late/ggr.scanmotifs.h5
#$SCRIPT_DIR/fig_2-d.0.plot.heatmaps.R $PREDICT_DIR/motifs.input_x_grad.dynamic.early/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.dynamic.mid/ggr.scanmotifs.h5 $PREDICT_DIR/motifs.input_x_grad.dynamic.late/ggr.scanmotifs.h5

# TF results from ENCODE data
#$SCRIPT_DIR/fig_S2-a.0.plot.encode_clf.TFs.R fig_S2-a $EVAL_DIR/encode-roadmap.basset.clf.testfold-*/by_task/*txt

# clf: trajectories
#$SCRIPT_DIR/fig_S2-b.0.plot.ggr_clf.traj.R fig_S2-b $EVAL_DIR/ggr.basset.clf.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.clf.pretrained.folds.testfold-*/by_task/*txt

# clf: histones
#$SCRIPT_DIR/fig_S2-b.1.plot.ggr_clf.histones.R fig_S2-b $EVAL_DIR/ggr.basset.clf.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.clf.pretrained.folds.testfold-*/by_task/*txt

# regr: histones
#$SCRIPT_DIR/fig_S2-c.0.plot.ggr_regr.histones.R fig_S2-c $EVAL_DIR/ggr.basset.regr.random_init.testfold-*/by_task/*txt $EVAL_DIR/ggr.basset.regr.pretrained.folds.testfold-*/by_task/*txt

# -----------
# C
# -----------

# MAIN: vignette
#$SCRIPT_DIR/choose_vignette.py
RNA_MAT=$GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz
#$GIT_DIR/figs/fig_1.datasets/plot.rna_bar.R $RNA_MAT ENSG00000189182 gene_vals.KRT77.pdf
#$GIT_DIR/figs/fig_1.datasets/plot.rna_bar.R $RNA_MAT ENSG00000073282 gene_vals.TP63.pdf
#$SCRIPT_DIR/plot_pwm.py $PWM_FILE TP53
#$SCRIPT_DIR/plot_pwm.py $PWM_FILE CEBP

# SUPPL: allelic imbalance ATAC analysis
#$SCRIPT_DIR/fig_3-b.3.prep.py
#$SCRIPT_DIR/fig_3-b.4.plot.R results.w_imptscore.txt results.sig_only.txt

# -----------
# D
# -----------

# Schematic of NeMo

# -----------
# E/G
# -----------

SUMMARY_DIR=$SIG_MOTIFS_DIR/summary
#$SCRIPT_DIR/fig_3-d.0.plot.motifs_and_tfs.R $SUMMARY_DIR/ggr.pwms_present_summary.txt $SUMMARY_DIR/ggr.pwms_patterns_summary.txt $SUMMARY_DIR/ggr.tfs_corr_summary.txt $SUMMARY_DIR/ggr.tfs_patterns_summary.txt

# -----------
# F
# -----------

TRONN_DIR=~/git/tronn

# MAIN: chipseq plots: TP63, KLF4, ZNF750
# DON'T FORGET UCSC TOOLS
MOTIFS=$SIG_MOTIFS_DIR/summary/ggr.pwms_present_summary.txt
MAPCHAIN_FILE=$ANNOT_DIR/hg19ToHg38.over.chain.gz
SCANMOTIFS=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/motifs.input_x_grad.lite/ggr.scanmotifs.h5
CISTROME_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/cistrome

# TP63
BIGWIG_FILE=$CISTROME_DIR/cistrome.44759.TP63.hg38.bw
#python $TRONN_DIR/scripts/ggr/ggr_agg_chipseq_signal.py --motif_list $MOTIFS --data_file $SCANMOTIFS --motif_string TP53 --bigwig_files $BIGWIG_FILE --mapchain $MAPCHAIN_FILE -o chipseq --prefix chipseq

BIGWIG_FILE=$CISTROME_DIR/cistrome.48800.KLF4.hg38.bw
#python $TRONN_DIR/scripts/ggr/ggr_agg_chipseq_signal.py --motif_list $MOTIFS --data_file $SCANMOTIFS --motif_string KLF12 --bigwig_files $BIGWIG_FILE --mapchain $MAPCHAIN_FILE -o chipseq --prefix chipseq

BIGWIG_FILE=$CISTROME_DIR/cistrome.48799.ZNF750.hg38.bw
#python $TRONN_DIR/scripts/ggr/ggr_agg_chipseq_signal.py --motif_list $MOTIFS --data_file $SCANMOTIFS --motif_string TFAP2A --bigwig_files $BIGWIG_FILE --mapchain $MAPCHAIN_FILE -o chipseq --prefix chipseq

# SUPPL: footprinting plots
# run two ways: 1) GC content matched, this gives the flanking accessibility diff between sites. 2) accessibility matched, this gives diff in footprint depth
D0_BAM="/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2019-06-04.bams_bp-resolution/primary_keratinocyte-d00.GGR.Stanford_Greenleaf.ATAC-seq.pooled.fixedtrim.PE2SE.nodup.bam"
D3_BAM="/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2019-06-04.bams_bp-resolution/primary_keratinocyte-d30.GGR.Stanford_Greenleaf.ATAC-seq.pooled.fixedtrim.PE2SE.nodup.bam"
D6_BAM="/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2019-06-04.bams_bp-resolution/primary_keratinocyte-d60.GGR.Stanford_Greenleaf.ATAC-seq.pooled.fixedtrim.PE2SE.nodup.bam"

# first pass to look at all of them
#python /users/dskim89/git/tronn/scripts/ggr/ggr_footprint_motifs.py --motif_list $MOTIFS --data_file $SCANMOTIFS --bam_files $D0_BAM $D3_BAM $D6_BAM -o footprints --prefix ggr --signal_matching_key sequence.active.gc_fract

#python /users/dskim89/git/tronn/scripts/ggr/ggr_footprint_motifs.py --motif_list $MOTIFS --data_file $SCANMOTIFS --bam_files $D0_BAM $D3_BAM $D6_BAM -o footprints --prefix ggr --filter_motifs HOXA1 CEBPD NFAT GRHL DLX3 --signal_matching_key sequence.active.gc_fract

# use if needed
##python /users/dskim89/git/tronn/scripts/ggr/ggr_footprint_motifs.py --motif_list $MOTIFS --data_file $SCANMOTIFS --bam_files $D0_BAM $D3_BAM $D6_BAM -o footprints --prefix ggr --filter_motifs NFKB NFAT GRHL --signal_matching_key ATAC_SIGNALS.NORM
