
# data
SCRIPT_DIR=~/git/ggr-project/figs/fig_4.homotypic
MOTIFS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.background/ggr.scanmotifs.h5
MOTIFS_EXTENDED_FILE=/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-17.footprinting/motifs.input_x_grad.lite/ggr.scanmotifs.h5
SIG_PWMS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_patterns_summary.txt

# genome: multiplicity
#$SCRIPT_DIR/fig_3-b.0.plot.counts.py $MOTIFS_FILE $SIG_PWMS_FILE

# genome: spacing
SIG_PWMS_FILE=ATAC_SIGNALS.NORM.pwms.keep.txt
#$SCRIPT_DIR/fig_3-c.0.analyze_spacing.py $MOTIFS_FILE $SIG_PWMS_FILE

# sim data
SIM_MULT_ALLNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.multiplicity.full_negatives
#SIM_MULT_ALLNEGS_DIR=/srv/scratch/dskim89/ggr/ggr.tronn.2019-08-30.motif_density/

# sims: multiplicity
OUT_DIR=sim.multiplicity.full_negatives
mkdir -p $OUT_DIR
$SCRIPT_DIR/analyze_density.py $SIM_MULT_ALLNEGS_DIR $OUT_DIR


# sims: spacing
