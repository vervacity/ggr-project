
# data
SCRIPT_DIR=~/git/ggr-project/figs/fig_4.homotypic
MOTIFS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.background/ggr.scanmotifs.h5
MOTIFS_EXTENDED_FILE=/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-17.footprinting/motifs.input_x_grad.lite/ggr.scanmotifs.h5
SIG_PWMS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_patterns_summary.txt

# genome: multiplicity (use all regions)
#$SCRIPT_DIR/fig_3-b.0.plot.counts.py $MOTIFS_FILE $SIG_PWMS_FILE

# genome: spacing (use all regions)
SIG_PWMS_FILE=ATAC_SIGNALS.NORM.pwms.keep.txt
#$SCRIPT_DIR/fig_3-c.0.analyze_spacing.py $MOTIFS_FILE $SIG_PWMS_FILE

# sims: multiplicity
SIM_MULT_ALLNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.multiplicity.full_negatives
OUT_DIR=sim.multiplicity.full_negatives
mkdir -p $OUT_DIR
#$SCRIPT_DIR/analyze_sims_multiplicity.py $SIM_MULT_ALLNEGS_DIR $OUT_DIR

SIM_MULT_DHSNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.multiplicity.dhs_negatives
OUT_DIR=sim.multiplicity.dhs_negatives
mkdir -p $OUT_DIR
#$SCRIPT_DIR/analyze_sims_multiplicity.py $SIM_MULT_DHSNEGS_DIR $OUT_DIR

# sims: spacing
SIM_SPACING_ALLNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.spacing.full_negatives
OUT_DIR=sim.spacing.full_negatives
mkdir -p $OUT_DIR
#$SCRIPT_DIR/analyze_sims_spacing.py $SIM_SPACING_ALLNEGS_DIR $OUT_DIR

SIM_SPACING_DHSNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.spacing.dhs_negatives
# TODO - run dhs negs version

# genome: multiplicity 1-3 motifs and gene sets (use dynamic regions)
SIG_PWMS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_patterns_summary.txt
TSS_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/annotations/ggr.rna.tss.expressed.bed.gz
BACKGROUND_GENES_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/data/ggr.rna.counts.pc.expressed.mat.txt.gz
EARLY_MOTIFS=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.early/ggr.scanmotifs.h5
MID_MOTIFS=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.mid/ggr.scanmotifs.h5
LATE_MOTIFS=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.late/ggr.scanmotifs.h5

OUT_DIR=multiplicity.gene_sets
mkdir -p $OUT_DIR
$SCRIPT_DIR/multiplicity_to_genes.py $OUT_DIR $SIG_PWMS_FILE $TSS_FILE $BACKGROUND_GENES_FILE $EARLY_MOTIFS $MID_MOTIFS $LATE_MOTIFS
