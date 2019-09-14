#!/bin/bash

# code
SCRIPT_DIR=~/git/ggr-project/figs/fig_4.homotypic

# dirs
GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a
SCANMOTIFS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05
INFER_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12

# files
SIG_PWMS_FILE=$INFER_DIR/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_patterns_summary.txt

# =================================================
# GENOME: HINTS TO MULTIPLICITY/SPACING
# =================================================

# data
UNFILT_SCAN_FILE=$SCANMOTIFS_DIR/motifs.input_x_grad.background/ggr.scanmotifs.h5

# multiplicity as seen in the genome (NN-active hits)
#$SCRIPT_DIR/analyze.genome.multiplicity.py $SCRIPT_DIR genome.multiplicity $UNFILT_SCAN_FILE $SIG_PWMS_FILE

# spacing as seen in the genome (NN-active hits)
#$SCRIPT_DIR/analyze.genome.spacing.py $SCRIPT_DIR genome.spacing $UNFILT_SCAN_FILE $SIG_PWMS_FILE

# =================================================
# SIMULATIONS: ISOLATING SINGLE MOTIF EFFECTS
# =================================================

# data
SIM_MULT_ALLNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.multiplicity.full_negatives
SIM_MULT_DHSNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.multiplicity.dhs_negatives
SIM_SPACING_ALLNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.spacing.full_negatives
SIM_SPACING_DHSNEGS_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/analysis_homotypic.2019-09-06/sims.spacing.dhs_negatives

# multiplicity as seen in simulations
#$SCRIPT_DIR/analyze.simulations.multiplicity.py $SCRIPT_DIR $SIM_MULT_ALLNEGS_DIR simulations.multiplicity.full_negs
#$SCRIPT_DIR/analyze.simulations.multiplicity.py $SCRIPT_DIR $SIM_MULT_DHSNEGS_DIR simulations.multiplicity.dhs_negs

# spacing as seen in simulations
$SCRIPT_DIR/analyze.simulations.spacing.py $SCRIPT_DIR $SIM_SPACING_ALLNEGS_DIR simulations.spacing.full_negs
#$SCRIPT_DIR/analyze.simulations.spacing.py $SCRIPT_DIR $SIM_SPACING_DHSNEGS_DIR simulations.spacing.dhs_negs

# =================================================
# BACK TO GENOME: USE MULTIPLICITY/SPACING FOR FUNCTION
# =================================================

# data
ATAC_MASTER_BED=$GGR_DIR/data/ggr.atac.idr.master.bed.gz
TSS_BED=$GGR_DIR/annotations/ggr.rna.tss.expressed.bed.gz
ATAC_SIGNAL_MAT=$GGR_DIR/data/ggr.atac.ends.counts.pooled.rlog.mat.txt.gz
RNA_SIGNAL_MAT=$GGR_DIR/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz

# set up links
OUT_DIR=linking
#$SCRIPT_DIR/setup_links.py $ATAC_MASTER_BED $TSS_BED $ATAC_SIGNAL_MAT $RNA_SIGNAL_MAT $OUT_DIR
LINKS_FILE=$OUT_DIR/links.bed.gz

# genome: multiplicity 1-3 motifs and gene sets (use dynamic regions)
BACKGROUND_GENES_FILE=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/data/ggr.rna.counts.pc.expressed.mat.txt.gz
EARLY_MOTIFS=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.early/ggr.scanmotifs.h5
MID_MOTIFS=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.mid/ggr.scanmotifs.h5
LATE_MOTIFS=/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.late/ggr.scanmotifs.h5

MULTIPLICITY_GENE_DIR=multiplicity.gene_sets
#$SCRIPT_DIR/multiplicity_to_genes.py $MULTIPLICITY_GENE_DIR $SIG_PWMS_FILE $LINKS_FILE $BACKGROUND_GENES_FILE $EARLY_MOTIFS $MID_MOTIFS $LATE_MOTIFS

SPACING_GENE_DIR=spacing.gene_sets
#$SCRIPT_DIR/spacing_to_genes.py $SPACING_GENE_DIR $SIG_PWMS_FILE $LINKS_FILE $BACKGROUND_GENES_FILE $EARLY_MOTIFS $MID_MOTIFS $LATE_MOTIFS

# TODO merge results
# build a manual file to aggregate information
OUT_DIR=enrichment.summaries
#$SCRIPT_DIR/summarize_enrichment_results.py $OUT_DIR $SIG_PWMS_FILE $MULTIPLICITY_GENE_DIR $SPACING_GENE_DIR
