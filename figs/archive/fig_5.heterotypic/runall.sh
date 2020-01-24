
MOTIFS_FILE=/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-17.footprinting/motifs.input_x_grad.lite/ggr.scanmotifs.h5
GRAMMAR_FILE=/srv/scratch/dskim89/ggr/inference.2019-03-12/dmim.shuffle/grammars.annotated.manual_filt.merged.final/grammars_summary.txt
SCRIPT_DIR=~/git/ggr-project/figs/fig_5.heterotypic

# make figs
$SCRIPT_DIR/fig_5-b.0.genomic_instances.R $MOTIFS_FILE $GRAMMAR_FILE
