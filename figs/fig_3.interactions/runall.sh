#!/bin/bash

SCRIPT_DIR=/users/dskim89/git/ggr-project/figs/fig_3.interactions
FIG2_SCRIPT_DIR=/users/dskim89/git/ggr-project/figs/fig_2.modelling
TRONN_R=/users/dskim89/git/tronn/R

GGR_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a
ANNOT_DIR=/mnt/lab_data3/dskim89/ggr/annotations

# -------------
# A
# -------------

# schematic

# -------------
# BCD
# -------------

# MAIN: pairwise interactions
WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle.complete
endogenous_file=$WORK_DIR/synergy_calcs.endog/interactions.summary.txt
simulations_file=$WORK_DIR/synergy_calcs.sims/interactions.summary.txt
$SCRIPT_DIR/fig_4-b.0.plot.interactions.R $endogenous_file $simulations_file

# SUPPL: pairwise combinations, using various signal sources
ordering_file=ordering.txt # from prev script run, apply ONLY to mut matrices
summary_files=/mnt/lab_data3/dskim89/ggr/nn/2020-01-13/grammars*/grammars.annotated/grammar_summary.txt
#$TRONN_R/plot.combinations.adjacencies.R fig_S3.unfiltered FALSE $summary_files
#$TRONN_R/plot.combinations.adjacencies.R fig_S3.go_filt TRUE $summary_files

# -------------
# E
# -------------

# data
synergy_dir=$WORK_DIR/synergy
#summary_file=$WORK_DIR/grammars.annotated.manual_filt.merged.final/grammars_summary.txt
summary_file=/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results/results/combinatorial_rules/summary.txt.gz # NOTE: utilizing only VALIDATED RULES for vignettes
#links_file=$GGR_DIR/results/linking/hichip/ggr.linking.ALL.overlap.interactions.txt.gz
links_file=$GGR_DIR/results/linking/proximity/ggr.linking.ALL.overlap.interactions.txt.gz
tss_file=$GGR_DIR/annotations/ggr.rna.tss.expressed.bed.gz
filter_genes_file=$GGR_DIR/annotations/hg19.genodermatoses.ensembl_ids.txt.gz
#filter_genes_file=$GGR_DIR/data/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat.txt.gz
chromsizes=$ANNOT_DIR/hg19.chrom.sizes
gene_id_mapping_file=$GGR_DIR/annotations/hg19.ensembl_geneids.pc.gencode19.mappings.mat.gz

pwm_file=$GGR_DIR/annotations/HOCOMOCOv11_core_pwms_HUMAN_mono.renamed.nonredundant.txt

# vignettes that are linked to genodermatoses genes
#rm -r vignettes.genodermatoses # REMOVE THIS LATER
#$SCRIPT_DIR/select_vignettes.py fig_3.vignettes $summary_file $synergy_dir $links_file $tss_file $filter_genes_file $chromsizes vignettes.genodermatoses
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file NFKB
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file CBFB
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file ATF4
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file GRHL
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file CEBP
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file TP53
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file ATF1
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file TAF1
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file TFAP2
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file NR2C1
#$FIG2_SCRIPT_DIR/plot_pwm.py $pwm_file MTF1

# vignettes with GWAS/GTEx
VARIANT_DIR=/mnt/lab_data/kundaje/users/dskim89/ggr/variants/validation.2019-04-27.processed

# gwas - note none found
filter_genes_file=$GGR_DIR/data/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat.txt.gz
#$SCRIPT_DIR/select_vignettes.GWAS_GTEx.py $VARIANT_DIR/gwas.skin.atac_filt.bed.gz vignettes.gwas $filter_genes_file $gene_id_mapping_file $synergy_dir $summary_file

# gtex - note 1 found
filter_genes_file=$GGR_DIR/data/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat.txt.gz
#$SCRIPT_DIR/select_vignettes.GWAS_GTEx.py $VARIANT_DIR/gtex.skin.atac_filt.bed.gz vignettes.gtex $filter_genes_file $gene_id_mapping_file $synergy_dir $summary_file


# SUPPL: predicted rule map
WORK_DIR=/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle/grammars.annotated.manual_filt.merged.final

# SUPPL: vignette of region with GWAS hit?


# just copy over
