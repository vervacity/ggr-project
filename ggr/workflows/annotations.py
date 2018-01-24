# Description: workflows for annotations

import os
import logging

from ggr.util.utils import run_shell_cmd

from ggr.analyses.annotations import get_proteincoding_gene_list_from_gtf
from ggr.analyses.annotations import get_proteincoding_tss_bed_from_gtf
from ggr.analyses.annotations import get_ensembl_to_geneid_mapping

from ggr.analyses.motifs import add_hocomoco_metadata
from ggr.analyses.motifs import reduce_pwm_redundancy


def setup_hocomoco_motifs_workflow(args, prefix):
    """Run workflow to set up hocomoco motif set
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)

    inputs = args.inputs["annot"][args.cluster]
    outputs = args.outputs["annotations"]
    out_dir = outputs["dir"]

    # assertions
    assert inputs.get("hocomoco_pwms") is not None
    assert inputs.get("hocomoco_metadata") is not None
    assert inputs.get("custom_pwms") is not None
    assert outputs.get("geneids.mappings.mat") is not None

    # -------------------------------------------------
    # ANALYSIS 0 - add metadata to get ensembl ids into PWM names
    # input: pwm file and metadata file
    # output: new pwm file with ensembl ids in name
    # -------------------------------------------------
    renamed_pwms_key = "pwms.renamed"
    outputs[renamed_pwms_key] = "{}/{}.renamed.txt".format(
        out_dir,
        os.path.basename(inputs["hocomoco_pwms"]).split(".txt")[0])
    if not os.path.isfile(outputs[renamed_pwms_key]):
        add_hocomoco_metadata(
            inputs["hocomoco_pwms"],
            outputs[renamed_pwms_key],
            inputs["hocomoco_metadata"],
            outputs["geneids.mappings.mat"])
        
    # -------------------------------------------------
    # ANALYSIS 1 - reduce redundancy using a hclust method (see RSAT matrix paper)
    # input: pwm file 
    # output: new pwm file with new metadata file
    # -------------------------------------------------
    pwm_dir = "{}/pwms".format(out_dir)
    run_shell_cmd("mkdir -p {}".format(pwm_dir))
    
    nonredundant_pwms_key = "{}.nonredundant".format(renamed_pwms_key)
    outputs[nonredundant_pwms_key] = "{}.nonredundant.txt".format(
        outputs[renamed_pwms_key].split(".txt")[0])

    nonredundant_pwms_metadata_key = "pwms.metadata.nonredundant".format(renamed_pwms_key)
    outputs[nonredundant_pwms_metadata_key] = "{}/{}.nonredundant.txt".format(
        out_dir, 
        os.path.basename(inputs["hocomoco_metadata"]).split(".tsv")[0])
    if not os.path.isfile(outputs[nonredundant_pwms_metadata_key]):
        reduce_pwm_redundancy(
            [(outputs[renamed_pwms_key], "log_likelihood"),
             (inputs["custom_pwms"], "probability")],
            outputs[nonredundant_pwms_key],
            outputs[nonredundant_pwms_metadata_key],
            tmp_prefix="{}/hocomoco".format(pwm_dir),
            num_threads=28)

    # later add in RNA expression information
    
    return args


def runall(args, prefix):
    """Run all workflows to generate important annotation files
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("MASTER_WORKFLOW: set up annotations")

    inputs = args.inputs["annot"][args.cluster]
    outputs = args.outputs["annotations"]
    
    # assertions
    assert inputs.get("gencode_proteincoding_gtf") is not None
    assert inputs.get("fantom5_transcriptionfactor_csv") is not None
    assert outputs.get("dir") is not None

    # set up dirs
    out_dir = outputs["dir"]
    run_shell_cmd("mkdir -p {}".format(out_dir))
    logging.debug("outputs going to {}".format(out_dir))
    
    # -------------------------------------------------
    # ANALYSIS 0 - get list of protein coding genes
    # input: GTF file
    # output: text file with protein coding gene ids
    # -------------------------------------------------
    logger.info("ANALYSIS: make protein coding list")
    geneids_pc_key = "geneids.pc.list"
    outputs[geneids_pc_key] = "{}/{}.ensembl_geneids.pc.gencode19.txt.gz".format(
        out_dir, prefix)
    if not os.path.isfile(outputs[geneids_pc_key]):
        get_proteincoding_gene_list_from_gtf(
            inputs["gencode_proteincoding_gtf"],
            outputs[geneids_pc_key])
        
    # -------------------------------------------------
    # ANALYSIS 1 - get list of TFs
    # input: transcription factor list (currently FANTOM5 list)
    # output: transcription factor list (converted to ensembl IDs)
    # -------------------------------------------------
    logger.info("ANALYSIS: get transcription factor list")
    geneids_tf_key = "geneids.tf.list"
    outputs[geneids_tf_key] = "{}/{}.ensembl_geneids.tf.fantom5.txt.gz".format(
        out_dir, prefix)
    if not os.path.isfile(outputs[geneids_tf_key]):
        convert_entrez_to_ensembl = "annot.entrez_to_ensembl.R {0} {1}".format(
            inputs["fantom5_transcriptionfactor_csv"],
            outputs[geneids_tf_key])
        run_shell_cmd(convert_entrez_to_ensembl)

    # -------------------------------------------------
    # ANALYSIS 2 - make a TSS file
    # input: GTF file
    # output: BED file
    # -------------------------------------------------
    logger.info("ANALYSIS: make a TSS file")
    tss_pc_key = "tss.pc.bed"
    outputs[tss_pc_key] = "{}/{}.tss.gencode19.bed.gz".format(
        out_dir, prefix)
    if not os.path.isfile(outputs[tss_pc_key]):
        get_proteincoding_tss_bed_from_gtf(
            inputs["gencode_proteincoding_gtf"],
            outputs[tss_pc_key])

    # -------------------------------------------------
    # ANALYSIS 3 - make an ensembl <-> HGNC mapping table
    # input: ensembl gene ids
    # output: conversion table
    # -------------------------------------------------
    logger.info("ANALYSIS: make a ensembl to HGNC mapping table")
    mapping_mat_key = "geneids.mappings.mat"
    outputs[mapping_mat_key] = "{}/{}.ensembl_geneids.pc.gencode19.mappings.mat.gz".format(
        out_dir, prefix)
    if not os.path.isfile(outputs[mapping_mat_key]):
        # make an analysis that produces this table
        get_ensembl_to_geneid_mapping(
            outputs[geneids_pc_key],
            outputs[mapping_mat_key])

    # -------------------------------------------------
    # ANALYSIS 4 - set up motif file with reduction
    # input: pwm file + metadata
    # output: reduced pwm file
    # -------------------------------------------------
    logger.info("ANALYSIS: generate a nonredundant PWM file")
    args = setup_hocomoco_motifs_workflow(args, prefix)

    logger.info("DONE")
    
    return args
