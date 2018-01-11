# Description: workflows for annotations

import os
import logging

from ggr.util.utils import run_shell_cmd

from ggr.analyses.annotations import get_proteincoding_gene_list_from_gtf
from ggr.analyses.annotations import get_proteincoding_tss_bed_from_gtf


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
    outputs[tss_pc_key] = "{}/hg19.tss.gencode19.bed.gz".format(
        out_dir, prefix)
    if not os.path.isfile(outputs[tss_pc_key]):
        get_proteincoding_tss_bed_from_gtf(
            inputs["gencode_proteincoding_gtf"],
            outputs[tss_pc_key])

    # add outputs and finish
    args.outputs["annotations"] = outputs
    logger.info("DONE")
    
    return args
