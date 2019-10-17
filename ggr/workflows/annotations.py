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
    assert inputs.get("gencode_all_gtf") is not None
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
            inputs["gencode_all_gtf"],
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
            inputs["gencode_all_gtf"],
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

    # -------------------------------------------------
    # ANALYSIS 5 - make dataset plot
    # input: dataset matrix
    # output: dataset summary plot
    # -------------------------------------------------
    logger.info("ANALYSIS: plot datasets")
    plot_file = "{}/dataset_summary.pdf".format(out_dir)
    if not os.path.isfile(plot_file):
        # plot
        plot_cmd = "plot.dataset_matrix.R {} {}".format(
            args.dataset_summary, plot_file)
        print plot_cmd
        os.system(plot_cmd)

    # -------------------------------------------------
    # ANALYSIS 6 - make interesting loci bed file
    # input: None
    # output: bed file
    # -------------------------------------------------
    logger.info("ANALYSIS: make interesting regions bed file")
    interesting_bed_file = "{}/interesting_loci.bed.gz".format(out_dir)
    if not os.path.isfile(interesting_bed_file):
        # EDC
        os.system("echo 'chr1\t151896923\t153687874' | gzip -c >> {}".format(interesting_bed_file))
        # KRT locus chr12
        os.system("echo 'chr12\t52490656\t53412661' | gzip -c >> {}".format(interesting_bed_file))
        # KRT locus chr17
        os.system("echo 'chr17\t38756691\t39834792' | gzip -c >> {}".format(interesting_bed_file))
        # MYC
        os.system("echo 'chr8\t128250000\t129250000' | gzip -c >> {}".format(interesting_bed_file))
        # TP63
        os.system("echo 'chr3\t189039644\t190039644' | gzip -c >> {}".format(interesting_bed_file))
        # ZNF750
        os.system("echo 'chr17\t80298130\t81163605' | gzip -c >> {}".format(interesting_bed_file))
        # KLF4
        os.system("echo 'chr9\t109752171\t110752171' | gzip -c >> {}".format(interesting_bed_file))
        # MAF
        os.system("echo 'chr16\t79134765\t80134765' | gzip -c >> {}".format(interesting_bed_file))
        # MAFB
        os.system("echo 'chr20\t38818007\t39818007' | gzip -c >> {}".format(interesting_bed_file))
        # OVOL1
        os.system("echo 'chr11\t65054469\t66054469' | gzip -c >> {}".format(interesting_bed_file))
        # GRHL3
        os.system("echo 'chr1\t24145515\t25145515' | gzip -c >> {}".format(interesting_bed_file))
        # IRF6
        os.system("echo 'chr1\t209479745\t210479745' | gzip -c >> {}".format(interesting_bed_file))
        # EHF
        os.system("echo 'chr11\t34142872\t35142872' | gzip -c >> {}".format(interesting_bed_file))
        # NFKB1
        os.system("echo 'chr4\t102917087\t103917087' | gzip -c >> {}".format(interesting_bed_file))
        # NFKB2
        os.system("echo 'chr10\t103661938\t104661938' | gzip -c >> {}".format(interesting_bed_file))
        # REL
        os.system("echo 'chr2\t60609403\t61609403' | gzip -c >> {}".format(interesting_bed_file))
        # RELA
        os.system("echo 'chr11\t64930583\t65930583' | gzip -c >> {}".format(interesting_bed_file))
        # RELB
        os.system("echo 'chr19\t45005173\t46005173' | gzip -c >> {}".format(interesting_bed_file))
        # RUNX1
        os.system("echo 'chr21\t35916378\t36916378' | gzip -c >> {}".format(interesting_bed_file))
        # RUNX3
        os.system("echo 'chr1\t24792697\t25792697' | gzip -c >> {}".format(interesting_bed_file))
        # TFAP2A
        os.system("echo 'chr6\t9903800\t10903800' | gzip -c >> {}".format(interesting_bed_file))
        # TFAP2B
        os.system("echo 'chr6\t50286761\t51286761' | gzip -c >> {}".format(interesting_bed_file))
        # CEBPA/G
        os.system("echo 'chr19\t33306665\t34306665' | gzip -c >> {}".format(interesting_bed_file))
        # CEBPB
        os.system("echo 'chr20\t48307410\t49307410' | gzip -c >> {}".format(interesting_bed_file))
        # CEBPD
        os.system("echo 'chr8\t48150442\t49150442' | gzip -c >> {}".format(interesting_bed_file))
        # CEBPE
        os.system("echo 'chr14\t23087544\t24087544' | gzip -c >> {}".format(interesting_bed_file))
        # CEBPZ
        os.system("echo 'chr2\t36958919\t37958919' | gzip -c >> {}".format(interesting_bed_file))
        
        
    logger.info("DONE")
    
    return args
