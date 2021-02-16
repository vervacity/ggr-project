

import os

# commands that use the tronn package



def tronn_intersect_with_expression_cmd(
        pvals_file,
        pwm_file,
        pwm_metadata_file,
        rna_file,
        out_dir,
        nn_motifs=False,
        min_cor=0.75): # adjust to 0.5 for homer runs
    """intersect sig motifs with pwms
    """
    # scanmotifs dir and files
    SCANMOTIFS_DIR = "/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze"
    foreground_main_groups = ["early", "mid", "late"]
    foreground_files = ["{}/motifs.input_x_grad.{}/ggr.scanmotifs.h5".format(SCANMOTIFS_DIR, val)
                        for val in foreground_main_groups]
    
    # params
    use_max_pwm_length = "--max_pwm_length 17 "
    out_dir = "{}/motifs.rna_filt".format(out_dir)

    # cmd
    cmd = "intersect_pwms_and_rna.py "
    cmd += "--dataset_files {} ".format(" ".join(foreground_files))
    cmd += "--pvals_file {} ".format(pvals_file)
    cmd += "--pwm_file {} ".format(pwm_file)
    cmd += "--pwm_metadata_file {} ".format(pwm_metadata_file)
    cmd += "--pwm_scores_key {} ".format("ATAC_SIGNALS.NORM")
    cmd += "--cor_thresh {} ".format(min_cor)
    cmd += "--other_targets ATAC_SIGNALS.NORM=0,2,3,4,5,6,9,10,12 logits=0,2,3,4,5,6,9,10,12 "
    cmd += "--rna_expression_file {} ".format(rna_file)
    cmd += use_max_pwm_length
    cmd += "-o {} ".format(out_dir)
    print cmd
    os.system(cmd)

    # also do the summary here too
    cmd = "ggr_summarize_motifs.py "
    cmd += "--data_file {}/pvals.rna_filt.corr_filt.h5 ".format(out_dir)
    cmd += "-o {}/summary ".format(out_dir)
    cmd += "--prefix ggr"
    print cmd
    os.system(cmd)
    


    
    return
