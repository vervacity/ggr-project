

import os
import re
import gzip

import numpy as np
import pandas as pd


def trim_and_attach_umi(fastq1, fastq2, fastq_out):
    """set up dependent, here the UMI is the first 15 bp of read1, but we want the
    20 bp after 17 bps in read 2
    """
    with gzip.open(fastq1, "r") as fp1:
        with gzip.open(fastq2, "r") as fp2:
            with gzip.open(fastq_out, "w") as out:
                
                reads_seen = 0
                while True:

                    # progress
                    if reads_seen % 1000000 == 0:
                        print reads_seen

                    # get a read pair
                    try:
                        read_1 = []
                        for i in range(4):
                            read_1.append(fp1.readline().strip())

                        if read_1[0] == "":
                            break

                        read_2 = []
                        for i in range(4):
                            read_2.append(fp2.readline().strip())

                        if read_2[0] == "":
                            break

                    except:
                        break
                    
                    # find the umi (read 1)
                    umi = read_1[1][:15]

                    # trim lines
                    read_2[0] = read_2[0].split()[0] + ":" + umi
                    read_2[1] = read_2[1][17:]
                    read_2[3] = read_2[3][17:]

                    # write out
                    out.write("{}\n".format("\n".join(read_2)))

                    # continue
                    reads_seen += 1
                    
    return


def dedup_and_count(sam_file, dedup_file):
    """go through sorted sam file and get deduped counts
    """
    umis = set()
    barcode = ""
    count = 0
    with open(sam_file, "r") as fp:
        with gzip.open(dedup_file, "w") as out:

            for line in fp:

                # don't read header
                if line.startswith("@"):
                    continue

                # read in important details
                fields = line.split("\t")
                line_barcode = fields[2]
                umi = fields[0].split(":")[-1]
                
                # check if current barcode or not
                if line_barcode != barcode:
                    # save out latest
                    if barcode != "":
                        result = "{}\t{}\n".format(barcode, count)
                        out.write(result)
                    
                    # and reset variables
                    umis = set([umi])
                    barcode = line_barcode
                    count = 1
                else:
                    # check if umi present or not
                    if umi in umis:
                        continue
                    # add to count
                    umis.add(umi)
                    count += 1
    
    return


def fastq_to_counts(fastq_file, ref_fasta_file, prefix):
    """take a fastq file through to deduped counts
    """
    
    # align
    out_sai = "{}.sai".format(prefix)
    align = "bwa aln -t 24 {} {} > {}".format(ref_fasta_file, fastq_file, out_sai)
    if not os.path.isfile(out_sai):
        print align
        os.system(align)

    out_sam = "{}.sam".format(prefix)
    samse = "bwa samse {} {} {} > {}".format(ref_fasta_file, out_sai, fastq_file, out_sam)
    if not os.path.isfile(out_sam):
        print samse
        os.system(samse)

    # filter/sort
    filter_sam = "{}.filt.sam".format(prefix)
    filter_fn = "samtools view -h -F 4 {} | samtools sort -@24 -o {} -".format(
        out_sam, filter_sam)
    if not os.path.isfile(filter_sam):
        print filter_fn
        os.system(filter_fn)

    # get counts
    counts_file = "{}.dedup.counts.txt.gz".format(prefix)
    if not os.path.isfile(counts_file):
        dedup_and_count(filter_sam, counts_file)

    return None


def counts_to_mat(counts_files, out_prefix, metadata_file=None):
    """read in all files and merge to matrix. produces:
    - count matrix (no metadata)
    - redistributed count  matrix (no metadata)
    - normalized signal matrix (no metadata)
    - normalized signal matrix (metadata)
    """
    mat_data = None
    for key in sorted(counts_files.keys()):
        print key
        count_data = pd.read_csv(counts_files[key], sep="\t", header=None, index_col=0)
        count_data.columns = [key]
        if mat_data is None:
            mat_data = count_data.copy()
        else:
            mat_data = mat_data.merge(
                count_data, how="outer",
                left_index=True, right_index=True)

    # count matrix
    count_mat_file = "{}.counts.raw.mat.txt.gz".format(out_prefix)
    mat_data.to_csv(count_mat_file, sep="\t", compression="gzip", header=True, index=True)

    # normalize: divide by read total (column norm)
    mat_data_norm = mat_data.fillna(0)
    mat_data_norm = mat_data_norm / np.sum(mat_data_norm, axis=0)

    # normalize: multiply by plasmid fract
    plasmid_fracts = mat_data_norm["plasmid"]
    mat_data_norm = mat_data_norm.drop("plasmid", axis=1)
    mat_data_norm = mat_data_norm.multiply(plasmid_fracts.values, axis=0)
    
    # renormalize: make probability (sum to 1) then multiply by total read count
    mat_data_norm = mat_data_norm / np.sum(mat_data_norm, axis=0)
    mat_data_norm = mat_data_norm.multiply(
        np.sum(mat_data.drop("plasmid", axis=1), axis=0), axis=1)
    
    norm_mat_file = "{}.counts.norm.mat.txt.gz".format(out_prefix)
    mat_data_norm.to_csv(norm_mat_file, sep="\t", compression="gzip", header=True, index=True)

    # and make pseudo TPM (log2)
    norm_signal_mat_file = "{}.signal.mat.txt.gz".format(out_prefix)
    signal_mat = mat_data_norm / np.sum(mat_data_norm, axis=0)
    signal_mat = signal_mat * 1000000
    signal_mat = np.log2(signal_mat)
    signal_mat[signal_mat < 0] = 0
    signal_mat.to_csv(
        norm_signal_mat_file, sep="\t", compression="gzip", header=True, index=True)
    
    # attach metadata
    if metadata_file is not None:
        metadata = pd.read_csv(metadata_file, header=0, sep="\t", index_col=0)
        signal_mat["unique_id"] = [
            re.sub("\\|.+", "", id_string) for id_string in signal_mat.index.values]
        mat_data_w_metadata = signal_mat.merge(metadata, on="unique_id")
        w_metadata_file = "{}.signal.w_metadata.mat.txt.gz".format(out_prefix)
        mat_data_w_metadata.to_csv(
            w_metadata_file, sep="\t",
            compression="gzip", header=True, index=False)
    
    return None


def test_grammar_trajectories(data_file, out_dir, drop_rep="b4"):
    """take in data file and look at values
    """
    # load data
    data = pd.read_csv(data_file, sep="\t", header=0)
    data["grammar"] = [re.sub("\\.\\d+$", "", name) for name in data["example_id"]]

    # filter? if so filter for just synergy set (remove distant and non-sig)
    data["synergy.max"] = data["synergy.scores.diff.sig.0"].apply(
        lambda x: np.max([int(val) for val in x.split(",")]))
    data = data[data["synergy.max"] != 0]
    
    # which headers to keep
    signal_headers = [name for name in data.columns if "_b" in name]
    if drop_rep is not None:
        signal_headers = [name for name in signal_headers if drop_rep not in name]
    metadata_headers = [
        "example_id",
        "combos"]
    
    # get hypotheses
    grammars = list(set(data["grammar"].values.tolist()))
    grammars = sorted([name for name in grammars if "ggr" in name])

    # go through hypotheses
    for grammar in grammars:
        print grammar
        grammar_prefix = "{}/{}".format(out_dir, grammar)
        
        # extract relevant data
        grammar_data = data[data["grammar"] == grammar]
        grammar_data = grammar_data[signal_headers+metadata_headers]

        # remove low signal barcodes
        signal_thresh = 0.5
        grammar_data = grammar_data.set_index(metadata_headers)
        grammar_data = grammar_data[grammar_data.max(axis=1) > signal_thresh]
        grammar_data = grammar_data.reset_index()
        
        # average the barcodes
        grammar_data_avg = grammar_data.groupby(metadata_headers).mean()

        # normalize by example?
        # hack for now
        grammar_data_avg["d3_b1"] = grammar_data_avg["d3_b1"] - grammar_data_avg["d0_b1"]
        grammar_data_avg["d6_b1"] = grammar_data_avg["d6_b1"] - grammar_data_avg["d0_b1"]
        grammar_data_avg["d0_b1"] = grammar_data_avg["d0_b1"] - grammar_data_avg["d0_b1"]

        grammar_data_avg["d3_b2"] = grammar_data_avg["d3_b2"] - grammar_data_avg["d0_b2"]
        grammar_data_avg["d6_b2"] = grammar_data_avg["d6_b2"] - grammar_data_avg["d0_b2"]
        grammar_data_avg["d0_b2"] = grammar_data_avg["d0_b2"] - grammar_data_avg["d0_b2"]

        grammar_data_avg["d3_b3"] = grammar_data_avg["d3_b3"] - grammar_data_avg["d0_b3"]
        grammar_data_avg["d6_b3"] = grammar_data_avg["d6_b3"] - grammar_data_avg["d0_b3"]
        grammar_data_avg["d0_b3"] = grammar_data_avg["d0_b3"] - grammar_data_avg["d0_b3"]

        # average across reps here
        days = ["d0", "d3", "d6"]
        for day in days:
            day_headers = [header for header in signal_headers if day in header]
            grammar_data_avg["{}_AVG".format(day)] = grammar_data_avg[day_headers].mean(axis=1)
        grammar_data_avg = grammar_data_avg.drop(signal_headers, axis=1)
        
        # reset index
        grammar_data_avg = grammar_data_avg.reset_index()
        
        # and only keep the 0,0
        grammar_data_avg = grammar_data_avg[grammar_data_avg["combos"] == "0,0"]
        grammar_data_avg = grammar_data_avg.drop("combos", axis=1)
        grammar_data_avg = grammar_data_avg.set_index("example_id")
        
        # save this out to plot in R
        plot_data_file = "{}.agg_data.txt.gz".format(grammar_prefix)
        grammar_data_avg.to_csv(
            plot_data_file, sep="\t", compression="gzip", header=True, index=True)

        # plot cmd
        plot_file = "{}.endogenous.traj.pdf".format(grammar_prefix)
        plot_cmd = "./plot.mpra.hypothesis.traj.R {} {}".format(
            plot_data_file, plot_file)
        print plot_cmd
        os.system(plot_cmd)

    return


def test_mutational_hypotheses(data_file, out_dir, drop_rep="b4"):
    """test mutational hypotheses
    """
    # load data, make grammar col
    data = pd.read_csv(data_file, sep="\t", header=0)
    data["grammar"] = [re.sub("\\.\\d+$", "", name) for name in data["example_id"]]
    
    # filter? if so filter for just synergy set (remove distant and non-sig)
    data["synergy.max"] = data["synergy.scores.diff.sig.0"].apply(
        lambda x: np.max([int(val) for val in x.split(",")]))
    data = data[data["synergy.max"] != 0]
    
    # which headers to keep
    signal_headers = [name for name in data.columns if "_b" in name]
    if drop_rep is not None:
        signal_headers = [name for name in signal_headers if drop_rep not in name]
    metadata_headers = [
        "example_id",
        "combos"]
    
    # get hypotheses
    grammars = list(set(data["grammar"].values.tolist()))
    grammars = sorted([name for name in grammars if "ggr" in name])

    # go through hypotheses
    for grammar in grammars:
        print grammar
        grammar_prefix = "{}/{}".format(out_dir, grammar)
        
        # extract relevant data
        grammar_data = data[data["grammar"] == grammar]
        grammar_data = grammar_data[signal_headers+metadata_headers]

        # remove low signal barcodes
        signal_thresh = 0.5
        grammar_data = grammar_data.set_index(metadata_headers)
        grammar_data = grammar_data[grammar_data.max(axis=1) > signal_thresh]
        grammar_data = grammar_data.reset_index()
        
        # average the barcodes
        grammar_data_avg = grammar_data.groupby(metadata_headers).mean()

        # TODO need to average across reps HERE
        days = ["d0", "d3", "d6"]
        for day in days:
            day_headers = [header for header in signal_headers if day in header]
            grammar_data_avg["{}_AVG".format(day)] = grammar_data_avg[day_headers].mean(axis=1)
        grammar_data_avg = grammar_data_avg.drop(signal_headers, axis=1)
        
        # reset index
        grammar_data_avg = grammar_data_avg.reset_index()
        
        # then per example need to normalize
        # per example subtract relative to background
        examples = sorted(list(set(grammar_data_avg["example_id"].values.tolist())))
        grammar_data_norm = []
        for example in examples:
            # get example data
            example_data = grammar_data_avg[grammar_data_avg["example_id"] == example]
            
            # check that all necessary examples exist
            if example_data[example_data["combos"] != "-1"].shape[0] != 4:
                continue
            
            # extract 1,1 (background) and subtract
            tmp_data = example_data.drop("example_id", axis=1)
            background_vals = tmp_data[tmp_data["combos"] == "1,1"].drop(
                "combos", axis=1)
            example_data = example_data.set_index(metadata_headers)
            example_data[:] = np.subtract(example_data.values, background_vals.values)
            example_data = example_data.reset_index()

            # append
            grammar_data_norm.append(example_data)

        if len(grammar_data_norm) == 0:
            continue

        # concat
        grammar_data_norm = pd.concat(grammar_data_norm, axis=0)
        
        # save this out to plot in R
        plot_data_file = "{}.agg_data.txt.gz".format(grammar_prefix)
        grammar_data_norm.to_csv(
            plot_data_file, sep="\t", compression="gzip", header=True, index=False)

        # plot cmd
        plot_prefix = "{}.endogenous.mutational".format(grammar_prefix)
        plot_cmd = "./plot.mpra.hypothesis.mut.R {} {}".format(
            plot_data_file, plot_prefix)
        print plot_cmd
        os.system(plot_cmd)
        
    return



def main():
    """analyze mpra test data
    """
    # inputs
    PLASMID_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-06-05.plasmid_lib_qc"
    DATA_DIR = "/srv/scratch/dskim89/ggr/ggr.mpra.2019-09-10.miseq_qc/raw"
    OUT_DIR = "/srv/scratch/dskim89/ggr/ggr.mpra.2019-09-10.miseq_qc/results"
    os.system("mkdir -p {}".format(OUT_DIR))

    ref_fasta_file = "/srv/scratch/dskim89/ggr/ggr.mpra.2019-06-05.plasmid_lib_qc/sequences.fa"
    metadata_file = "/srv/scratch/dskim89/ggr/inference.2019-03-12/mpra/mpra.seqs.run-0.txt"
    
    # make sure to load bwa 0.7.10
    
    # ===========
    # GET DNA PLASMID COUNTS
    # ===========
    
    # analyze the plasmid library (to get normalization fractions)
    # need to look at r1 and r2 to get the UMIs to remove PCR duplicates
    plasmid_fastq1 = "{}/raw/DK-C2-plasmidlib_S1_L001_R1_001.fastq.gz".format(PLASMID_DIR)
    plasmid_fastq2 = "{}/raw/DK-C2-plasmidlib_S1_L001_R2_001.fastq.gz".format(PLASMID_DIR)
    plasmid_trim_fastq = "{}/plasmid.trim.fastq.gz".format(PLASMID_DIR)

    # trim and get umis
    if not os.path.isfile(plasmid_trim_fastq):
        trim_and_attach_umi(plasmid_fastq1, plasmid_fastq2, plasmid_trim_fastq)

    # get counts
    plasmid_counts = "{}/plasmid.dedup.counts.txt.gz".format(PLASMID_DIR)
    if not os.path.isfile(plasmid_counts):
        fastq_to_counts(plasmid_trim_fastq, ref_fasta_file, "{}/plasmid".format(PLASMID_DIR))
        
    # ===========
    # GET EXPRESSION COUNTS
    # ===========
    
    # args
    mixed_fastq_file = "{}/Undetermined_S0_L001_R1_001.fastq.gz".format(DATA_DIR)
    demux_fastqs = {
        "d0_b1": "{}/1-D0_S1_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d3_b1": "{}/1-D3_S2_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d6_b1": "{}/1-D6_S3_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d0_b2": "{}/2-D0_S4_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d3_b2": "{}/2-D3_S5_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d6_b2": "{}/2-D6_S6_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d0_b3": "{}/3-D0_S7_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d3_b3": "{}/3-D3_S8_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d6_b3": "{}/3-D6_S9_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d0_b4": "{}/4-D0_S10_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d3_b4": "{}/4-D3_S11_L001_R1_001.fastq.gz".format(DATA_DIR),
        "d6_b4": "{}/4-D6_S12_L001_R1_001.fastq.gz".format(DATA_DIR)}
    
    # chop the fastq file AND demultiplex AND attach UMIs to ID
    print "demultiplex, trim, attach UMIs to ID"
    trimmed_fastqs = {}
    for key in demux_fastqs.keys():
        trimmed_fastq = "{}/{}.trim.fastq.gz".format(OUT_DIR, key)
        trimmed_fastqs[key] = trimmed_fastq
        
    reads_seen = 0
    skipped = 0
    if not os.path.isfile(trimmed_fastq):

        # load in demultiplexed IDs
        print "loading demultiplexed IDs"
        read_ids = {}
        for key in sorted(demux_fastqs.keys()):
            fastq_file = demux_fastqs[key]
            print fastq_file
            line_num = 0
            with gzip.open(fastq_file, "r") as fp:
                for line in fp:
                    if line_num % 4000000 == 0:
                        print "reads:", line_num / 4
                    if line_num % 4 == 0:
                        read_id = line.strip().split()[0]
                        read_ids[read_id] = key
                    line_num += 1

        # read in non demultiplexed file and split out
        with gzip.open(mixed_fastq_file, "r") as fp:
            while True:
                if reads_seen % 1000000 == 0:
                    print reads_seen
                
                # get a read (4 lines)
                read = []
                for i in range(4):
                    read.append(fp.readline().strip())

                    if read[0] == "":
                        break
                    
                if read[0] == "":
                    break
                    
                # keep id
                read_id = read[0].split()[0]
                
                # attach UMI to read name
                new_id = read[0].split()[0] + ":" + read[0][-12:]
                read[0] = new_id
                
                # trim lines
                read[1] = read[1].strip()[0:20]
                read[3] = read[3].strip()[0:20]

                try:
                    fastq_key = read_ids[read_id]
                    with gzip.open(trimmed_fastqs[fastq_key], "a") as out:
                        out.write("{}\n".format("\n".join(read)))
                except:
                    skipped += 1
                    pass

                reads_seen += 1

        print "reads seen", reads_seen
        print "skipped", skipped

    # get counts
    counts_files = {"plasmid": plasmid_counts}
    for key in sorted(trimmed_fastqs.keys()):
        out_prefix = "{}/{}".format(OUT_DIR, key)
        fastq_to_counts(trimmed_fastqs[key], ref_fasta_file, out_prefix)
        counts_files[key] = "{}.dedup.counts.txt.gz".format(out_prefix)

    # and then merge all counts into a matrix file
    out_prefix = "{}/ggr".format(OUT_DIR)
    count_mat_file = "{}.counts.raw.mat.txt.gz".format(out_prefix)
    signal_mat_file = "{}.signal.mat.txt.gz".format(out_prefix)
    if not os.path.isfile(signal_mat_file):
        counts_to_mat(counts_files, out_prefix, metadata_file=metadata_file)
        
    # run QC
    qc_cmd = "./run_qc.R {} {}".format(count_mat_file, signal_mat_file)
    print qc_cmd
    #os.system(qc_cmd)

    # signal file w metadata
    signal_mat_file = "{}.signal.w_metadata.mat.txt.gz".format(out_prefix)
    
    # now try check hypotheses:
    # endogenous seq - trajectories
    traj_dir = "{}/hypotheses.trajectories".format(OUT_DIR)
    os.system("mkdir -p {}".format(traj_dir))
    #test_grammar_trajectories(signal_mat_file, traj_dir)
    
    # mutational analyses
    mut_dir = "{}/hypotheses.mutational".format(OUT_DIR)
    os.system("mkdir -p {}".format(mut_dir))
    test_mutational_hypotheses(signal_mat_file, mut_dir)
    
    return

main()
