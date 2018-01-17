# description: code for working with motifs

import os
import gzip

import numpy as np
import pandas as pd

from scipy.stats import pearsonr

from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

from multiprocessing import Pool


def read_pwm_file(pwm_file, as_dict=False):
    """Extracts motifs into PWM class format
    """
    # option to set up as dict or list
    if as_dict:
        pwms = {}
    else:
        pwms = []

    # open motif file and read
    with open(pwm_file) as fp:
        line = fp.readline().strip()
        while True:
            if line == '':
                break
            header = line.strip('>').strip()
            weights = []
            
            while True:
                line = fp.readline()
                if line == '' or line[0] == '>': break
                weights.append(map(float, line.split()))

            pwm = PWM(np.array(weights).transpose(1,0), header)

            # store into dict or list
            if as_dict:
                pwms[header] = pwm
            else:
                pwms.append(pwm)
                
    return pwms



class PWM(object):
    """PWM class for PWM operations
    """
    
    def __init__(self, weights, name=None, threshold=None):
        self.weights = weights
        self.name = name
        self.threshold = threshold

        
    def normalize(self, style="gaussian", in_place=True):
        """Normalize pwm
        """
        if style == "gaussian":
            mean = np.mean(self.weights)
            std = np.std(self.weights)
            normalized_weights = (self.weights - mean) / std
        elif style == "probabilities":
            col_sums = self.weights.sum(axis=0)
            normalized_pwm_tmp = self.weights / np.amax(col_sums[np.newaxis,:])
            normalized_weights = np.nan_to_num(normalized_pwm_tmp)
        elif style == "log_odds":
            print "Not yet implemented"
        else:
            print "Style not recognized"
            
        if in_place:
            self.weights = normalized_weights
            return self
        else:
            new_pwm = PWM(normalized_weights, "{}.norm".format(self.name))
            return new_pwm
        
        
    def xcor(self, pwm, normalize=True):
        """Compute xcor score with other motif, return score and offset relative to first pwm
        """
        if normalize:
            pwm1_norm = self.normalize(in_place=False)
            pwm2_norm = pwm.normalize(in_place=False)
        else:
            pwm1_norm = pwm1
            pwm2_norm = pwm2

        # calculate xcor
        xcor_vals = correlate2d(pwm1_norm.weights, pwm2_norm.weights, mode='same')
        xcor_norm = xcor_vals / (pwm1_norm.weights.shape[0]*pwm1_norm.weights.shape[1])
        score = np.max(xcor_norm[1,:])
        offset = np.argmax(xcor_norm[1,:]) - int(math.ceil(pwm2_norm.weights.shape[1] / 2.) - 1)

        return score, offset


    def pearson_xcor(self, pwm, ncor=False):
        """Calculate pearson across offsets, return best score
        and best position
        """
        # get total offsets
        offset_total = self.weights.shape[1] + pwm.weights.shape[1] - 1
                
        # set up values
        max_cor_val = -1
        best_offset = 0

        for i in xrange(offset_total):
            self_padded_weights, other_padded_weights = self.pad_by_offset(pwm, i)
            start_idx, stop_idx = PWM(np.maximum(self_padded_weights,other_padded_weights)).chomp_points()

            self_padded_weights_chomped = self_padded_weights[:,start_idx:stop_idx]
            other_padded_weights_chomped = other_padded_weights[:,start_idx:stop_idx]
            
            # take both and calculate
            cor_val, pval = pearsonr(
                self_padded_weights_chomped.flatten(),
                other_padded_weights_chomped.flatten())

            # normalization (RSAT)
            if ncor:
                width_norm_val = (
                    self.weights.shape[1] + pwm.weights.shape[1] - self_padded_weights_chomped.shape[1]) / float(
                        self_padded_weights_chomped.shape[1])
                cor_val = cor_val * width_norm_val
                
            if cor_val > max_cor_val:
                max_cor_val = cor_val
                best_offset = i

        return max_cor_val, best_offset

    
    def rsat_cor(self, pwm, ncor=False, offset=None):
        """Calculate a pearson correlation across all positions
        """
        # tODO - dont really need this
        # return the pearson
        val, offset = self.pearson_xcor(pwm, ncor=ncor)
        
        return val

    
    def get_probs(self, count_factor=500, epsilon=0.01):
        """Take weights and convert to a PFM
        """
        #pseudo_counts = count_factor * np.exp(self.weights) + epsilon

        probs = np.exp(self.weights) / np.sum(np.exp(self.weights), axis=0)

        return probs

    
    def get_ic(self):
        """Get information content per each position per base pair
        """
        probs = self.get_probs()
        ic = 2 + np.sum(probs * np.log2(probs), axis=0)
        
        return ic
    
    
    def chomp_points(self, ic_thresh=0.4):
        """Remove leading/trailing Ns. In place, but also outputs self
        """
        ic = self.get_ic()
        
        # find starting point
        # iterate through positions unti you find the last
        # position before a high IC position
        start_idx = 0
        while start_idx < self.weights.shape[1]:
            # calculate IC of position
            if ic[start_idx] > ic_thresh:
                break
            start_idx += 1
        if start_idx == self.weights.shape[1]:
            start_idx = 0

        # find stop point
        stop_idx = self.weights.shape[1] - 1
        while stop_idx > 0:
            # calculate IC of position
            if ic[stop_idx] > ic_thresh:
                break
            stop_idx -= 1
        if stop_idx == 0:
            stop_idx = self.weights.shape[1]

        return start_idx, stop_idx + 1

    
    def chomp(self, ic_thresh=0.4):
        """Remove leading/trailing Ns. In place, but also outputs self
        """
        start_idx, stop_idx = self.chomp_points(ic_thresh=ic_thresh)

        # chomp
        self.weights = self.weights[:,start_idx:stop_idx+1]
        
        return self


    def pad_weights(self, start_pad, end_pad, in_place=False):
        """Pad weights with start_pad bp in the front
        and end_pad bp in the back
        """
        padded_weights = np.concatenate(
            (np.zeros((4, start_pad)),
             self.weights,
             np.zeros((4, end_pad))),
            axis=1)
        
        return padded_weights


    def pad_by_offset(self, pwm, offset, chomp=False):
        """Pads self and other pwm to be same length
        """
        total_length = self.weights.shape[1] + 2*(pwm.weights.shape[1] - 1) #-offset
        
        # self pwm
        front_pad = pwm.weights.shape[1] - 1
        end_pad = total_length - (front_pad + self.weights.shape[1])
        self_padded_weights = self.pad_weights(front_pad, end_pad)

        # other pwm
        front_pad = offset
        end_pad = total_length - (front_pad + pwm.weights.shape[1])
        other_padded_weights = pwm.pad_weights(front_pad, end_pad)

        return self_padded_weights, other_padded_weights

    
    def merge(
            self,
            pwm,
            offset,
            weights=(1.0, 1.0),
            ic_thresh=0.4,
            background_freq=0.25,
            new_name=None,
            normalize=False):
        """Merge in another PWM and output a new PWM
        """
        self_padded_weights, other_padded_weights = self.pad_by_offset(pwm, offset)

        self_padded_probs = np.exp(self_padded_weights) / np.sum(np.exp(self_padded_weights), axis=0)
        other_padded_probs = np.exp(other_padded_weights) / np.sum(np.exp(other_padded_weights), axis=0)
        
        # merge
        # merging by first moving back to prob space and then
        # returning to log space
        weighted_summed_probs = weights[0] * self_padded_probs + weights[1] * other_padded_probs
        weight_sum = weights[0] + weights[1]
        new_pwm = PWM(
            np.log2(
                weighted_summed_probs / (weight_sum * background_freq)),
            name=new_name)

        # chomp
        new_pwm.chomp(ic_thresh=ic_thresh)

        # normalize if desired
        if normalize:
            new_pwm.normalize()

        import ipdb
        ipdb.set_trace()
            
        return new_pwm

    
    def to_motif_file(
            self,
            motif_file,
            motif_format="homer",
            pseudo_counts=500):
        """Write PWM out to file
        """
        # TODO allow various formats of output
        # such as transfac, homer, etc
        with open(motif_file, 'a') as fp:
            if motif_format == "homer":
                fp.write('>{}\n'.format(self.name))
                for i in range(self.weights.shape[1]):
                    vals = self.weights[:,i].tolist()
                    val_strings = [str(val) for val in vals]
                    fp.write('{}\n'.format('\t'.join(val_strings)))
            elif motif_format == "transfac":
                # TODO does not include consensus letter at the moment
                fp.write('ID {}\n'.format(self.name))
                fp.write('BF Homo_sapiens\n')
                fp.write('P0\tA\tC\tG\tT\n')
                for i in range(self.weights.shape[1]):
                    exp_vals = np.exp(self.weights[:,i])
                    vals = pseudo_counts * (exp_vals / np.sum(exp_vals))
                    val_strings = [str(val) for val in vals.tolist()]
                    fp.write("{num:02d}\t{}\n".format("\t".join(val_strings), num=i+1))
                fp.write("XX\n")
                fp.write("//\n")
        
        return None

    
    def plot(self, out_file, tmp_dir="."):
        """Plot out PWM to visualize
        """
        # save out in transfac format
        tmp_out_file = "{}/motif.{}.vals.transfac.tmp".format(
            tmp_dir, self.name.strip().split("_")[0])
        self.to_motif_file(tmp_out_file, motif_format="transfac")

        # and call weblogo
        weblogo_cmd = (
            "weblogo "
            "-X NO --errorbars NO --fineprint \"\" "
            "-C \"#CB2026\" A A "
            "-C \"#34459C\" C C "
            "-C \"#FBB116\" G G "
            "-C \"#0C8040\" T T "
            "-f {0} "
            "-D transfac "
            "-F pdf "
            "-o {1}").format(
                tmp_out_file, out_file)
        print weblogo_cmd
        os.system(weblogo_cmd)

        # and remove tmp file
        os.system("rm {}".format(tmp_out_file))
        
        return None


def correlate_pwm_pair(input_list):
    """get cor and ncor for pwm1 and pwm2
    Set up this way because multiprocessing pool only takes 1
    input
    """
    i = input_list[0]
    j = input_list[1]
    pwm1 = input_list[2]
    pwm2 = input_list[3]
    
    motif_cor = pwm1.rsat_cor(pwm2)
    motif_ncor = pwm1.rsat_cor(pwm2, ncor=True)

    return i, j, motif_cor, motif_ncor


def correlate_pwms(
        pwms,
        cor_thresh=0.6,
        ncor_thresh=0.4,
        num_threads=24):
    """Correlate PWMS
    """
    # set up
    num_pwms = len(pwms)
    cor_mat = np.zeros((num_pwms, num_pwms))
    ncor_mat = np.zeros((num_pwms, num_pwms))

    pool = Pool(processes=num_threads)
    pool_inputs = []
    # for each pair of motifs, get correlation information
    for i in xrange(num_pwms):
        for j in xrange(num_pwms):

            # only calculate upper triangle
            if i > j:
                continue

            pwm_i = pwms[i]
            pwm_j = pwms[j]
            
            pool_inputs.append((i, j, pwm_i, pwm_j))

    # run multiprocessing
    pool_outputs = pool.map(correlate_pwm_pair, pool_inputs)

    for i, j, motif_cor, motif_ncor in pool_outputs:
        # if passes cutoffs, save out to matrix
        if (motif_cor >= cor_thresh) and (motif_ncor >= ncor_thresh):
            cor_mat[i,j] = motif_cor
            ncor_mat[i,j] = motif_ncor        

    # and reflect over the triangle
    lower_triangle_indices = np.tril_indices(cor_mat.shape[0], -1)
    cor_mat[lower_triangle_indices] = cor_mat.T[lower_triangle_indices]
    ncor_mat[lower_triangle_indices] = ncor_mat.T[lower_triangle_indices]

    # multiply each by the other to double threshold
    cor_present = (cor_mat > 0).astype(float)
    ncor_present = (ncor_mat > 0).astype(float)

    # and mask
    cor_filt_mat = cor_mat * ncor_present
    ncor_filt_mat = ncor_mat * cor_present

    return cor_filt_mat, ncor_filt_mat


def correlate_pwms_old(
        pwms,
        cor_thresh=0.6,
        ncor_thresh=0.4):
    """Correlate PWMS
    """
    # set up
    pwms_ids = [pwm.name for pwm in pwms]
    num_pwms = len(pwms)
    cor_mat = np.zeros((num_pwms, num_pwms))
    ncor_mat = np.zeros((num_pwms, num_pwms))

    # for each pair of motifs, get correlation information
    for i in xrange(num_pwms):
        if i % 10 == 0:
            print "Finished {} rows...".format(i)
        for j in xrange(num_pwms):

            # only calculate upper triangle
            if i > j:
                continue

            pwm_i = pwms[i]
            pwm_j = pwms[j]

            # TODO(dk): remove gaps when comparing
            # for now don't ungap
            motif_cor = pwm_i.rsat_cor(pwm_j)
            motif_ncor = pwm_i.rsat_cor(pwm_j, ncor=True)

            # if passes cutoffs, save out to matrix
            if (motif_cor >= cor_thresh) and (motif_ncor >= ncor_thresh):
                cor_mat[i,j] = motif_cor
                ncor_mat[i,j] = motif_ncor

    # and reflect over the triangle
    lower_triangle_indices = np.tril_indices(cor_mat.shape[0], -1)
    cor_mat[lower_triangle_indices] = cor_mat.T[lower_triangle_indices]
    ncor_mat[lower_triangle_indices] = ncor_mat.T[lower_triangle_indices]

    # multiply each by the other to double threshold
    cor_present = (cor_mat > 0).astype(float)
    ncor_present = (ncor_mat > 0).astype(float)

    # and mask
    cor_filt_mat = cor_mat * ncor_present
    ncor_filt_mat = ncor_mat * cor_present

    # pandas and save out
    cor_df = pd.DataFrame(cor_filt_mat, index=pwms_ids, columns=pwms_ids)
    cor_df.to_csv(cor_mat_file, sep="\t")
    ncor_df = pd.DataFrame(ncor_filt_mat, index=pwms_ids, columns=pwms_ids)
    cor_df.to_csv(ncor_mat_file, sep="\t")


    return cor_filt_mat, ncor_filt_mat


def hagglom_pwms(
        cor_mat_file,
        pwm_dict,
        ic_thresh=0.4,
        cor_thresh=0.8,
        ncor_thresh=0.65):
    """hAgglom on the PWMs to reduce redundancy
    """
    # read in table
    cor_df = pd.read_table(cor_mat_file, index_col=0)

    # set up pwm lists
    # set up (PWM, weight)
    hclust_pwms = [(pwm_dict[key], 1.0) for key in cor_df.columns.tolist()]
    non_redundant_pwms = []

    # hierarchically cluster
    hclust = linkage(squareform(1 - cor_df.as_matrix()), method="ward")

    # keep a list of pwms in hclust, when things get merged add to end
    # (to match the scipy hclust structure)
    # put a none if not merging
    # if the motif did not successfully merge with its partner, pull out
    # it and its partner. if there was a successful merge, keep in there
    for i in xrange(hclust.shape[0]):
        idx1, idx2, dist, cluster_size = hclust[i,:]

        # check if indices are None
        pwm1, pwm1_weight = hclust_pwms[int(idx1)]
        pwm2, pwm2_weight = hclust_pwms[int(idx2)]

        if (pwm1 is None) and (pwm2 is None):
            hclust_pwms.append((None, None))
            continue
        elif (pwm1 is None):
            # save out PWM 2
            print "saving out {}".format(pwm2.name)
            non_redundant_pwms.append(pwm2)
            hclust_pwms.append((None, None))
            continue
        elif (pwm2 is None):
            # save out PWM1
            print "saving out {}".format(pwm1.name)
            non_redundant_pwms.append(pwm1)
            hclust_pwms.append((None, None))
            continue

        # try check
        try:
            cor_val, offset = pwm1.pearson_xcor(pwm2, ncor=True)
            ncor_val, offset = pwm1.pearson_xcor(pwm2, ncor=False)
        except:
            import ipdb
            ipdb.set_trace()

        if (cor_val > cor_thresh) and (ncor_val >= ncor_thresh):
            # TODO check ncor vals
            
            # store new merged pwm
            name = "{};{}".format(pwm1.name, pwm2.name)
            print name, cor_val, ncor_val
            new_pwm = pwm1.merge(
                pwm2,
                offset,
                weights=(pwm1_weight, pwm2_weight),
                new_name=name)
            hclust_pwms.append((new_pwm, pwm1_weight+pwm2_weight))
        else:
            print "saving out {}".format(pwm1.name)
            print "saving out {}".format(pwm2.name)
            non_redundant_pwms.append(pwm1)
            non_redundant_pwms.append(pwm2)
            hclust_pwms.append((None, None))

    return non_redundant_pwms
    

def reduce_pwm_redundancy(
        pwm_file,
        out_pwm_file,
        tmp_prefix="motif",
        ic_thresh=0.4,
        cor_thresh=0.8,
        ncor_thresh=0.65,
        num_threads=28):
    """Take in a PWM file, reduce redundancy, and
    output a reduced PWM file

    Note that RSAT stringent thresholds were ncor 0.65, cor 0.8
    """
    # read in pwm file
    pwms = read_pwm_file(pwm_file, as_dict=False)
    pwm_dict = read_pwm_file(pwm_file, as_dict=True)
    num_pwms = len(pwms)

    # trim pwms
    pwms = [pwm.chomp(ic_thresh=ic_thresh) for pwm in pwms]
    for key in pwm_dict.keys():
        pwm_dict[key] = pwm_dict[key].chomp(ic_thresh=ic_thresh)
    pwms_ids = [pwm.name for pwm in pwms]
    
    # correlate pwms - uses multiprocessing
    cor_mat_file = "{}.cor.motifs.mat.txt".format(tmp_prefix)
    ncor_mat_file = "{}.ncor.motifs.mat.txt".format(tmp_prefix)

    if True:
        cor_filt_mat, ncor_filt_mat = correlate_pwms(
            pwms,
            cor_thresh=cor_thresh,
            ncor_thresh=ncor_thresh,
            num_threads=num_threads)
        
        # pandas and save out
        cor_df = pd.DataFrame(cor_filt_mat, index=pwms_ids, columns=pwms_ids)
        cor_df.to_csv(cor_mat_file, sep="\t")
        ncor_df = pd.DataFrame(ncor_filt_mat, index=pwms_ids, columns=pwms_ids)
        cor_df.to_csv(ncor_mat_file, sep="\t")

    # TODO(dk) plot here (put in R script when stable)
    
    # read in matrix to save time
    non_redundant_pwms = hagglom_pwms(
        ncor_mat_file,
        pwm_dict,
        ic_thresh=ic_thresh,
        cor_thresh=cor_thresh,
        ncor_thresh=ncor_thresh)
    
    # save out reduced list to file
    # TODO - save cluster groups to metadata file and make new IDs for file
    
    for pwm in non_redundant_pwms:
        pwm.to_motif_file(out_pwm_file)
        # and plot
        # keep name simple for now
        plot_file = "{}.pwm.plot.pdf".format(pwm.name.strip().split("_")[0])
        pwm.plot(plot_file)

    return


def visualize_pwms(
        pwm_file,
        tmp_file="pwm_array.tmp"):
    """Visualize pwms: do this in R
    """
    # for each PWM, make a tmp array file
    # pass to R to visualize PWM
    line_num = 0
    with open(pwm_file, "r") as fp:
        for line in fp:
            if line.starts_with(">"):
                if line_num != 0:
                    # visualize previous array
                    pass
                # start a new file
                pass
            line_num += 1

    # and remove the tmp file
    os.system("rm {}".format(tmp_file))
    
    return


def add_hocomoco_metadata(
        pwm_file,
        out_pwm_file,
        metadata_file,
        conversion_file):
    """Given a motif file, adjust the names of the 
    motifs by info in the metadata
    """

    # read in the metadata file (want model to entrez)
    metadata = {}
    with open(metadata_file, "r") as fp:
        for line in fp:
            if line.startswith("Model"):
                continue
            fields = line.strip().split("\t")
            model = fields[0]
            entrez_id = fields[16]
            metadata[model] = entrez_id
    
    # read in conversion file (want entrez to ensembl and HGNC)
    entrez_to_ids = {}
    hgnc_to_ids = {}
    with gzip.open(conversion_file, "r") as fp:
        for line in fp:
            if line.startswith("ensembl"):
                continue
            fields = line.strip().split("\t")
            if len(fields[2]) == 0:
                continue
            if len(fields[1]) == 0:
                continue
            entrez_to_ids[fields[2]] = (fields[0], fields[1])
            hgnc_to_ids[fields[1]] = (fields[0], fields[1])
            
    # and then go through PWM file and adjust the names
    with open(pwm_file, "r") as fp:
        with open(out_pwm_file, "w") as out:
            for line in fp:
                if line.startswith(">"):
                    model = line.strip()[1:]
                    try:
                        ensembl_id, hgnc_id = entrez_to_ids[metadata[model]]
                    except:
                        try:
                            # try by hgnc
                            hgnc_id = model.split("_")[0]
                            ensembl_id, hgnc_id = hgnc_to_ids[hgnc_id]
                        except:
                            hgnc_id = model.split("_")[0]
                            ensembl_id = "UNK"
                    new_name = "{}_H11MO.{}.{}".format(hgnc_id, ensembl_id, line.strip().split(".")[-1])
                    out.write(">{}\n".format(new_name))
                else:
                    out.write(line)
                    
    return None


# HOCOMOCO v11
pwm_file = "/mnt/lab_data/kundaje/users/dskim89/annotations/hocomoco/v11/HOCOMOCOv11_core_pwms_HUMAN_mono.txt"
metadata_file = "/mnt/lab_data/kundaje/users/dskim89/annotations/hocomoco/v11/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv"
conversion_file = "/srv/scratch/shared/indra/dskim89/ggr/integrative/v0.2.4/annotations/hg19.ensembl_geneids.pc.gencode19.mappings.mat.gz"
adjusted_pwm_file = "HOCOMOCOv11_core_pwms_HUMAN_mono.renamed.txt"
reduced_pwm_file = "HOCOMOCOv11_core_pwms_HUMAN_mono.renamed.reduced.txt"

# test
add_hocomoco_metadata(pwm_file, adjusted_pwm_file, metadata_file, conversion_file)
reduce_pwm_redundancy(
    adjusted_pwm_file,
    reduced_pwm_file,
    num_threads=28)

# TODO write a pwm file vis function

# TODO write a renaming function? redundancy file - new ids?


