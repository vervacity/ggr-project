
import os
import glob
import json
import subprocess

import pandas as pd


def get_md5sum(filename):
    """get md5sum
    """
    md5sum_val = subprocess.check_output(["md5sum", filename])
    md5sum_val = md5sum_val.split()[0]
    
    return md5sum_val


def extract_names_and_md5sums(model_dir, model_json, data_dict):
    """extract data
    """
    # get model json
    with open(model_json, "r" ) as fp:
        model = json.load(fp)

    # get prefix 
    model_ckpt_prefix = "{}/{}".format(
        model_dir,
        "/".join(model["checkpoint"].split("/")[-3:]))

    # model name
    model_name = model["checkpoint"].split("/")[-3]
    data_dict["names"].append(model_name)
    
    # model data
    model_data = "{}.data-00000-of-00001".format(model_ckpt_prefix)
    model_data_md5 = get_md5sum(model_data)
    data_dict["args_data_url"].append(model_data)
    data_dict["args_data_md5"].append(model_data_md5)

    # model index
    model_index = "{}.index".format(model_ckpt_prefix)
    model_index_md5 = get_md5sum(model_index)
    data_dict["args_index_url"].append(model_index)
    data_dict["args_index_md5"].append(model_index_md5)
    
    # model meta
    model_meta = "{}.meta".format(model_ckpt_prefix)
    model_meta_md5 = get_md5sum(model_meta)
    data_dict["args_meta_url"].append(model_meta)
    data_dict["args_meta_md5"].append(model_meta_md5)

    return 


def main():
    """upload models to kipoi
    """
    # inputs
    MODEL_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/models.2018-12-03"

    # data dict
    data_dict = {
        "names": [],
        "args_data_url": [],
        "args_data_md5": [],
        "args_index_url": [],
        "args_index_md5": [],
        "args_meta_url": [],
        "args_meta_md5": []
    }
    
    # encode models
    model_jsons = sorted(
        glob.glob("{}/encode-roadmap.basset.clf.testfold-*/model.json".format(MODEL_DIR)))
    for model_json in model_jsons:
        extract_names_and_md5sums(MODEL_DIR, model_json, data_dict)

    # GGR classification models
    model_jsons = sorted(
        glob.glob("{}/ggr.basset.clf.pretrained.folds.testfold-*/model.json".format(MODEL_DIR)))
    for model_json in model_jsons:
        print model_json
        extract_names_and_md5sums(MODEL_DIR, model_json, data_dict)
        
    # GGR regr models
    model_jsons = sorted(
        glob.glob("{}/ggr.basset.regr.pretrained.folds.testfold-*/model.json".format(MODEL_DIR)))
    for model_json in model_jsons:
        extract_names_and_md5sums(MODEL_DIR, model_json, data_dict)
        
    # save out to table
    data = pd.DataFrame(data_dict)
    print data
    headers_ordered = [
        "names",
        "args_meta_url",
        "args_meta_md5",
        "args_index_url",
        "args_index_md5",
        "args_data_url",
        "args_data_md5"
        ]
    data = data[headers_ordered]
    data.to_csv(
        "models_orig.tsv", sep="\t", header=True, index=False)
        
    return


main()
