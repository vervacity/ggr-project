
import os
import sys
import json
import requests


def _upload_file(
        filename,
        full_path_filename,
        bucket_url,
        params):
    """ we pass the file object (fp) directly to the request as the data to be uploaded
    the target URL is a combination of the buckets link with the desired filename separated by a slash
    """
    with open(full_path_filename, "rb") as fp:
        r = requests.put(
            "{}/{}".format(bucket_url, filename),
            data=fp,
            # No headers included in the request, since it's a raw byte request
            params=params)
        
    return



def main():
    """upload script
    """
    # USER INPUT
    ACCESS_TOKEN = sys.argv[1] # always read in access token, do not save to script!!
    deposition_id = None # either none or a numerical ID

    # file(s) - give with full paths
    UPLOAD_DIR = "/srv/scratch/dskim89/ggr/ggr.paper.2020-09-27.datasets_compressed/ggr/submit_models"
    filenames = [
        "{}/nn.ggr.models.basset.clf.tar.gz".format(UPLOAD_DIR),
        "{}/nn.ggr.models.basset.regr.tar.gz".format(UPLOAD_DIR)
    ]

    # deposition metadata - fill this out as needed
    metadata = {
        "metadata": {
            "title": "Convolutional Neural Net (CNN) models for epigenomic landscapes in epithelial differentiation - Basset architecture, classification and regression",
            "upload_type": "dataset",
            "description": "Deep learning models trained on epigenomic landscapes in keratinocyte differentiation",
            "creators": [
                {"name": "Kim, Daniel Sunwook",
                 "affiliation": "Stanford School of Medicine"},
                {"name": "Kundaje, Anshul",
                 "affiliation": "Stanford School of Medicine"}
            ]
        }
    }    

    # END USER INPUT
    
    
    # set up
    if deposition_id is None:
        new_deposition = True
    else:
        new_deposition = False
    
    # API set up
    headers = {"Content-Type": "application/json"}
    params = {'access_token': ACCESS_TOKEN}

    # get the deposition: either make new, or grab old one
    if new_deposition:
        r = requests.post(
            'https://zenodo.org/api/deposit/depositions',
            params=params,
            json={},
            # Headers are not necessary here since "requests" automatically
            # adds "Content-Type: application/json", because we're using
            # the "json=" keyword argument
            # headers=headers, 
            headers=headers)
        deposition_id = r.json()["id"]

        # add deposition metadata
        r = requests.put(
            "https://zenodo.org/api/deposit/depositions/{}".format(deposition_id),
            params={'access_token': ACCESS_TOKEN},
            data=json.dumps(metadata),
            headers=headers)
        
    else:
        r = requests.get(
            "https://zenodo.org/api/deposit/depositions/{}".format(deposition_id),
                         params=params)
    #print r.status_code
    print r.json()
    bucket_url = r.json()["links"]["bucket"]
    
    # upload files
    for filename in filenames:
        full_path_filename = filename
        base_filename = os.path.basename(filename)
        _upload_file(base_filename, full_path_filename, bucket_url, params)

    
    return

main()
