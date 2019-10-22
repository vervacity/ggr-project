"""description: h5 util functions
"""

import h5py

import numpy as np


def load_data_from_multiple_h5_files(h5_files, key, example_indices=None):
    """convenience wrapper
    example indices is a list of arrays
    """
    key_data = []
    for h5_idx in range(h5_files):
        h5_file = h5_files[h5_idx]
        h5_indices = example_indices[h5_idx]
        
        with h5py.File(h5_file, "r") as hf:
            if example_indices is not None:
                key_data.append(hf[key][h5_indices])
            else:
                key_data.append(hf[key][:])
    key_data = np.concatenate(key_data, axis=0)
    
    return key_data
