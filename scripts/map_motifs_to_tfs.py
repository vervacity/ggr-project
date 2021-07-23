
import re
import sys

import pandas as pd


def main():
    """
    """
    
    motif_ordering_file = sys.argv[1]
    hocomoco_file = sys.argv[2]
    tf_patterns_file = sys.argv[3]
    out_file = sys.argv[4]

    motifs_ordered = pd.read_csv(motif_ordering_file, sep="\t")
    tfs_non_ordered = pd.read_csv(tf_patterns_file, sep="\t", index_col=0)

    # make mapping from hocomoco annotation file
    hocomoco = pd.read_csv(hocomoco_file, sep="\t")
    hocomoco["motif"] = [re.sub("HCLUST-\d+_", "", motif) for motif in hocomoco["hclust_model_name"].values]
    hocomoco["motif"] = [re.sub(".UNK.0.A", "", motif) for motif in hocomoco["motif"].values]
    hocomoco = hocomoco[["motif", "expressed_hgnc"]]
    hocomoco = hocomoco[~hocomoco["expressed_hgnc"].isna()]
    hocomoco["expressed_hgnc"] = hocomoco["expressed_hgnc"].str.split(';').tolist()
    if False:
        tmp = hocomoco.apply(
            lambda x: pd.Series(x['expressed_hgnc']),axis=1).stack().reset_index(level=1, drop=True)
        tmp.name = "expressed_hgnc"
        hocomoco = hocomoco.drop("expressed_hgnc", axis=1).join(tmp)
    hocomoco = hocomoco.set_index("motif")

    # go through ordered motifs and save out TFs in order of motifs
    ordered_tfs = []
    for motif in motifs_ordered["x"].values:
        expressed_tfs = hocomoco.loc[motif, "expressed_hgnc"]
        for tf in expressed_tfs:
            if tf in tfs_non_ordered.index:
                ordered_tfs.append(tf)

    ordered_tfs = pd.DataFrame({"ordered": ordered_tfs})
    ordered_tfs.to_csv(out_file, sep="\t")
    
    return

main()
