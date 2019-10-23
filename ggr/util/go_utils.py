"""description: helpful fns for GO terms
"""

import pandas as pd

class GoParams(object):
    """clean up and othe fns
    """

    REMOVE_SUBSTRINGS = [
        "anatomical",
        "ameboidal",
        "animal organ",
        "multicellular organism",
        "cellular developmental",
        "tube",
        "regulation of",
        "embryonic",
        "cardiovascular",
        "angiogenesis",
        "blood vessel",
        "vasculature",
        "immune",
        "defense",
        "signaling",
        "response to",
        "movement of"]

    REMOVE_EXACT_STRINGS = [
        "system process",
        "system development",
        "developmental process",
        "tissue development"]

    GOOD_GO_TERMS = [
        "stem cell differentiation",
        "hemidesmosome",
        "hair",
        "cell migration",
        "skin",
        "keratinocyte",
        "cell cycle",
        "epiderm",
        "cell junction",
        "cell proliferation",
        "adhesion",
        "lipase activity",
        "fatty acid",
        "sphingolipid",
        "glycerolipid"]


def is_enriched(go_file, filter_good_terms=False):
    """clean up terms and check if any terms in desired list
    """
    # pull in go terms
    go_terms = pd.read_csv(go_file, sep="\t")
    go_terms = go_terms[go_terms["domain"] == "BP"]
    go_terms = go_terms["term.name"].values.tolist()

    # filter out bad terms
    keep_terms = []
    for func_term in go_terms:
        keep = True
        for bad_term_str in GoParams.REMOVE_SUBSTRINGS:
            if bad_term_str in func_term:
                keep = False
        if func_term in GoParams.REMOVE_EXACT_STRINGS:
            keep = False
            
        if keep:
            keep_terms.append(func_term)

    # if filter good terms, do this here
            
    # if any remain, then enriched
    if len(keep_terms) > 0:
        enriched = True
    else:
        enriched = False
        
    return enriched
    


