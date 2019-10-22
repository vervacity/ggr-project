"""description: helpful fns for GO terms
"""


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
    

    


