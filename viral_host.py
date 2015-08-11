#!/usr/bin/python

def extract_putative_host(db_names, nuc_matches):
    """ voting scheme nt > bacteria > viral
    Args: list of db names
          list of hits corresponding to db names
    Return:
        putatitive host annotation
    """
    length_db = len(db_names)
    host = "NA"
    assert length_db == len(nuc_matches)
    db2match = {}
    for index in range(0, length_db):
        db2match[db_names[index]] = nuc_matches[index]
    if ("viral" in db2match) and (db2match["viral"] != "NA"):
        host = db2match["viral"]
    elif ("bacteria" in db2match) and (db2match["bacteria"] != "NA"):
        host = db2match["bacteria"]
    elif ("nt" in db2match) and (db2match["nt"] != "NA"):
        host = db2match["nt"]
    return host
            

