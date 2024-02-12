import numpy as np
import pandas as pd
import awkward as ak
import uproot
from tqdm import tqdm

def read_log(path, tag, header):
    """
    Reads an input log file and extracts lines with the specified tag
    into a Pandas DataFrame.

    Parameters
    ----------
    path: str
        The full path of the input log file.
    tag: str
        The identifier that tags relevant lines in the log file.
    header: list[str]
        The list of column names for the CSV file.

    Returns
    -------
    df: Pandas.DataFrame
        The DataFrame containing the requested information.
    """
    input_file = open(path)
    lines = input_file.readlines()
    selected = [x.strip('\n').split(',')[1:] for x in lines if tag in x]
    selected = [x if x[-1] != '' else x[:-1] for x in selected]
    df = pd.DataFrame(selected, columns=header)
    for k in header:
        df[k] = pd.to_numeric(df[k], errors='coerce', downcast='float')
        if df[k].apply(float.is_integer).all():
            df[k] = df[k].astype(int)
    return df

def extract_weights(path, selected, widx):
    """
    Extracts weights from the specified CAF file for the requested
    systematic source/parameter. 

    Parameters
    ----------
    path: str
        The full path of the input CAF file.
    selected: pandas.DataFrame
        The DataFrame containing information about the selected
        neutrinos (run, subrun, event, nu_id).
    widx: int
        The index of the systematic source/parameter.

    Returns
    -------
    weights: numpy.array
        The weights for each selected candidate for the specified
        systematic parameter with shape ()
    """
    rf = uproot.open(path)['recTree']
    index_keys = {'rec.hdr.run': 'run', 'rec.hdr.subrun': 'subrun',
                  'rec.hdr.evt': 'event', 'rec.mc.nu.index': 'nu_id'}
    caf = rf.arrays(list(index_keys.keys()), library='ak')
    caf = ak.to_dataframe(caf).rename(columns=index_keys)
    caf['i'] = [x[0] for x in caf.index]

    common = caf.merge(selected, on=list(index_keys.values()), how='inner')
    begin = rf['rec.mc.nu.wgt.univ..idx'].array()[0][widx]
    end = begin + rf['rec.mc.nu.wgt.univ..length'].array()[0][widx]
    l = int(rf['rec.mc.nu.wgt.univ..totarraysize'].array()[0] / len(rf['rec.mc.nu.wgt..length'].array()[0]))
    universes = rf['rec.mc.nu.wgt.univ'].array()
    weights = np.array([universes[ni][l*n + begin : l*n + end] for ni, n in common[['i', 'nu_id']].to_numpy()])

    return weights