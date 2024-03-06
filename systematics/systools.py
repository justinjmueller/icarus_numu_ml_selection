import numpy as np
import pandas as pd
import awkward as ak
import uproot
from tqdm import tqdm

def calc_fractional_error(cov, cv):
    """
    Normalize each element of the covariance matrix by the product of
    the CV vector elements.

    Parameters
    ----------
    cov: numpy.array
        The covariance matrix to be normalized.
    cv: numpy.array
        The CV counts per bin to be used in the normalization.
    
    Returns
    -------
    numpy.array
        The fractional covariance matrix.
    """
    norm = np.outer(cv, cv)        
    return np.divide(cov, norm, where=norm!=0)

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
    df = pd.DataFrame(selected, columns=header[:len(selected[0])])
    for k in header[:len(df.columns)]:
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
        The full path to the input CAF file.
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

    cols = list(index_keys.values())
    common = caf[~caf.duplicated(subset=cols, keep='first')].merge(selected[selected['nu_id'] != -1], on=cols, how='right')
    begin = rf['rec.mc.nu.wgt.univ..idx'].array()[0][widx]
    end = begin + rf['rec.mc.nu.wgt.univ..length'].array()[0][widx]
    l = int(rf['rec.mc.nu.wgt.univ..totarraysize'].array()[0] / len(rf['rec.mc.nu.wgt..length'].array()[0]))
    weights = np.empty((len(common), end - begin))
    ioffset, eoffset = 0, 0
    for wgt in rf['rec.mc.nu.wgt.univ'].iterate(step_size='1 GB'):
        wgt = wgt['rec.mc.nu.wgt.univ']
        batch = common[['i', 'nu_id']].to_numpy()
        ext = np.array([wgt[ni-ioffset][l*n + begin : l*n + end] for ni, n in batch if ni-ioffset >= 0 and ni-ioffset < len(wgt)])
        weights[eoffset:eoffset+ext.shape[0], :] = ext
        ioffset += len(wgt)
        eoffset += ext.shape[0]

    return weights

def calc_multisim_covariance(sys, caf, header, var, bins):
    """
    Calculates the covariance matrix of the binned reconstructed variable for
    the selected candidates using the systematic weights.

    Parameters
    ----------
    sys: dict
        The dictionary containing the configuration details for the systematic
        variation. 
    caf: str
        The full path to the input CAF file.
    header: list[str]
        The list of column names for the selected interactions in the
        input log file.
    var: str
        The name of the reconstructed variable.
    bins: list[float]
        The number of bins, lower edge, and upper edge of the variable.

    Returns
    -------
    cov: numpy.ndarray
        The covariance matrix with shape (bins,bins) calculated with respect
        to the selected candidates and the weights.
    cv: numpy.ndarray
        The central value of the selected candidates per bin.
    """
    # Load the selected interactions
    selected = read_log(sys['cv_log'], f'SELECTED_{sys["channel"].upper()}', header)

    # Calculate the bin edges and the bin index for each selected
    # interaction.
    _, bin_edges = np.histogram(selected[var], bins=int(bins[0]), range=bins[1:])
    nbins = len(bin_edges) - 1
    selected['bidx'] = np.digitize(selected[var], bin_edges) - 1

    # Extract the weights for the systematic parameter.
    weights = extract_weights(caf, selected, sys['index'])
    is_cosmic = selected['nu_id'] == -1
    ensemble = list()
    cv = list()

    # Loop over the bins and store 1) the total count per bin for each
    # universe (ensemble has shape (nbins, nuniverses)) and 2) the
    # central value for each bin (cv has shape (nbins,)). Interactions
    # with cosmic origin are treated as having unit weight across all
    # universes.
    for b in range(int(bins[0])):
        mask = selected.bidx[~is_cosmic] == b
        e = np.sum(weights[mask, :], axis=0)
        cosmics = np.sum(selected.bidx[is_cosmic] == b)
        if cosmics > 0:
            e += np.repeat(cosmics, e.shape[0])
        ensemble.append(e)
        cv.append(np.sum(selected.bidx == b))

    # Calculate the covariance matrix of the ensemble about the central
    # value points.
    ensemble = np.stack(ensemble, axis=0)
    cv = np.stack(cv, axis=0)
    cov = np.cov(np.subtract(cv[:, np.newaxis], ensemble))
    return cov, cv

def load_detector_variation(header, sys):
    """
    Loads the signal interactions and selected candidates common to the
    CV sample and the systematic variation sample.

    Parameters
    ----------
    header: list[str]
        The list of column names for the selected interactions in the
        input log file.
    sys: dict
        The dictionary containing the configuration details for the systematic
        variation.
    
    Returns
    -------
    common_events: pandas.DataFrame
        The set of events common to both samples.
    cv_selected: pandas.DataFrame
        The selected candidates for the CV sample.
    sys_selected: pandas.DataFrame
        The selected candidates for the systematic variation sample.
    """
    cv_events = read_log(sys['cv_log'], 'EVENT', header)
    sys_events = read_log(sys['sys_log'], 'EVENT', header)
    cv_selected = read_log(sys['cv_log'], f'SELECTED_{sys["channel"].upper()}', header)
    sys_selected = read_log(sys['sys_log'], f'SELECTED_{sys["channel"].upper()}', header)
    common_events = cv_events.merge(sys_events, how='inner', on=['run', 'subrun', 'event'])
    return common_events, cv_selected, sys_selected

def bootstrap_iterate(common, selected, nbins, stats_limit=1):
    """
    Performs one bootstrap iteration of the set of signal interactions
    common to all samples and produces a count of selected events for
    each population.

    Parameters
    ----------
    common: pandas.DataFrame
        The set of signal interactions common to all samples.
    selected: list[pandas.DataFrame]
        The selected signal candidates for each sample.
    nbin: int
        The number of bins used in the bootstrap.
    stats_limit: float
        The fraction of the common signal interactions to be used in the bootstrap.

    Returns
    -------
    bins: numpy.array
        The number of selected interactions in each bin.
     """
    choice = np.random.choice(np.arange(int(stats_limit * len(common)), dtype=int), size=int(stats_limit * len(common)), replace=True)
    merged = [s.merge(common.iloc[choice], how='inner', on=['run', 'subrun', 'event']) for s in selected]
    bins = np.zeros(shape=(len(merged), nbins), dtype=np.int16)
    for si, s in enumerate(merged):
        bins[si, :] = np.array([np.sum(s['bidx'] == b) for b in range(nbins)])
    return bins

def random_cholesky_variable(cov):
    """
    Generates a random vector using the Cholesky decomposition of
    the give covariance matrix.

    Parameters
    ----------
    cov: numpy.array
        The covariance matrix to be used.

    Returns
    -------
    numpy.array
        A random vector drawn from the given covariance matrix.
    """
    L = np.linalg.cholesky(cov)
    return L @ np.random.normal(size=cov.shape[0])

def generate_detector_universe(cov, nominal):
    """
    Generates a random detector universe using the Cholesky decomposition
    of the covariance matrix and the nominal difference vector.

    Parameters
    ----------
    cov: numpy.array
        The covariance matrix to be used.
    nominal: numpy.array
        The nominal difference vector to be used.
    
    Returns
    -------
    numpy.array
        A random detector universe drawn from the given covariance matrix.
    """
    perturbed = nominal + random_cholesky_variable(cov)
    return np.random.normal() * perturbed

def calc_detector_covariance(sys, header, var, bins):
    """
    Calculates the covariance matrix for the detector variation using
    the CV and systematic variation samples.

    Parameters
    ----------
    sys: dict
        The dictionary containing the configuration details for the systematic
        variation.
    header: list[str]
        The list of column names for the selected interactions in the
        input log file.
    var: str
        The name of the reconstructed variable.
    bins: list[float]
        The number of bins, lower edge, and upper edge of the variable.

    Returns
    -------
    cv: numpy.array
        The counts of selected interactions per bin.
    vnominal: numpy.array
        The nominal difference vector.
    rmatrix: numpy.array
        The covariance matrix of the nominal difference vector.
    dmatrix: numpy.array
        The covariance matrix of the detector response.
    vratio: numpy.array
        The bootstrapped ratio of the systematic variation to the CV sample.
    cratio: numpy.array
        The covariance matrix of the bootstrapped ratio.
    """
    common, cv_selected, sys_selected = load_detector_variation(header, sys)

    # Calculate the bin edges and the bin index for each interaction.
    _, bin_edges = np.histogram(cv_selected[var], bins=int(bins[0]), range=bins[1:])
    nbins = len(bin_edges) - 1
    cv_selected['bidx'] = np.digitize(cv_selected[var], bin_edges) - 1
    sys_selected['bidx'] = np.digitize(sys_selected[var], bin_edges) - 1

    # Bootstrap the signal events to characterize the covariance of the bins.
    nboots = sys['nboots']
    bins = np.zeros((2, nbins, nboots), dtype=np.int16)
    for i in range(nboots):
        bins[:, :, i] = bootstrap_iterate(common, [cv_selected, sys_selected], nbins)

    # Calculate V_nominal and the associated covariance matrix (M_R).
    vnominal = np.mean(bins[1,:,:] - bins[0,:,:], axis=1)
    rmatrix = np.cov(bins[1,:,:] - bins[0,:,:] - vnominal[:,np.newaxis])

    # Calculate the ratio and the associated covariance matrix.
    vratio = np.ones((nbins, nboots), dtype=np.float64)
    np.divide(bins[1,:,:], bins[0,:,:], where=bins[0,:,:]!=0, out=vratio)
    vratio = np.mean(vratio, axis=1)
    cratio = np.cov(np.divide(bins[1,:,:], bins[0,:,:], where=bins[0,:,:]!=0) - vratio[:,np.newaxis])

    # Mask bins which will cause a singular response matrix.
    singular_mask = (vnominal != 0)

    # Calculate the detector response matrix (M_D).
    nuniverses = sys['nuniverses']
    universes = np.zeros((np.sum(singular_mask), nuniverses), dtype=np.float64)
    for n in range(nuniverses):
        universes[:, n] = generate_detector_universe(rmatrix[singular_mask, :][:, singular_mask], vnominal[singular_mask])
    dmatrix = np.zeros((nbins, nbins), dtype=np.float64)
    cov = np.cov(universes)
    for bi, b in enumerate(np.arange(nbins)[singular_mask]):
        dmatrix[b, singular_mask] = cov[bi, :]
    cv = np.array([np.sum(cv_selected['bidx'] == b) for b in range(nbins)])

    return cv, vnominal, rmatrix, dmatrix, vratio, cratio

def calc_statistical_covariance(sys, header, var, bins):
    """
    Calculate the statistical covariance matrix for the selected
    interactions.

    Parameters
    ----------
    sys: dict
        The dictionary containing the configuration details for the
        uncertainty source.
    header: list[str]
        The list of column names for the selected interactions in the
        input log file.
    var: str
        The name of the reconstructed variable.
    bins: list[float]
        The number of bins, lower edge, and upper edge of the variable.

    Returns
    -------
    cov: numpy.array
        The statistical covariance matrix.
    """
    selected = read_log(sys['cv_log'], f'SELECTED_{sys["channel"].upper()}', header)

    # Calculate the bin edges and the bin index for each selected
    # interaction.
    _, bin_edges = np.histogram(selected[var], bins=int(bins[0]), range=bins[1:])
    nbins = len(bin_edges) - 1
    selected['bidx'] = np.digitize(selected[var], bin_edges) - 1

    # Calculate the statistical covariance matrix by constructing a
    # diagonal matrix with the counts of selected interactions per bin.
    cov = np.diag([np.sum(selected['bidx'] == b) for b in range(nbins)])
    
    return cov