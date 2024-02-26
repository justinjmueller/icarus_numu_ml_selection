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

def calc_covariance(selected, weights, bins):
    """
    Calculates the covariance matrix of the binned reconstructed variable for
    the selected candidates using the systematic weights.

    Parameters
    ----------
    selected: pandas.DataFrame
        The DataFrame containing the tabular data about the selected candidates
        and their binning.
    weights: numpy.ndarray
        The array of weights for each of the systematic universes.
    bins: int
        The number of bins.

    Returns
    -------
    cov: numpy.ndarray
        The covariance matrix with shape (bins,bins) calculated with respect
        to the selected candidates and the weights.
    """
    is_cosmic = selected['nu_id'] == -1
    ensemble = list()
    cv = list()
    for b in range(bins):
        mask = selected.bidx[~is_cosmic] == b
        e = np.sum(weights[mask, :], axis=0)
        cosmics = np.sum(selected.bidx[is_cosmic] == b)
        if cosmics > 0:
            e += np.repeat(cosmics, e.shape[0])
        ensemble.append(e)
        cv.append(np.sum(selected.bidx == b))
    ensemble = np.stack(ensemble, axis=0)
    cv = np.stack(cv, axis=0)
    cov = np.cov(np.subtract(cv[:, np.newaxis], ensemble))
    return cov

def load_detector_variation(sys):
    """
    Loads the signal interactions and selected candidates common to the
    CV sample and the systematic variation sample.

    Parameters
    ----------
    sys: dict
        The dictionary containing the configuration details for the systematic
        variation
    
    Returns
    -------
    common_signal: pandas.DataFrame
        The set of signal interactions common to both samples.
    cv_selected: pandas.DataFrame
        The selected candidates for the CV sample.
    sys_selected: pandas.DataFrame
        The selected candidates for the systematic variation sample.
    """
    header = ['run', 'subrun', 'event', 'nu_id', 'image_id', 'id', 'category', 'category_topology',
              'category_mode', 'visible_energy', 'leading_muon_ke', 'leading_proton_ke', 'leading_muon_pt',
              'leading_proton_pt', 'interaction_pt', 'leading_muon_cosine_theta_xz',
              'leading_proton_cosine_theta_xz', 'cosine_opening_angle', 'cosine_opening_angle_transverse']
    cv_signal = read_log(sys['cv_log'], 'NEUTRINO', header)
    sys_signal = read_log(sys['sys_log'], 'NEUTRINO', header)
    cv_selected = read_log(sys['cv_log'], f'SELECTED_{sys["channel"].upper()}', header).merge(cv_signal, how='inner', on=['run', 'subrun', 'event', 'nu_id'])
    sys_selected = read_log(sys['sys_log'], f'SELECTED_{sys["channel"].upper()}', header).merge(sys_signal, how='inner', on=['run', 'subrun', 'event', 'nu_id'])
    common_signal = cv_signal.merge(sys_signal, how='inner', on=['run', 'subrun', 'event', 'nu_id'])
    return common_signal, cv_selected, sys_selected

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
    merged = [s.merge(common.iloc[choice], how='inner', on=['run', 'subrun', 'event', 'nu_id']) for s in selected]
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

def calc_detector_response(sys, var, bins):
    """
    Calculates the covariance matrix for the detector variation using
    the CV and systematic variation samples.

    Parameters
    ----------
    sys: dict
        The dictionary containing the configuration details for the systematic
        variation.
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
    """
    common, cv_selected, sys_selected = load_detector_variation(sys)

    # Calculate the bin edges and the bin index for each interaction.
    _, bin_edges = np.histogram(cv_selected[var], bins=bins[0], range=bins[1:])
    nbins = len(bin_edges) - 1
    cv_selected['bidx'] = np.digitize(cv_selected[var], bin_edges) - 1
    sys_selected['bidx'] = np.digitize(sys_selected[var], bin_edges) - 1

    # Bootstrap the signal events to characterize the covariance of the bins.
    nboots = sys['nboots']
    bins = np.zeros((2, nbins, nboots), dtype=np.int16)
    for i in range(nboots):
        bins[:, :, i] = bootstrap_iterate(common, [cv_selected, sys_selected], nbins, stats_limit=1.0)

    # Calculate V_nominal and the associated covariance matrix (M_R).
    vnominal = np.mean(bins[1,:,:] - bins[0,:,:], axis=1)
    rmatrix = np.cov(bins[1,:,:] - bins[0,:,:] - vnominal[:,np.newaxis])

    # Mask bins which will cause a singular response matrix.
    rank = np.linalg.matrix_rank(rmatrix)
    singular_mask = (vnominal != 0)

    # Calculate the detector response matrix (M_D).
    nuniverses = sys['nuniverses']
    universes = np.zeros((rank, nuniverses), dtype=np.float64)
    for n in range(nuniverses):
        universes[:, n] = generate_detector_universe(rmatrix[singular_mask, :][:, singular_mask], vnominal[singular_mask])
    dmatrix = np.zeros((nbins, nbins), dtype=np.float64)
    cov = np.cov(universes)
    for bi, b in enumerate(np.arange(nbins)[singular_mask]):
        dmatrix[b, singular_mask] = cov[bi, :]
    cv = np.array([np.sum(cv_selected['bidx'] == b) for b in range(nbins)])

    return cv, vnominal, rmatrix, dmatrix