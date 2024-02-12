import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
import argparse
import uproot
import toml
from tqdm import tqdm
from systools import read_log, extract_weights

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

def main(configuration, log, caf, output):
    header = [  'run', 'subrun', 'event', 'nu_id', 'category', 'category_topology', 'category_mode',
                'visible_energy', 'leading_muon_ke', 'leading_proton_ke', 'leading_muon_pt',
                'leading_proton_pt', 'interaction_pt']
    selected = read_log(log, 'SELECTED_1MU1P', header)

    cfg = dict()
    if isinstance(configuration[0], list):
        configuration = [c for y in configuration for c in y]
    for c in configuration:
        with open(c, 'r') as f:
            cfg = {**cfg, **toml.load(f)}
    variables = cfg['variables']
    cfg.pop('variables')
    allsys = {f'{ki}_{vj}': kj for ki, vi in cfg.items() for kj, vj in vi.items()}

    for var_key, var_bins in variables.items():
        _, e = np.histogram(selected[var_key], bins=var_bins[0], range=var_bins[1:])
        selected['bidx'] = np.digitize(selected[var_key], e) - 1

        covariances = dict()
        pbar = tqdm(allsys.items())
        for k, v in pbar:
            pbar.set_description(f'Systematic: "{k}"')
            weights = extract_weights(caf, selected, int(v))
            covariances[k] = calc_covariance(selected, weights, var_bins[0])

        total = np.zeros((variables[var_key][0], variables[var_key][0]))
        genie = np.zeros((variables[var_key][0], variables[var_key][0]))
        flux = np.zeros((variables[var_key][0], variables[var_key][0]))
        for k, v in covariances.items():
            total += v
            if 'genie' in k:
                genie += v
            elif 'flux' in k:
                flux += v

        stat = np.zeros(shape=(variables[var_key][0], variables[var_key][0]))
        for b in range(variables[var_key][0]):
            stat[b,b] = np.sum(selected['bidx'] == b)
        
        covariances['total'] = total + stat
        covariances['TotalNoStat'] = total
        covariances['genie'] = genie
        covariances['flux'] = flux
        covariances['statistical'] = stat
        np.savez(f'{output}_{var_key}.npz', **covariances)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log', default='../output.log')
    parser.add_argument('-w', '--weights')
    parser.add_argument('-c', '--config', action='append', nargs='*')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.config, args.log, args.weights, args.output)