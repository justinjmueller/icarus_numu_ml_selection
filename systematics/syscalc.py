import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
import argparse
import uproot
import toml
from tqdm import tqdm
from systools import read_log, extract_weights, calc_covariance, calc_detector_response

def main(configuration, log, caf, output):
    cfg = dict()
    if isinstance(configuration[0], list):
        configuration = [c for y in configuration for c in y]
    for c in configuration:
        with open(c, 'r') as f:
            cfg = {**cfg, **toml.load(f)}
    variables = cfg['general']['variables']

    header = ['run', 'subrun', 'event', 'nu_id'] + cfg['general']['columns']
    selected = read_log(log, 'SELECTED_1MU1P', header)

    covariances = dict()
    for var_key, var_bins in variables.items():
        total = np.zeros((var_bins[0], var_bins[0]))
        total_no_stat = np.zeros((var_bins[0], var_bins[0]))
        _, e = np.histogram(selected[var_key], bins=var_bins[0], range=var_bins[1:])
        selected['bidx'] = np.digitize(selected[var_key], e) - 1
        for name, sys in cfg['sys'].items():
            if sys['type'] == 'multisim':
                systot = np.zeros((var_bins[0], var_bins[0]))
                pbar = tqdm([(k, v) for k,v in sys.items() if k != 'type'])
                for k, v in pbar:
                    pbar.set_description(f'Systematic: "{v}"')
                    weights = extract_weights(caf, selected, int(k))
                    covariances[f'{name}_{var_key}_{v}'] = calc_covariance(selected, weights, var_bins[0])
                    systot += covariances[f'{name}_{var_key}_{v}']
                covariances[f'{name}_{var_key}'] = systot
                total += covariances[f'{name}_{var_key}']
                total_no_stat += covariances[f'{name}_{var_key}']
            elif sys['type'] == 'stat':
                covariances[f'statistical_{var_key}'] = np.zeros((var_bins[0], var_bins[0]))
                for b in range(variables[var_key][0]):
                    covariances[f'statistical_{var_key}'][b,b] = np.sum(selected['bidx'] == b)    
                total += covariances[f'statistical_{var_key}']
            elif sys['type'] == 'detector':
                c, v, r, d = calc_detector_response(sys, var_key, var_bins)
                covariances[f'{name}_{var_key}'] = d

    np.savez(f'{output}covariances.npz', **covariances)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log', default='../output.log')
    parser.add_argument('-w', '--weights')
    parser.add_argument('-c', '--config', action='append', nargs='*')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.config, args.log, args.weights, args.output)