import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
import argparse
import uproot
import toml
from tqdm import tqdm
from systools import read_log, extract_weights

class PlotDescription:
    """
    Class container for plot configuration details.

    Attributes
    ----------
    name: str
        The name assigned to the plot which is also used as the name of
        the saved image file.
    output: str
        The full path specifying the save location for the plot.
    """
    def __init__(self, k, cfg, output):
        plt.style.use('../plotting/plot_style.mplstyle')
        self.name = k
        if output[-1] != '/':
            output = output + '/'
        self.save = f'{output}{self.name}.png'
        for k, v in cfg.items():
            setattr(self, k, v)

def plot_histogram_with_systematics(desc, selected, cov_path):
    """
    Plots a 1D histogram with the specified plot configurations for the
    input selected candidates.

    Parameters
    ----------
    desc: PlotDescription
        Contains the plot configuration details/settings from the input
        TOML configuration file.
    selected: pandas.DataFrame
        A DataFrame containing basic information about the selected
        neutrinos (run, subrun, event, nu_id) along with any binned
        variables (e.g. reconstructed energy).
    cov_path: str
        The full path to the numpy zip file containing the covariance
        matrix.

    Returns
    -------
    None.
    """
    figure = plt.figure(figsize=(8,6))
    ax = figure.add_subplot()

    contents = list()
    centers = list()
    labels = list()
    for m in desc.merge[::-1]:
        mask = np.isin(selected[desc.categorical_var], m)
        c, e = np.histogram(selected[desc.var][mask], bins=desc.bins[0], range=desc.bins[1:])
        contents.append(c)
        centers.append((e[1:] + e[:-1]) / 2.0)
        labels.append(desc.categories[m[0]])

    covariance = np.load(cov_path)
    error = np.sqrt(np.diagonal(covariance[desc.sys_var]))

    colors = [f'C{i}' for i in desc.colors][::-1]
    ax.hist(centers, weights=contents, bins=desc.bins[0], range=desc.bins[1:], label=labels, color=colors, histtype='barstacked')
    ax.errorbar(centers[0], np.sum(contents, axis=0), yerr=error, fmt='o', c='black')
    ax.set_xlim(*desc.bins[1:])
    ax.set_ylim(*desc.ylim)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1])
    ax.set_xlabel(desc.xlabel)
    ax.set_ylabel('Entries')
    figure.suptitle(f'Selected 1$\mu$1p Candidates\nError Source: {desc.sys_var.capitalize()}')
    figure.savefig(desc.save)
    plt.close(figure)

def main(configuration, plot, log, cov, output):
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
    plot_keys = list(cfg.keys()) if plot is None else plot
    
    pbar = tqdm(plot_keys)
    for k in pbar:
        pbar.set_description(f'Making plot "{k}"')
        desc = PlotDescription(k, cfg[k], output)
        plot_histogram_with_systematics(desc, selected, f'{cov}_{desc.var}.npz')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log', default='../output.log')
    parser.add_argument('-e', '--covariance')
    parser.add_argument('-c', '--config', action='append', nargs='*')
    parser.add_argument('-p', '--plot', action='append')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.config, args.plot, args.log, args.covariance, args.output)