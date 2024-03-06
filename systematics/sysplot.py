import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch
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

def add_error_boxes(ax, x, y, xerr, yerr, **kwargs):
    """
    Adds error boxes to the input axis.

    Parameters
    ----------
    ax: matplotlib.axes.Axes
        The axis to which the error boxes are to be added.
    x: numpy.array
        The x-coordinates of the error boxes.
    y: numpy.array
        The y-coordinates of the error boxes.
    xerr: numpy.array
        The x-error values of the error boxes.
    yerr: numpy.array
        The y-error values of the error boxes.
    kwargs: dict
        Keyword arguments to be passed to the errorbar function.

    Returns
    -------
    None.
    """
    boxes = [Rectangle((x[i] - xerr[i], y[i] - yerr[i]), 2 * np.abs(xerr[i]), 2 * yerr[i]) for i in range(len(x))]
    pc = PatchCollection(boxes, **kwargs)
    ax.add_collection(pc)
    return boxes[0]

def plot_histogram(desc, selected, cov):
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
    cov: str
        The full path to the numpy zip file containing the covariance
        matrix.

    Returns
    -------
    figure: matplotlib.figure.Figure
        The figure object containing the histogram plot.
    """
    if desc.multiplot != 'none':
        figure = plt.figure(figsize=(8,8))
        gspec = figure.add_gridspec(2, 1, height_ratios=[4, 1])
        rax = figure.add_subplot(gspec[1])
        ax = figure.add_subplot(gspec[0], sharex=rax)
    else:
        figure = plt.figure(figsize=(8,6))
        ax = figure.add_subplot()

    contents = list()
    centers = list()
    labels = list()
    width = list()
    for m in desc.merge[::-1]:
        mask = np.isin(selected[desc.categorical_var], m)
        c, e = np.histogram(selected[desc.var][mask], bins=desc.bins[0], range=desc.bins[1:])
        contents.append(c)
        centers.append((e[1:] + e[:-1]) / 2.0)
        width.append(np.diff(e))
        labels.append(desc.categories[m[0]])

    colors = [f'C{i}' for i in desc.colors][::-1]
    ax.hist(centers, weights=contents, bins=desc.bins[0], range=desc.bins[1:], label=labels, color=colors, histtype='barstacked')
    h, l = ax.get_legend_handles_labels()
    show_percentage = hasattr(desc, 'show_percentage') and desc.show_percentage
    if show_percentage:
        l = [f'{l} ({np.sum(contents[li]):.0f}, {np.sum(contents[li]) / np.sum(contents):.02%})'for li, l in enumerate(l)]
    else:
        l = [f'{l} ({np.sum(contents[li]):.0f})'for li, l in enumerate(l)]
    h = h[::-1]
    l = l[::-1]

    covariance = np.load(cov)
    ccycle = iter(['gray', 'rebeccapurple'])
    hcycle = iter(['///', ''])
    for syskey, sysname in desc.systematics.items():
        error = np.sqrt(np.diagonal(covariance[f'{syskey}_{desc.var}']))
        c = next(ccycle)
        a = next(hcycle)
        add_error_boxes(ax, centers[0], np.sum(contents, axis=0), width[0] / 2, error, facecolor=c, edgecolor='none', alpha=0.5, hatch=a)
        h.append(Patch(facecolor=c, edgecolor=c, alpha=0.5, hatch=a))
        l.append(sysname)
    
    ax.legend(h, l)

    if desc.multiplot == 'ratio':
        vnominal = covariance[f'{syskey}_{desc.var}_vnominal']
        singular_mask = vnominal != 0
        chi2 = vnominal[singular_mask] @ np.linalg.inv(covariance[f'{syskey}_{desc.var}_rmatrix'][singular_mask,:][:,singular_mask]) @ vnominal[singular_mask].transpose()
        #print(chi2, 1-stats.chi2.cdf(chi2, np.sum(singular_mask)))

        sys = list(desc.systematics.items())
        ratio = covariance[f'{sys[0][0]}_{desc.var}_ratio']
        undefined = (ratio == 1)
        error = np.sqrt(np.diagonal(covariance[f'{sys[0][0]}_{desc.var}_cratio']))
        rax.errorbar(x=centers[0][~undefined], y=ratio[~undefined], yerr=error[~undefined], c='black', fmt='o')
        rax.axhline(1, c='black', ls='--')
        rax.set_ylim(0.8, 1.2)
        rax.set_xlim(*desc.bins[1:])
        rax.set_xlabel(desc.xlabel)
        rax.set_ylabel('Sys. / CV')
        plt.setp(ax.get_xticklabels(), visible=False)
    elif desc.multiplot == 'error':
        error = np.sqrt(np.diagonal(covariance[f'{syskey}_{desc.var}']))
        total = np.sum(contents, axis=0)
        rax.errorbar(centers[0], 100*np.divide(error, total, where=total!=0), xerr=width[0] / 2, yerr=0, fmt='_', c='black')
        rax.set_xlim(*desc.bins[1:])
        rax.set_xlabel(desc.xlabel)
        rax.set_ylim(0, 50)
        rax.set_ylabel('Rel. Error [%]')
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        ax.set_xlim(*desc.bins[1:])
        ax.set_xlabel(desc.xlabel)

    ax.set_ylim(*desc.ylim)
    ax.set_ylabel(desc.ylabel)
    ax.set_title(desc.title)
    figure.savefig(desc.save)

    return figure

def main(configuration, plots, log, cov_path, output):
    """
    Main function for making systematic-related plots.

    Parameters
    ----------
    configuration: str
        The full path to the TOML configuration file containing the
        plot settings.
    plots: list[str]
        The list of the plots to be made. If None, all plots in the
        configuration file are made.
    log: str
        The full path to the log file containing the selected events.
    cov_path: str
        The full path to the numpy zip file containing the covariance
        matrix.
    output: str
        The full path specifying the save location for the plot.
    """
    # Read the configuration file and extract the plot settings
    cfg = toml.load(configuration)
    header = cfg['metadata']['columns'] + list(cfg['variables'].keys())
    all = {f'{plot_name}_{var_name}': dict(**plot_cfg, **var_cfg) for plot_name, plot_cfg in cfg['plots'].items() for var_name, var_cfg in cfg['variables'].items()}
    if plots != None:
        all = {k: v for k, v in all.items() if k in plots}

    # Loop over the configured plots and make them
    pbar = tqdm(all.items())
    for k, v in pbar:
        pbar.set_description(f'{"Making plot `" + k + "`":^45}')
        selected = read_log(log, f'SELECTED_{v["channel"].upper()}', header)
        cov = f'{cov_path}covariances_{v["channel"]}.npz'
        v['ylim'] = v['ylim'][v['channel']]
        plot_histogram(PlotDescription(k, v, output), selected, cov)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log', default='../output.log')
    parser.add_argument('-e', '--covariance')
    parser.add_argument('-c', '--config')
    parser.add_argument('-p', '--plot', action='append')
    parser.add_argument('-o', '--output', default='./')
    args = parser.parse_args()
    main(args.config, args.plot, args.log, args.covariance, args.output)