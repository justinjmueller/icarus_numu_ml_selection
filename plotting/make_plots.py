import yaml
import uproot
import sys
from plot_funcs import PlotDescription, plot_histogram_1d, plot_histogram_2d, plot_confusion, plot_flow

import uproot
import sys
import argparse
from plot_funcs import PlotDescription, plot_histogram_1d, plot_histogram_2d, plot_confusion, plot_flow
import toml
from tqdm import tqdm
import time

def main(spectra_file, configuration, plot, output):
    """
    Main function that handles the plotting of spectra. Configuration
    of plots is done using an external TOML file (or multiple).

    Parameters
    ----------
    spectra_file: str
        The full path to the ROOT file containing the input spectra.
    configuration: list[str]
        The list of full paths to the TOML files containing the
        configuration details for the plots.
    plot: list[str]
        The list of plots within the configuration files to produce. If
        None, then all plots will be produced.
    output: str
        The full path to the plot output directory.

    Returns
    -------
    None.
    """
    cfg = dict()
    for c in configuration:
        with open(c, 'r') as f:
            cfg = cfg | toml.load(f)
    plot_keys = list(cfg.keys()) if plot is None else plot
    for k in (pbar := tqdm(plot_keys)):
        pbar.set_description(f'Making plot "{k}"')
        desc = PlotDescription(k, cfg[k], output)
        time.sleep(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--spectra', default='../spectra_combined.root')
    parser.add_argument('-c', '--config', action='append')
    parser.add_argument('-p', '--plot', action='append')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.spectra, args.config, args.plot, args.output)
