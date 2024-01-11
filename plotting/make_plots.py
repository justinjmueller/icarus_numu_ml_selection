import yaml
import uproot
import sys
from plot_funcs import PlotDescription, plot_histogram_1d, plot_histogram_2d, plot_confusion, plot_flow

def main(args):
  cfg_file = open('plot_configurations.yml', 'r')
  cfg = yaml.safe_load(cfg_file)
  rf = uproot.open('../spectra.root')

  for h in cfg:
    desc = PlotDescription(h)
    desc.save = f'plots/{args[0] + "/" if args[0] != "" else ""}{desc.save}'
    if desc.type == 'hist1d':
      plot_histogram_1d(rf, desc)
    elif desc.type == 'hist2d':
      plot_histogram_2d(rf, desc)
    elif desc.type == 'confusion':
      plot_confusion(rf, desc)
    elif desc.type == 'flow':
      plot_flow(rf, desc)

if __name__ == '__main__':
  if len(sys.argv) == 2:
    main(sys.argv[1:])
  else:
    main(['',])
