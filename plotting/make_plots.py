import yaml
import uproot
from plot_funcs import PlotDescription, plot_histogram_1d, plot_confusion, plot_flow

def main():
    cfg = """
    - type: hist1d
      name: particles
      var: [sCountParticles, sCountPrimaries]
      xlabel: Particles
      ylabel: Entries (Norm. Event)
      title: Interaction Multiplicity
      save: plots/count_particles.png
      plot_kwargs:
        histtype: barstacked
    - type: confusion
      name: pid
      var: sPID_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: PID Efficiency
      title: Particle PID Confusion
      entries: [Photon, Electron, Muon, Pion, Proton]
      save: plots/pid_confusion.png
    - type: flow
      pops: [1mu1p, OtherNu, Cosmic]
      labels: [1$\mu$1p, Other $\\nu$, Cosmic]
      cuts: [NoCut, FVCut, FVConCut, FVConTopCut, AllCut]
      clabels: [No Cut, Fiducial\nVolume, Containment, Topological, Flash Time]
      title: 1$\mu$1p Selection Statistics
      span: [-4,3]
      save: plots/flow_1mu1p.png
    """
    cfg = yaml.safe_load(cfg)

    rf = uproot.open('../spectra.root')

    for h in cfg:
        desc = PlotDescription(h)
        if desc.type == 'hist1d':
            plot_histogram_1d(rf, desc)
        elif desc.type == 'confusion':
            plot_confusion(rf, desc)
        elif desc.type == 'flow':
            plot_flow(rf, desc)

if __name__ == '__main__':
    main()
