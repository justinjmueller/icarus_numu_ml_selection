import yaml
import uproot
from plot_funcs import PlotDescription, plot_histogram_1d, plot_histogram_2d, plot_confusion, plot_flow

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
    - type: hist2d
      name: lowx
      var: sLowX
      xlabel: True Low X-Coordinate [cm]
      ylabel: Reco Low X-Coordinate [cm]
      title: Particle Lowest X-Coordinate
      save: plots/lowx.png
    - type: confusion
      name: primary
      var: sPrimary_confusion
      xlabel: True Particle Category
      ylabel: Reconstructed Particle Category
      clabel: Primary Efficiency
      title: Particle Primary/Non-primary Confusion
      entries: [Non-primary, Primary]
      save: plots/primary_confusion.png
    - type: confusion
      name: pid
      var: sPID_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: PID Efficiency
      title: Particle PID Confusion
      entries: [Photon, Electron, Muon, Pion, Proton]
      save: plots/pid_confusion.png
    - type: confusion
      name: primary_pid
      var: sPrimaryPID_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: Efficiency
      title: Particle Primary+PID Confusion
      entries: [sPhoton,sElectron,sMuon,sPion,sProton,pPhoton,pElectron,pMuon,pPion,pProton]
      save: plots/primary_pid_confusion.png
    - type: confusion
      name: primary_neutrino
      var: sPrimary_Neutrino_confusion
      xlabel: True Particle Category
      ylabel: Reconstructed Particle Category
      clabel: Primary Efficiency
      title: Particle Primary/Non-primary Confusion (Neutrinos)
      entries: [Non-primary, Primary]
      save: plots/primary_neutrino_confusion.png
    - type: confusion
      name: pid
      var: sPID_Neutrino_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: PID Efficiency
      title: Particle PID Confusion (Neutrinos)
      entries: [Photon, Electron, Muon, Pion, Proton]
      save: plots/pid_neutrino_confusion.png
    - type: confusion
      name: primary_pid_neutrino
      var: sPrimaryPID_Neutrino_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: Efficiency
      title: Particle Primary+PID Confusion (Neutrinos)
      entries: [nPhoton,nElectron,nMuon,nPion,nProton,pPhoton,pElectron,pMuon,pPion,pProton]
      save: plots/primary_pid_neutrino_confusion.png
    - type: confusion
      name: primary_cosmic
      var: sPrimary_Cosmic_confusion
      xlabel: True Particle Category
      ylabel: Reconstructed Particle Category
      clabel: Primary Efficiency
      title: Particle Primary/Non-primary Confusion (Cosmics)
      entries: [Non-primary, Primary]
      save: plots/primary_cosmic_confusion.png
    - type: confusion
      name: pid_cosmic
      var: sPID_Cosmic_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: PID Efficiency
      title: Particle PID Confusion (Cosmics)
      entries: [Photon, Electron, Muon, Pion, Proton]
      save: plots/pid_cosmic_confusion.png
    - type: confusion
      name: primary_pid_cosmic
      var: sPrimaryPID_Cosmic_confusion
      xlabel: True Particle Type
      ylabel: Reconstructed Particle Type
      clabel: Efficiency
      title: Particle Primary+PID Confusion (Cosmics)
      entries: [nPhoton,nElectron,nMuon,nPion,nProton,pPhoton,pElectron,pMuon,pPion,pProton]
      save: plots/primary_pid_cosmic_confusion.png
    - type: flow
      direction: ttp
      pops: [0, 1, 2]
      labels: [1$\mu$1p, Other $\\nu$, Cosmic]
      cuts: [NoCut, FVCut, FVConCut, FVConTopCut, AllCut]
      clabels: [No Cut, Fiducial\nVolume, Containment, Topological, Flash Time]
      title: 1$\mu$1p Selection Statistics
      span: [-4,3]
      save: plots/flow_1mu1p_efficiency.png
    - type: flow
      direction: ptt
      pops: [0, 1, 2]
      labels: [1$\mu$1p, Other $\\nu$, Cosmic]
      cuts: [NoCut, FVCut, FVConCut, FVConTopCut, AllCut]
      clabels: [No Cut, Fiducial\nVolume, Containment, Topological, Flash Time]
      title: 1$\mu$1p Selection Statistics
      span: [-4,3]
      save: plots/flow_1mu1p_purity.png
    """
    cfg = yaml.safe_load(cfg)

    rf = uproot.open('../spectra.root')

    for h in cfg:
        desc = PlotDescription(h)
        if desc.type == 'hist1d':
            plot_histogram_1d(rf, desc)
        elif desc.type == 'hist2d':
            plot_histogram_2d(rf, desc)
        elif desc.type == 'confusion':
            plot_confusion(rf, desc)
        elif desc.type == 'flow':
            plot_flow(rf, desc)

if __name__ == '__main__':
    main()
