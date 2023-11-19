# Introduction
This repository contains the code necessary for constructing a muon neutrino selection at the ICARUS experiment using the outputs of a machine learning reconstruction (see [lartpc_mlreco3d](https://github.com/DeepLearnPhysics/lartpc_mlreco3d)). The selection is implemented using the CAFAna framework (see [sbnana](https://github.com/SBNSoftware/sbnana)), which takes its name from the Common Analysis Files (CAFs) that contain the input analysis-level informatino necessary for developing a selection.

# CAFs

## CAF Format
The CAF format employs a hierarchical organization of C++ objects stored in ROOT files (traditionally using `.caf.root` as an extension). At a basic level, the hierarchy follows something similar to the below:
* `StandardRecord` - The top-level object representing a single "spill" or event.
    * `SRSlice` (slc) - A single reco interaction or "slice" as defined by Pandora (vector).
    * `SRInteractionDLP` (dlp) - The reco interactions as defined by the ML reconstruction (vector).
        * `SRParticleDLP` (particles) - The reco particles comprising the interaction (vector).
    * `SRInteractionTruthDLP` (dlp_true) - The true interactions as defined by the ML reconstruction (vector).
        * `SRParticleTruthDLP` (particles) - The true particles comprising the interaction (vector).

The CAF format contains many nested layers, and so is often considered to be a bit unwieldy or slow to use directly - a cost associated with navigating the complex structure. For this reason, a "flattened" version of the files is often used instead (traditionally using `.flat.root` as an extension) which broadcasts all branches to match the deepest level. This greatly simplifies the navigation and results in a significant speed up for any framework using them as input (e.g. CAFAna). The flattening is performed by an executable that ships with `sbnana` called `flatten_caf`.

## Generating CAFs
The analysis-level output of the machine learning reconstruction is stored in the HDF5 format. The advantage of the HDF5 format is its portability and the "self-describing" nature of the dataset format. The disadvantage is that it requires a little bit of work to get the analysis outputs into a CAF file. This functionality is implemented in [sbn_ml_cafmaker](https://github.com/justinjmueller/sbn_ml_cafmaker) and won't be described in great detail here. At a basic level, it reads the HDF5 input and organizes the truth and reco information in the new branches `dlp_true` and `dlp` within the `StandardRecord`. The resulting CAFs have been verified to work with CAFAna directly (albeit with some reduced functionality) and are able to be flattened using the `flatten_caf` executable. 