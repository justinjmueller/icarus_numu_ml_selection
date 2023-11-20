/**
 * @file analysis.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author justin.mueller@colostate.edu
*/

#include "include/analysis.h"
#include "include/container.h"
#include "sbnana/CAFAna/Core/Binning.h"

using namespace ana;

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/
void analysis()
{
    SpecContainer spectra("/exp/icarus/app/users/mueller/sbn_ml_cafmaker/sbn_ml_cafmaker/build/bnb_nucosmics.flat.root", "spectra.root");

    spectra.add_spectrum("sCountParticles", Binning::Simple(20, 0, 20), kCountParticles);
    spectra.add_spectrum("sCountPrimaries", Binning::Simple(20, 0, 20), kCountPrimaries);
    spectra.add_spectrum("sCountParticlesTruth", Binning::Simple(20, 0, 20), kCountParticlesTruth);
    spectra.add_spectrum("sCountPrimariesTruth", Binning::Simple(20, 0, 20), kCountPrimariesTruth);

    spectra.add_spectrum("sCount_1mu1p_NoCut", Binning::Simple(1, 0, 2), kCount_1mu1p_NoCut);
    spectra.add_spectrum("sCount_1mu1p_FVCut", Binning::Simple(1, 0, 2), kCount_1mu1p_FVCut);
    spectra.add_spectrum("sCount_1mu1p_FVConCut", Binning::Simple(1, 0, 2), kCount_1mu1p_FVConCut);
    spectra.add_spectrum("sCount_1mu1p_FVConTopCut", Binning::Simple(1, 0, 2), kCount_1mu1p_FVConTopCut);
    spectra.add_spectrum("sCount_1mu1p_AllCut", Binning::Simple(1, 0, 2), kCount_1mu1p_AllCut);

    spectra.add_spectrum("sCount_OtherNu_NoCut", Binning::Simple(1, 0, 2), kCount_OtherNu_NoCut);
    spectra.add_spectrum("sCount_OtherNu_FVCut", Binning::Simple(1, 0, 2), kCount_OtherNu_FVCut);
    spectra.add_spectrum("sCount_OtherNu_FVConCut", Binning::Simple(1, 0, 2), kCount_OtherNu_FVConCut);
    spectra.add_spectrum("sCount_OtherNu_FVConTopCut", Binning::Simple(1, 0, 2), kCount_OtherNu_FVConTopCut);
    spectra.add_spectrum("sCount_OtherNu_AllCut", Binning::Simple(1, 0, 2), kCount_OtherNu_AllCut);

    spectra.add_spectrum("sCount_Cosmic_NoCut", Binning::Simple(1, 0, 2), kCount_Cosmic_NoCut);
    spectra.add_spectrum("sCount_Cosmic_FVCut", Binning::Simple(1, 0, 2), kCount_Cosmic_FVCut);
    spectra.add_spectrum("sCount_Cosmic_FVConCut", Binning::Simple(1, 0, 2), kCount_Cosmic_FVConCut);
    spectra.add_spectrum("sCount_Cosmic_FVConTopCut", Binning::Simple(1, 0, 2), kCount_Cosmic_FVConTopCut);
    spectra.add_spectrum("sCount_Cosmic_AllCut", Binning::Simple(1, 0, 2), kCount_Cosmic_AllCut);

    spectra.add_spectrum("sPID_muon", Binning::Simple(5, 0, 5), kPID_muon);

    spectra.run();
}