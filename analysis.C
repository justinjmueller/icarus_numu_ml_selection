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

    /**
     * Spectra (1D) for interactions.
    */
    spectra.add_spectrum1d("sCountParticles", Binning::Simple(20, 0, 20), kCountParticles);
    spectra.add_spectrum1d("sCountPrimaries", Binning::Simple(20, 0, 20), kCountPrimaries);
    spectra.add_spectrum1d("sCountParticlesTruth", Binning::Simple(20, 0, 20), kCountParticlesTruth);
    spectra.add_spectrum1d("sCountPrimariesTruth", Binning::Simple(20, 0, 20), kCountPrimariesTruth);

    /**
     * Spectra (1D) for counting selection statistics by category (efficiency).
    */
    spectra.add_spectrum1d("sCountTTP_1mu1p_NoCut", Binning::Simple(1, 0, 2), kCountTTP_1mu1p_NoCut);
    spectra.add_spectrum1d("sCountTTP_1mu1p_FVCut", Binning::Simple(1, 0, 2), kCountTTP_1mu1p_FVCut);
    spectra.add_spectrum1d("sCountTTP_1mu1p_FVConCut", Binning::Simple(1, 0, 2), kCountTTP_1mu1p_FVConCut);
    spectra.add_spectrum1d("sCountTTP_1mu1p_FVConTopCut", Binning::Simple(1, 0, 2), kCountTTP_1mu1p_FVConTopCut);
    spectra.add_spectrum1d("sCountTTP_1mu1p_AllCut", Binning::Simple(1, 0, 2), kCountTTP_1mu1p_AllCut);

    spectra.add_spectrum1d("sCountTTP_OtherNu_NoCut", Binning::Simple(1, 0, 2), kCountTTP_OtherNu_NoCut);
    spectra.add_spectrum1d("sCountTTP_OtherNu_FVCut", Binning::Simple(1, 0, 2), kCountTTP_OtherNu_FVCut);
    spectra.add_spectrum1d("sCountTTP_OtherNu_FVConCut", Binning::Simple(1, 0, 2), kCountTTP_OtherNu_FVConCut);
    spectra.add_spectrum1d("sCountTTP_OtherNu_FVConTopCut", Binning::Simple(1, 0, 2), kCountTTP_OtherNu_FVConTopCut);
    spectra.add_spectrum1d("sCountTTP_OtherNu_AllCut", Binning::Simple(1, 0, 2), kCountTTP_OtherNu_AllCut);

    spectra.add_spectrum1d("sCountTTP_Cosmic_NoCut", Binning::Simple(1, 0, 2), kCountTTP_Cosmic_NoCut);
    spectra.add_spectrum1d("sCountTTP_Cosmic_FVCut", Binning::Simple(1, 0, 2), kCountTTP_Cosmic_FVCut);
    spectra.add_spectrum1d("sCountTTP_Cosmic_FVConCut", Binning::Simple(1, 0, 2), kCountTTP_Cosmic_FVConCut);
    spectra.add_spectrum1d("sCountTTP_Cosmic_FVConTopCut", Binning::Simple(1, 0, 2), kCountTTP_Cosmic_FVConTopCut);
    spectra.add_spectrum1d("sCountTTP_Cosmic_AllCut", Binning::Simple(1, 0, 2), kCountTTP_Cosmic_AllCut);

    /**
     * Spectra (1D) for counting selection statistics by category (purity).
    */
    spectra.add_spectrum1d("sCountPTT_1mu1p_NoCut", Binning::Simple(1, 0, 2), kCountPTT_1mu1p_NoCut);
    spectra.add_spectrum1d("sCountPTT_1mu1p_FVCut", Binning::Simple(1, 0, 2), kCountPTT_1mu1p_FVCut);
    spectra.add_spectrum1d("sCountPTT_1mu1p_FVConCut", Binning::Simple(1, 0, 2), kCountPTT_1mu1p_FVConCut);
    spectra.add_spectrum1d("sCountPTT_1mu1p_FVConTopCut", Binning::Simple(1, 0, 2), kCountPTT_1mu1p_FVConTopCut);
    spectra.add_spectrum1d("sCountPTT_1mu1p_AllCut", Binning::Simple(1, 0, 2), kCountPTT_1mu1p_AllCut);

    spectra.add_spectrum1d("sCountPTT_OtherNu_NoCut", Binning::Simple(1, 0, 2), kCountPTT_OtherNu_NoCut);
    spectra.add_spectrum1d("sCountPTT_OtherNu_FVCut", Binning::Simple(1, 0, 2), kCountPTT_OtherNu_FVCut);
    spectra.add_spectrum1d("sCountPTT_OtherNu_FVConCut", Binning::Simple(1, 0, 2), kCountPTT_OtherNu_FVConCut);
    spectra.add_spectrum1d("sCountPTT_OtherNu_FVConTopCut", Binning::Simple(1, 0, 2), kCountPTT_OtherNu_FVConTopCut);
    spectra.add_spectrum1d("sCountPTT_OtherNu_AllCut", Binning::Simple(1, 0, 2), kCountPTT_OtherNu_AllCut);

    spectra.add_spectrum1d("sCountPTT_Cosmic_NoCut", Binning::Simple(1, 0, 2), kCountPTT_Cosmic_NoCut);
    spectra.add_spectrum1d("sCountPTT_Cosmic_FVCut", Binning::Simple(1, 0, 2), kCountPTT_Cosmic_FVCut);
    spectra.add_spectrum1d("sCountPTT_Cosmic_FVConCut", Binning::Simple(1, 0, 2), kCountPTT_Cosmic_FVConCut);
    spectra.add_spectrum1d("sCountPTT_Cosmic_FVConTopCut", Binning::Simple(1, 0, 2), kCountPTT_Cosmic_FVConTopCut);
    spectra.add_spectrum1d("sCountPTT_Cosmic_AllCut", Binning::Simple(1, 0, 2), kCountPTT_Cosmic_AllCut);

    /**
     * Spectra (1D) for particles.
    */
    spectra.add_spectrum1d("sCSDA", Binning::Simple(50, 0, 2500), kCSDA);
    spectra.add_spectrum1d("sCSDATruth", Binning::Simple(50, 0, 2500), kCSDATruth);

    /**
     * Spectra (1D) for matched (truth-to-predicted) particles.
    */
    spectra.add_spectrum1d("sPID_muon", Binning::Simple(5, 0, 5), kPID_muon);

    /**
     * Spectra (2D) for matched (truth-to-predicted) particles.
    */
    spectra.add_spectrum2d("sPID_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth, kPID);

    spectra.run();
}