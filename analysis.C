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
     * Spectra (2D) for counting selection statistics by interaction categorization (efficiency).
    */
    spectra.add_spectrum2d("sCountTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_NoCut, kCountTTP_NoCut);
    spectra.add_spectrum2d("sCountTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVCut, kCountTTP_FVCut);
    spectra.add_spectrum2d("sCountTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConCut, kCountTTP_FVConCut);
    spectra.add_spectrum2d("sCountTTP_FVConTopCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTopCut, kCountTTP_FVConTopCut);
    spectra.add_spectrum2d("sCountTTP_AllCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_AllCut, kCountTTP_AllCut);

    /**
     * Spectra (2D) for visible energy.
    */
    spectra.add_spectrum2d("sVisibleEnergyTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(10, 0, 2500), kCategoryTTP_NoCut, kVisibleEnergyTTP_NoCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(10, 0, 2500), kCategoryTTP_FVCut, kVisibleEnergyTTP_FVCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(10, 0, 2500), kCategoryTTP_FVConCut, kVisibleEnergyTTP_FVConCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTopCut", Binning::Simple(10, 0, 10), Binning::Simple(10, 0, 2500), kCategoryTTP_FVConTopCut, kVisibleEnergyTTP_FVConTopCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_AllCut", Binning::Simple(10, 0, 10), Binning::Simple(10, 0, 2500), kCategoryTTP_AllCut, kVisibleEnergyTTP_AllCut);

    /**
     * Spectra (2D) for counting selection statistics by interaction categorization (purity).
    */
    spectra.add_spectrum2d("sCountPTT_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_NoCut, kCountPTT_NoCut);
    spectra.add_spectrum2d("sCountPTT_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVCut, kCountPTT_FVCut);
    spectra.add_spectrum2d("sCountPTT_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConCut, kCountPTT_FVConCut);
    spectra.add_spectrum2d("sCountPTT_FVConTopCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTopCut, kCountPTT_FVConTopCut);
    spectra.add_spectrum2d("sCountPTT_AllCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_AllCut, kCountPTT_AllCut);

    /**
     * Spectra (1D) for particles.
    */
    spectra.add_spectrum1d("sCSDA", Binning::Simple(50, 0, 2500), kCSDA);
    spectra.add_spectrum1d("sCSDATruth", Binning::Simple(50, 0, 2500), kCSDATruth);

    /**
     * Spectra (2D) for matched (truth-to-predicted) particles.
    */
    spectra.add_spectrum2d("sPrimary_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth, kPrimary);
    spectra.add_spectrum2d("sPID_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth, kPID);
    spectra.add_spectrum2d("sPrimaryPID_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth, kPrimaryPID);

    spectra.run();
}