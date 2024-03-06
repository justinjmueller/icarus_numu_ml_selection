/**
 * @file analysis.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author justin.mueller@colostate.edu
*/

#include "include/analysis.h"
#include "include/container.h"
#include "include/csv_maker.h"
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
    /**
     * 1. BNB neutrino (full flux) + out-of-time cosmics (v09_63_01).
     * 2. BNB in-time cosmics + out-of-time cosmics (v09_63_01).
    */
    SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/bnb_nucosmics_v6.flat.root", "spectra_nucosmics.root", 1.253e19, 2.5e20);
    //SpecContainer spectra("/exp/icarus/data/users/mueller/mlcafs/bnb_intime.flat.root", "spectra_intime.root", 9070*2.05e14, 2.5e20);

    /**
     * 3. BNB neutrino (full flux) + out-of-time cosmics *     Central Value    * (v09_82_02_01).
     * 4. BNB neutrino (full flux) + out-of-time cosmics * Coherent Noise +4.5% * (v09_82_02_01).
     * 5. BNB neutrino (full flux) + out-of-time cosmics *  Elli. Recombination * (v09_82_02_01).
     * 6. BNB neutrino (full flux) + out-of-time cosmics * Untuned Signal Shape * (v09_82_02_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_cv_v2.flat.root", "spectra_cv.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_tpcnoise_coh_p1_v2.flat.root", "spectra_tpcnoise_coh_p1.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_recombination.flat.root", "spectra_recombination.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_untunedsigshape.flat.root", "spectra_untunedsigshape.root", -1, 2.5e20);
    
    /**
     * 7. BNB neutrino-only (full flux)  *     Central Value    * (v09_82_02_01).
     * 8. BNB neutrino-only (full flux)  * Coherent Noise +4.5% * (v09_82_02_01).
     * 9. BNB neutrino-only (full flux)  *  Elli. Recombination * (v09_82_02_01).
     * 10. BNB neutrino-only (full flux) * Untuned Signal Shape * (v09_82_02_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.flat.root", "spectra_cv.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cohnoise.flat.root", "spectra_tpcnoise_coh_p1.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_intnoise.flat.root", "spectra_intnoise.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_recombination.flat.root", "spectra_recombination.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_sigshape.flat.root", "spectra_untunedsigshape.root", -1, 2.5e20);
    
    /**
     * 11. MPV/MPR sample (v09_63_00).
    */
    //SpecContainer spectra("/exp/icarus/data/users/mueller/mlcafs/mpv_mpr.flat.root", "spectra_mpvmpr.root", 1e19, 2.5e20);

    /**
     * Spectra (1D) for interactions.
    */
    spectra.add_spectrum1d("sCountParticles", Binning::Simple(20, 0, 20), kCountParticles);
    spectra.add_spectrum1d("sCountPrimaries", Binning::Simple(20, 0, 20), kCountPrimaries);
    spectra.add_spectrum1d("sCountParticlesTruth", Binning::Simple(20, 0, 20), kCountParticlesTruth);
    spectra.add_spectrum1d("sCountPrimariesTruth", Binning::Simple(20, 0, 20), kCountPrimariesTruth);
    spectra.add_spectrum1d("sEnergy_1mu1p_signal_bias", Binning::Simple(50,-1,1), kEnergy_1mu1p_signal_bias);
    spectra.add_spectrum1d("sEnergy_1mu1p_othernu_bias", Binning::Simple(50,-1,1), kEnergy_1mu1p_othernu_bias);
    spectra.add_spectrum1d("sEnergy_1mu1p_cosmic_bias", Binning::Simple(50,-1,1), kEnergy_1mu1p_cosmic_bias);
    spectra.add_spectrum1d("sEnergy_1muNp_1p_signal_bias", Binning::Simple(50,-1,1), kEnergy_1muNp_1p_signal_bias);
    spectra.add_spectrum1d("sEnergy_1muNp_Np_signal_bias", Binning::Simple(50,-1,1), kEnergy_1muNp_Np_signal_bias);
    spectra.add_spectrum1d("sEnergy_1muNp_othernu_bias", Binning::Simple(50,-1,1), kEnergy_1muNp_othernu_bias);
    spectra.add_spectrum1d("sEnergy_1muNp_cosmic_bias", Binning::Simple(50,-1,1), kEnergy_1muNp_cosmic_bias);

    spectra.add_spectrum1d("sNuEnergy_1mu1p_signal_bias", Binning::Simple(50,-1,1), kNuEnergy_1mu1p_signal_bias);
    spectra.add_spectrum1d("sNuEnergy_1mu1p_othernu_bias", Binning::Simple(50,-1,1), kNuEnergy_1mu1p_othernu_bias);
    spectra.add_spectrum1d("sNuEnergy_1muNp_1p_signal_bias", Binning::Simple(50,-1,1), kNuEnergy_1muNp_1p_signal_bias);
    spectra.add_spectrum1d("sNuEnergy_1muNp_Np_signal_bias", Binning::Simple(50,-1,1), kNuEnergy_1muNp_Np_signal_bias);
    spectra.add_spectrum1d("sNuEnergy_1muNp_othernu_bias", Binning::Simple(50,-1,1), kNuEnergy_1muNp_othernu_bias);

    /**
     * Spectra (2D) for counting selection statistics by interaction categorization (efficiency).
    */
    spectra.add_spectrum2d("sCountTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_NoCut, kCountTTP_NoCut);
    spectra.add_spectrum2d("sCountTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVCut, kCountTTP_FVCut);
    spectra.add_spectrum2d("sCountTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConCut, kCountTTP_FVConCut);
    spectra.add_spectrum2d("sCountTTP_FVConTop1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1mu1pCut, kCountTTP_FVConTop1mu1pCut);
    spectra.add_spectrum2d("sCountTTP_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1mu1pCut, kCountTTP_All1mu1pCut);
    spectra.add_spectrum2d("sCountTTP_FVConTop1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1muNpCut, kCountTTP_FVConTop1muNpCut);
    spectra.add_spectrum2d("sCountTTP_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1muNpCut, kCountTTP_All1muNpCut);
    spectra.add_spectrum2d("sCountTTP_FVConTop1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1muXCut, kCountTTP_FVConTop1muXCut);
    spectra.add_spectrum2d("sCountTTP_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1muXCut, kCountTTP_All1muXCut);

    /**
     * Spectra (2D) for counting selection statistics by interaction categorization (purity).
    */
    spectra.add_spectrum2d("sCountPTT_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_NoCut, kCountPTT_NoCut);
    spectra.add_spectrum2d("sCountPTT_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVCut, kCountPTT_FVCut);
    spectra.add_spectrum2d("sCountPTT_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConCut, kCountPTT_FVConCut);
    spectra.add_spectrum2d("sCountPTT_FVConTop1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1mu1pCut, kCountPTT_FVConTop1mu1pCut);
    spectra.add_spectrum2d("sCountPTT_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1mu1pCut, kCountPTT_All1mu1pCut);
    spectra.add_spectrum2d("sCountPTT_FVConTop1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1muNpCut, kCountPTT_FVConTop1muNpCut);
    spectra.add_spectrum2d("sCountPTT_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1muNpCut, kCountPTT_All1muNpCut);
    spectra.add_spectrum2d("sCountPTT_FVConTop1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1muXCut, kCountPTT_FVConTop1muXCut);
    spectra.add_spectrum2d("sCountPTT_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1muXCut, kCountPTT_All1muXCut);

    /**
     * Spectra (2D) for visible energy.
    */
    spectra.add_spectrum2d("sVisibleEnergyTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_NoCut, kVisibleEnergyTTP_NoCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVCut, kVisibleEnergyTTP_FVCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConCut, kVisibleEnergyTTP_FVConCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTop1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConTop1mu1pCut, kVisibleEnergyTTP_FVConTop1mu1pCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1muNpCut, kVisibleEnergyTTP_All1muNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTop1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConTop1muNpCut, kVisibleEnergyTTP_FVConTop1muNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1muNpCut, kVisibleEnergyTTP_All1muNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTop1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConTop1muXCut, kVisibleEnergyTTP_FVConTop1muXCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1muXCut, kVisibleEnergyTTP_All1muXCut);

    /**
     * Spectra (2D) for flash time.
    */
    spectra.add_spectrum2d("sFlashTime_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_NoCut, kFlashTimePTT_NoCut);
    spectra.add_spectrum2d("sFlashTime_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVCut, kFlashTimePTT_FVCut);
    spectra.add_spectrum2d("sFlashTime_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConCut, kFlashTimePTT_FVConCut);
    spectra.add_spectrum2d("sFlashTime_FVConTop1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConTop1mu1pCut, kFlashTimePTT_FVConTop1mu1pCut);
    spectra.add_spectrum2d("sFlashTime_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_All1mu1pCut, kFlashTimePTT_All1mu1pCut);
    spectra.add_spectrum2d("sFlashTime_FVConTop1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConTop1muNpCut, kFlashTimePTT_FVConTop1muNpCut);
    spectra.add_spectrum2d("sFlashTime_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_All1muNpCut, kFlashTimePTT_All1muNpCut);
    spectra.add_spectrum2d("sFlashTime_FVConTop1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConTop1muXCut, kFlashTimePTT_FVConTop1muXCut);
    spectra.add_spectrum2d("sFlashTime_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_All1muXCut, kFlashTimePTT_All1muXCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_NoCut, kFlashTimePTT_NoCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_FVCut, kFlashTimePTT_FVCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_FVConCut, kFlashTimePTT_FVConCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConTop1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_FVConTop1mu1pCut, kFlashTimePTT_FVConTop1mu1pCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_All1mu1pCut, kFlashTimePTT_All1mu1pCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConTop1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_FVConTop1muNpCut, kFlashTimePTT_FVConTop1muNpCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_All1muNpCut, kFlashTimePTT_All1muNpCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConTop1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_FVConTop1muXCut, kFlashTimePTT_FVConTop1muXCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 4), kCategoryPTT_All1muXCut, kFlashTimePTT_All1muXCut);

    /**
     * Spectra (2D) for (stacked) reconstructed quantities.
    */
    spectra.add_spectrum2d("sFlashTimePTT_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 5.6), kCategoryTopologyPTT_NoCut, kFlashTimePTT_NoCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1mu1pCut, kVisibleEnergyPTT_All1mu1pCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_InteractionMode_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModePTT_All1mu1pCut, kVisibleEnergyPTT_All1mu1pCut);
    spectra.add_spectrum2d("sFlashTimePTT_Topology_All1mu1pCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 5.6), kCategoryTopologyPTT_All1mu1pCut, kFlashTimePTT_All1mu1pCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1muNpCut, kVisibleEnergyPTT_All1muNpCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_InteractionMode_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModePTT_All1muNpCut, kVisibleEnergyPTT_All1muNpCut);
    spectra.add_spectrum2d("sFlashTimePTT_Topology_All1muNpCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 5.6), kCategoryTopologyPTT_All1muNpCut, kFlashTimePTT_All1muNpCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryPTT_All1muXCut, kVisibleEnergyPTT_All1muXCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1muXCut, kVisibleEnergyPTT_All1muXCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_InteractionMode_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModePTT_All1muXCut, kVisibleEnergyPTT_All1muXCut);
    spectra.add_spectrum2d("sFlashTimePTT_Topology_All1muXCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 5.6), kCategoryTopologyPTT_All1muXCut, kFlashTimePTT_All1muXCut);

    /**
     * Spectra (2D) for particles.
    */
    spectra.add_spectrum2d("sCSDA_muon", Binning::Simple(50, 0, 1000), Binning::Simple(50, 0, 1000), kCSDATruth_muon, kCSDA_muon);
    spectra.add_spectrum2d("sCSDA_muon2muon", Binning::Simple(50, 0, 1000), Binning::Simple(50, 0, 1000), kCSDATruth_muon, kCSDA_muon2muon);
    spectra.add_spectrum2d("sCSDA_muon_bias2d", Binning::Simple(10, 0, 1000), Binning::Simple(250,-0.25,0.25), kCSDATruth_muon, kCSDA_muon_bias);
    spectra.add_spectrum1d("sCSDA_muon_bias", Binning::Simple(75,-1,1), kCSDA_muon_bias);
    spectra.add_spectrum1d("sCSDA_noncc_muon_bias", Binning::Simple(75,-1,1), kCSDA_noncc_muon_bias);
    spectra.add_spectrum1d("sCSDA_wellreco_muon_bias", Binning::Simple(75,-1,1), kCSDA_wellreco_muon_bias);
    spectra.add_spectrum1d("sCCOverlap", Binning::Simple(50, 0, 1), kCCOverlap);
    spectra.add_spectrum1d("sNonCCOverlap", Binning::Simple(50, 0, 1), kNonCCOverlap);

    /**
     * Spectra (2D) for matched (truth-to-predicted) particles.
    */
    spectra.add_spectrum2d("sPrimary_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth, kPrimary);
    spectra.add_spectrum2d("sPID_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth, kPID);
    spectra.add_spectrum2d("sPrimaryPID_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth, kPrimaryPID);
    spectra.add_spectrum2d("sPrimary_Neutrino_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth_Neutrino, kPrimary_Neutrino);
    spectra.add_spectrum2d("sPID_Neutrino_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth_Neutrino, kPID_Neutrino);
    spectra.add_spectrum2d("sPrimaryPID_Neutrino_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth_Neutrino, kPrimaryPID_Neutrino);
    spectra.add_spectrum2d("sPrimary_Cosmic_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth_Cosmic, kPrimary_Cosmic);
    spectra.add_spectrum2d("sPID_Cosmic_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth_Cosmic, kPID_Cosmic);
    spectra.add_spectrum2d("sPrimaryPID_Cosmic_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth_Cosmic, kPrimaryPID_Cosmic);
    spectra.add_spectrum2d("sLowX", Binning::Simple(100,-400,400), Binning::Simple(100,-400,400), kLowX, kLowXTruth);

    spectra.add_spectrum2d("sPrimaryWellReco_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryWellRecoTruth, kPrimaryWellReco);
    spectra.add_spectrum2d("sPIDWellReco_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDWellRecoTruth, kPIDWellReco);
    spectra.add_spectrum2d("sPrimaryPIDWellReco_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDWellRecoTruth, kPrimaryPIDWellReco);

    spectra.add_spectrum2d("sPrimaryWellReco_Neutrino_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryWellRecoTruth_Neutrino, kPrimaryWellReco_Neutrino);
    spectra.add_spectrum2d("sPIDWellReco_Neutrino_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDWellRecoTruth_Neutrino, kPIDWellReco_Neutrino);
    spectra.add_spectrum2d("sPrimaryPIDWellReco_Neutrino_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDWellRecoTruth_Neutrino, kPrimaryPIDWellReco_Neutrino);

    /**
     * Spectra (2D) for correlating truth quantities.
    */
    spectra.add_spectrum2d("sScatteringProtonOverlap", Binning::Simple(50, 0.25, 1), Binning::Simple(25, 0, 1), kProtonScattering, kLeadingProtonOverlap);

    /**
     * Dummy spectra for dumping particle-level information to a CSV log file.
    */
    spectra.add_spectrum1d("sSelected", Binning::Simple(1, 0, 2), kSelected);
    spectra.add_spectrum1d("sSignal", Binning::Simple(1, 0, 2), kSignal);

    spectra.run();
}