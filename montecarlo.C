/**
 * @file montecarlo.C
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

void montecarlo()
{
    /**
     * Available samples for the analysis:
     * 
     * Central Value
     * 1. BNB neutrino (full flux) + out-of-time cosmics (v09_84_00_01)
     * 
     * Systematics
     * 2. BNB neutrino (full flux) no cosmics (v09_89_01_01) [CV]
     * 3. BNB neutrino (full flux) no cosmics (v09_89_01_01) [TPC Untuned Signal Shape]
     * 4. BNB neutrino (full flux) no cosmics (v09_89_01_01) [TPC Ind2 Opaque]
     * 5. BNB neutrino (full flux) no cosmics (v09_89_01_01) [TPC Ind2 Transparent]
     * 6. BNB neutrino (full flux) no cosmics (v09_89_01_01) [TPC Ind1 Increase Gain]
     * 7. BNB neutrino (full flux) no cosmics (v09_89_01_01) [PMT Decreased QE]
    */
    SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/mc_run2_new_weights/flat/*.flat.root", "spectra_mc_v09_84_01_01r3.root", -1, 2.68171e20);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_cv_v09_89_01_01r3.flat.root", "spectra_cv_v09_89_01_01r3.root", -1, -1);

    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_cv_v09_89_01_01r3.flat.root", "spectra_cv_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpcuntunedsigshape_v09_89_01_01.flat.root", "spectra_tpcuntunedsigshape_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpcind2opaque_v09_89_01_01.flat.root", "spectra_tpcind2opaque_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpcind2transparent_v09_89_01_01.flat.root", "spectra_tpcind2transparent_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpcind1increasegain_v09_89_01_01.flat.root", "spectra_tpcind1increasegain_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpcind1decreasegain_v09_89_01_01.flat.root", "spectra_tpcind1decreasegain_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_pmtdecreasedqe2_v09_89_01_01.flat.root", "spectra_pmtdecreasedqe_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_ellipsoidalrecomb_v09_89_01_01.flat.root", "spectra_ellipsoidalrecomb_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpccohnoisep1sigma_v09_89_01_01.flat.root", "spectra_tpccohnoisep1sigma_v09_89_01_01r3.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/neutrino2024/systematics/sample_tpcintnoisep1sigma_v09_89_01_01.flat.root", "spectra_tpcintnoisep1sigma_v09_89_01_01r3.root", -1, -1);

    spectra.add_spectrum1d("sInfoVar", Binning::Simple(1, 0, 2), kInfoVar);

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

    spectra.run();
}