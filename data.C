/**
 * @file data.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection on
 * data events.
 * @author justin.mueller@colostate.edu
*/

//#include "include/analysis.h"
//#include "include/datalog.h"
#include "include/definitions.h"
#include "include/cuts.h"
#include "include/variables.h"
#include "include/numu_variables.h"
#include "include/container.h"
#include "sbnana/CAFAna/Core/Binning.h"

#include <fstream>

using namespace ana;

#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

std::ofstream output("output.log");

/**
 * Writes reconstructed variables selected interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param j the reco interaction (selected).
 * @return None.
*/
void write_reco(const caf::SRSpillProxy* sr, const caf::SRInteractionDLPProxy& j)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
            << CSV(vars::image_id(j)) << CSV(vars::id(j))
            << CSV(vars::leading_muon_ke(j))
            << CSV(vars::leading_proton_ke(j))
            << CSV(vars::visible_energy(j))
            << CSV(vars::leading_muon_pt(j))
            << CSV(vars::leading_proton_pt(j))
            << CSV(vars::muon_polar_angle(j))
            << CSV(vars::muon_azimuthal_angle(j))
            << CSV(vars::opening_angle(j))
            << CSV(vars::interaction_pt(j))
            << CSV(vars::phiT(j))
            << CSV(vars::alphaT(j))
            << CSV(vars::muon_softmax(j))
            << CSV(vars::proton_softmax(j))
            << CSV(cuts::all_1mu1p_data_cut(j))
            << CSV(cuts::all_1muNp_data_cut(j))
            << CSV(cuts::all_1muX_data_cut(j))
            << std::endl;
}

/**
 * Writes the reconstructed variables for the selected interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A dummy vector of doubles.
*/
const SpillMultiVar kDataInfo([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over reconstructed interactions and log interaction-level
     * information. 
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::all_1muX_data_cut(i) || cuts::all_1muNp_data_cut(i) || cuts::all_1mu1p_data_cut(i))
        {
            OUT(output,"DATA");
            write_reco(sr, i);
        }
    }
    return std::vector<double>{1};
});

const SpillMultiVar kHandscanInfo([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over reconstructed interactions and log interaction-level
     * information. No truth information can be used.
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::topological_1muNp_cut(i))
        {
            size_t leading_muon(vars::leading_particle_index(i, 2));
            size_t leading_proton(vars::leading_particle_index(i, 4));
            OUT(output,"INTERACTION")   << CSV(sr->hdr.run) << CSV(sr->hdr.evt)
                                        << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                        << CSV(vars::cryostat(i)) << CSV(i.is_fiducial)
                                        << CSV(i.is_contained) << CSV(cuts::topology(i))
                                        << CSV(cuts::flash_cut_data(i))
                                        << CSV(i.vertex[0]) << CSV(i.vertex[1]) << CSV(i.vertex[2])
                                        << CSV(i.particles[leading_muon].length)
                                        << CSV(vars::leading_muon_ke(i))
                                        << CSV(i.particles[leading_proton].length)
                                        << CSV(vars::leading_proton_ke(i))
                                        << CSV(vars::flash_time(i))
                                        << CSV(i.particles[leading_muon].end_point[0])
                                        << CSV(i.particles[leading_muon].end_point[1])
                                        << CSV(i.particles[leading_muon].end_point[2])
                                        << CSV(i.particles[leading_muon].start_dir[0])
                                        << CSV(i.particles[leading_muon].start_dir[1])
                                        << CSV(i.particles[leading_muon].start_dir[2])
                                        << CSV(i.particles[leading_proton].end_point[0])
                                        << CSV(i.particles[leading_proton].end_point[1])
                                        << CSV(i.particles[leading_proton].end_point[2])
                                        << CSV(i.particles[leading_proton].start_dir[0])
                                        << CSV(i.particles[leading_proton].start_dir[1])
                                        << CSV(i.particles[leading_proton].start_dir[2])
                                        << std::endl;
        }
    }

    return std::vector<double>{1};
});

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/
void data()
{
    RECO_SIGNAL_VAR(kVisibleEnergy, vars::visible_energy);
    RECO_SIGNAL_VAR(kLeadingMuonKE, vars::leading_muon_ke);
    RECO_SIGNAL_VAR(kLeadingProtonKE, vars::leading_proton_ke);
    RECO_SIGNAL_VAR(kLeadingMuonPT, vars::leading_muon_pt);
    RECO_SIGNAL_VAR(kLeadingProtonPT, vars::leading_proton_pt);
    RECO_SIGNAL_VAR(kInteractionPT, vars::interaction_pt);
    RECO_SIGNAL_VAR(kLeadingMuonCosineThetaXZ, vars::leading_muon_cosine_theta_xz);
    RECO_SIGNAL_VAR(kLeadingProtonCosineThetaXZ, vars::leading_proton_cosine_theta_xz);
    RECO_SIGNAL_VAR(kCosineOpeningAngle, vars::cosine_opening_angle);
    RECO_SIGNAL_VAR(kCosineOpeningAngleTransverse, vars::cosine_opening_angle_transverse);
    RECO_SIGNAL_VAR(kLeadingMuonSoftmax, vars::leading_muon_softmax);
    RECO_SIGNAL_VAR(kLeadingProtonSoftmax, vars::leading_proton_softmax);

    VARDLP_RECO(kFlashTime, vars::flash_time, cuts::fiducial_containment_topological_1muNp_cut);

    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/run9435/*.flat.root", "spectra_run9435.root", -1, -1);
    SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2/flat/*.flat.root", "spectra_data.root", -1, -1);

    spectra.add_spectrum1d("sFlashTime", Binning::Simple(100, -10, 10), kFlashTime);

    spectra.add_spectrum1d("sVisibleEnergy_1mu1p", Binning::Simple(25, 0, 3000), kVisibleEnergy_1mu1p);
    spectra.add_spectrum1d("sVisibleEnergy_1muNp", Binning::Simple(25, 0, 3000), kVisibleEnergy_1muNp);
    spectra.add_spectrum1d("sVisibleEnergy_1muX", Binning::Simple(25, 0, 3000), kVisibleEnergy_1muX);
    spectra.add_spectrum1d("sLeadingMuonKE_1mu1p", Binning::Simple(25, 0, 2000), kLeadingMuonKE_1mu1p);
    spectra.add_spectrum1d("sLeadingMuonKE_1muNp", Binning::Simple(25, 0, 2000), kLeadingMuonKE_1muNp);
    spectra.add_spectrum1d("sLeadingMuonKE_1muX", Binning::Simple(25, 0, 2000), kLeadingMuonKE_1muX);
    spectra.add_spectrum1d("sLeadingProtonKE_1mu1p", Binning::Simple(25, 0, 600), kLeadingProtonKE_1mu1p);
    spectra.add_spectrum1d("sLeadingProtonKE_1muNp", Binning::Simple(25, 0, 600), kLeadingProtonKE_1muNp);
    //spectra.add_spectrum1d("sLeadingProtonKE_1muX", Binning::Simple(25, 0, 600), kLeadingProtonKE_1muX);
    spectra.add_spectrum1d("sLeadingMuonPT_1mu1p", Binning::Simple(25, 0, 750), kLeadingMuonPT_1mu1p);
    spectra.add_spectrum1d("sLeadingMuonPT_1muNp", Binning::Simple(25, 0, 750), kLeadingMuonPT_1muNp);
    spectra.add_spectrum1d("sLeadingMuonPT_1muX", Binning::Simple(25, 0, 750), kLeadingMuonPT_1muX);
    spectra.add_spectrum1d("sLeadingProtonPT_1mu1p", Binning::Simple(25, 0, 750), kLeadingProtonPT_1mu1p);
    spectra.add_spectrum1d("sLeadingProtonPT_1muNp", Binning::Simple(25, 0, 750), kLeadingProtonPT_1muNp);
    //spectra.add_spectrum1d("sLeadingProtonPT_1muX", Binning::Simple(25, 0, 750), kLeadingProtonPT_1muX);
    spectra.add_spectrum1d("sInteractionPT_1mu1p", Binning::Simple(25, 0, 1200), kInteractionPT_1mu1p);
    spectra.add_spectrum1d("sInteractionPT_1muNp", Binning::Simple(25, 0, 1200), kInteractionPT_1muNp);
    spectra.add_spectrum1d("sInteractionPT_1muX", Binning::Simple(25, 0, 1200), kInteractionPT_1muX);
    spectra.add_spectrum1d("sLeadingMuonCosineThetaXZ_1mu1p", Binning::Simple(25, -1, 1), kLeadingMuonCosineThetaXZ_1mu1p);
    spectra.add_spectrum1d("sLeadingMuonCosineThetaXZ_1muNp", Binning::Simple(25, -1, 1), kLeadingMuonCosineThetaXZ_1muNp);
    spectra.add_spectrum1d("sLeadingMuonCosineThetaXZ_1muX", Binning::Simple(25, -1, 1), kLeadingMuonCosineThetaXZ_1muX);
    spectra.add_spectrum1d("sLeadingProtonCosineThetaXZ_1mu1p", Binning::Simple(25, -1, 1), kLeadingProtonCosineThetaXZ_1mu1p);
    spectra.add_spectrum1d("sLeadingProtonCosineThetaXZ_1muNp", Binning::Simple(25, -1, 1), kLeadingProtonCosineThetaXZ_1muNp);
    //spectra.add_spectrum1d("sLeadingProtonCosineThetaXZ_1muX", Binning::Simple(25, -1, 1), kLeadingProtonCosineThetaXZ_1muX);
    spectra.add_spectrum1d("sCosineOpeningAngle_1mu1p", Binning::Simple(25, -1, 1), kCosineOpeningAngle_1mu1p);
    spectra.add_spectrum1d("sCosineOpeningAngle_1muNp", Binning::Simple(25, -1, 1), kCosineOpeningAngle_1muNp);
    //spectra.add_spectrum1d("sCosineOpeningAngle_1muX", Binning::Simple(25, -1, 1), kCosineOpeningAngle_1muX);
    spectra.add_spectrum1d("sCosineOpeningAngleTransverse_1mu1p", Binning::Simple(25, -1, 1), kCosineOpeningAngleTransverse_1mu1p);
    spectra.add_spectrum1d("sCosineOpeningAngleTransverse_1muNp", Binning::Simple(25, -1, 1), kCosineOpeningAngleTransverse_1muNp);
    //spectra.add_spectrum1d("sCosineOpeningAngleTransverse_1muX", Binning::Simple(25, -1, 1), kCosineOpeningAngleTransverse_1muX);
    spectra.add_spectrum1d("sLeadingMuonSoftmax_1mu1p", Binning::Simple(25, 0, 1), kLeadingMuonSoftmax_1mu1p);
    spectra.add_spectrum1d("sLeadingMuonSoftmax_1muNp", Binning::Simple(25, 0, 1), kLeadingMuonSoftmax_1muNp);
    spectra.add_spectrum1d("sLeadingMuonSoftmax_1muX", Binning::Simple(25, 0, 1), kLeadingMuonSoftmax_1muX);
    spectra.add_spectrum1d("sLeadingProtonSoftmax_1mu1p", Binning::Simple(25, 0.8, 1), kLeadingProtonSoftmax_1mu1p);
    spectra.add_spectrum1d("sLeadingProtonSoftmax_1muNp", Binning::Simple(25, 0.8, 1), kLeadingProtonSoftmax_1muNp);
    //spectra.add_spectrum1d("sLeadingProtonSoftmax_1muX", Binning::Simple(25, 0.8, 1), kLeadingProtonSoftmax_1muX);

    spectra.add_spectrum1d("sDataInfo", Binning::Simple(1, 0, 2), kDataInfo);

    spectra.run();
}