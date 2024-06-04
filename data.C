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

std::ofstream output("output_data_crtpmt.log");
std::ofstream output_evt("output_evt.log");

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
            << CSV(cuts::crtpmt_veto_data(sr))
            << CSV(j.volume_id)
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

    output_evt  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun) << std::endl;

    return std::vector<double>{1};
});

/**
 * Enumerates the cut that each interaction passes.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A vector of doubles containing the enumeration of cuts passed
 * by each interaction.
*/
const SpillMultiVar kOffbeam1mu1pCut([](const caf::SRSpillProxy* sr)
{
    std::vector<double> cut_vector;
    for(auto const & i : sr->dlp)
    {
        // No cut: 0, fiducial: 1, contained: 2, topological: 3, flash: 4
        double cut(0);
        if(cuts::fiducial_cut(i))
            cut = 1;
        if(cuts::fiducial_containment_cut(i))
            cut = 2;
        if(cuts::fiducial_containment_topological_1mu1p_cut(i))
            cut = 3;
        if(cuts::all_1mu1p_data_cut(i))
            cut = 4;
        cut_vector.push_back(cut);
    }
    return cut_vector;
});

/**
 * Enumerates the cut that each interaction passes.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A vector of doubles containing the enumeration of cuts passed
 * by each interaction.
*/
const SpillMultiVar kOffbeam1muNpCut([](const caf::SRSpillProxy* sr)
{
    std::vector<double> cut_vector;
    for(auto const & i : sr->dlp)
    {
        // No cut: 0, fiducial: 1, contained: 2, topological: 3, flash: 4
        double cut(0);
        if(cuts::fiducial_cut(i))
            cut = 1;
        if(cuts::fiducial_containment_cut(i))
            cut = 2;
        if(cuts::fiducial_containment_topological_1muNp_cut(i))
            cut = 3;
        if(cuts::all_1muNp_data_cut(i))
            cut = 4;
        cut_vector.push_back(cut);
    }
    return cut_vector;
});

/**
 * Enumerates the cut that each interaction passes.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A vector of doubles containing the enumeration of cuts passed
 * by each interaction.
*/
const SpillMultiVar kOffbeam1muXCut([](const caf::SRSpillProxy* sr)
{
    std::vector<double> cut_vector;
    for(auto const & i : sr->dlp)
    {
        // No cut: 0, fiducial: 1, contained: 2, topological: 3, flash: 4
        double cut(0);
        if(cuts::fiducial_cut(i))
            cut = 1;
        if(cuts::fiducial_containment_cut(i))
            cut = 2;
        if(cuts::fiducial_containment_topological_1muX_cut(i))
            cut = 3;
        if(cuts::all_1muX_data_cut(i))
            cut = 4;
        cut_vector.push_back(cut);
    }
    return cut_vector;
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

    //SpecContainer spectra("/pnfs/icarus/persistent/users/mueller/run9435_new_weights/*.flat.root", "spectra_run9435.root", -1, -1);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2_actual_new_weights/offbeam/hdf5/*.flat.root", "spectra_data_offbeam.root", -1, 266267);
    SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2_actual_new_weights/onbeam/flat/*.flat.root", "spectra_data_new_weights.root", -1, -1);

    spectra.add_spectrum1d("sDataInfo", Binning::Simple(1, 0, 2), kDataInfo);
    spectra.add_spectrum1d("sOffbeam1mu1pCut", Binning::Simple(5, 0, 5), kOffbeam1mu1pCut);
    spectra.add_spectrum1d("sOffbeam1muNpCut", Binning::Simple(5, 0, 5), kOffbeam1muNpCut);
    spectra.add_spectrum1d("sOffbeam1muXCut", Binning::Simple(5, 0, 5), kOffbeam1muXCut);
    //spectra.add_spectrum1d("sHandscanInfo", Binning::Simple(1, 0, 2), kHandscanInfo);

    spectra.run();
}
