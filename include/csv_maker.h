/**
 * @file csv_maker.h
 * @brief Header file defining a dummy SpillMultiVar for dumping particle info.
 * @author justin.mueller@colostate.edu
*/
#ifndef CSV_MAKER_H
#define CSV_MAKER_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "cuts.h"
#include "variables.h"
#include "numu_variables.h"

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

std::ofstream output("output_mc_crtpmt.log");
//std::ofstream output("output_tpcuntunedsigshape.log");

#define GUARD(VAL) std::isinf(VAL) ? -9999 : VAL
#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

/**
 * Writes information about a failed containment cut.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the truth interaction (signal)
 * @return None.
*/
void write_file_info(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
            << CSV(i.nu_id) << CSV(vars::image_id(i)) << CSV(vars::id(i))
            << CSV(std::string(sr->hdr.sourceName))
            << std::endl;
}

/**
 * Writes reconstructed variables (truth and reco) for selected/signal
 * interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the truth interaction (signal)
 * @param j the reco interaction (selected).
 * @return None.
*/
void write_pair(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i, const caf::SRInteractionDLPProxy& j)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
            //<< CSV(i.nu_energy_init + i.nu_position[2]) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
            << CSV(i.nu_id) << CSV(vars::image_id(i)) << CSV(vars::id(i))
            << CSV(sr->hdr.triggerinfo.global_trigger_det_time)
            << CSV(vars::category(i))
            << CSV(vars::category_topology(i))
            << CSV(vars::category_interaction_mode(i))
            << CSV(vars::leading_muon_ke(i))
            << CSV(vars::leading_muon_ke(j))
            << CSV(vars::leading_proton_ke(i))
            << CSV(vars::leading_proton_ke(j))
            << CSV(vars::visible_energy(i))
            << CSV(vars::visible_energy(j))
            << CSV(vars::leading_muon_pt(i))
            << CSV(vars::leading_muon_pt(j))
            << CSV(vars::leading_proton_pt(i))
            << CSV(vars::leading_proton_pt(j))
            << CSV(vars::muon_polar_angle(i))
            << CSV(vars::muon_polar_angle(j))
            << CSV(vars::muon_azimuthal_angle(i))
            << CSV(vars::muon_azimuthal_angle(j))
            << CSV(vars::opening_angle(i))
            << CSV(vars::opening_angle(j))
            << CSV(vars::interaction_pt(i))
            << CSV(vars::interaction_pt(j))
            << CSV(vars::phiT(i)) << CSV(vars::phiT(j))
            << CSV(vars::alphaT(i)) << CSV(vars::alphaT(j))
            << CSV(vars::muon_softmax(j)) << CSV(vars::proton_softmax(j))
            << CSV(cuts::all_1mu1p_cut(j))
            << CSV(cuts::all_1muNp_cut(j))
            << CSV(cuts::all_1muX_cut(j))
            << CSV(cuts::crtpmt_veto(sr))
            << CSV(j.volume_id)
            << std::endl;
}

const SpillMultiVar kInfoVar([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over truth interactions for efficiency metrics and for signal-level
     * variables of interest.
    */
    for(auto const & i : sr->dlp_true)
    {
        if(cuts::neutrino(i))
        {
            int category(vars::category(i));
            if(category % 2 == 0 && category < 5)
            {
                if(cuts::matched(i))
                {
                    OUT(output, "SIGNAL");
                    const auto & r = sr->dlp[i.match[0]];
                    write_pair(sr, i, r);

                    if(cuts::fiducial_cut(r) && !cuts::containment_cut(r))
                    {
                        OUT(output, "CONTAINMENT");
                        write_file_info(sr, i);
                    }
                }
            }
        }
    }

    /**
     * Loop over reconstructed interactions for purity metrics and for
     * reconstructed variables of interest.
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::all_1muX_cut(i) || cuts::all_1muNp_cut(i) || cuts::all_1mu1p_cut(i))
        {
            if(cuts::matched(i))
            {
                const auto & t = sr->dlp_true[i.match[0]];
                OUT(output, "SELECTED");
                write_pair(sr, t, i);
            }
        }
    }

    return std::vector<double>{1};
});

#endif