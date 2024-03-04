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

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

std::ofstream output("output.log");

#define GUARD(VAL) std::isinf(VAL) ? -9999 : VAL
#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

/**
 * Writes information about the event to an output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return None.
*/
void write_event(const caf::SRSpillProxy* sr)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
            << CSV(sr->hdr.evt) << std::endl;
}

/**
 * Writes information about signal interactions to an output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the signal interaction.
 * @return None.
*/
void write_signal(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
            << CSV(sr->hdr.evt) << CSV(i.nu_id) << std::endl;
}

/**
 * Writes information about selected interactions to an output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the truth interaction matched to by the selected interaction.
 * @param j the selected interaction.
 * @return None.
*/
void write_selected(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i, const caf::SRInteractionDLPProxy& j)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
            << CSV(sr->hdr.evt) << CSV(i.nu_id)
            << CSV(vars::image_id(i)) << CSV(vars::id(i))
            << CSV(vars::category(i))
            << CSV(vars::category_topology(i))
            << CSV(vars::category_interaction_mode(i))
            << CSV(vars::visible_energy(j))
            << CSV(vars::leading_muon_ke(j))
            << CSV(vars::leading_proton_ke(j))
            << CSV(vars::leading_muon_pt(j))
            << CSV(vars::leading_proton_pt(j))
            << CSV(vars::interaction_pt(j))
            << CSV(vars::leading_muon_cosine_theta_xz(j))
            << CSV(vars::leading_proton_cosine_theta_xz(j))
            << CSV(vars::cosine_opening_angle(j))
            << CSV(vars::cosine_opening_angle_transverse(j))
            << CSV(vars::leading_muon_softmax(j))
            << CSV(vars::leading_proton_softmax(j))
            << std::endl;
}

/**
 * Writes empty output for each selected interaction with no truth match.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return None.
*/
void write_none(const caf::SRSpillProxy* sr)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
            << CSV(sr->hdr.evt) << CSV(-1) << CSV(-1)
            << CSV(-1) << CSV(-1) << CSV(-1)
            << CSV(-1) << CSV(-1) << CSV(-1)
            << CSV(-1) << CSV(-1) << CSV(-1)
            << CSV(-1) << CSV(-1) << CSV(-1)
            << CSV(-1) << CSV(-1) << std::endl;
}

/**
 * Writes information about selected non-signal interactions to an output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the selected interaction.
 * @return None.
*/
void write_mistake(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i)
{
    const auto & p = i.particle_counts;
    const auto & pr = cuts::count_primaries(i);
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
            << CSV(sr->hdr.evt) << CSV(i.nu_id)
            << CSV(vars::image_id(i)) << CSV(vars::id(i))
            << CSV(p[0]) << CSV(p[1]) << CSV(p[2])
            << CSV(p[3]) << CSV(p[4]) << CSV(pr[0])
            << CSV(pr[1]) << CSV(pr[2]) << CSV(pr[3])
            << CSV(pr[4]) << std::endl;
}

/**
 * "Dummy" SpillMultiVar for  writing information about each signal (using
 * truth information) into a CSV log file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return a vector with a single dummy entry.
*/
const SpillMultiVar kSignal([](const caf::SRSpillProxy* sr)
{
    for(auto const & i : sr->dlp_true)
    {
        /**
         * General neutrinos
        */
        if(cuts::neutrino(i))
        {
            OUT(output, "NEUTRINO");
            write_signal(sr, i);
        }
        /**
         * Signal: 1mu1p
        */
        if(cuts::signal_1mu1p(i) && cuts::fiducial_containment_cut(i))
        {
            OUT(output, "SIGNAL_1MU1P");
            write_signal(sr, i);
        }
        /**
         * Signal: 1muNp
        */
        if(cuts::signal_1muNp(i) && cuts::fiducial_containment_cut(i))
        {
            OUT(output, "SIGNAL_1MUNP");
            write_signal(sr, i);
        }
        /**
         * Signal: 1muX
        */
        if(cuts::signal_1muX(i) && cuts::fiducial_containment_cut(i))
        {
            OUT(output, "SIGNAL_1MUX");
            write_signal(sr, i);
        }
        /**
         * Mistake: 1mu1p
        */
        if(cuts::matched(i) && !cuts::signal_1mu1p(i) && cuts::all_1mu1p_cut(sr->dlp[i.match[0]]))
        {
            OUT(output, "MISTAKE_1MU1P");
            write_mistake(sr, i);
        }
        /**
         * Mistake: 1muNp
        */
        if(cuts::matched(i) && !cuts::signal_1muNp(i) && cuts::all_1muNp_cut(sr->dlp[i.match[0]]))
        {
            OUT(output, "MISTAKE_1MUNP");
            write_mistake(sr, i);
        }
        /**
         * Mistake: 1muX
        */
        if(cuts::matched(i) && !cuts::signal_1muX(i) && cuts::all_1muX_cut(sr->dlp[i.match[0]]))
        {
            OUT(output, "MISTAKE_1MUX");
            write_mistake(sr, i);
        }
    }
    return std::vector<double>{1};
});

/**
 * "Dummy" SpillMultiVar for writing information about each selected signal
 * candidate into a CSV log file. 
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return a vector with a single dummy entry.
*/
const SpillMultiVar kSelected([](const caf::SRSpillProxy* sr)
{
    OUT(output,"EVENT");
    write_event(sr);
    for(auto const & i : sr->dlp)
    {
        /**
         * Selected: 1mu1p
        */
        if(cuts::all_1mu1p_cut(i))
        {
            if(cuts::matched(i))
            {
                const auto & t = sr->dlp_true[i.match[0]];
                OUT(output,"SELECTED_1MU1P");
                write_selected(sr, t, i);
            }
            else
            {
                OUT(output,"SELECTED_NONE");
                write_none(sr);
            }
        }
        /**
         * Selected: 1muNp
        */
        if(cuts::all_1muNp_cut(i))
        {
            if(cuts::matched(i))
            {
                const auto & t = sr->dlp_true[i.match[0]];
                OUT(output,"SELECTED_1MUNP");
                write_selected(sr, t, i);
            }
            else
            {
                OUT(output,"SELECTED_NONE");
                write_none(sr);
            }
        }
        /**
         * Selected: 1muX
        */
        if(cuts::all_1muX_cut(i))
        {
            if(cuts::matched(i))
            {
                const auto & t = sr->dlp_true[i.match[0]];
                OUT(output,"SELECTED_1MUX");
                write_selected(sr, t, i);
            }
            else
            {
                OUT(output,"SELECTED_NONE");
                write_none(sr);
            }
        }
    }
    return std::vector<double>{1};
});

const SpillMultiVar kDataLogger([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over reconstructed interactions and log interaction-level
     * information. No truth information can be used.
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::topological_1muNp_cut(i))
        {
            size_t leading_muon(0), leading_proton(0), index(0);
            double leading_muon_ke(0), leading_proton_ke(0);
            for(auto & p : i.particles)
            {
                if(p.pid == 2 && p.csda_ke > leading_muon_ke)
                {
                    leading_muon = index;
                    leading_muon_ke = p.csda_ke;
                }
                else if(p.pid == 4 && p.csda_ke > leading_proton_ke)
                {
                    leading_proton = index;
                    leading_proton_ke = p.csda_ke;
                }
                ++index;
            }
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

#endif