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
 * "Dummy" SpillMultiVar for writing information about each selected 1mu1p
 * candidate into a CSV log file. 
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return a vector with a single dummy entry.
*/
const SpillMultiVar kSelected1mu1p([](const caf::SRSpillProxy* sr)
{
    for(auto const & i : sr->dlp)
    {
        if(cuts::all_1mu1p_cut(i))
        {
            if(cuts::matched(i))
            {
                const auto & t = sr->dlp_true[i.match[0]];
                OUT(output,"SELECTED_1MU1P")    << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
                                                << CSV(sr->hdr.evt) << CSV(t.nu_id)
                                                << CSV(vars::category(t))
                                                << CSV(vars::category_topology(t))
                                                << CSV(vars::category_interaction_mode(t))
                                                << CSV(vars::visible_energy(i))
                                                << CSV(vars::leading_muon_ke(i))
                                                << CSV(vars::leading_proton_ke(i))
                                                << CSV(vars::leading_muon_pt(i))
                                                << CSV(vars::leading_proton_pt(i))
                                                << CSV(vars::interaction_pt(i))
                                                << std::endl;
            }
            else
            {
                OUT(output,"SELECTED_NONE") << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
                                            << CSV(sr->hdr.evt) << CSV(-1) << CSV(-1)
                                            << CSV(-1) << CSV(-1) << CSV(-1)
                                            << CSV(-1) << CSV(-1) << CSV(-1)
                                            << CSV(-1) << CSV(-1) << std::endl;
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