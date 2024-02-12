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
 * "Dummy" SpillMultiVar for writing weight information about each neutrino
 * into an output CSV log file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return a vector with a single dummy entry.
*/
const SpillMultiVar kReweight([](const caf::SRSpillProxy* sr)
{
    for(auto const & i : sr->dlp_true)
    {
        if(cuts::neutrino(i))
        {
            OUT(output,"REWEIGHT")  << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
                                    << CSV(sr->hdr.evt) << CSV(i.nu_id);
            for(size_t j(0); j < sr->mc.nu[i.nu_id].wgt[52].univ.size(); ++j)
                output << CSV(sr->mc.nu[i.nu_id].wgt[52].univ[j]);
            output << std::endl;
        }
    }
    return std::vector<double>{1};
});

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

const SpillMultiVar kCSVLogger([](const caf::SRSpillProxy* sr)
{
    /**
     * Build a map of reco_particles and true_particles to make reco->true and
     * true->reco matching of particles possible.
    */
    std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles;
    std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleTruthDLP> *> true_particles;
    for(auto const& i : sr->dlp)
    {
        for(auto const& p : i.particles)
            reco_particles.insert(std::make_pair(p.id, &p));
    }
    for(auto const& i : sr->dlp_true)
    {
        for(auto const & p : i.particles)
            true_particles.insert(std::make_pair(p.id, &p));
    }

    /**
     * Loop over true interactions and log the relevant interaction-level
     * information. 
    */
    for(auto const& i : sr->dlp_true)
    {
        /**
         * Loop over true particles within the interaction and log the relevant
         * particle-level information.
        */
        for(auto const & p : i.particles)
        {
            if(cuts::neutrino(i) && cuts::matched_muon(p) && cuts::muon(*reco_particles[p.match[0]]))
            {
                auto rp = reco_particles[p.match[0]];
                OUT(output,"PARTICLE")  << CSV(vars::image_id(i)) << CSV(vars::id(i)) << CSV("2")
                                        << CSV(vars::ke_init(p)) << CSV(vars::energy_deposit(p))
                                        << CSV(vars::csda_ke(*rp)) << std::endl;
            }
            else if(cuts::neutrino(i) && cuts::matched_proton(p) && cuts::proton(*reco_particles[p.match[0]]))
            {
                auto rp = reco_particles[p.match[0]];
                OUT(output,"PARTICLE")  << CSV(vars::image_id(i)) << CSV(vars::id(i)) << CSV("4")
                                        << CSV(vars::ke_init(p)) << CSV(vars::energy_deposit(p))
                                        << CSV(vars::csda_ke(*rp)) << std::endl;
            }
        }

        /**
         * Log information about neutrino-genic true interactions.
        */
        if(cuts::neutrino(i))
            OUT(output,"NEUTRINO")  << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                    << CSV(i.nu_interaction_mode) << CSV(i.truth_vertex[0])
                                    << CSV(i.truth_vertex[1]) << CSV(i.truth_vertex[2])
                                    << CSV(i.nu_id) << std::endl;
        
        /**
         * Log information about signal interactions.
        */
        if(cuts::signal_1mu1p(i) && cuts::fiducial_containment_cut(i))
        {
            bool selected_fv(cuts::matched(i) && cuts::fiducial_cut(sr->dlp[i.match[0]]));
            bool selected_cn(cuts::matched(i) && cuts::containment_cut(sr->dlp[i.match[0]]));
            OUT(output,"SIG1MU1P")  << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                    << CSV(vars::leading_muon_ke(i))
                                    << CSV(vars::leading_proton_ke(i))
                                    << CSV(selected_fv) << CSV(selected_cn)
                                    << std::endl;

            /**
             * Log information about particles belonging to the signal
             * interaction.
            */
            for(auto const & p : i.particles)
                OUT(output,"SIGPART")   << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                        << CSV(p.pid) << CSV(p.is_primary)
                                        << CSV(p.is_contained) << CSV(p.size)
                                        << CSV(p.energy_init) << CSV(p.start_point[0])
                                        << CSV(p.start_point[1]) << CSV(p.start_point[2])
                                        << CSV(p.end_point[0]) << CSV(p.end_point[1])
                                        << CSV(p.end_point[2]) << std::endl;
            
            if(cuts::matched(i))
            {
                /**
                 * Log information about interactions which have been matched
                 * to by true signal interactions.
                */
                auto & mi = sr->dlp[i.match[0]];
                OUT(output,"SIGRECO")   << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                        << CSV(vars::image_id(mi)) << CSV(vars::id(mi))
                                        << CSV(mi.is_fiducial) << CSV(mi.is_contained)
                                        << std::endl;

                /**
                 * Log information about particles in the interactions that have
                 * been matched to by true signal interactions.
                */
                for(auto const & p : mi.particles)
                {
                    float t0(-9999);
                    if(cuts::matched(p))
                        t0 = true_particles[p.match[0]]->t;
                    OUT(output,"SIGPART")   << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                            << CSV(p.pid) << CSV(p.is_primary)
                                            << CSV(p.is_contained) << CSV(p.size)
                                            << CSV(p.csda_ke) << CSV(p.start_point[0])
                                            << CSV(p.start_point[0]) << CSV(p.start_point[1])
                                            << CSV(p.start_point[2]) << CSV(p.end_point[1])
                                            << CSV(p.end_point[2]) << CSV(t0) << std::endl;
                }
            }
        }
    }

    /**
     * Loop over reco interactions and log the relevant interaction-level
     * information. 
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::fiducial_containment_cut(i))
        {
            /**
             * 
            */
            std::vector<uint32_t> counts = cuts::count_primaries(i);
            double truth_category(-1);
            int64_t true_interaction_id(-1), nu_id(-1);
            if(cuts::matched(i))
            {
                truth_category = vars::category(sr->dlp_true[i.match[0]]);
                true_interaction_id = sr->dlp_true[i.match[0]].id;
                nu_id = sr->dlp_true[i.match[0]].nu_id;
            }
            OUT(output,"INTIME")    << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                    << CSV((counts[0])) << CSV((counts[1]))
                                    << CSV((counts[2])) << CSV((counts[3]))
                                    << CSV((counts[4])) << CSV(cuts::valid_flashmatch(i))
                                    << CSV(cuts::flash_cut(i)) << CSV(truth_category)
                                    << CSV(true_interaction_id) << CSV(nu_id) << std::endl;
            for(auto const & p : i.particles)
            {
                float t0(-10000);
                int64_t parent_pdg_code(0);
                int64_t interaction_id(-1);
                if(cuts::matched(p))
                {
                    t0 = true_particles[p.match[0]]->t;
                    parent_pdg_code = true_particles[p.match[0]]->parent_pdg_code;
                    interaction_id = true_particles[p.match[0]]->interaction_id;
                }
                OUT(output,"PINT")  << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                    << CSV(p.pid) << CSV(p.is_primary)
                                    << CSV(p.is_contained) << CSV(p.size)
                                    << CSV(p.csda_ke) << CSV(p.start_point[0])
                                    << CSV(p.start_point[1]) << CSV(p.start_point[2])
                                    << CSV(p.end_point[0]) << CSV(p.end_point[1])
                                    << CSV(p.end_point[2]) << CSV(t0)
                                    << CSV(parent_pdg_code) << CSV(interaction_id)
                                    << std::endl;
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