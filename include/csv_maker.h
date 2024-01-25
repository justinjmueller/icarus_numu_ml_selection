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

const SpillMultiVar kParticleCSV([](const caf::SRSpillProxy* sr)
{
    std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles;
    for(auto const& i : sr->dlp)
    {
        for(auto const& p : i.particles)
            reco_particles.insert(std::make_pair(p.id, &p));
    }

    std::vector<double> var{1};
    for(auto const& i : sr->dlp_true)
    {
        for(auto const & p : i.particles)
        {

            if(cuts::neutrino(i) && cuts::matched_muon(p) && cuts::muon(*reco_particles[p.match[0]]))
            {
                auto rp = reco_particles[p.match[0]];
                output << std::fixed << "PARTICLE," << vars::image_id(i) << "," << vars::id(i) << "," << "2," << vars::ke_init(p) << "," << vars::energy_deposit(p) << "," << vars::csda_ke(*rp) << std::endl;
            }
            else if(cuts::neutrino(i) && cuts::matched_proton(p) && cuts::proton(*reco_particles[p.match[0]]))
            {
                auto rp = reco_particles[p.match[0]];
                output << std::fixed << "PARTICLE," << vars::image_id(i) << "," << vars::id(i) << "," << "4," << vars::ke_init(p) << "," << vars::energy_deposit(p) << "," << vars::csda_ke(*rp) << std::endl;
            }
        }
        if(cuts::signal_1mu1p(i) && cuts::fiducial_containment_cut(i))
        {
            bool selected_fv(cuts::matched(i) && cuts::fiducial_cut(sr->dlp[i.match[0]]));
            bool selected_cn(cuts::matched(i) && cuts::containment_cut(sr->dlp[i.match[0]]));
            output << std::fixed << "SIG1MU1P," << vars::image_id(i) << "," << vars::id(i) << "," << vars::leading_muon_ke(i) << "," << vars::leading_proton_ke(i) << "," << selected_fv << "," << selected_cn << std::endl;
            for(auto const & p : i.particles)
                output << std::fixed << "SIGPART," << vars::image_id(i) << "," << vars::id(i) << "," << p.pid << "," << p.is_primary << "," << p.is_contained << "," << p.num_voxels << "," << p.energy_init << "," << p.start_point[0] << "," << p.end_point[0] << std::endl;
        }
    }
    return var;
});

#endif