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

#include "cuts.h"
#include "variables.h"

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

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
                std::cout << "2," << vars::ke_init(p) << "," << vars::csda_ke(*rp) << std::endl;
            }
            else if(cuts::neutrino(i) && cuts::matched_proton(p) && cuts::proton(*reco_particles[p.match[0]]))
            {
                auto rp = reco_particles[p.match[0]];
                std::cout << "4," << vars::ke_init(p) << "," << vars::csda_ke(*rp) << std::endl;
            }
        }
    }
    return var;
});

#endif