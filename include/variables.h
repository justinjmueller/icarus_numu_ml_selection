/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include <algorithm>

namespace vars
{

    /**
     * Variable for counting interactions/particles.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return 1.0 (always).
     */
    template<class T>
        double count(const T & obj) { return 1.0; }

    /**
     * Variable for image_id (unique identifier for the event).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the image_id of the interaction/particle.
    */
    template<class T>
        double image_id(const T & obj) { return obj.image_id; }

    /**
     * Variable for id (unique identifier for the object).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the id of the interaction/particle.
    */
    template<class T>
        double id(const T & obj) { return obj.id; }

    /**
     * Variable for enumerating interaction categories. This is a basic
     * categorization using only signal, neutrino background, and cosmic
     * background as the three categories.
     * 0: 1mu1p (contained and fiducial), 1: 1mu1p (not contained or fiducial)
     * 2: 1muNp (N > 1, contained and fiducial), 3: 1muNp (N > 1, not contained
     * or fiducial), 4: Other nu, 5: cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category(const T & interaction)
        {
            double cat(5);
            if(cuts::signal_1mu1p(interaction))
            {
                if(cuts::fiducial_containment_cut(interaction)) cat = 0;
                else cat = 1;
            }
            else if(cuts::signal_1muNp(interaction))
            {
                if(cuts::fiducial_containment_cut(interaction)) cat = 2;
                else cat = 3;
            }
            else if(cuts::other_nu_1muNp(interaction)) cat = 4;
            return cat;
        }

    /**
     * Variable for enumerating interaction categories. This classifies the
     * interactions based on the visible final states.
     * 0: 1mu1p, 1: 1mu0h, 2: 1muNp (N>1), 3: 1mu1p1pi, 4: nu_mu CC Other, 5: NC, 6: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_topology(const T & interaction)
        {
            uint16_t cat(6);
            if(interaction.is_neutrino)
            {
                uint32_t counts[5] = {};
                for(auto &p : interaction.particles)
                {
                    if(p.is_primary)
                    {
                        double energy(p.csda_ke);
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                            energy = p.energy_deposit;
                        if((p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 50))
                            ++counts[p.pid];
                    }
                }
                if(counts[0] == 0 && counts[1] == 0 && counts[2] == 1)
                {
                    if(counts[3] == 0 && counts[4] == 1) cat = 0;
                    else if(counts[3] == 0 && counts[4] == 0) cat = 1;
                    else if(counts[3] == 0 && counts[4] > 1) cat = 2;
                    else if(counts[3] == 1 && counts[4] == 1) cat = 3;
                    else if(interaction.nu_current_type == 0) cat = 4;
                }
                else if(interaction.nu_current_type == 0) cat = 4;
                else if(interaction.nu_current_type == 1) cat = 5;
            }
            return cat;
        }

    /**
     * Variable for enumerating interaction categories. This categorization
     * uses the interaction type (generator truth) classify the interactions
     * 0: nu_mu CC QE, 1: nu_mu CC Res, 2: nu_mu CC MEC, 3: nu_mu CC DIS, 4: nu_mu CC Coh, 5: nu_e CC, 6: NC, 7: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_interaction_mode(const T & interaction)
        {
            double cat(7);

            if(interaction.is_neutrino)
            {
                if(interaction.nu_current_type == 0)
                {
                    if(abs(interaction.nu_pdg_code) == 14)
                    {
                        if(interaction.nu_interaction_mode == 0) cat = 0;
                        else if(interaction.nu_interaction_mode == 1) cat = 1;
                        else if(interaction.nu_interaction_mode == 10) cat = 2;
                        else if(interaction.nu_interaction_mode == 2) cat = 3;
                        else if(interaction.nu_interaction_mode == 3) cat = 4;
                        else cat = 8;
                    }
                    else cat = 5;
                }
                else cat = 6;
            }

            return cat;
        }

    /**
     * Variable for counting particles in interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the number of particles in the interaction.
     */
    template<class T>
        double count_particles(const T & interaction) { return interaction.num_particles; }

    /**
     * Variable for counting primaries in interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the number of primaries in the interaction.
     */
    template<class T>
        double count_primaries(const T & interaction) { return interaction.num_primaries; }
    
    /**
     * Variable for total visible energy of interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the total visible energy of the interaction.
    */
    template<class T>
        double visible_energy(const T & interaction)
        {
            double energy(0);
            for(const auto & p : interaction.particles)
            {
                if(p.is_primary)
                {
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        energy += p.energy_deposit;
                    }
                    else
                    {
                        if(p.pid < 2) energy += p.calo_ke;
                        else energy += p.csda_ke;
                    }
                    if(p.pid == 2) energy += MUON_MASS;
                    else if(p.pid == 3) energy += PION_MASS;
                }
            }
            return energy;
        }
    /**
     * Variable for energy of the neutrino primary of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the neutrino energy.
    */
    template<class T>
        double neutrino_energy(const T & interaction) { return 1000*interaction.nu_energy_init; }

    /**
     * Variable for matched interaction flash time.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the matched flash time of the interaction.
    */
    template<class T>
        double flash_time(const T & interaction)
        {
            if(!cuts::valid_flashmatch(interaction))
                return -100000.0;
            else
                return interaction.flash_time;
        }

    /**
     * Variable for particle primary categorizations.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the primary/non-primary designation of the particle.
    */
    template<class T>
        double primary(const T & particle) { return particle.is_primary ? 1 : 0; }

    /**
     * Variable for particle PID.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID of the particle.
    */
    template<class T>
        double pid(const T & particle) { return particle.pid; }

    /**
     * Variable for particle PID + primary.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID+primary information for the particle.
    */
    template<class T>
        double primary_pid(const T & particle) { return particle.pid + (particle.is_primary ? 5 : 0); }

    /**
     * Variable for particle csda_ke.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle.
    */
    template<class T>
        double csda_ke(const T & particle) { return particle.csda_ke; }

    /**
     * Variable for particle csda_ke (muons only).
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle (if a muon).
    */
    template<class T>
        double csda_ke_muon(const T & particle) { return (cuts::muon(particle)) ? csda_ke(particle) : -1; }

    /**
     * Variable for true particle energy deposited.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the energy deposited by the particle (truth only).
    */
    template<class T>
        double energy_deposit(const T & particle) { return particle.energy_deposit; }

    /**
     * Variable for true particle energy starting kinetic energy.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the starting kinetic energy of the particle.
    */
    template<class T>
        double ke_init(const T & particle)
        {
            double energy(particle.energy_init);
            switch (particle.pid)
            {
            case 2:
                energy -= MUON_MASS;
                break;
            case 3:
                energy -= PION_MASS;
                break;
            case 4:
                energy -= PROTON_MASS;
                break;
            default:
                break;
            }
            return energy;
        }
    
    /**
     * Variable for particle overlap (IoU) of best match.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the particle overlap of the best match (IoU).
    */
    template<class T>
        double overlap(const T & particle) { return particle.match.size() > 0 ? (double)particle.match_overlap[0] : 0.0; }

    /**
     * Variable for the particles lowest x-coordinate.
     * @tparam T the type of particle (true or reco)
     * @param particle to apply the variable on.
     * @return the lowest x-coordinate of the particle start/end points.
    */
    template<class T>
        double lowx(const T & particle)
        {
            if(std::isnan(particle.start_point[0]) || std::isnan(particle.end_point[0]) || std::isinf(particle.start_point[0]) || std::isinf(particle.end_point[0]))
                return -100000.;
            else
                return std::min(particle.start_point[0], particle.end_point[0]);
        }
}

#endif